// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// Phase Field Fracture (PFF) main driver
// Uses the staggered solver implementation.
//
// This example reproduces the Mode 2 brittle fracture tutorial from Raccoon:
// /Users/chunhuizhao/projects/farms_cdms/raccoon/tutorials/mode2_brittle_fracture/
//
// Problem description:
// - Domain: [0,1] x [-0.5,0.5] with horizontal pre-crack at y=0, x < 0.5
// - Top boundary (y=0.5): u_x = t (shear), u_y = 0 (fixed)
// - Bottom boundary (y=-0.5): u_x = 0 (fixed), u_y = 0 (fixed)
// - Pre-crack: represented by initial damage field
// - Local refinement: 2 levels in oriented box along expected crack path
//
// Material parameters (Raccoon defaults):
// - E = 2.1e5, nu = 0.3
// - Gc = 2.7, l = 0.02
// - eta = 1e-6, p = 2 (AT2 model)

#include "mfem.hpp"
#include "pff_solver.hpp"

#include <memory>
#include <cmath>
#include <fstream>

using namespace mfem;
using namespace mfem::pff;

/** @brief Prescribed displacement coefficient for Mode 2 shear loading.
 *
 * Applies a time-dependent horizontal displacement ONLY on the top boundary.
 * Bottom boundary and elsewhere: u_x = 0, u_y = 0
 * Top boundary: u_x = amplitude * t, u_y = 0
 *
 * Uses a relative tolerance based on domain size for robust detection.
 */
class Mode2ShearLoad : public VectorCoefficient
{
   real_t amplitude_;
   real_t t_;
   real_t y_top_;   // y-coordinate of top boundary
   real_t y_bot_;   // y-coordinate of bottom boundary
public:
   Mode2ShearLoad(int dim, real_t amplitude, real_t y_top = 0.5, real_t y_bot = -0.5)
      : VectorCoefficient(dim), amplitude_(amplitude), t_(0.0),
        y_top_(y_top), y_bot_(y_bot) {}

   void SetTime(real_t t) override { t_ = t; }

   void Eval(Vector &V, ElementTransformation &T,
             const IntegrationPoint &ip) override
   {
      V.SetSize(vdim);
      V = 0.0;  // Default: zero displacement

      // Get physical coordinates
      Vector x;
      T.Transform(ip, x);

      // Use relative tolerance (1% of domain height)
      real_t domain_height = y_top_ - y_bot_;
      real_t tol = 0.01 * domain_height;

      // Only apply shear displacement at TOP boundary (y close to y_top_)
      if (x(1) > y_top_ - tol)
      {
         V(0) = amplitude_ * t_;  // u_x = amplitude * t (shear)
         V(1) = 0.0;              // u_y = 0 (fixed vertical)
      }
      // Bottom boundary (y close to y_bot_): u_x = 0, u_y = 0 (already set)
   }

   void SetAmplitude(real_t a) { amplitude_ = a; }
   void SetTopY(real_t y) { y_top_ = y; }
   void SetBotY(real_t y) { y_bot_ = y; }
};

/** @brief Check if a point is inside an oriented box.
 *
 * The oriented box is defined by:
 * - center: center point of the box
 * - length_dir: direction of the long axis (will be normalized)
 * - width_dir: direction of the short axis (will be normalized)
 * - half_length: half the length along length_dir
 * - half_width: half the width along width_dir
 *
 * @param[in] pt Point to check (2D)
 * @param[in] center Center of the box
 * @param[in] length_dir Direction of long axis
 * @param[in] width_dir Direction of short axis
 * @param[in] half_length Half-length
 * @param[in] half_width Half-width
 * @return true if point is inside the box
 */
bool PointInOrientedBox(const Vector &pt, const Vector &center,
                        const Vector &length_dir, const Vector &width_dir,
                        real_t half_length, real_t half_width)
{
   // Vector from center to point
   Vector diff(2);
   diff(0) = pt(0) - center(0);
   diff(1) = pt(1) - center(1);

   // Normalize directions
   real_t len_norm = std::sqrt(length_dir(0)*length_dir(0) + length_dir(1)*length_dir(1));
   real_t wid_norm = std::sqrt(width_dir(0)*width_dir(0) + width_dir(1)*width_dir(1));

   Vector len_hat(2), wid_hat(2);
   len_hat(0) = length_dir(0) / len_norm;
   len_hat(1) = length_dir(1) / len_norm;
   wid_hat(0) = width_dir(0) / wid_norm;
   wid_hat(1) = width_dir(1) / wid_norm;

   // Project onto axes
   real_t proj_len = diff(0)*len_hat(0) + diff(1)*len_hat(1);
   real_t proj_wid = diff(0)*wid_hat(0) + diff(1)*wid_hat(1);

   // Check if inside
   return (std::abs(proj_len) <= half_length && std::abs(proj_wid) <= half_width);
}

/** @brief Mark elements for refinement in an oriented box region.
 *
 * Matches Raccoon's OrientedBoxMarker:
 * - center = '0.65 -0.25 0'
 * - length = 0.8, width = 0.2
 * - length_direction = '1 -1.5 0'
 * - width_direction = '1.5 1 0'
 *
 * @param[in] mesh The mesh to mark
 * @param[out] refinement_marker Array to mark elements (1 = refine)
 */
void MarkElementsForRefinement(Mesh &mesh, Array<int> &refinement_marker)
{
   // Raccoon OrientedBoxMarker parameters
   Vector center(2);
   center(0) = 0.65;
   center(1) = -0.25;

   Vector length_dir(2), width_dir(2);
   length_dir(0) = 1.0;
   length_dir(1) = -1.5;
   width_dir(0) = 1.5;
   width_dir(1) = 1.0;

   real_t half_length = 0.8 / 2.0;  // length = 0.8
   real_t half_width = 0.2 / 2.0;   // width = 0.2

   refinement_marker.SetSize(mesh.GetNE());
   refinement_marker = 0;

   Vector centroid(2);
   for (int i = 0; i < mesh.GetNE(); i++)
   {
      // Get element centroid
      mesh.GetElementCenter(i, centroid);

      // Mark if centroid is inside the oriented box
      if (PointInOrientedBox(centroid, center, length_dir, width_dir,
                             half_length, half_width))
      {
         refinement_marker[i] = 1;
      }
   }
}

/** @brief Apply local refinement matching Raccoon's adaptive mesh.
 *
 * Performs initial_steps levels of refinement in the oriented box region.
 * Uses triangle mesh for proper conforming refinement.
 *
 * @param[in,out] mesh Mesh to refine
 * @param[in] initial_steps Number of refinement passes (Raccoon uses 2)
 */
void ApplyLocalRefinement(Mesh &mesh, int initial_steps)
{
   for (int step = 0; step < initial_steps; step++)
   {
      Array<int> refinement_marker;
      MarkElementsForRefinement(mesh, refinement_marker);

      // Count elements to refine
      int count = 0;
      for (int i = 0; i < refinement_marker.Size(); i++)
      {
         if (refinement_marker[i]) { count++; }
      }

      if (count == 0) { break; }

      // Build refinement array
      Array<Refinement> refs;
      for (int i = 0; i < refinement_marker.Size(); i++)
      {
         if (refinement_marker[i])
         {
            refs.Append(Refinement(i));  // Default refinement type
         }
      }

      // Apply refinement
      mesh.GeneralRefinement(refs, -1, true);
   }
}

/** @brief Create mesh for Mode 2 fracture with STITCHED geometry (like RACCOON).
 *
 * Creates a mesh on [0, Lx] x [-Ly/2, Ly/2] with a physical crack (gap)
 * at y=0 for x < crack_tip. This is done by creating two half-meshes and
 * only connecting them where x >= crack_tip.
 *
 * This matches RACCOON's StitchedMeshGenerator approach where:
 * - Top half and bottom half are separate meshes
 * - They're only stitched at x >= crack_tip along y=0
 * - The crack region (x < crack_tip) has no connectivity
 *
 * @param[in] nx Number of elements in x (30)
 * @param[in] ny_half Number of elements in y per half (15)
 * @param[in] Lx Domain length (1.0)
 * @param[in] Ly_half Domain half-height (0.5)
 * @param[in] crack_tip X-coordinate where stitching begins (0.5)
 * @return Mesh pointer with physical crack gap
 */
Mesh* CreateMode2MeshStitched(int nx, int ny_half, real_t Lx, real_t Ly_half,
                               real_t crack_tip)
{
   // Calculate grid spacing
   real_t dx = Lx / nx;
   real_t dy = Ly_half / ny_half;

   // Find the x-index where crack tip is (rounded to nearest node)
   int crack_idx = static_cast<int>(std::round(crack_tip / dx));

   // Total vertices: two halves, with shared vertices only at x >= crack_tip
   // Top half: (nx+1) * (ny_half+1) vertices
   // Bottom half: (nx+1) * (ny_half+1) vertices
   // Shared at y=0 for x >= crack_tip: (nx+1 - crack_idx) vertices
   int nv_half = (nx + 1) * (ny_half + 1);
   int nv_shared = (nx + 1) - crack_idx;  // Vertices shared at y=0
   int nv_total = 2 * nv_half - nv_shared;

   // Total elements: nx * ny_half per half = nx * 2*ny_half total
   int ne_total = nx * 2 * ny_half;

   // Boundary edges: top, bottom, left, right + crack faces
   // - Top: nx edges at y=Ly_half (attribute 3)
   // - Bottom: nx edges at y=-Ly_half (attribute 1)
   // - Left: 2*ny_half edges at x=0 (attribute 4)
   // - Right: 2*ny_half edges at x=Lx (attribute 2)
   // - Crack faces: crack_idx edges on top face, crack_idx on bottom face (attribute 5)
   int nb_total = 2*nx + 4*ny_half + 2*crack_idx;

   // Create mesh with boundary
   Mesh *mesh = new Mesh(2, nv_total, ne_total, nb_total);

   // Add vertices
   // Top half vertices: y from 0 to Ly_half
   // Indices: 0 to nv_half-1
   for (int j = 0; j <= ny_half; j++)
   {
      for (int i = 0; i <= nx; i++)
      {
         real_t x = i * dx;
         real_t y = j * dy;  // y from 0 to Ly_half
         mesh->AddVertex(x, y);
      }
   }

   // Bottom half vertices: y from -Ly_half to 0
   // For x < crack_tip at y=0: add separate vertices (crack gap)
   // For x >= crack_tip at y=0: use shared vertices from top half
   // Indices: nv_half to nv_total-1
   for (int j = 0; j <= ny_half; j++)
   {
      for (int i = 0; i <= nx; i++)
      {
         real_t x = i * dx;
         real_t y_bottom = -Ly_half + j * dy;  // y from -Ly_half to 0

         // At j=ny_half, y=0 (interface)
         if (j == ny_half)
         {
            // Only add vertex if x < crack_tip (unstitched region)
            // Otherwise, this vertex is shared with top half
            if (i < crack_idx)
            {
               mesh->AddVertex(x, 0.0);
            }
            // Else: shared vertex, don't add
         }
         else
         {
            mesh->AddVertex(x, y_bottom);
         }
      }
   }

   // Helper to get vertex index
   // Top half: vertex at (i,j) where j=0 is y=0, j=ny_half is y=Ly_half
   auto top_vertex = [nx](int i, int j) {
      return j * (nx + 1) + i;
   };

   // Bottom half: vertex at (i,j) where j=0 is y=-Ly_half, j=ny_half is y=0
   auto bottom_vertex = [nx, ny_half, nv_half, crack_idx, &top_vertex](int i, int j) {
      if (j == ny_half)  // At y=0 interface
      {
         if (i < crack_idx)
         {
            // Unstitched: use bottom's own vertex
            // These are added after all non-interface vertices
            return nv_half + ny_half * (nx + 1) + i;
         }
         else
         {
            // Stitched: use top's vertex at y=0 (j=0 in top half)
            return top_vertex(i, 0);
         }
      }
      else
      {
         // Non-interface: sequential after top half
         return nv_half + j * (nx + 1) + i;
      }
   };

   // Add elements
   // Top half elements
   for (int j = 0; j < ny_half; j++)
   {
      for (int i = 0; i < nx; i++)
      {
         // Quad with vertices ordered counterclockwise
         int v0 = top_vertex(i, j);
         int v1 = top_vertex(i + 1, j);
         int v2 = top_vertex(i + 1, j + 1);
         int v3 = top_vertex(i, j + 1);
         mesh->AddQuad(v0, v1, v2, v3, 1);  // Attribute 1
      }
   }

   // Bottom half elements
   for (int j = 0; j < ny_half; j++)
   {
      for (int i = 0; i < nx; i++)
      {
         int v0 = bottom_vertex(i, j);
         int v1 = bottom_vertex(i + 1, j);
         int v2 = bottom_vertex(i + 1, j + 1);
         int v3 = bottom_vertex(i, j + 1);
         mesh->AddQuad(v0, v1, v2, v3, 1);  // Attribute 1
      }
   }

   // Add boundary elements
   // Boundary attributes:
   // 1 = bottom (y = -Ly_half)
   // 2 = right (x = Lx)
   // 3 = top (y = Ly_half)
   // 4 = left (x = 0)
   // 5 = crack faces (internal at y=0 for x < crack_tip)

   // Bottom boundary (y = -Ly_half)
   for (int i = 0; i < nx; i++)
   {
      int v0 = bottom_vertex(i, 0);
      int v1 = bottom_vertex(i + 1, 0);
      mesh->AddBdrSegment(v0, v1, 1);
   }

   // Top boundary (y = Ly_half)
   for (int i = 0; i < nx; i++)
   {
      int v0 = top_vertex(i, ny_half);
      int v1 = top_vertex(i + 1, ny_half);
      mesh->AddBdrSegment(v0, v1, 3);
   }

   // Left boundary (x = 0): top half then bottom half
   for (int j = 0; j < ny_half; j++)
   {
      int v0 = top_vertex(0, j);
      int v1 = top_vertex(0, j + 1);
      mesh->AddBdrSegment(v0, v1, 4);
   }
   for (int j = 0; j < ny_half; j++)
   {
      int v0 = bottom_vertex(0, j);
      int v1 = bottom_vertex(0, j + 1);
      mesh->AddBdrSegment(v0, v1, 4);
   }

   // Right boundary (x = Lx): top half then bottom half
   for (int j = 0; j < ny_half; j++)
   {
      int v0 = top_vertex(nx, j);
      int v1 = top_vertex(nx, j + 1);
      mesh->AddBdrSegment(v0, v1, 2);
   }
   for (int j = 0; j < ny_half; j++)
   {
      int v0 = bottom_vertex(nx, j);
      int v1 = bottom_vertex(nx, j + 1);
      mesh->AddBdrSegment(v0, v1, 2);
   }

   // Crack faces at y=0 for x < crack_tip (attribute 5)
   // Top crack face (bottom edge of top half elements at y=0)
   for (int i = 0; i < crack_idx; i++)
   {
      int v0 = top_vertex(i, 0);
      int v1 = top_vertex(i + 1, 0);
      mesh->AddBdrSegment(v0, v1, 5);
   }
   // Bottom crack face (top edge of bottom half elements at y=0)
   for (int i = 0; i < crack_idx; i++)
   {
      int v0 = bottom_vertex(i, ny_half);
      int v1 = bottom_vertex(i + 1, ny_half);
      mesh->AddBdrSegment(v0, v1, 5);
   }

   // Finalize mesh
   mesh->FinalizeQuadMesh(true);

   return mesh;
}

/** @brief Create simple continuous mesh for Mode 2 fracture.
 *
 * Creates a mesh on [0, Lx] x [-Ly/2, Ly/2].
 * Matches Raccoon: 30x15 per half = 30x30 total.
 * NOTE: This creates a CONTINUOUS mesh without physical crack gap.
 *
 * @param[in] nx Number of elements in x (30)
 * @param[in] ny_half Number of elements in y per half (15)
 * @param[in] Lx Domain length (1.0)
 * @param[in] Ly_half Domain half-height (0.5)
 * @return Mesh pointer
 */
Mesh* CreateMode2Mesh(int nx, int ny_half, real_t Lx, real_t Ly_half)
{
   // Create mesh: [0, Lx] x [0, 2*Ly_half]
   Mesh *mesh = new Mesh(Mesh::MakeCartesian2D(
      nx, 2*ny_half, Element::QUADRILATERAL, true, Lx, 2*Ly_half));

   // Shift y from [0, 2*Ly_half] to [-Ly_half, Ly_half]
   for (int i = 0; i < mesh->GetNV(); i++)
   {
      real_t *v = mesh->GetVertex(i);
      v[1] -= Ly_half;
   }

   return mesh;
}

/** @brief Initialize damage field with pre-existing crack.
 *
 * Sets damage along the crack region (y ≈ 0, x < crack_tip).
 * Uses a smooth profile consistent with phase field regularization.
 *
 * The crack is a horizontal line at y=0 from x=0 to x=crack_tip.
 * Damage profile: d = exp(-(y/l)^2) for x < crack_tip, with smooth tip.
 *
 * @param[in,out] d Damage GridFunction to initialize
 * @param[in] crack_tip X-coordinate of crack tip
 * @param[in] l Phase field regularization length
 */
void InitializeCrack(ParGridFunction &d, real_t crack_tip, real_t l)
{
   // Crack coefficient: d = 1 at crack centerline (y=0), decaying over length l
   class CrackCoefficient : public Coefficient
   {
      real_t crack_tip_, l_;
   public:
      CrackCoefficient(real_t ct, real_t ll) : crack_tip_(ct), l_(ll) {}
      real_t Eval(ElementTransformation &T, const IntegrationPoint &ip) override
      {
         Vector x;
         T.Transform(ip, x);

         // Only apply damage for x < crack_tip (the crack region)
         // Beyond crack tip, no damage
         if (x(0) > crack_tip_)
         {
            return 0.0;
         }

         // Distance from crack line (y=0)
         real_t dist_y = std::abs(x(1));

         // Gaussian profile perpendicular to crack: d = exp(-(y/l)^2)
         // This gives d=1 at y=0, decaying to ~0 at |y| > 2*l
         real_t d_val = std::exp(-dist_y * dist_y / (l_ * l_));

         // Smooth transition near crack tip
         if (x(0) > crack_tip_ - 2.0 * l_)
         {
            // Smooth ramp from 1 to 0 as x approaches crack_tip
            real_t x_factor = 0.5 * (1.0 - std::tanh(3.0 * (x(0) - crack_tip_ + l_) / l_));
            d_val *= x_factor;
         }

         return d_val;
      }
   };

   CrackCoefficient crack_coeff(crack_tip, l);
   d.ProjectCoefficient(crack_coeff);

   // Clamp to [0, 1]
   for (int i = 0; i < d.Size(); i++)
   {
      d(i) = std::max(0.0, std::min(1.0, d(i)));
   }
}

int main(int argc, char *argv[])
{
   // Initialize MPI
   Mpi::Init(argc, argv);
   int myid = Mpi::WorldRank();

   // Default parameters matching Raccoon mode2_brittle_fracture tutorial
   const char *mesh_file = "";
   int order = 1;
   int ref_levels = 0;       // Uniform refinement levels
   int local_ref_steps = 2;  // Local refinement steps (Raccoon: initial_steps = 2)

   // Material parameters (Raccoon defaults)
   real_t E = 2.1e5;       // Young's modulus
   real_t nu = 0.3;        // Poisson's ratio
   real_t l = 0.02;        // Regularization length
   real_t Gc = 2.7;        // Fracture toughness

   // Mesh parameters matching Raccoon: 30x15 per half
   int nx = 30;            // Elements in x-direction
   int ny_half = 15;       // Elements in y-direction per half (total = 2*ny_half)
   real_t Lx = 1.0;        // Domain length [0, 1]
   real_t Ly_half = 0.5;   // Domain half-height [-0.5, 0.5]
   real_t crack_tip = 0.5; // Pre-crack tip x-coordinate

   // Time/load parameters (Raccoon defaults for Mode 2)
   real_t t_final = 2e-2;  // 20 ms (Raccoon: end_time = 2e-2)
   real_t dt = 2e-5;       // Time step (Raccoon: dt = 2e-5)
   real_t load_amp = 1;  // u_x = load_amp * t (shear)

   // Solver parameters (Raccoon defaults)
   int max_fp_iter = 20;
   real_t fp_rel_tol = 1e-8;
   real_t fp_abs_tol = 1e-10;

   // Output parameters
   int vis_steps = 10;
   int print_level = 1;
   bool output_csv = true;
   bool init_crack = true;
   bool local_refine = true;  // Apply local refinement
   bool stitched_mesh = false; // Use stitched mesh (physical crack gap) like RACCOON - DISABLED for now

   // Parse command line options
   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
   args.AddOption(&order, "-o", "--order", "Finite element order (default 1).");
   args.AddOption(&ref_levels, "-r", "--refine", "Uniform mesh refinements.");
   args.AddOption(&local_ref_steps, "-lr", "--local-refine",
                  "Local refinement steps near crack (Raccoon: 2).");
   args.AddOption(&E, "--E", "--youngs", "Young's modulus.");
   args.AddOption(&nu, "--nu", "--poisson", "Poisson's ratio.");
   args.AddOption(&l, "--l", "--length", "Regularization length.");
   args.AddOption(&Gc, "--Gc", "--fracture-toughness", "Fracture toughness.");
   args.AddOption(&nx, "--nx", "--num-x", "Number of elements in x (30).");
   args.AddOption(&ny_half, "--ny", "--num-y-half",
                  "Number of elements in y per half (15).");
   args.AddOption(&crack_tip, "--crack", "--crack-tip",
                  "Pre-crack tip x-coordinate.");
   args.AddOption(&t_final, "-tf", "--t-final", "Final time.");
   args.AddOption(&dt, "-dt", "--time-step", "Time step.");
   args.AddOption(&load_amp, "-la", "--load-amplitude",
                  "Shear load amplitude (displacement per unit time).");
   args.AddOption(&max_fp_iter, "--fp-max", "--fp-max-iter",
                  "Max staggered iterations.");
   args.AddOption(&fp_rel_tol, "--fp-rtol", "--fp-rel-tol",
                  "Staggered relative tolerance.");
   args.AddOption(&fp_abs_tol, "--fp-atol", "--fp-abs-tol",
                  "Staggered absolute tolerance.");
   args.AddOption(&vis_steps, "-vs", "--vis-steps", "Visualization steps.");
   args.AddOption(&print_level, "-pl", "--print-level",
                  "Print level (0=silent, 1=summary, 2=verbose).");
   args.AddOption(&output_csv, "-csv", "--output-csv", "-no-csv", "--no-csv",
                  "Output force-displacement CSV file.");
   args.AddOption(&init_crack, "-ic", "--init-crack", "-no-ic", "--no-init-crack",
                  "Initialize damage with pre-existing crack.");
   args.AddOption(&local_refine, "-lref", "--local-refine",
                  "-no-lref", "--no-local-refine",
                  "Apply local refinement near crack path.");
   args.AddOption(&stitched_mesh, "-stitch", "--stitched-mesh",
                  "-no-stitch", "--no-stitched-mesh",
                  "Use stitched mesh with physical crack gap (like RACCOON).");

   args.Parse();
   if (!args.Good())
   {
      if (myid == 0) { args.PrintUsage(std::cout); }
      return 1;
   }
   if (myid == 0) { args.PrintOptions(std::cout); }

   // Create mesh
   std::unique_ptr<Mesh> serial_mesh;
   if (mesh_file && mesh_file[0] != '\0')
   {
      serial_mesh = std::make_unique<Mesh>(mesh_file, 1, 1);
   }
   else if (stitched_mesh)
   {
      // Create STITCHED mesh with physical crack gap (like RACCOON)
      // This creates a mesh where nodes along y=0 for x < crack_tip are NOT connected
      if (myid == 0)
      {
         mfem::out << "\nUsing STITCHED mesh with physical crack gap at y=0, x < "
                   << crack_tip << "\n";
      }
      serial_mesh.reset(CreateMode2MeshStitched(nx, ny_half, Lx, Ly_half, crack_tip));
      // With stitched mesh, we don't need to initialize damage at the crack
      // The crack is already represented by the mesh geometry
      init_crack = false;
   }
   else
   {
      // Create continuous mesh (requires damage initialization for crack)
      serial_mesh.reset(CreateMode2Mesh(nx, ny_half, Lx, Ly_half));
   }

   // Uniform refinement (before local refinement)
   for (int lev = 0; lev < ref_levels; lev++)
   {
      serial_mesh->UniformRefinement();
   }

   // Local refinement in oriented box (matching Raccoon's OrientedBoxMarker)
   if (local_refine && local_ref_steps > 0)
   {
      if (myid == 0)
      {
         mfem::out << "\nApplying " << local_ref_steps
                   << " local refinement steps in crack region...\n";
         mfem::out << "  Initial elements: " << serial_mesh->GetNE() << "\n";
      }
      ApplyLocalRefinement(*serial_mesh, local_ref_steps);
      if (myid == 0)
      {
         mfem::out << "  Final elements: " << serial_mesh->GetNE() << "\n";
      }
   }

   ParMesh pmesh(MPI_COMM_WORLD, *serial_mesh);
   serial_mesh.reset();
   int dim = pmesh.Dimension();

   // Debug: check boundary attributes
   if (myid == 0)
   {
      mfem::out << "\nBoundary attributes on mesh:\n";
      mfem::out << "  Max bdr attribute: " << pmesh.bdr_attributes.Max() << "\n";
      // Count elements per attribute
      Array<int> bdr_count(pmesh.bdr_attributes.Max());
      bdr_count = 0;
      for (int be = 0; be < pmesh.GetNBE(); be++)
      {
         int attr = pmesh.GetBdrAttribute(be);
         if (attr >= 1 && attr <= bdr_count.Size())
         {
            bdr_count[attr - 1]++;
         }
      }
      for (int i = 0; i < bdr_count.Size(); i++)
      {
         mfem::out << "  Attribute " << (i + 1) << ": " << bdr_count[i] << " elements\n";
      }
   }

   // Material parameters with explicit Gc
   PFFMaterialParameters mat(E, nu, l, Gc, 1e-6, 2.0, true);
   if (myid == 0)
   {
      mfem::out << "\n";
      mat.Print(mfem::out);
   }

   // Create solver
   PFFSolver solver(pmesh, mat, order);
   solver.SetPrintLevel(print_level);
   solver.SetMaxIterations(max_fp_iter);
   solver.SetRelativeTolerance(fp_rel_tol);
   solver.SetAbsoluteTolerance(fp_abs_tol);

   // Configure Mode 2 (shear) fracture boundary conditions
   solver.ConfigureMode2FractureBCs();

   // Prescribed shear displacement - only apply at TOP boundary (y = Ly_half)
   // Bottom boundary (y = -Ly_half) gets zero displacement
   Mode2ShearLoad load(dim, load_amp, Ly_half, -Ly_half);
   solver.SetPrescribedDisplacement(load);

   // Initialize pre-existing crack in damage field
   if (init_crack)
   {
      InitializeCrack(solver.GetDamage(), crack_tip, l);
      real_t d_max_local = solver.GetDamage().Max();
      real_t d_max;
      MPI_Allreduce(&d_max_local, &d_max, 1, MPI_DOUBLE, MPI_MAX, pmesh.GetComm());
      if (myid == 0)
      {
         mfem::out << "Pre-existing crack initialized (d_max = " << d_max << ")\n";
      }

      // Debug: Check initial damage at right boundary
      const ParGridFunction &d = solver.GetDamage();
      const ParFiniteElementSpace &d_fes = solver.GetDamageFES();
      int d_ndof = d_fes.GetNDofs();

      real_t corner_tol = 0.1;
      real_t d_right_max = 0.0;
      real_t right_x = 0, right_y = 0;

      // For scalar FE space, get node coordinates from mesh vertices
      // For order 1, DOFs are at vertices
      for (int i = 0; i < d_ndof && i < pmesh.GetNV(); i++)
      {
         const real_t *vc = pmesh.GetVertex(i);
         real_t x = vc[0];
         real_t y = vc[1];
         real_t d_val = d(i);

         // Check right boundary region
         if (x > Lx - corner_tol && d_val > d_right_max)
         {
            d_right_max = d_val;
            right_x = x;
            right_y = y;
         }
      }
      if (myid == 0)
      {
         mfem::out << "Initial damage at right boundary: d=" << d_right_max
                   << " at (" << right_x << ", " << right_y << ")\n";
      }
   }

   // ParaView output - include strain energy for debugging
   ParaViewDataCollection pv_dc("pff_mode2", &pmesh);
   pv_dc.SetDataFormat(VTKFormat::BINARY);
   pv_dc.SetHighOrderOutput(true);
   pv_dc.SetLevelsOfDetail(order);
   pv_dc.RegisterField("displacement", &solver.GetDisplacement());
   pv_dc.RegisterField("damage", &solver.GetDamage());

   // CSV output for force-displacement curve
   std::ofstream csv_file;
   if (output_csv && myid == 0)
   {
      csv_file.open("pff_mode2_force_displacement.csv");
      csv_file << "time,displacement,reaction_force,max_damage\n";
   }

   // Time stepping
   real_t t = 0.0;
   int step = 0;
   int total_steps = static_cast<int>(std::ceil(t_final / dt));

   // Call global functions on all processes before printing (MPI collective ops)
   long long global_ne = pmesh.GetGlobalNE();
   HYPRE_BigInt global_u_dofs = solver.GetDisplacementFES().GlobalTrueVSize();

   if (myid == 0)
   {
      mfem::out << "\n=== Starting Mode 2 (shear) fracture simulation ===\n"
                << "  Domain: [0, " << Lx << "] x [" << -Ly_half << ", " << Ly_half << "]\n"
                << "  Base mesh: " << nx << " x " << 2*ny_half << " elements\n"
                << "  Local refinement: " << local_ref_steps << " steps\n"
                << "  Total mesh elements: " << global_ne << "\n"
                << "  Displacement DOFs: " << global_u_dofs << "\n"
                << "  Pre-crack: y=0, x < " << crack_tip << "\n"
                << "  Loading: shear u_x = " << load_amp << " * t on top\n"
                << "  t_final = " << t_final << ", dt = " << dt << "\n"
                << "  Total steps: " << total_steps << "\n\n";
   }

   // Initial output
   pv_dc.SetCycle(0);
   pv_dc.SetTime(0.0);
   pv_dc.Save();

   while (t < t_final - 1e-14)
   {
      step++;

      // Update load time
      load.SetTime(t + dt);

      if (myid == 0 && print_level > 0)
      {
         mfem::out << "--- Step " << step << "/" << total_steps
                   << ", t = " << t + dt
                   << ", u_x(top) = " << load_amp * (t + dt) << " ---\n";
      }

      // Solve time step
      int fp_iters = solver.SolveTimeStep(dt);

      // Get damage statistics
      const ParGridFunction &d = solver.GetDamage();
      real_t d_max_local = d.Max();
      real_t d_max;
      MPI_Allreduce(&d_max_local, &d_max, 1, MPI_DOUBLE, MPI_MAX,
                    pmesh.GetComm());

      // Get strain energy statistics
      const ParGridFunction &psi = solver.GetStrainEnergy();
      real_t psi_max_local = psi.Max();
      real_t psi_max;
      MPI_Allreduce(&psi_max_local, &psi_max, 1, MPI_DOUBLE, MPI_MAX,
                    pmesh.GetComm());

      const ParGridFunction &H = solver.GetStrainEnergyHistory();
      real_t H_max_local = H.Max();
      real_t H_max;
      MPI_Allreduce(&H_max_local, &H_max, 1, MPI_DOUBLE, MPI_MAX,
                    pmesh.GetComm());

      // Critical strain energy: Gc / (c0 * l)
      real_t psi_crit = mat.Gc / (mat.c0 * mat.l);
      if (myid == 0 && print_level > 0)
      {
         mfem::out << "  Strain energy: psi_max = " << psi_max
                   << ", H_max = " << H_max
                   << ", psi_crit = " << psi_crit << "\n";
      }

      // Debug: check displacement at top and bottom boundaries
      // Use coordinate projection approach for robustness
      const ParGridFunction &u = solver.GetDisplacement();
      const ParFiniteElementSpace &u_fes = solver.GetDisplacementFES();
      int ndof = u_fes.GetNDofs();  // Scalar DOFs per component

      // Project coordinates to get node positions
      ParGridFunction coords(const_cast<ParFiniteElementSpace*>(&u_fes));
      VectorFunctionCoefficient coord_coeff(dim,
         [](const Vector &x, Vector &v) { v = x; });
      coords.ProjectCoefficient(coord_coeff);

      // Debug: Check damage at corners
      const ParFiniteElementSpace &d_fes = solver.GetDamageFES();
      int d_ndof = d_fes.GetNDofs();

      real_t corner_tol = 0.05;  // Tolerance for corner detection
      real_t d_corner_max = 0.0;
      real_t corner_x = 0, corner_y = 0;

      // For scalar FE space, get node coordinates from mesh vertices
      for (int i = 0; i < d_ndof && i < pmesh.GetNV(); i++)
      {
         const real_t *vc = pmesh.GetVertex(i);
         real_t x = vc[0];
         real_t y = vc[1];
         real_t d_val = d(i);

         // Check if near right boundary (x close to 1)
         if (x > Lx - corner_tol)
         {
            if (d_val > d_corner_max)
            {
               d_corner_max = d_val;
               corner_x = x;
               corner_y = y;
            }
         }
      }

      // Print corner damage info
      if (myid == 0 && print_level > 0 && d_corner_max > 0.01)
      {
         mfem::out << "  WARNING: Damage at right boundary: d=" << d_corner_max
                   << " at (" << corner_x << ", " << corner_y << ")\n";
      }

      // Use tight tolerance to only get actual boundary nodes
      real_t bdr_tol = 1e-8;

      real_t u_top_x_max = -1e20, u_top_x_min = 1e20;
      real_t u_bot_x_max = -1e20, u_bot_x_min = 1e20;
      real_t u_top_y_max = -1e20, u_top_y_min = 1e20;
      real_t u_bot_y_max = -1e20, u_bot_y_min = 1e20;
      for (int i = 0; i < ndof; i++)
      {
         // byNODES ordering: y-coordinate at node i is at index ndof + i
         real_t y = coords(ndof + i);
         // x-component of displacement at node i is at index i (byNODES)
         real_t ux = u(i);
         // y-component of displacement at node i is at index ndof + i (byNODES)
         real_t uy = u(ndof + i);

         if (std::abs(y - Ly_half) < bdr_tol)  // Top boundary
         {
            u_top_x_max = std::max(u_top_x_max, ux);
            u_top_x_min = std::min(u_top_x_min, ux);
            u_top_y_max = std::max(u_top_y_max, uy);
            u_top_y_min = std::min(u_top_y_min, uy);
         }
         if (std::abs(y + Ly_half) < bdr_tol)  // Bottom boundary
         {
            u_bot_x_max = std::max(u_bot_x_max, ux);
            u_bot_x_min = std::min(u_bot_x_min, ux);
            u_bot_y_max = std::max(u_bot_y_max, uy);
            u_bot_y_min = std::min(u_bot_y_min, uy);
         }
      }

      // Also get global min/max of displacement field
      real_t u_max_local = u.Max();
      real_t u_min_local = u.Min();
      real_t u_max_global, u_min_global;
      MPI_Allreduce(&u_max_local, &u_max_global, 1, MPI_DOUBLE, MPI_MAX, pmesh.GetComm());
      MPI_Allreduce(&u_min_local, &u_min_global, 1, MPI_DOUBLE, MPI_MIN, pmesh.GetComm());

      if (myid == 0 && print_level > 0)
      {
         mfem::out << "  Staggered iters: " << fp_iters
                   << ", d_max = " << d_max << "\n";
         mfem::out << "  u_x at top: [" << u_top_x_min << ", " << u_top_x_max << "]"
                   << ", expected: " << load_amp * (t + dt) << "\n";
         mfem::out << "  u_y at top: [" << u_top_y_min << ", " << u_top_y_max << "]"
                   << ", expected: 0\n";
         mfem::out << "  u_x at bot: [" << u_bot_x_min << ", " << u_bot_x_max << "]"
                   << ", expected: 0\n";
         mfem::out << "  u_y at bot: [" << u_bot_y_min << ", " << u_bot_y_max << "]"
                   << ", expected: 0\n";
         mfem::out << "  Displacement field: min=" << u_min_global
                   << ", max=" << u_max_global << "\n";
      }

      // Reaction force placeholder
      real_t reaction_force = 0.0;

      // CSV output
      if (output_csv && myid == 0)
      {
         csv_file << t + dt << ","
                  << load_amp * (t + dt) << ","
                  << reaction_force << ","
                  << d_max << "\n";
         csv_file.flush();
      }

      // Save visualization
      if (vis_steps > 0 && (step % vis_steps == 0))
      {
         pv_dc.SetCycle(step);
         pv_dc.SetTime(t + dt);
         pv_dc.Save();

         if (myid == 0 && print_level > 0)
         {
            mfem::out << "  Saved visualization at step " << step << "\n";
         }
      }

      // Advance time
      solver.AdvanceTime(dt);
      t += dt;
   }

   // Final output
   pv_dc.SetCycle(step);
   pv_dc.SetTime(t);
   pv_dc.Save();

   if (output_csv && myid == 0)
   {
      csv_file.close();
   }

   if (myid == 0)
   {
      mfem::out << "\n=== Mode 2 simulation completed ===\n"
                << "  Total steps: " << step << "\n"
                << "  Final time: " << t << "\n"
                << "  Output saved to: pff_mode2/\n";
      if (output_csv)
      {
         mfem::out << "  CSV file: pff_mode2_force_displacement.csv\n";
      }
   }

   return 0;
}
