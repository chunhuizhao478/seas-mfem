// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// Phase Field Fracture (PFF) main driver - Mode 2 Pure Shear with Gmsh mesh
//
// This example reproduces the Mode 2 brittle fracture tutorial from Raccoon,
// but uses a Gmsh-generated mesh instead of procedural mesh generation.
//
// Problem description:
// - Domain: [0,1] x [-0.5,0.5] with horizontal pre-crack at y=0, x < 0.5
// - Top boundary (y=0.5): u_x = t (shear), u_y = 0 (fixed)
// - Bottom boundary (y=-0.5): u_x = 0 (fixed), u_y = 0 (fixed)
// - Pre-crack: represented by mesh geometry (stitched mesh) or initial damage
// - Local refinement: done in Gmsh geo file to match expected crack path
//
// Material parameters (Raccoon defaults):
// - E = 2.1e5, nu = 0.3
// - Gc = 2.7, l = 0.02
// - eta = 1e-6, p = 2 (AT2 model)
//
// Gmsh mesh boundary attributes (from mesh_quad.geo):
// - 1 = Bottom (y = -0.5)
// - 2 = Right (x = 1)
// - 3 = Top (y = 0.5)
// - 4 = Left (x = 0)
// - 5 = Notch faces (upper and lower faces of the physical notch)

#include "mfem.hpp"
#include "pff_solver.hpp"

#include <memory>
#include <cmath>
#include <fstream>
#include <iomanip>

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

/** @brief Check if the mesh has a physical notch/crack geometry.
 *
 * A mesh with a physical notch has boundary elements marking the notch faces,
 * which means the crack is represented by mesh geometry rather than damage field.
 * This is detected by checking if notch boundary attributes exist.
 *
 * @param[in] mesh The mesh to check
 * @param[in] notch_attr Boundary attribute for notch faces (default: 5)
 * @return true if notch boundary elements are found
 */
bool HasPhysicalNotch(Mesh &mesh, int notch_attr = 5)
{
   int max_attr = mesh.bdr_attributes.Max();
   if (max_attr < notch_attr) { return false; }

   // Count boundary elements with notch attribute
   int notch_count = 0;
   for (int be = 0; be < mesh.GetNBE(); be++)
   {
      int attr = mesh.GetBdrAttribute(be);
      if (attr == notch_attr)
      {
         notch_count++;
      }
   }

   return notch_count > 0;
}

int main(int argc, char *argv[])
{
   // Initialize MPI
   Mpi::Init(argc, argv);
   int myid = Mpi::WorldRank();

   // Default parameters matching Raccoon mode2_brittle_fracture tutorial
   const char *mesh_file = "../mesh/mesh_quad.msh";  // Default Gmsh mesh file
   int order = 1;
   int ref_levels = 0;       // Uniform refinement levels (after loading mesh)

   // Material parameters (Raccoon defaults)
   real_t E = 2.1e5;       // Young's modulus
   real_t nu = 0.3;        // Poisson's ratio
   real_t l = 0.02;        // Regularization length
   real_t Gc = 2.7;        // Fracture toughness

   // Domain parameters (should match Gmsh mesh)
   real_t Ly_half = 0.5;   // Domain half-height [-0.5, 0.5]
   real_t crack_tip = 0.5; // Pre-crack tip x-coordinate

   // Time/load parameters (Raccoon defaults for Mode 2)
   real_t t_final = 2e-2;  // 20 ms (Raccoon: end_time = 2e-2)
   real_t dt = 2e-5;       // Time step (Raccoon: dt = 2e-5)
   real_t load_amp = 1;    // u_x = load_amp * t (shear)

   // Solver parameters (Raccoon defaults)
   int max_fp_iter = 10;
   real_t fp_rel_tol = 1e-6;
   real_t fp_abs_tol = 1e-8;
   bool accept_on_max_fp_iter = false;  // Like MOOSE's accept_on_max_fixed_point_iteration

   // Adaptive time stepping parameters (like MOOSE's FarmsIterationAdaptiveDT)
   real_t dt_cutback_factor = 0.5;      // Factor to reduce dt on failure
   real_t dt_growth_factor = 1.25;      // Factor to grow dt on success
   real_t dt_min = 1e-10;               // Minimum allowed dt
   real_t dt_max = 0.0;                 // Maximum allowed dt (0 = use initial dt)
   int optimal_fp_iters = 5;            // Target number of fixed-point iterations

   // Output parameters
   int vis_steps = 10;
   int print_level = 1;
   bool output_csv = true;
   bool init_crack = true;  // Initialize damage with pre-existing crack (for continuous mesh)

   // Parse command line options
   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh", "Gmsh mesh file to use.");
   args.AddOption(&order, "-o", "--order", "Finite element order (default 1).");
   args.AddOption(&ref_levels, "-r", "--refine", "Uniform mesh refinements after loading.");
   args.AddOption(&E, "--E", "--youngs", "Young's modulus.");
   args.AddOption(&nu, "--nu", "--poisson", "Poisson's ratio.");
   args.AddOption(&l, "--l", "--length", "Regularization length.");
   args.AddOption(&Gc, "--Gc", "--fracture-toughness", "Fracture toughness.");
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
   args.AddOption(&accept_on_max_fp_iter, "--fp-accept", "--accept-on-max-fp",
                  "-no-fp-accept", "--no-accept-on-max-fp",
                  "Accept solution when max fixed-point iterations reached.");
   args.AddOption(&dt_cutback_factor, "--dt-cut", "--dt-cutback",
                  "Time step cutback factor on failure.");
   args.AddOption(&dt_growth_factor, "--dt-grow", "--dt-growth",
                  "Time step growth factor on success.");
   args.AddOption(&dt_min, "--dt-min", "--min-time-step",
                  "Minimum allowed time step.");
   args.AddOption(&dt_max, "--dt-max", "--max-time-step",
                  "Maximum allowed time step (0 = use initial dt).");
   args.AddOption(&optimal_fp_iters, "--fp-optimal", "--optimal-fp-iters",
                  "Target number of fixed-point iterations for adaptive dt.");
   args.AddOption(&vis_steps, "-vs", "--vis-steps", "Visualization steps.");
   args.AddOption(&print_level, "-pl", "--print-level",
                  "Print level (0=silent, 1=summary, 2=verbose).");
   args.AddOption(&output_csv, "-csv", "--output-csv", "-no-csv", "--no-csv",
                  "Output force-displacement CSV file.");
   args.AddOption(&init_crack, "-ic", "--init-crack", "-no-ic", "--no-init-crack",
                  "Initialize damage with pre-existing crack (for continuous mesh).");

   args.Parse();
   if (!args.Good())
   {
      if (myid == 0) { args.PrintUsage(std::cout); }
      return 1;
   }
   if (myid == 0) { args.PrintOptions(std::cout); }

   // Load Gmsh mesh
   if (myid == 0)
   {
      mfem::out << "\nLoading Gmsh mesh from: " << mesh_file << "\n";
   }

   std::unique_ptr<Mesh> serial_mesh;
   try
   {
      // Load mesh - MFEM automatically handles Gmsh format
      serial_mesh = std::make_unique<Mesh>(mesh_file, 1, 1);
   }
   catch (const std::exception &e)
   {
      if (myid == 0)
      {
         mfem::err << "Error loading mesh: " << e.what() << "\n";
         mfem::err << "Make sure to generate the mesh first:\n";
         mfem::err << "  cd ../mesh && gmsh -2 mesh_quad.geo -o mesh_quad.msh\n";
      }
      return 1;
   }

   int dim = serial_mesh->Dimension();
   if (myid == 0)
   {
      mfem::out << "  Mesh dimension: " << dim << "\n";
      mfem::out << "  Initial elements: " << serial_mesh->GetNE() << "\n";
      mfem::out << "  Initial vertices: " << serial_mesh->GetNV() << "\n";
   }

   // Check if mesh has physical notch
   bool has_physical_notch = HasPhysicalNotch(*serial_mesh);
   if (myid == 0)
   {
      if (has_physical_notch)
      {
         mfem::out << "  Mesh has PHYSICAL NOTCH (crack represented by mesh geometry)\n";
         mfem::out << "  Damage initialization will be skipped\n";
      }
      else
      {
         mfem::out << "  Mesh is CONTINUOUS (crack represented by damage field)\n";
      }
   }

   // If mesh has physical notch, don't initialize damage
   if (has_physical_notch)
   {
      init_crack = false;
   }

   // Uniform refinement (if requested)
   for (int lev = 0; lev < ref_levels; lev++)
   {
      serial_mesh->UniformRefinement();
   }

   if (ref_levels > 0 && myid == 0)
   {
      mfem::out << "  After " << ref_levels << " refinements: "
                << serial_mesh->GetNE() << " elements\n";
   }

   // Create parallel mesh
   ParMesh pmesh(MPI_COMM_WORLD, *serial_mesh);
   serial_mesh.reset();

   // Get actual domain bounds from mesh
   Vector bb_min, bb_max;
   pmesh.GetBoundingBox(bb_min, bb_max);
   real_t y_top = bb_max(1);
   real_t y_bot = bb_min(1);
   Ly_half = (y_top - y_bot) / 2.0;

   if (myid == 0)
   {
      mfem::out << "  Domain: [" << bb_min(0) << ", " << bb_max(0) << "] x ["
                << bb_min(1) << ", " << bb_max(1) << "]\n";
   }

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
   solver.SetAcceptOnMaxFixedPointIteration(accept_on_max_fp_iter);

   // Set dt_max if not specified
   if (dt_max <= 0.0) { dt_max = dt; }

   // Configure Mode 2 (shear) fracture boundary conditions
   solver.ConfigureMode2FractureBCs();

   // Prescribed shear displacement - only apply at TOP boundary (y = y_top)
   // Bottom boundary (y = y_bot) gets zero displacement
   Mode2ShearLoad load(dim, load_amp, y_top, y_bot);
   solver.SetPrescribedDisplacement(load);

   // Initialize pre-existing crack in damage field (for continuous mesh only)
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
   }

   // ParaView output
   ParaViewDataCollection pv_dc("pff_mode2_gmsh", &pmesh);
   pv_dc.SetDataFormat(VTKFormat::BINARY);
   pv_dc.SetHighOrderOutput(true);
   pv_dc.SetLevelsOfDetail(order);
   pv_dc.RegisterField("displacement", &solver.GetDisplacement());
   pv_dc.RegisterField("damage", &solver.GetDamage());

   // CSV output for force-displacement curve
   std::ofstream csv_file;
   if (output_csv && myid == 0)
   {
      csv_file.open("pff_mode2_gmsh_force_displacement.csv");
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
                << "  Domain: [" << bb_min(0) << ", " << bb_max(0) << "] x ["
                << bb_min(1) << ", " << bb_max(1) << "]\n"
                << "  Mesh type: Gmsh triangles\n"
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

   // Adaptive time stepping state
   real_t current_dt = dt;
   int failed_steps = 0;
   int consecutive_successes = 0;

   while (t < t_final - 1e-14)
   {
      step++;

      // Ensure we don't overshoot t_final
      if (t + current_dt > t_final)
      {
         current_dt = t_final - t;
      }

      // Update load time
      load.SetTime(t + current_dt);

      if (myid == 0 && print_level > 0)
      {
         mfem::out << "\nTime Step " << step << "/" << total_steps
                   << ", t = " << std::scientific << std::setprecision(4)
                   << t + current_dt << std::defaultfloat
                   << ", dt = " << std::scientific << current_dt << std::defaultfloat
                   << ", u_x(top) = " << load_amp * (t + current_dt) << "\n";
      }

      // Solve time step (returns convergence status, like MOOSE Picard iteration)
      int fp_iters = 0;
      bool converged = solver.SolveTimeStep(current_dt, &fp_iters);

      // Handle convergence failure with adaptive time stepping
      // (Like MOOSE's FarmsIterationAdaptiveDT with cutback_factor_at_failure)
      if (!converged && !accept_on_max_fp_iter)
      {
         failed_steps++;
         consecutive_successes = 0;

         // Cut back time step
         real_t new_dt = current_dt * dt_cutback_factor;

         if (new_dt < dt_min)
         {
            if (myid == 0)
            {
               mfem::err << "ERROR: Time step " << new_dt
                         << " below minimum " << dt_min << "\n";
               mfem::err << "Simulation failed at t = " << t << "\n";
            }
            MPI_Abort(MPI_COMM_WORLD, 1);
            return 1;
         }

         if (myid == 0 && print_level >= 0)
         {
            mfem::out << " Step FAILED - cutting dt from " << std::scientific
                      << current_dt << " to " << new_dt << std::defaultfloat << "\n";
         }

         current_dt = new_dt;
         step--;  // Redo this step
         continue;
      }

      // Get damage statistics
      const ParGridFunction &d = solver.GetDamage();
      real_t d_max_local = d.Max();
      real_t d_max;
      MPI_Allreduce(&d_max_local, &d_max, 1, MPI_DOUBLE, MPI_MAX,
                    pmesh.GetComm());

      // Print summary after solve
      if (myid == 0 && print_level > 0)
      {
         mfem::out << " Solve converged in " << fp_iters << " iterations, d_max = "
                   << std::scientific << std::setprecision(4) << d_max
                   << std::defaultfloat << "\n";
      }

      // Verbose output (print_level > 1): strain energy and boundary checks
      if (print_level > 1)
      {
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

         real_t psi_crit = mat.Gc / (mat.c0 * mat.l);
         if (myid == 0)
         {
            mfem::out << "  Strain energy: psi_max = " << psi_max
                      << ", H_max = " << H_max
                      << ", psi_crit = " << psi_crit << "\n";
         }

         // Debug: check displacement at boundaries
         const ParGridFunction &u = solver.GetDisplacement();
         const ParFiniteElementSpace &u_fes = solver.GetDisplacementFES();
         int ndof = u_fes.GetNDofs();

         ParGridFunction coords(const_cast<ParFiniteElementSpace*>(&u_fes));
         VectorFunctionCoefficient coord_coeff(dim,
            [](const Vector &x, Vector &v) { v = x; });
         coords.ProjectCoefficient(coord_coeff);

         real_t bdr_tol = 1e-8;
         real_t u_top_x_max = -1e20, u_top_x_min = 1e20;
         real_t u_bot_x_max = -1e20, u_bot_x_min = 1e20;

         for (int i = 0; i < ndof; i++)
         {
            real_t y = coords(ndof + i);
            real_t ux = u(i);

            if (std::abs(y - y_top) < bdr_tol)
            {
               u_top_x_max = std::max(u_top_x_max, ux);
               u_top_x_min = std::min(u_top_x_min, ux);
            }
            if (std::abs(y - y_bot) < bdr_tol)
            {
               u_bot_x_max = std::max(u_bot_x_max, ux);
               u_bot_x_min = std::min(u_bot_x_min, ux);
            }
         }

         if (myid == 0)
         {
            mfem::out << "  u_x at top: [" << u_top_x_min << ", " << u_top_x_max << "]"
                      << ", expected: " << load_amp * (t + current_dt) << "\n";
            mfem::out << "  u_x at bot: [" << u_bot_x_min << ", " << u_bot_x_max << "]"
                      << ", expected: 0\n";
         }
      }

      // Compute reaction force (shear force Fx on top boundary)
      real_t reaction_force = solver.ComputeShearReactionForce();

      if (myid == 0 && print_level > 0)
      {
         mfem::out << " Reaction force Fx(top) = " << std::scientific
                   << std::setprecision(6) << reaction_force
                   << std::defaultfloat << "\n";
      }

      // CSV output
      if (output_csv && myid == 0)
      {
         csv_file << t + current_dt << ","
                  << load_amp * (t + current_dt) << ","
                  << reaction_force << ","
                  << d_max << "\n";
         csv_file.flush();
      }

      // Save visualization
      if (vis_steps > 0 && (step % vis_steps == 0))
      {
         pv_dc.SetCycle(step);
         pv_dc.SetTime(t + current_dt);
         pv_dc.Save();

         if (myid == 0 && print_level > 0)
         {
            mfem::out << "  Saved visualization at step " << step << "\n";
         }
      }

      // Advance time
      solver.AdvanceTime(current_dt);
      t += current_dt;

      // Adaptive time stepping: grow dt if converged quickly
      // (Like MOOSE's FarmsIterationAdaptiveDT with optimal_iterations)
      consecutive_successes++;
      if (fp_iters <= optimal_fp_iters && consecutive_successes >= 2)
      {
         real_t new_dt = std::min(current_dt * dt_growth_factor, dt_max);
         if (new_dt > current_dt && myid == 0 && print_level > 0)
         {
            mfem::out << " Growing dt from " << std::scientific << current_dt
                      << " to " << new_dt << std::defaultfloat << "\n";
         }
         current_dt = new_dt;
      }
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
                << "  Failed steps: " << failed_steps << "\n"
                << "  Final time: " << t << "\n"
                << "  Final dt: " << current_dt << "\n"
                << "  Output saved to: pff_mode2_gmsh/\n";
      if (output_csv)
      {
         mfem::out << "  CSV file: pff_mode2_gmsh_force_displacement.csv\n";
      }
   }

   return 0;
}
