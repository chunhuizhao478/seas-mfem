// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// fault_frame_field.hpp — 2D in-plane design rev 3, phase A: ONE fault
// orientation system for the quasi-dynamic and the dynamic operators.
//
// The dynamic solver's orientation machinery (fault_orientation.hpp: the
// census J1, propagation O1-O3, one pin per patch O4/O5, the side test S3,
// partition-invariant keys P1; fault_frame.hpp: the F1/F5/F6 frame) runs
// ONCE per run on the fault faces, and once more, with the same pin, on the
// faces that carry a prescribed slip (interior Dirichlet faces: BP5's creep
// plane, BP3's creep segment).  Every ElasticityDomainOperator of the run
// (the QD phase's operator, and the Nitsche operator's elastic and
// Kelvin-Voigt copies) and the Nitsche operator itself read their
// orientation from this one object; the first operator that needs it builds
// it from its own face lists, the others receive it and check that their
// face lists are the same.
//
// What each consumer takes from it:
//   QD operator  the sign of every fault-face normal (the reference normal
//                handed to FaultBasis::ComputeOrientedFrame becomes the
//                face's oriented normal n), the up vector, and the sign of
//                every prescribed-slip face (ComputeSkeletonDirichletSign).
//                The frame arithmetic itself (per quadrature point and per
//                DOF, Tandem's sign-baked storage) is unchanged, so on any
//                mesh where n has the sign of the old constant reference
//                normal (0, -1, 0) and up = (0, 0, 1), every stored number
//                is bit for bit what it was before phase A.
//   Nitsche      the fault census itself (NitscheWaveOperator::Orientation)
//                and the N-09 frame of every fault face, computed by the
//                same calls the operator made before phase A.
//
// sigma_f (per fault face) = +1 when the operator's Elem1 is the S3 plus
// element (the oriented normal n points away from Elem1), -1 otherwise.
// The constructor asserts sigma_f = sign(n_raw . n), n_raw = CalcOrtho at
// the face centre (MFEM's Elem1-outward normal).  The QD kernel triad is then
// sigma_f (n, d, s) up to rounding, and the dynamic operator's eps_f equals
// sigma_f.
//
// Seam rule: where a fault face and a prescribed-slip face share an edge and
// no other face of either set does, the two censuses must orient the two
// faces consistently (O3: they run the shared edge in opposite directions).
// The field refuses the run otherwise, so the two censuses give the signs
// that one census over the union of the two sets would give.  A fault face
// and a prescribed-slip face meeting at a junction edge (three or more faces
// of the two sets) is refused too (junction_policy "refuse").
//
// 2D (phase B, Chunhui Zhao, 2026-10-25): on a triangle mesh the faces are
// segments and FaultOrientation runs its edge mode (tips and junctions are
// vertices).  The frame is the padded F1 rule of derivation 10.2, the same
// call as in 3D: up must lie in the mesh plane (BP3: up = [0, 1], the slip in
// the dip slot) or be +-e_z (TPV205-2D: up = [0, 0, 1], the slip in the
// strike slot), so that one tangential slot is an exact zero; the field
// refuses any other up, and a ref_normal or pin point off the plane.  The
// seam rule's shared entity is the vertex where a fault segment meets a
// prescribed-slip segment (BP3: x_d = W_f).
//
// Vertex normals (vertex-normal plan rev1, phase P1; Chunhui Zhao,
// 2026-11-01: "send this document to 2d thread and start implementation and
// test"; /mnt/project-files/2d-inplane/nonplanar-literature/
// vertex-normal-plan-rev1.pdf): FaultFrameConfig::averaged_vertex_normals
// ([fault_geometry] vertex_normals = "averaged"; default "face", off).  On a
// 2-D mesh every fault vertex v shared by exactly two fault faces A, B with
// turn angle dtheta_v = angle(n_A, n_B) <= vertex_normal_max_turn_deg gets
//    n_v = (n_A + n_B) / |n_A + n_B|                                  (1)
// (n_v = n_A exactly when n_A = n_B component by component, so a straight
// fault keeps its face normals bit for bit).  A vertex above the threshold
// (a kink), a tip (one face) and so the end of the fault at a prescribed-slip
// seam (the creep faces are not in the fault census) keep each face's own
// normal; three or more fault faces at a vertex are refused (junction).  The
// faces come from every rank (gathered, de-duplicated by their
// partition-invariant vertex keys), so n_v is the same on every partition.
// FrameAt(fi, x) is the frame at a point x of face fi:
//    n(xi) = ((1 - xi) nhat_0 + xi nhat_1) / |...|,  then the F1 rule  (2)
// with xi in [0, 1] along the chord from the face's vertex 0 to vertex 1
// (key order; xi within 1e-12 of 0 or 1 is snapped, so both faces at a
// vertex use exactly n_v there), nhat_i = n_v at a smooth vertex, else the
// face normal.  A face whose nhat_0 and nhat_1 both equal its normal returns
// Frame(fi) unchanged, and so does every face with the option off.  P1 only
// builds and reports these; the operators read them from phase P2/P3.  3-D
// meshes are refused with the option on.
//
// Scope (Chunhui Zhao, 2026-10-24): tetrahedral meshes, planar triangular
// faces (the scope of FaultOrientation); 2D triangle meshes from phase B.

#ifndef MFEM_SEAS_FAULT_FRAME_FIELD_HPP
#define MFEM_SEAS_FAULT_FRAME_FIELD_HPP

#include "mfem.hpp"
#include "fault_orientation.hpp"
#include "fault_frame.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <memory>
#include <ostream>
#include <set>
#include <vector>

namespace mfem
{
namespace seas
{

#ifdef MFEM_USE_MPI

/// The deck input of the orientation system ([fault_geometry]).
struct FaultFrameConfig
{
   FaultOrientationConfig orientation;   ///< ref_normal, pin, margins, fault attrs
   real_t up[3] = {0.0, 0.0, 1.0};       ///< F1 up vector, both regimes
   /// [fault_geometry] vertex_normals: false = "face" (default), true =
   /// "averaged" (the file comment; 2-D only)
   bool averaged_vertex_normals = false;
   /// [fault_geometry] vertex_normal_max_turn_deg: the turn above which a
   /// vertex is a kink and keeps the face normals (> 0)
   real_t vertex_normal_max_turn_deg = 5.0;
};

/// True when two configurations ask for the same orientation (the fields the
/// census and the frame read; the attribute sets only name persistent ids).
inline bool SameFrameConfig(const FaultFrameConfig &a, const FaultFrameConfig &b)
{
   const FaultOrientationConfig &p = a.orientation, &q = b.orientation;
   for (int i = 0; i < 3; i++)
   {
      if (a.up[i] != b.up[i] || p.ref_normal[i] != q.ref_normal[i]
          || p.pin_point[i] != q.pin_point[i]) { return false; }
   }
   return p.pin == q.pin && p.c_min == q.c_min && p.side_floor_frac == q.side_floor_frac
          && p.planarity_rel_tol == q.planarity_rel_tol
          && a.averaged_vertex_normals == b.averaged_vertex_normals
          && a.vertex_normal_max_turn_deg == b.vertex_normal_max_turn_deg;
}

struct FaultFrameReport
{
   FaultOrientationReport fault;   ///< the fault census
   FaultOrientationReport creep;   ///< the prescribed-slip census (faces = 0: none)
   int seam_edges = 0;             ///< edges shared by one fault and one creep face
   int frame_fallbacks = 0;        ///< fault faces whose frame used F6 (F7)
   // averaged vertex normals (all 0 with the option off); counts over the
   // fault's distinct vertices, the same on every rank
   int smooth_vertices = 0;        ///< two faces, turn <= the threshold: n_v by (1)
   int kink_vertices = 0;          ///< two faces, turn above it: face normals kept
   int tip_vertices = 0;           ///< one face: its normal kept
   real_t max_turn_smooth_deg = 0; ///< largest turn among the smooth vertices
   real_t max_turn_kink_deg = 0;   ///< largest turn among the kinks
};

class FaultFrameField
{
public:
   /// fault_int / creep_int: local interior face indices; fault_sh /
   /// creep_sh: shared face indices (ParMesh::GetSharedFace numbering), in
   /// the operator's order (interior faces first, then shared).  Collective.
   FaultFrameField(ParMesh &pm, const Array<int> &fault_int, const Array<int> &fault_sh,
                   const Array<int> &creep_int, const Array<int> &creep_sh,
                   const std::set<int> &creep_attrs, const FaultFrameConfig &cfg)
      : cfg_(cfg)
   {
      fault_int_.assign(fault_int.begin(), fault_int.end());
      fault_sh_.assign(fault_sh.begin(), fault_sh.end());
      creep_int_.assign(creep_int.begin(), creep_int.end());
      creep_sh_.assign(creep_sh.begin(), creep_sh.end());
      MPI_Comm comm = pm.GetComm();
      dim_ = pm.Dimension();
      if (dim_ == 2) { Check2D_(cfg); }

      // The fault census: the call NitscheWaveOperator::SetupFault_ made
      // before phase A (same faces, same order, same configuration).
      fault_ = std::make_unique<FaultOrientation>(pm, fault_int, fault_sh, cfg.orientation);
      rep_.fault = fault_->GetReport();

      // The prescribed-slip census: same pin, its own attributes.
      int nc = creep_int.Size() + creep_sh.Size();
      MPI_Allreduce(MPI_IN_PLACE, &nc, 1, MPI_INT, MPI_SUM, comm);
      creep_slot_.assign(pm.GetNumFaces(), -1);
      if (nc > 0)
      {
         FaultOrientationConfig cc = cfg.orientation;
         cc.fault_attrs = creep_attrs;
         creep_ = std::make_unique<FaultOrientation>(pm, creep_int, creep_sh, cc);
         rep_.creep = creep_->GetReport();
         for (int ci = 0; ci < creep_->NumLocalFaces(); ci++)
         {
            const int lf = (ci < creep_int.Size())
                           ? creep_int[ci] : pm.GetSharedFace(creep_sh[ci - creep_int.Size()]);
            creep_slot_[lf] = ci;
         }
      }
      else { rep_.creep.faces = 0; }

      // Per fault face: the N-09 frame (as SetupFault_ computed it) and sigma_f.
      const int nf = fault_->NumLocalFaces();
      const int nint = fault_int.Size();
      frame_.resize(nf);
      sigma_.resize(nf);
      const real_t eps_deg = fault_frame::DegeneracyThreshold(1.0);
      int fallbacks = 0;
      for (int fi = 0; fi < nf; fi++)
      {
         fault_frame::FaceVertices fv;
         fv.nv = fault_->FaceVertices();
         fault_->FaceVertexData(fi, fv.key, fv.x);
         const bool ok = fault_frame::ComputeFrame(fault_->Normal(fi), cfg.up, eps_deg, fv,
                                                   frame_[fi]);
         MFEM_VERIFY(ok, "FaultFrameField: degenerate frame (F6(b)) on local fault face "
                     << fi);
         if (frame_[fi].fallback) { fallbacks++; }

         const bool shared = (fi >= nint);
         FaceElementTransformations *Tr = shared
            ? pm.GetSharedFaceTransformations(fault_sh[fi - nint])
            : pm.GetInteriorFaceTransformations(fault_int[fi]);
         MFEM_VERIFY(Tr, "FaultFrameField: no face transformation for local fault face " << fi);
         sigma_[fi] = (fault_->PlusElement(fi) == Tr->Elem1No) ? 1.0 : -1.0;
         const IntegrationPoint &ipc = Geometries.GetCenter(Tr->GetGeometryType());
         Tr->SetAllIntPoints(&ipc);
         Vector nr2(dim_);
         CalcOrtho(Tr->Jacobian(), nr2);
         real_t nr[3] = {0.0, 0.0, 0.0};
         for (int d = 0; d < dim_; d++) { nr[d] = nr2(d); }
         const real_t c = fault_frame::Dot(nr, fault_->Normal(fi));
         MFEM_VERIFY((c > 0.0) == (sigma_[fi] > 0.0),
                     "FaultFrameField: on local fault face " << fi << " the S3 plus element ("
                     << fault_->PlusElement(fi) << ", Elem1 = " << Tr->Elem1No << ") and the "
                     "Elem1-outward normal (n_raw . n = " << c << ") disagree");
      }
      MPI_Allreduce(&fallbacks, &rep_.frame_fallbacks, 1, MPI_INT, MPI_SUM, comm);

      SeamCheck_(pm);
      VertexNormals_(pm);
   }

   const FaultFrameConfig &Config() const { return cfg_; }
   const real_t *Up() const { return cfg_.up; }
   /// Mesh dimension (2 or 3).
   int Dimension() const { return dim_; }
   /// 2D: the kernel slot that carries the in-plane slip (derivation 10.2):
   /// 0 (dip) when up lies in the mesh plane (BP3: a vertical section), 1
   /// (strike) when up = +-e_z (TPV205-2D: map view).  The other slot is an
   /// exact zero (Check2D_ refuses any other up).
   int InPlaneSlot() const
   {
      MFEM_VERIFY(dim_ == 2, "FaultFrameField::InPlaneSlot: a 2-D mesh only");
      return (cfg_.up[2] != 0.0) ? 1 : 0;
   }

   // ---- fault faces, operator order fi (interior faces first, then shared) ----
   const FaultOrientation &FaultCensus() const { return *fault_; }
   int NumFaultFaces() const { return fault_->NumLocalFaces(); }
   /// Oriented unit normal n of fault face fi (points away from the + side).
   const real_t *Normal(int fi) const { return fault_->Normal(fi); }
   /// N-09 frame (n, d, s) of fault face fi.
   const fault_frame::Frame &Frame(int fi) const { return frame_[fi]; }
   /// +1 when Elem1 of fault face fi is the + element, else -1.
   real_t Sigma(int fi) const { return sigma_[fi]; }

   // ---- vertex normals (file comment; P1 of the vertex-normal plan) ----
   bool AveragedVertexNormals() const { return cfg_.averaged_vertex_normals; }
   /// Vertex a (0 or 1, key order) of 2-D fault face fi: its coordinates,
   /// the normal nhat_a that (2) uses there (n_v at a smooth vertex, else
   /// the face normal; always the face normal with the option off), its kind
   /// (0 tip, 1 smooth, 2 kink; 1 also for the faces of a straight line) and
   /// its turn angle [deg] (0 at a tip).
   const real_t *VertexX(int fi, int a) const { return vx_[fi][a].data(); }
   const real_t *VertexNormal(int fi, int a) const { return vn_[fi][a].data(); }
   int VertexKind(int fi, int a) const { return vkind_[fi][a]; }
   real_t VertexTurnDeg(int fi, int a) const { return vturn_[fi][a]; }
   /// The frame (n, d, s) at the point x of fault face fi, by (2); Frame(fi)
   /// itself with the option off or where both nhat equal the face normal.
   void FrameAt(int fi, const real_t x[3], fault_frame::Frame &fr) const
   {
      const real_t *n0 = vn_[fi][0].data(), *n1 = vn_[fi][1].data();
      const real_t *nf = Normal(fi);
      if (!cfg_.averaged_vertex_normals
          || (n0[0] == nf[0] && n0[1] == nf[1] && n0[2] == nf[2]
              && n1[0] == nf[0] && n1[1] == nf[1] && n1[2] == nf[2]))
      {
         fr = frame_[fi];
         return;
      }
      const real_t *a = vx_[fi][0].data(), *b = vx_[fi][1].data();
      real_t ab2 = 0.0, t = 0.0;
      for (int d = 0; d < 3; d++)
      {
         const real_t ab = b[d] - a[d];
         ab2 += ab * ab;
         t += (x[d] - a[d]) * ab;
      }
      real_t xi = t / ab2;
      if (xi < 1e-12) { xi = 0.0; }
      if (xi > 1.0 - 1e-12) { xi = 1.0; }
      real_t n[3], nl = 0.0;
      for (int d = 0; d < 3; d++)
      {
         n[d] = (1.0 - xi) * n0[d] + xi * n1[d];
         nl += n[d] * n[d];
      }
      nl = std::sqrt(nl);
      for (int d = 0; d < 3; d++) { n[d] /= nl; }
      fault_frame::FaceVertices fv;
      fv.nv = fault_->FaceVertices();
      fault_->FaceVertexData(fi, fv.key, fv.x);
      const bool ok = fault_frame::ComputeFrame(n, cfg_.up, fault_frame::DegeneracyThreshold(1.0),
                                                fv, fr);
      MFEM_VERIFY(ok, "FaultFrameField::FrameAt: degenerate frame (F6(b)) on local fault face "
                  << fi);
   }

   // ---- prescribed-slip (creep) faces, by local face index ----
   bool HasCreep() const { return creep_ != nullptr; }
   bool IsCreepFace(int local_face) const
   {
      return local_face >= 0 && local_face < static_cast<int>(creep_slot_.size())
             && creep_slot_[local_face] >= 0;
   }
   /// Oriented unit normal of the prescribed-slip face with local index
   /// `local_face`; refuses a face that is not one.
   const real_t *CreepNormal(int local_face) const
   {
      MFEM_VERIFY(IsCreepFace(local_face), "FaultFrameField: local face " << local_face
                  << " is not a prescribed-slip (interior Dirichlet) face of the census");
      return creep_->Normal(creep_slot_[local_face]);
   }
   const FaultOrientation *CreepCensus() const { return creep_.get(); }

   /// An operator handed this field checks that its own face lists are the
   /// ones the field was built from.
   bool SameFaces(const Array<int> &fault_int, const Array<int> &fault_sh,
                  const Array<int> &creep_int, const Array<int> &creep_sh) const
   {
      auto same = [](const std::vector<int> &a, const Array<int> &b) {
         if (static_cast<int>(a.size()) != b.Size()) { return false; }
         for (int i = 0; i < b.Size(); i++) { if (a[i] != b[i]) { return false; } }
         return true;
      };
      return same(fault_int_, fault_int) && same(fault_sh_, fault_sh)
             && same(creep_int_, creep_int) && same(creep_sh_, creep_sh);
   }

   const FaultFrameReport &Report() const { return rep_; }

   void Print(std::ostream &os) const
   {
      const real_t *r = cfg_.orientation.ref_normal;
      os << "  FaultFrameField (one orientation for QD and dynamic): ref_normal ("
         << r[0] << ", " << r[1] << ", " << r[2] << "), up (" << cfg_.up[0] << ", "
         << cfg_.up[1] << ", " << cfg_.up[2] << ")\n"
         << "    fault: " << rep_.fault.faces << " faces, " << rep_.fault.components
         << " patch(es), " << rep_.fault.junction_edges
         << ((dim_ == 2) ? " junction vertices, " : " junction edges, ");
      if (dim_ == 2) { os << rep_.fault.tip_edges << " tips, "; }
      os << "flipped vs raw "
         << rep_.fault.flipped_vs_raw << ", pin margin " << rep_.fault.min_pin_margin
         << ", side margin " << rep_.fault.min_side_margin << ", F6 frames "
         << rep_.frame_fallbacks << "\n"
         << "    prescribed slip: " << rep_.creep.faces << " faces";
      if (rep_.creep.faces > 0)
      {
         os << ", " << rep_.creep.components << " patch(es), flipped vs raw "
            << rep_.creep.flipped_vs_raw << ", pin margin " << rep_.creep.min_pin_margin
            << ", seam edges with the fault " << rep_.seam_edges << " (all consistent)";
      }
      os << "\n";
      if (cfg_.averaged_vertex_normals)
      {
         os << "    vertex normals averaged (turn <= " << cfg_.vertex_normal_max_turn_deg
            << " deg): " << rep_.smooth_vertices << " smooth (largest turn "
            << rep_.max_turn_smooth_deg << " deg), " << rep_.kink_vertices
            << " kinks kept (largest turn " << rep_.max_turn_kink_deg << " deg), "
            << rep_.tip_vertices << " tips\n";
      }
   }

private:
   // The seam rule (see the file comment), on gathered, de-duplicated faces:
   // the same data and the same verdict on every rank (P1).
   void SeamCheck_(ParMesh &pm)
   {
      if (!creep_) { return; }
      MPI_Comm comm = pm.GetComm();
      // record: type (0 fault, 1 creep), 3 x 3 coordinates, oriented normal
      const int R = 13;
      std::vector<double> loc;
      auto add = [&](const FaultOrientation &O, int type) {
         for (int i = 0; i < O.NumLocalFaces(); i++)
         {
            long long key[3];
            real_t x[3][3];
            O.FaceVertexData(i, key, x);
            loc.push_back(type);
            for (int a = 0; a < 3; a++)
               for (int d = 0; d < 3; d++) { loc.push_back(x[a][d]); }
            for (int d = 0; d < 3; d++) { loc.push_back(O.Normal(i)[d]); }
         }
      };
      add(*fault_, 0);
      add(*creep_, 1);
      int nl = static_cast<int>(loc.size()), np = 0;
      MPI_Comm_size(comm, &np);
      std::vector<int> cnt(np), dsp(np, 0);
      MPI_Allgather(&nl, 1, MPI_INT, cnt.data(), 1, MPI_INT, comm);
      for (int r = 1; r < np; r++) { dsp[r] = dsp[r - 1] + cnt[r - 1]; }
      std::vector<double> all(dsp[np - 1] + cnt[np - 1]);
      MPI_Allgatherv(loc.data(), nl, MPI_DOUBLE, all.data(), cnt.data(), dsp.data(),
                     MPI_DOUBLE, comm);

      using P3 = std::array<real_t, 3>;
      struct SF { int type; P3 x[3]; real_t n[3]; };
      const int nv = fault_->FaceVertices();          // 2D: x[2] = 0 is padding
      std::map<std::array<real_t, 10>, int> seen;   // (type, sorted vertices)
      std::vector<SF> F;
      for (size_t o = 0; o + R <= all.size(); o += R)
      {
         SF f;
         f.type = static_cast<int>(all[o]);
         for (int a = 0; a < 3; a++)
            for (int d = 0; d < 3; d++) { f.x[a][d] = all[o + 1 + 3 * a + d]; }
         for (int d = 0; d < 3; d++) { f.n[d] = all[o + 10 + d]; }
         std::sort(f.x, f.x + nv);
         std::array<real_t, 10> k;
         k[0] = f.type;
         for (int a = 0; a < 3; a++)
            for (int d = 0; d < 3; d++) { k[1 + 3 * a + d] = f.x[a][d]; }
         if (seen.emplace(k, static_cast<int>(F.size())).second) { F.push_back(f); }
      }
      // edges by their (sorted) vertex pair -> (face, third vertex); 2D: the
      // shared entity is a vertex a (key (a, 0)) and the other vertex of the
      // segment plays the third vertex's part
      std::map<std::array<real_t, 6>, std::vector<std::pair<int, int>>> E;
      for (int g = 0; g < static_cast<int>(F.size()); g++)
      {
         if (nv == 2)
         {
            for (int i = 0; i < 2; i++)
            {
               const P3 &a = F[g].x[i];
               E[{a[0], a[1], a[2], 0.0, 0.0, 0.0}].push_back({g, 1 - i});
            }
            continue;
         }
         const int pairs[3][3] = {{0, 1, 2}, {0, 2, 1}, {1, 2, 0}};
         for (const auto &pr : pairs)
         {
            const P3 &a = F[g].x[pr[0]], &b = F[g].x[pr[1]];
            E[{a[0], a[1], a[2], b[0], b[1], b[2]}].push_back({g, pr[2]});
         }
      }
      int seam = 0, bad = 0, junction = 0;
      for (const auto &e : E)
      {
         int nf = 0, nc = 0;
         for (const auto &fr : e.second) { (F[fr.first].type == 0) ? nf++ : nc++; }
         if (nf == 0 || nc == 0) { continue; }
         if (nf + nc > 2) { junction++; continue; }
         seam++;
         // direction in which face g runs the edge a -> b under its normal
         // (2D: b - a = e_z, the segment extruded along +z as in its census)
         const real_t a[3] = {e.first[0], e.first[1], e.first[2]};
         const real_t b[3] = {e.first[3], e.first[4], e.first[5]};
         auto dir = [&](const std::pair<int, int> &fr) {
            const SF &f = F[fr.first];
            real_t ab[3], ar[3], c[3];
            for (int d = 0; d < 3; d++)
            {
               ab[d] = (nv == 2) ? ((d == 2) ? 1.0 : 0.0) : b[d] - a[d];
               ar[d] = f.x[fr.second][d] - a[d];
            }
            fault_frame::Cross(ab, ar, c);
            return fault_frame::Dot(c, f.n);
         };
         if (dir(e.second[0]) * dir(e.second[1]) >= 0.0) { bad++; }
      }
      MFEM_VERIFY(junction == 0, "FaultFrameField: " << junction << " edges (2D: vertices) "
                  "where fault and "
                  "prescribed-slip faces meet three or more at a time (a junction between "
                  "the fault and the creep surface; junction_policy \"refuse\")");
      MFEM_VERIFY(bad == 0, "FaultFrameField (seam rule): on " << bad << " of " << seam
                  << " edges (2D: vertices) shared by a fault face and a prescribed-slip "
                  "face the two "
                  "censuses orient the faces inconsistently; the one pin "
                  "([fault_geometry] ref_normal) does not give the fault and the creep "
                  "surface the same side. Refusing the run.");
      rep_.seam_edges = seam;
   }

   // The vertex normals (file comment): every rank gathers every fault face
   // (vertex keys, oriented normal), de-duplicates, and decides each vertex
   // from the same data, so the result is the same on every partition.
   void VertexNormals_(ParMesh &pm)
   {
      const int nf = NumFaultFaces();
      vx_.resize(nf);
      vn_.resize(nf);
      vkind_.assign(nf, {{0, 0}});
      vturn_.assign(nf, {{0.0, 0.0}});
      for (int fi = 0; fi < nf; fi++)
      {
         long long key[3];
         real_t x[3][3];
         fault_->FaceVertexData(fi, key, x);
         for (int a = 0; a < 2; a++)
            for (int d = 0; d < 3; d++)
            {
               vx_[fi][a][d] = x[a][d];
               vn_[fi][a][d] = Normal(fi)[d];
            }
      }
      if (!cfg_.averaged_vertex_normals) { return; }
      MFEM_VERIFY(dim_ == 2, "FaultFrameField: [fault_geometry] vertex_normals = \"averaged\" "
                  "is implemented for 2-D meshes only (vertex-normal plan rev1, D3); this mesh "
                  "is " << dim_ << "-D");
      MFEM_VERIFY(cfg_.vertex_normal_max_turn_deg > 0.0, "FaultFrameField: "
                  "vertex_normal_max_turn_deg = " << cfg_.vertex_normal_max_turn_deg
                  << " must be > 0");
      MPI_Comm comm = pm.GetComm();
      // record: the two vertex keys, the oriented normal
      const int R = 5;
      std::vector<double> loc;
      for (int fi = 0; fi < nf; fi++)
      {
         long long key[3];
         real_t x[3][3];
         fault_->FaceVertexData(fi, key, x);
         loc.push_back(static_cast<double>(key[0]));
         loc.push_back(static_cast<double>(key[1]));
         for (int d = 0; d < 3; d++) { loc.push_back(Normal(fi)[d]); }
      }
      int nl = static_cast<int>(loc.size()), np = 0;
      MPI_Comm_size(comm, &np);
      std::vector<int> cnt(np), dsp(np, 0);
      MPI_Allgather(&nl, 1, MPI_INT, cnt.data(), 1, MPI_INT, comm);
      for (int r = 1; r < np; r++) { dsp[r] = dsp[r - 1] + cnt[r - 1]; }
      std::vector<double> all(dsp[np - 1] + cnt[np - 1]);
      MPI_Allgatherv(loc.data(), nl, MPI_DOUBLE, all.data(), cnt.data(), dsp.data(),
                     MPI_DOUBLE, comm);
      using N3 = std::array<real_t, 3>;
      std::map<std::pair<long long, long long>, N3> faces;   // shared faces: twice
      for (size_t o = 0; o + R <= all.size(); o += R)
      {
         const long long k0 = std::llround(all[o]), k1 = std::llround(all[o + 1]);
         faces.emplace(std::make_pair(std::min(k0, k1), std::max(k0, k1)),
                       N3{all[o + 2], all[o + 3], all[o + 4]});
      }
      std::map<long long, std::vector<N3>> at;
      for (const auto &f : faces)
      {
         at[f.first.first].push_back(f.second);
         at[f.first.second].push_back(f.second);
      }
      struct V { int kind = 0; N3 n{{0.0, 0.0, 0.0}}; real_t turn = 0.0; };
      std::map<long long, V> vert;
      int junction = 0;
      const real_t rad = M_PI / 180.0;
      for (const auto &kv : at)
      {
         V v;
         const std::vector<N3> &ns = kv.second;
         if (ns.size() >= 3) { junction++; continue; }
         if (ns.size() == 1) { v.kind = 0; rep_.tip_vertices++; vert[kv.first] = v; continue; }
         const N3 &A = ns[0], &B = ns[1];
         real_t c[3];
         fault_frame::Cross(A.data(), B.data(), c);
         v.turn = std::atan2(std::sqrt(fault_frame::Dot(c, c)),
                             fault_frame::Dot(A.data(), B.data())) / rad;
         if (A == B) { v.kind = 1; v.n = A; v.turn = 0.0; }
         else if (v.turn <= cfg_.vertex_normal_max_turn_deg)
         {
            v.kind = 1;
            real_t s = 0.0;
            for (int d = 0; d < 3; d++) { v.n[d] = A[d] + B[d]; s += v.n[d] * v.n[d]; }
            s = std::sqrt(s);
            for (int d = 0; d < 3; d++) { v.n[d] /= s; }
         }
         else { v.kind = 2; }
         if (v.kind == 1)
         {
            rep_.smooth_vertices++;
            rep_.max_turn_smooth_deg = std::max(rep_.max_turn_smooth_deg, v.turn);
         }
         else
         {
            rep_.kink_vertices++;
            rep_.max_turn_kink_deg = std::max(rep_.max_turn_kink_deg, v.turn);
         }
         vert[kv.first] = v;
      }
      MFEM_VERIFY(junction == 0, "FaultFrameField: " << junction << " fault vertices are "
                  "shared by three or more fault faces (a junction): [fault_geometry] "
                  "vertex_normals = \"averaged\" refuses junctions (junction_policy \"refuse\")");
      for (int fi = 0; fi < nf; fi++)
      {
         long long key[3];
         real_t x[3][3];
         fault_->FaceVertexData(fi, key, x);
         for (int a = 0; a < 2; a++)
         {
            const V &v = vert.at(key[a]);
            vkind_[fi][a] = v.kind;
            vturn_[fi][a] = v.turn;
            if (v.kind == 1)
            {
               for (int d = 0; d < 3; d++) { vn_[fi][a][d] = v.n[d]; }
            }
         }
      }
   }

   // 2D: the inputs that must lie in the mesh plane (derivation 10.2).  up in
   // the plane gives s = -+e_z (slip in the dip slot), up = +-e_z gives d =
   // -+e_z (slip in the strike slot); any other up mixes the slots.
   static void Check2D_(const FaultFrameConfig &cfg)
   {
      const real_t *u = cfg.up, *r = cfg.orientation.ref_normal;
      MFEM_VERIFY((u[2] == 0.0 && (u[0] != 0.0 || u[1] != 0.0))
                  || (u[0] == 0.0 && u[1] == 0.0 && u[2] != 0.0),
                  "FaultFrameField (2D): up = (" << u[0] << ", " << u[1] << ", " << u[2]
                  << ") must lie in the mesh plane ([fault_geometry] up = [ux, uy]) or be "
                  "+-e_z (up = [0, 0, 1]), so that one tangential slot is an exact zero");
      MFEM_VERIFY(r[2] == 0.0, "FaultFrameField (2D): ref_normal = (" << r[0] << ", "
                  << r[1] << ", " << r[2] << ") must lie in the mesh plane");
      const real_t *q = cfg.orientation.pin_point;
      MFEM_VERIFY(cfg.orientation.pin != FaultOrientationConfig::Pin::Point || q[2] == 0.0,
                  "FaultFrameField (2D): the pin point must lie in the mesh plane (z = "
                  << q[2] << ")");
   }

   int dim_ = 3;
   FaultFrameConfig cfg_;
   std::vector<int> fault_int_, fault_sh_, creep_int_, creep_sh_;
   std::unique_ptr<FaultOrientation> fault_, creep_;
   std::vector<int> creep_slot_;               ///< local face -> creep census index
   std::vector<fault_frame::Frame> frame_;
   std::vector<real_t> sigma_;
   // vertex normals, per fault face and vertex (key order); file comment
   std::vector<std::array<std::array<real_t, 3>, 2>> vx_, vn_;
   std::vector<std::array<int, 2>> vkind_;
   std::vector<std::array<real_t, 2>> vturn_;
   FaultFrameReport rep_;
};

#endif // MFEM_USE_MPI

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_FAULT_FRAME_FIELD_HPP
