# Phase 7.2 of miniapps/seas/document/io_dev/PLAN_bulk_compression_and_size_estimator_2026-05-09.md.
#
# Per-driver registered-field tables for the size estimator.  The
# `FieldSchema` dataclass mirrors the (name, dof_basis, n_components,
# is_fault, is_volume, is_static) tuple of a `pv_out->RegisterDomainField(...)`
# or `pv_dc_->RegisterField(...)` call in the C++ code.
#
# Drift detection (`get_driver_schema(driver)` is called from the
# estimator AND from `test_estimate_output_size.py::test_schemas_match
# _driver_register_calls`).  See plan §Phase 7.2 R-308 for the regex
# AND the failure-message format.
#
# LIMITATION (R-308): the drift regex matches only string-literal field
# names.  Indirect calls like `RegisterField(get_name(...), ...)` are
# NOT matched and require manual schema entry.

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List

from _io_size_mesh import MeshSummary


@dataclass
class FieldSchema:
    name: str
    dof_basis: str         # "L2_p0" / "L2_p0_face" / "L2_pK_byNODES" / "H1_pK"
    n_components: int = 1
    sizeof_dtype: int = 8  # double
    is_fault: bool = False
    is_volume: bool = True
    is_static: bool = False  # only emitted at first save


# -----------------------------------------------------------------------------
# Per-driver schemas
# -----------------------------------------------------------------------------

# 12 standard fault fields + 5 _k4 diagnostic fields, all L2_p0 on the
# fault-face submesh.  param_a / param_Dc / fault_x2 / fault_x3 are
# emitted once at the first save (is_static=True) per Phase 2b §2.
_FAULT_FIELDS_BP5 = [
    FieldSchema("slip_dip", "L2_p0_face",
                is_fault=True, is_volume=False),
    FieldSchema("slip_strike", "L2_p0_face",
                is_fault=True, is_volume=False),
    FieldSchema("slip_rate_dip", "L2_p0_face",
                is_fault=True, is_volume=False),
    FieldSchema("slip_rate_strike", "L2_p0_face",
                is_fault=True, is_volume=False),
    FieldSchema("traction_dip", "L2_p0_face",
                is_fault=True, is_volume=False),
    FieldSchema("traction_strike", "L2_p0_face",
                is_fault=True, is_volume=False),
    FieldSchema("state_variable", "L2_p0_face",
                is_fault=True, is_volume=False),
    FieldSchema("normal_stress", "L2_p0_face",
                is_fault=True, is_volume=False),
    FieldSchema("param_a", "L2_p0_face",
                is_fault=True, is_volume=False, is_static=True),
    FieldSchema("param_Dc", "L2_p0_face",
                is_fault=True, is_volume=False, is_static=True),
    FieldSchema("fault_x2", "L2_p0_face",
                is_fault=True, is_volume=False, is_static=True),
    FieldSchema("fault_x3", "L2_p0_face",
                is_fault=True, is_volume=False, is_static=True),
]

# `_k4` diagnostics — only present in TPV* runs that supply stage-4
# buffers.  Estimator counts them when the user requests it.
_K4_FAULT_FIELDS = [
    FieldSchema("slip_rate_dip_k4", "L2_p0_face",
                is_fault=True, is_volume=False),
    FieldSchema("slip_rate_strike_k4", "L2_p0_face",
                is_fault=True, is_volume=False),
    FieldSchema("traction_dip_k4", "L2_p0_face",
                is_fault=True, is_volume=False),
    FieldSchema("traction_strike_k4", "L2_p0_face",
                is_fault=True, is_volume=False),
    FieldSchema("normal_stress_k4", "L2_p0_face",
                is_fault=True, is_volume=False),
]


TPV102_FIELDS = [
    # Volume side — primary collection (velocity + rank).
    FieldSchema("velocity", "L2_pK_byNODES", n_components=3),
    FieldSchema("mpi_rank", "L2_p0", is_static=True),
    # Bulk-side stress fields (registered on pv_bulk_out for TPV*).
    FieldSchema("sigma_yy", "L2_pK_byNODES"),
    FieldSchema("sigma_xy", "L2_pK_byNODES"),
    FieldSchema("sigma_xz", "L2_pK_byNODES"),
    # Fault side (Phase 2b InitFaultOutputBP5 in paraview_output.hpp).
    *_FAULT_FIELDS_BP5,
    *_K4_FAULT_FIELDS,
]

TPV104_FIELDS = list(TPV102_FIELDS)
TPV205_FIELDS = list(TPV102_FIELDS)

# BP5 driver: displacement on a vector DG FE space (vdim=3) — see
# `miniapps/seas/domain/elasticity_operator.hpp:fec_` (DG_FECollection).
# DOFs = n_elements * per-element-nodes (no continuity across faces),
# matching `L2_pK_byNODES`.  Fault fields are the same as the rest.
BP5_FIELDS = [
    FieldSchema("displacement", "L2_pK_byNODES", n_components=3),
    FieldSchema("mpi_rank", "L2_p0", is_static=True),
    *_FAULT_FIELDS_BP5,  # _k4 only when has_k4 is true; estimator
                          # doesn't add them by default for BP5.
]


_DRIVER_TO_SCHEMA: Dict[str, List[FieldSchema]] = {
    "tpv102": TPV102_FIELDS,
    "tpv104": TPV104_FIELDS,
    "tpv205": TPV205_FIELDS,
    "bp5": BP5_FIELDS,
}


def get_driver_schema(driver_name: str) -> List[FieldSchema]:
    """Return the registered-field list for one of `tpv102`, `tpv104`,
    `tpv205`, `bp5`.  Raises `ValueError` listing supported drivers
    on unknown name."""
    if driver_name not in _DRIVER_TO_SCHEMA:
        raise ValueError(
            f"unknown driver {driver_name!r}; supported: "
            f"{sorted(_DRIVER_TO_SCHEMA.keys())}")
    return list(_DRIVER_TO_SCHEMA[driver_name])


# -----------------------------------------------------------------------------
# DOF count formulas
# -----------------------------------------------------------------------------

def _binom(n: int, k: int) -> int:
    """Binomial coefficient (n choose k).  Uses math.comb on Python 3.8+
    via integer math; written explicitly so the script works on any
    standard library."""
    if k < 0 or k > n:
        return 0
    if hasattr(math, "comb"):
        return math.comb(n, k)
    num = 1
    for i in range(k):
        num = num * (n - i) // (i + 1)
    return num


def estimate_dofs(field: FieldSchema,
                  mesh: MeshSummary,
                  order: int) -> int:
    """Return an estimated DOF count for `field` on `mesh` at the given
    polynomial `order`.

    DOF count formulas (per plan §Phase 7.2 step 2):
      - L2_p0 element-wise:  n_dofs = n_elements
      - L2_p0_face (fault):  n_dofs = n_fault_faces
      - L2_pK_byNODES:       n_dofs = n_elements × ((K+1)(K+2)(K+3)/6)
      - H1_pK (tet K=1):     n_dofs ~= n_vertices
        (K>=2):              n_dofs ~= n_elements × C(K+3, 3) / shared_factor

    The `n_components` multiplier is NOT applied here; the caller
    multiplies separately so component handling stays explicit.
    """
    if field.dof_basis == "L2_p0":
        return mesh.n_elements
    if field.dof_basis == "L2_p0_face":
        return mesh.n_fault_faces
    if field.dof_basis == "L2_pK_byNODES":
        # Tet/hex node count for polynomial order K:
        # tet has (K+1)(K+2)(K+3)/6 nodes; hex has (K+1)^3.
        if mesh.element_type == "tet":
            per_elem = (order + 1) * (order + 2) * (order + 3) // 6
        elif mesh.element_type == "hex":
            per_elem = (order + 1) ** 3
        else:
            # Fallback: tet formula.
            per_elem = (order + 1) * (order + 2) * (order + 3) // 6
        return mesh.n_elements * per_elem
    if field.dof_basis == "H1_pK":
        if order == 1:
            return mesh.n_vertices
        # Approximate H1 DOF count for higher orders.  For tets at
        # K=2..4 the shared factor (= number of elements sharing each
        # vertex / edge / face DOF) averages ~4 in tetrahedral meshes.
        # This is documented as "approximation" in plan §Phase 7.2.
        if mesh.element_type == "tet":
            per_elem_hi = _binom(order + 3, 3)
            return max(mesh.n_vertices,
                       mesh.n_elements * per_elem_hi // 4)
        # Hex H1 closed-form is also approximate; use the tet recipe as
        # an upper bound.
        per_elem_hi = (order + 1) ** 3
        return max(mesh.n_vertices,
                   mesh.n_elements * per_elem_hi // 8)
    raise ValueError(
        f"unknown dof_basis {field.dof_basis!r} for field {field.name!r}")
