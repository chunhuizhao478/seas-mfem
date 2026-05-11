"""data_projection — offline preprocessor for the data-projection feature.

See miniapps/seas/document/features_dev/data_projection_feature_plan_v2.md
and data_projection_schema_v1.md for the full specification.
"""

from .raw_readers import (
    RawSliceStack,
    RawGrid,
    read_cvmh_ascii,
    read_csv_grid,
)
from .crs import geographic_to_utm11n, utm11n_to_geographic
from .sidecar import write_sidecar, read_sidecar_attrs, SCHEMA_VERSION
from .bbox_check import (
    mesh_msh_bbox,
    grid_bbox,
    grid_contains_bbox,
    assert_grid_contains_mesh,
    BBoxContainmentError,
)

__all__ = [
    "RawSliceStack",
    "RawGrid",
    "read_cvmh_ascii",
    "read_csv_grid",
    "geographic_to_utm11n",
    "utm11n_to_geographic",
    "write_sidecar",
    "read_sidecar_attrs",
    "SCHEMA_VERSION",
    "mesh_msh_bbox",
    "grid_bbox",
    "grid_contains_bbox",
    "assert_grid_contains_mesh",
    "BBoxContainmentError",
]
