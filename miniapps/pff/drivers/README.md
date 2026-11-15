# drivers (phase 2)

`pff_driver.cpp -c cases/<model>/<case>.toml [-o key=value ...]`: reads the case, the mesh
(Gmsh 2.2, names resolved by `io/gmsh_physical_names.hpp`), builds the operators and steps to the end time.
