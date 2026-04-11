#include "mfem.hpp"
#include "../../config/seas_config.hpp"
#include <iostream>
int main() {
   std::cout << "before config" << std::endl;
   mfem::seas::SEASConfig config;
   std::cout << "config created" << std::endl;
   std::cout << "face_basis=" << config.solver.face_basis_type << std::endl;
   return 0;
}
