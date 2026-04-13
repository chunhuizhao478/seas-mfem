#include <iostream>
#include "../../config/seas_config.hpp"
using namespace mfem::seas;
int main() {
   std::cout << "A" << std::endl;
   SEASConfig config;
   std::cout << "B face=" << config.solver.face_basis_type << std::endl;
   std::cout << "C dg=" << config.solver.dg_method << std::endl;
   std::cout << "done" << std::endl;
   return 0;
}
