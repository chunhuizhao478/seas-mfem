#include "mfem.hpp"
#include "../../config/seas_config.hpp"
#include "test_macros.hpp"
#include <iostream>
using namespace mfem; using namespace mfem::seas;
int main() {
   std::cout << "start" << std::endl;
   SEASConfig config;
   std::cout << "face=" << config.solver.face_basis_type << std::endl;
   TEST_ASSERT(config.solver.face_basis_type == BasisType::GaussLobatto, "GL");
   TEST_PRINT_RESULTS();
   return 0;
}
