#include "mfem.hpp"
#include "../../common/mpi_context.hpp"
#include "../../config/seas_config.hpp"
#include "../../config/seas_config_parser.hpp"
#include "../../config/bp5_params.hpp"
#include <cstdio>
using namespace mfem; using namespace mfem::seas;
int main(int argc, char *argv[]) {
   MPIContext mpi(&argc, &argv);
   printf("start\n");
   SEASConfig config;
   printf("config ok\n");
   return 0;
}
