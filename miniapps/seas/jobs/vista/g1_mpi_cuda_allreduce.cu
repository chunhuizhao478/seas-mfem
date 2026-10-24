// Device-buffer MPI_Allreduce test (GPU plan G1 acceptance): every rank
// fills a GPU array with rank + 1 and sums it in place with MPI; each entry
// must then equal size (size + 1) / 2.  Needs a CUDA-aware MPI.
//
//   nvcc -ccbin mpicxx -arch=sm_90 -o g1_mpi_cuda_allreduce g1_mpi_cuda_allreduce.cu
#include <mpi.h>
#include <cuda_runtime.h>
#include <cstdio>

__global__ void Fill(double *x, int n, double v)
{
   const int i = blockIdx.x * blockDim.x + threadIdx.x;
   if (i < n) { x[i] = v; }
}

int main(int argc, char **argv)
{
   MPI_Init(&argc, &argv);
   int rank, size, node_rank, ndev;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   MPI_Comm node;
   MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL,
                       &node);
   MPI_Comm_rank(node, &node_rank);
   MPI_Comm_free(&node);
   cudaGetDeviceCount(&ndev);
   cudaSetDevice(node_rank % ndev);

   const int n = 1 << 20;
   double *d = nullptr;
   cudaMalloc(&d, n * sizeof(double));
   Fill<<<(n + 255) / 256, 256>>>(d, n, rank + 1.0);
   cudaDeviceSynchronize();
   MPI_Allreduce(MPI_IN_PLACE, d, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

   double ends[2];
   cudaMemcpy(&ends[0], d, sizeof(double), cudaMemcpyDeviceToHost);
   cudaMemcpy(&ends[1], d + n - 1, sizeof(double), cudaMemcpyDeviceToHost);
   const double expect = 0.5 * size * (size + 1);
   int ok = (ends[0] == expect && ends[1] == expect) ? 1 : 0, all_ok = 0;
   MPI_Allreduce(&ok, &all_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
   if (rank == 0)
   {
      std::printf("device-buffer MPI_Allreduce on %d ranks: %s (got %g, %g;"
                  " expected %g)\n", size, all_ok ? "PASS" : "FAIL",
                  ends[0], ends[1], expect);
   }
   cudaFree(d);
   MPI_Finalize();
   return all_ok ? 0 : 1;
}
