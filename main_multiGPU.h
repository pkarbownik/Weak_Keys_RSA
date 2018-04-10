#ifndef MAIN_MULTIGPU_H
#define MAIN_MULTIGPU_H

#include "cuda_bignum.h"

typedef struct
{
    //Host-side input data
    int dataN;

    //Partial sum for this GPU
    U_BN   *h_A, *h_B, *h_C;

    //Device buffers
    U_BN   *device_U_BN_A, *device_U_BN_B, *device_U_BN_C;

    //Reduction copied back from GPU
    float *h_Sum_from_device;

    //Stream for asynchronous command execution
    cudaStream_t stream;

} TGPUplan;

extern "C"
void launch_reduceKernel(float *d_Result, float *d_Input, int N, int BLOCK_N, int THREAD_N, cudaStream_t &s);

#endif // MAIN_MULTIGPU_H
