#ifndef _DEVICE_CUDA_BIGNUM_H_
#define _DEVICE_CUDA_BIGNUM_H_

#include "test.h"
#include "cuda_bignum.h"
#include "files_manager.h"
#include <time.h>



__host__ __device__ int cu_dev_bn_ucmp(const U_BN *a, const U_BN *b);
__host__ __device__ long cu_dev_long_abs(long number);

__host__ __device__ int cu_dev_bn_usub(const U_BN *a, const U_BN *b, U_BN *r);


__host__ __device__ int cu_dev_bn_rshift1(U_BN *a);

__host__ __device__ int cu_dev_bn_lshift(U_BN *a, unsigned n);


__host__ __device__ U_BN *cu_dev_binary_gcd(U_BN *a, U_BN *b);


__host__ __device__ U_BN *cu_dev_fast_binary_euclid(U_BN *a, U_BN *b);

__host__ __device__ U_BN *cu_dev_classic_euclid(U_BN *a, U_BN *b);

void CPU_computation(unsigned number_of_keys, unsigned key_size, char *keys_directory);

__global__ void orgEuclideanKernel(U_BN *A, U_BN *B, U_BN *C, unsigned n);

__global__ void binEuclideanKernel(U_BN *A, U_BN *B, U_BN *C, unsigned n);

__global__ void fastBinaryKernel(U_BN *A, U_BN *B, U_BN *C, unsigned n);

#endif // #ifndef _DEVICE_CUDA_BIGNUM_H_