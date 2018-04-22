#include <stdio.h>
#include <openssl/rsa.h>
#include <openssl/pem.h>
#include <openssl/bn.h>
#include <string.h>
#include <time.h>
#include "cuda_runtime.h"
#include "devices.cuh"

// helper functions and utilities to work with CUDA
#include "helper_functions.h"
#include "helper_cuda.h"
#include "cuPrintf.cu"
#include <assert.h>

#define BIGNUM_SIZE 128*10
#define abs(x) ((x)<0 ? -(x) : (x))


typedef struct cu_bignum_st {
    BN_ULONG *d;                /* Pointer to an array of 'BN_BITS2' bit
                                 * chunks. */
    int top;                    /* Index of last used d +1. */
    /* The next are internal book keeping for bn_expand. */
    int dmax;                   /* Size of the d array. */
    int neg;                    /* one if the number is negative */
    int flags;
} CU_BIGNUM;

#ifndef MAX
#define MAX(a,b) (a > b ? a : b)
#endif

//The macro CUPRINTF is defined for architectures
//with different compute capabilities.
#if __CUDA_ARCH__ < 200     //Compute capability 1.x architectures
#define CUPRINTF cuPrintf
#else                       //Compute capability 2.x architectures
#define CUPRINTF(fmt, ...) printf("[%d, %d]:\t" fmt, \
                                  blockIdx.y*gridDim.x+blockIdx.x,\
                                  threadIdx.z*blockDim.x*blockDim.y+threadIdx.y*blockDim.x+threadIdx.x,\
                                  __VA_ARGS__)
#endif

/* This is used both by bn_expand2() and bn_dup_expand() */
/* The caller MUST check that words > b->dmax before calling this */

__device__ CU_BIGNUM *cu_bn_copy(CU_BIGNUM *a, const CU_BIGNUM *b)
{
    int i;
    BN_ULONG *A;
    const BN_ULONG *B;

    bn_check_top(b);

    if (a == b)
        return (a);
    /*if (bn_wexpand(a, b->top) == NULL)
        return (NULL);*/

#if 1
    A = a->d;
    B = b->d;
    for (i = b->top >> 2; i > 0; i--, A += 4, B += 4) {
        BN_ULONG a0, a1, a2, a3;
        a0 = B[0];
        a1 = B[1];
        a2 = B[2];
        a3 = B[3];
        A[0] = a0;
        A[1] = a1;
        A[2] = a2;
        A[3] = a3;
    }
    /* ultrix cc workaround, see comments in bn_expand_internal */
    switch (b->top & 3) {
    case 3:
        A[2] = B[2];
    case 2:
        A[1] = B[1];
    case 1:
        A[0] = B[0];
    case 0:;
    }
#else
    memcpy(a->d, b->d, sizeof(b->d[0]) * b->top);
#endif

    a->top = b->top;
    a->neg = b->neg;
    bn_check_top(a);
    return (a);
}

__device__ int cu_comparision(CU_BIGNUM* A, CU_BIGNUM* B){
    int i;
    int gt, lt;
    BN_ULONG t1, t2;

    if ((A == NULL) || (B == NULL)) {
        if (A != NULL)
            return (-1);
        else if (B != NULL)
            return (1);
        else
            return (0);
    }

    bn_check_top(A);
    bn_check_top(B);

    if (A->neg != B->neg) {
        if (A->neg)
            return (-1);
        else
            return (1);
    }
    if (A->neg == 0) {
        gt = 1;
        lt = -1;
    } else {
        gt = -1;
        lt = 1;
    }

    if (A->top > B->top)
        return (gt);
    if (A->top < B->top)
        return (lt);
    for (i = A->top - 1; i >= 0; i--) {
        t1 = A->d[i];
        t2 = B->d[i];
        if (t1 > t2)
            return (gt);
        if (t1 < t2)
            return (lt);
    }
    return (0);
}



__device__ int cu_subtraction(CU_BIGNUM* A, CU_BIGNUM* B, CU_BIGNUM* C){
    int max, min, dif;
    register BN_ULONG t1, t2, *ap, *bp, *rp;
    int i, carry;
    bn_check_top(A);
    bn_check_top(B);

    max = A->top;
    min = B->top;
    dif = max - min;
    if (dif < 0) {              /* unsigned subtraction of B from A, A must be larger than B. */
        return (0);
    }

   /* if (bn_wexpand(C, max) == NULL)
        return (0);*/

    ap = A->d;
    bp = B->d;
    rp = C->d;

    carry = 0;
    for (i = min; i != 0; i--) {
        t1 = *(ap++);
        t2 = *(bp++);
        if (carry) {
            carry = (t1 <= t2);
            t1 = (t1 - t2 - 1) & BN_MASK2;
        } else {
            carry = (t1 < t2);
            t1 = (t1 - t2) & BN_MASK2;
        }
        *(rp++) = t1 & BN_MASK2;
    }

    if (carry) {                /* subtracted */
        if (!dif)
            /* error: a < b */
            return 0;
        while (dif) {
            dif--;
            t1 = *(ap++);
            t2 = (t1 - 1) & BN_MASK2;
            *(rp++) = t2;
            if (t1)
                break;
        }
    }

    if (rp != ap) {
        for (;;) {
            if (!dif--)
                break;
            rp[0] = ap[0];
            if (!dif--)
                break;
            rp[1] = ap[1];
            if (!dif--)
                break;
            rp[2] = ap[2];
            if (!dif--)
                break;
            rp[3] = ap[3];
            rp += 4;
            ap += 4;
        }
    }

    C->top = max;
    C->neg = 0;
    bn_correct_top(C);
    return (1);
}

__device__ int cu_GCD(CU_BIGNUM* A, CU_BIGNUM* B, CU_BIGNUM* C){
	bn_check_top(A);
    bn_check_top(B);
	if(BN_is_zero(A)) 
		if (cu_bn_copy(C, B) == NULL)
			return 0;
		else
			return 1;
	if(BN_is_zero(B)) 
		if (cu_bn_copy(C, A) == NULL)
			return 0;
		else
			return 1;
	while(!BN_is_zero(B)){
		if(cu_comparision(A, B)>0){
			cu_subtraction(A,B,C);
			if (cu_bn_copy(A, C) == NULL)
				return 0;
		} else {
			cu_subtraction(B,A,C);
			if (cu_bn_copy(B, C) == NULL)
				return 0;
		}
	}

}



// GCD KERNEL 
__global__ void kernel(CU_BIGNUM* A, CU_BIGNUM* B, CU_BIGNUM* C){
	int id = threadIdx.x + 512 * blockDim.x;
	int tid = threadIdx.x;
	__syncthreads();
	if(id > 0){
		//cu_GCD(A[tid], B[tid], C[id]);
	}
	//id = A->top;
	CUPRINTF("kernel comparision: %d\n",C[1].d[0]=A[0].d[1]);
}



int main(int argc, char **argv)
{
	cudaEvent_t start_cuda, stop_cuda;
	cudaEventCreate(&start_cuda);
	cudaEventCreate(&stop_cuda);
	cudaEventRecord(start_cuda, 0);
    cudaError_t err = cudaSuccess;
	int devID;
    cudaDeviceProp props;

    // This will pick the best possible CUDA capable device
    devID = findCudaDevice(argc, (const char **)argv);

    //Get GPU information
    checkCudaErrors(cudaGetDevice(&devID));
    checkCudaErrors(cudaGetDeviceProperties(&props, devID));
    printf("Device %d: \"%s\" with Compute %d.%d capability\n", devID, props.name, props.major, props.minor);

    //Architectures with compute capability 1.x, function
    //cuPrintf() is used. Otherwise, function printf() is called.
    bool use_cuPrintf = (props.major < 2);

    if (use_cuPrintf)
    {
        //Initializaton, allocate buffers on both host
        //and device for data to be printed.
        cudaPrintfInit();
        printf("cuPrintf() is called. Output:\n\n");
    }
    //Architecture with compute capability 2.x, function
    //printf() is called.
    else
    {
        printf("printf() is called. Output:\n\n");
    }
    // Allocate the device input vector B
    CU_BIGNUM  *A =  (CU_BIGNUM*)malloc(BIGNUM_SIZE);
    CU_BIGNUM  *B =  (CU_BIGNUM*)malloc(BIGNUM_SIZE);
    CU_BIGNUM  *C =  (CU_BIGNUM*)malloc(BIGNUM_SIZE);
    CU_BIGNUM  *d_A =  NULL;
    CU_BIGNUM  *d_B =  NULL;
    CU_BIGNUM  *d_C =  NULL;	
	printf("lengft of A: %d\n", BN_dec2bn((BIGNUM**)&A, "1058184458759909426852660465011110"));
	printf("lengft of B: %d\n", BN_dec2bn((BIGNUM**)&B, "58599201514707945705960384879780048245246610"));
    printf("sizeof of CU_BIGNUM: %d\n", BIGNUM_SIZE);
	//printf("Result is %s\n", BN_bn2dec(A));
	//printf("Result is %s\n", BN_bn2dec(B));
	printf("A[0]->d[0] is %lu\n", A[0].d[0]);
	printf("B[0]->d[0] is %lu\n", B[0].d[0]);
	printf("A[0]->d[1] is %lu\n", A[0].d[1]);
	printf("B[0]->d[1] is %lu\n", B[0].d[1]);
	printf("A[0]->d[2] is %lu\n", A[0].d[2]);
	printf("B[0]->d[2] is %lu\n", B[0].d[2]);
	printf("A[0]->d[3] is %lu\n", A[0].d[3]);
	printf("B[0]->d[3] is %lu\n", B[0].d[3]);
    if(cudaMalloc(&d_A, BIGNUM_SIZE ) != cudaSuccess )
		 printf("cudaMalloc error!\n");

    if(cudaMalloc(&d_B, BIGNUM_SIZE ) != cudaSuccess )
		 printf("cudaMalloc error!\n");

    if(cudaMalloc(&d_C, BIGNUM_SIZE ) != cudaSuccess )
		 printf("cudaMalloc error!\n");

  	cudaMemcpy( d_A, A, BIGNUM_SIZE, cudaMemcpyHostToDevice );
  	cudaMemcpy( d_B, B, BIGNUM_SIZE, cudaMemcpyHostToDevice );
	kernel<<<1, 1>>>(d_A, d_B, d_C);
	cudaEventSynchronize(stop_cuda);
    err = cudaGetLastError();

    if (err != cudaSuccess){
        fprintf(stderr, "Failed to launch kernel (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
	cudaMemcpy(C , d_C, BIGNUM_SIZE, cudaMemcpyDeviceToHost);
	//printf("after C[1]->d[0] is %d\n", C[1]->d[0]);



    if (use_cuPrintf)
    {
        //Dump current contents of output buffer to standard
        //output, and origin (block id and thread id) of each line
        //of output is enabled(true).
        cudaPrintfDisplay(stdout, true);

        //Free allocated buffers by cudaPrintfInit().
        cudaPrintfEnd();
    }


    if (cudaFree(d_A) != cudaSuccess){
        fprintf(stderr, "Failed to free device vector A (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    if (cudaFree(d_B) != cudaSuccess){
        fprintf(stderr, "Failed to free device vector B (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    if (cudaFree(d_C) != cudaSuccess){
        fprintf(stderr, "Failed to free device vector C (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    if (cudaDeviceReset() != cudaSuccess){
        fprintf(stderr, "Failed to deinitialize the device! error=%s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    free(A);
    free(B);
    free(C);
    printf("Done\n");
	return 0;
}

