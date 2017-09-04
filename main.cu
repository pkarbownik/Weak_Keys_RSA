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

int N  = 100;	
EVP_PKEY* pPubKey  = NULL;
FILE*     pemFile   = NULL;
BIGNUM* modulus = NULL;
BIGNUM* GCD_result = NULL;
RSA* rsa = NULL;
unsigned int i, j;
char *nr = NULL;

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


void init_variables(){
	pPubKey = EVP_PKEY_new();
	modulus = BN_new();
	GCD_result = BN_new();
}

void free_variables(){
	BN_free(modulus);
	BN_free(GCD_result);
	EVP_PKEY_free(pPubKey);
    fclose(pemFile);
}

__device__ float multiplyByTwo(float number)
{
    return number * 2.0f;
}

// GCD KERNEL 
__global__ void kernel(BIGNUM** A, BIGNUM** B, BIGNUM** C){
	CUPRINTF("Value is: %f\n", multiplyByTwo(12.67));
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
    printf("Device %d: \"%s\" with Compute %d.%d capability\n",
           devID, props.name, props.major, props.minor);

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
    BIGNUM  **A =  (BIGNUM**)malloc(10*sizeof(BIGNUM));
    BIGNUM  **B =  (BIGNUM**)malloc(10*sizeof(BIGNUM));
    BIGNUM  **C =  (BIGNUM**)malloc(10*sizeof(BIGNUM));
    BIGNUM  **d_A = NULL;
    BIGNUM  **d_B = NULL;
    BIGNUM  **d_C = NULL;
	printf("lengft of A: %d\n", BN_dec2bn(&A[1], "12345676876543567845673456"));
	printf("lengft of B: %d\n", BN_dec2bn(&B[1], "4294967295"));
    printf("sizeof of BIGNUM: %d\n",  sizeof(BIGNUM));
	printf("Result is %s\n", BN_bn2dec(A[1]));
	printf("Result is %s\n", BN_bn2dec(B[1]));
	printf("B[1]->d[0] is %u\n", B[1]->d[0]);
	//printf("C[1]->d[0] is %d\n", C[1]->d[0]);

    if(cudaMalloc( (void**) &d_A, 10*sizeof(BIGNUM) ) != cudaSuccess )
		 printf("cudaMalloc error!\n");

    if(cudaMalloc( (void**) &d_B, 10*sizeof(BIGNUM) ) != cudaSuccess )
		 printf("cudaMalloc error!\n");

    if(cudaMalloc( (void**) &d_C, 10*sizeof(BIGNUM) ) != cudaSuccess )
		 printf("cudaMalloc error!\n");

  	cudaMemcpy( d_A, A, 10*sizeof(BIGNUM), cudaMemcpyHostToDevice );
  	cudaMemcpy( d_B, B, 10*sizeof(BIGNUM), cudaMemcpyHostToDevice );
	kernel<<<1, 5>>>(d_A, d_B, d_C);
	cudaEventSynchronize(stop_cuda);
    err = cudaGetLastError();

    if (err != cudaSuccess){
        fprintf(stderr, "Failed to launch vectorAdd kernel (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
	cudaMemcpy(C , d_C, 10*sizeof(BIGNUM), cudaMemcpyDeviceToHost);
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

