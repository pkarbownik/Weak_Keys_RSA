#include "cuPrintf.cu"
#include "test.h"
#include "cuda_bignum.h"
//#include <openssl/bn.h>
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


__global__ void testKernel(VQ_VECTOR *X, int N, unsigned *out){
    int i= blockIdx.x*blockDim.x + threadIdx.x;
    /*for(int k=0; k<X[i].top; k++)
        cuPrintf("testKernel entrance by the global threadIdx= %d value: %u\n", i ,X[i].d[k]);*/
}



//extern __global__ void testKernel(VQ_VECTOR *X, int N);
//__device__ int x;


int main(void){

    unit_test();

    int L = 128, //.Data length
        N = 3;

    VQ_VECTOR   *A,
                *device_VQ_VECTOR;

    cudaError_t cudaStatus;

    A =   (VQ_VECTOR*)malloc(N*sizeof(VQ_VECTOR));
    for(int i=0; i<N; i++){
        VQ_VECTOR a;
        a.d = (unsigned*)malloc(L*sizeof(unsigned));;
        a.top =   L;
        for(int j=0; j<L; j++)
            a.d[j]=0;

        A[i] = a;
    }

    /*for (int i=0;i<N;i++){
        cu_BN_dec2bn(&A[i], "1844657685765856823456789023456789336478689676456476786468976475687647658767864576475744073709551617");
    }*/
    L=A[0].top;


    //Allocate and Copy data from A to device_VQ_VECTORon the GPU memory

    cudaDeviceReset();
    cudaStatus = cudaMalloc((void**)&device_VQ_VECTOR, N*sizeof(VQ_VECTOR));    
    cudaStatus = cudaMemcpy(device_VQ_VECTOR, A, N*sizeof(VQ_VECTOR), cudaMemcpyHostToDevice);
    unsigned *out;
    for(int i = 0; i != N; ++i) {
        /* can't access device_VQ_VECTOR[i].d directly from host-side,
         * working around it with proxy variable */
    	    //Prinf of all the elements of A

        /*printf("\nA[%d]={", i);
        for(int j=0; j<L; j++)
            printf("%u ",A[i].d[j]);
        printf("}\n");
    	printf("\n\n");*/
        cudaMalloc(&out, L*sizeof(unsigned));
        cudaMemcpy(out, A[i].d, L*sizeof(unsigned),
                cudaMemcpyHostToDevice);
        cudaMemcpy(&device_VQ_VECTOR[i].d, &out, sizeof(void*),
                cudaMemcpyHostToDevice);

        // will re-allocate later, for simplicity sake
        free(A[i].d);
    }

    cudaPrintfInit();
    testKernel<<<1,N>>>(device_VQ_VECTOR, N, out);//to test and see on a sigle thread
    cudaPrintfDisplay(stdout, true);
    cudaPrintfEnd();
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "\n testKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
        return 1;
    }
    cudaStatus = cudaMemcpy(A, device_VQ_VECTOR, N*sizeof(VQ_VECTOR), cudaMemcpyDeviceToHost);
    for(int i = 0; i != N; ++i) {
        // allocate array, copy data
        unsigned *array = (unsigned*)malloc(L*sizeof(unsigned));
        cudaMemcpy(array, A[i].d, L*sizeof(unsigned),
                cudaMemcpyDeviceToHost);

        // assign new array to A[i]
        A[i].d = array;
    }
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "\n testKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
        return 1;
    }
    cudaFree(device_VQ_VECTOR);

    // don't forget to free A and all its Data

    return 0;
}