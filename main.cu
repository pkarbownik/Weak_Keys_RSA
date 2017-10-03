#include <stdio.h>
#include <openssl/rsa.h>
#include <openssl/pem.h>
#include <string.h>
#include <time.h>
#include "cuda_runtime.h"
#include "cuPrintf.cu"
#include <openssl/bn.h>
#include <ctype.h>



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







struct   __Q_VECTOR__{
    unsigned long* d;       
    int     top;    
};

typedef struct __Q_VECTOR__     VQ_VECTOR;

__global__ void testKernel(VQ_VECTOR *X, int N){
    int i= blockIdx.x*blockDim.x + threadIdx.x;
    //cuPrintf("\n testKernel entrance by the global threadIdx= %d\n", i);
    for(int k=0; k<X[i].top; k++)
        cuPrintf("testKernel entrance by the global threadIdx= %d value: %lu\n", i ,X[i].d[k]);
    //cuPrintf("\n");
}

# define cu_bn_correct_top(a) \
        { \
        BN_ULONG *ftl; \
        int tmp_top = (a)->top; \
        if (tmp_top > 0) \
                { \
                for (ftl= &((a)->d[tmp_top-1]); tmp_top > 0; tmp_top--) \
                        if (*(ftl--)) break; \
                (a)->top = tmp_top; \
                } \
        }

VQ_VECTOR *cu_BN_new(void)
{
    VQ_VECTOR *ret=NULL;
    ret->top = 0;
    ret->d = NULL;
    bn_check_top(ret);
    return (ret);
}

void cu_BN_free(VQ_VECTOR *a)
{
    if (a == NULL)
        return;
    bn_check_top(a);
    a->d = NULL;
}

int cu_BN_set_word(VQ_VECTOR *a, BN_ULONG w)
{
    bn_check_top(a);
    a->d[0] = w;
    a->top = (w ? 1 : 0);
    bn_check_top(a);
    return (1);
}
#  define cu_BN_zero(a)      (cu_BN_set_word((a),0))

int cu_BN_mul_word(VQ_VECTOR *a, BN_ULONG w)
{
    BN_ULONG ll;

    bn_check_top(a);
    w &= BN_MASK2;
    if (a->top) {
        if (w == 0)
            cu_BN_zero(a);
        else {
            ll = bn_mul_words(a->d, a->d, a->top, w);
            if (ll) {
                a->d[a->top++] = ll;
            }
        }
    }
    bn_check_top(a);
    return (1);
}

int cu_BN_add_word(VQ_VECTOR *a, BN_ULONG w)
{
    BN_ULONG l;
    int i;

    bn_check_top(a);
    w &= BN_MASK2;

    /* degenerate case: w is zero */
    if (!w)
        return 1;
    /* degenerate case: a is zero */
    if (BN_is_zero(a))
        return cu_BN_set_word(a, w);

    for (i = 0; w != 0 && i < a->top; i++) {
        a->d[i] = l = (a->d[i] + w) & BN_MASK2;
        w = (w > l) ? 1 : 0;
    }
    if (w && i == a->top) {
        a->top++;
        a->d[i] = w;
    }
    bn_check_top(a);
    return (1);
}

int cu_BN_dec2bn(VQ_VECTOR *bn, const char *a)
{
    VQ_VECTOR *ret = NULL;
    unsigned long l = 0;
    int i, j;
    int num;

    if ((a == NULL) || (*a == '\0'))
        return (0);

    for (i = 0; i <= (INT_MAX/4) && isdigit((unsigned char)a[i]); i++)
        continue;

    if (i > INT_MAX/4)
        goto err;

    num = i;
    if (bn == NULL)
        return (num);

    /*
     * a is the start of the digits, and it is 'i' long. We chop it into
     * BN_DEC_NUM digits at a time
     */
    if (bn == NULL) {
        if ((ret = cu_BN_new()) == NULL)
            return (0);
    } else {
        ret = bn;
        cu_BN_zero(ret);
    }

    j = BN_DEC_NUM - (i % BN_DEC_NUM);
    if (j == BN_DEC_NUM)
        j = 0;
    l = 0;
    while (--i >= 0) {
        l *= 10;
        l += *a - '0';
        a++;
        if (++j == BN_DEC_NUM) {
            cu_BN_mul_word(ret, BN_DEC_CONV);
            cu_BN_add_word(ret, l);
            l = 0;
            j = 0;
        }
    }

    cu_bn_correct_top(ret);
    bn = ret;
    bn_check_top(ret);
    /* Don't set the negative flag if it's zero. */
    return (num);
 err:
    if (bn == NULL)
        cu_BN_free(ret);
    return (0);
}



int main(void){
    int L = 128, //.Data length
        N = 1;

    VQ_VECTOR   *A,
                *device_VQ_VECTOR;

    cudaError_t cudaStatus;

    A =   (VQ_VECTOR*)malloc(N*sizeof(VQ_VECTOR));
    for(int i=0; i<N; i++){
        VQ_VECTOR a;
        a.d = (unsigned long*)malloc(L*sizeof(unsigned long));;
        a.top =   L;
        for(int j=0; j<L; j++)
            a.d[j]=0;

        A[i] = a;
    }

    cu_BN_dec2bn(&A[0], "437768685634765");
    L=A[0].top;
    //Prinf of all the elements of A
    for(int i=0; i<N; i++){
        printf("\nA[%d]={", i);
        for(int j=0; j<L; j++)
            printf("%lu ",A[i].d[j]);
        printf("}\n");
    }

    printf("\n\n");
    //I Allocate and Copy data from A to device_VQ_VECTORon the GPU memory

    cudaDeviceReset();
    cudaStatus = cudaMalloc((void**)&device_VQ_VECTOR, N*sizeof(VQ_VECTOR));    
    cudaStatus = cudaMemcpy(device_VQ_VECTOR, A, N*sizeof(VQ_VECTOR), cudaMemcpyHostToDevice);

    for(int i = 0; i != N; ++i) {
        /* can't access device_VQ_VECTOR[i].d directly from host-side,
         * working around it with proxy variable */
        unsigned long *out;
        cudaMalloc(&out, L*sizeof(unsigned long));
        cudaMemcpy(out, A[i].d, L*sizeof(unsigned long),
                cudaMemcpyHostToDevice);
        cudaMemcpy(&device_VQ_VECTOR[i].d, &out, sizeof(void*),
                cudaMemcpyHostToDevice);

        // will re-allocate later, for simplicity sake
        free(A[i].d);
    }

    cudaPrintfInit();
    testKernel<<<1,N>>>(device_VQ_VECTOR, N);//to test and see on a sigle thread
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
        unsigned long *array = (unsigned long*)malloc(L*sizeof(unsigned long));
        cudaMemcpy(array, A[i].d, L*sizeof(unsigned long),
                cudaMemcpyDeviceToHost);

        // assign new array to A[i]
        A[i].d = array;
    }
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "\n testKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
        return 1;
    }
/*  for(int i=0; i<2; i++){
        printf("\nA[%d]={", i);
        for(int j=0; j<L; j++)
            printf("%.3f",A[i].Data[j]);
        printf("}\n");
    }*/
    cudaFree(device_VQ_VECTOR);

    // don't forget to free A and all its Data

    return 0;
}