#include "cuPrintf.cu"
#include "test.h"
#include "cuda_bignum.h"
#include "files_manager.h"
//#include <openssl/bn.h>

#if __CUDA_ARCH__ < 200     //Compute capability 1.x architectures
#define CUPRINTF cuPrintf
#else                       //Compute capability 2.x architectures
#define CUPRINTF(fmt, ...) printf("[%d, %d]:\t" fmt, \
                                  blockIdx.y*gridDim.x+blockIdx.x,\
                                  threadIdx.z*blockDim.x*blockDim.y+threadIdx.y*blockDim.x+threadIdx.x,\
                                  __VA_ARGS__)
#endif

#define MAIN_COMPUTATIONS 1


__device__ int cu_dev_BN_ucmp(const U_BN *a, const U_BN *b){

    int i;
    unsigned t1, t2, *ap, *bp;

    i = a->top - b->top;
    if (i != 0)
        return (i);
    ap = a->d;
    bp = b->d;
    for (i = a->top - 1; i >= 0; i--) {
        t1 = ap[i];
        t2 = bp[i];
        if (t1 != t2)
            return ((t1 > t2) ? 1 : -1);
    }
    return (0);

}

__device__ long cu_dev_long_abs(long number){

    if(number<0)
        return -number;
    else
        return number;

}

__device__ int cu_dev_bn_usub(const U_BN *a, const U_BN *b, U_BN *r){

    unsigned max, min, dif;
    register unsigned t1, t2, *ap, *bp, *rp;
    int i, carry;

    if(NULL == a || NULL == b || NULL == r)
        return 0;

    if(NULL == a->d || NULL == b->d || NULL == r->d)
        return 0;

    max = a->top;
    min = b->top;
    dif = cu_dev_long_abs(max - min);

    ap = a->d;
    bp = b->d;
    rp = r->d;

#if 1
    carry = 0;
    for (i = min; i != 0; i--) {
        t1 = *(ap++);
        t2 = *(bp++);
        if (carry) {
            carry = (t1 <= t2);
            t1 = (t1 - t2 - 1);// & BN_MASK2;
        } else {
            carry = (t1 < t2);
            t1 = (t1 - t2);// & BN_MASK2;
        }
# if defined(IRIX_CC_BUG) && !defined(LINT)
        dummy = t1;
# endif
        *(rp++) = t1;// & BN_MASK2;
    }
#else
    carry = bn_sub_words(rp, ap, bp, min);
    ap += min;
    bp += min;
    rp += min;
#endif
    if (carry) {     
        if (!dif)

            return 0;
        while (dif) {
            dif--;
            t1 = *(ap++);
            t2 = (t1 - 1);// & BN_MASK2;
            *(rp++) = t2;
            if (t1)
                break;
        }
    }
#if 0
    memcpy(rp, ap, sizeof(*rp) * (max - i));
#else
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
#endif

    r->top = max;
    cu_bn_correct_top(r);
    return (1);

}


__device__ int cu_dev_bn_rshift1(U_BN *a){


    if(NULL == a)
        return 0;

    if(NULL == a->d)
        return 0;

    if (CU_BN_is_zero(a))
        return 0;

    unsigned *ap, *rp , t, c;
    int i, j;

    i = a->top;
    ap = a->d;

    j = i - (ap[i - 1] == 1);

    rp = a->d;
    t = ap[--i];
    c = (t & 1) ? CU_BN_TBIT : 0;
    if (t >>= 1)
        rp[i] = t;
    while (i > 0) {
        t = ap[--i];
        rp[i] = ((t >> 1) & CU_BN_MASK2) | c;
        c = (t & 1) ? CU_BN_TBIT : 0;
    }
    a->top = j;
    return (1);

}

__device__ int cu_dev_bn_lshift(U_BN *a, unsigned n){

    if(NULL == a)
        return 0;

    if(NULL == a->d)
        return 0;

    if (CU_BN_is_zero(a))
        return 0;

    if (0 == n)
        return 0;

    unsigned nw = 0, lb, rb, l;
    int i;
    unsigned nwb = 0, c = 0;

    nw = (n / CU_BN_BITS2);
    lb = (n % CU_BN_BITS2);
    rb = (CU_BN_BITS2 - lb);

    l=a->d[a->top-1];
    if( (l >> rb) > 0 ) nwb = 1;
    if(nw || nwb){
        //a->d = (unsigned*)realloc(a->d, (a->top + nw + nwb)*sizeof(unsigned)) ;
        a->d[a->top]=0;
        //memset((a->d+a->top-1), 0, (nw + nwb));
        //memset(a->d, 0, (nw + nwb)*sizeof(unsigned));
    }

    if (lb == 0 && nw != 0 ){
        for (i = a->top - 1; i >= 0; i--){
            a->d[nw + i] = a->d[i];
        }
    } else {
        for (i = 0; i < (a->top + nw + nwb); i++) {
            l = a->d[i];
            a->d[i] = (l << lb) | c;
            c = (l >> rb);

        }

    }
    a->top += (nw + nwb);
    return (1);

}


__device__ U_BN *cu_dev_euclid(U_BN *a, U_BN *b){
    U_BN *t = NULL;
    unsigned shifts = 0;
    while (!CU_BN_is_zero(b)) {
        if (cu_BN_is_odd(a)) {
            if (cu_BN_is_odd(b)) {
                cu_dev_bn_usub(a, b, a);
                cu_dev_bn_rshift1(a);
                if (cu_dev_BN_ucmp(a, b) < 0) {
                    t = a;
                    a = b;
                    b = t;
                }
            } else {      
                cu_dev_bn_rshift1(b);
                if (cu_dev_BN_ucmp(a, b) < 0) {
                    t = a;
                    a = b;
                    b = t;
                }
            }
        } else {              
            if (cu_BN_is_odd(b)) {
                cu_dev_bn_rshift1(a);
                if (cu_dev_BN_ucmp(a, b) < 0) {
                    t = a;
                    a = b;
                    b = t;
                }
            } else {      
                cu_dev_bn_rshift1(a);
                cu_dev_bn_rshift1(b);
                shifts++;
            }
        }
    }

    if (shifts) {
        cu_dev_bn_lshift(a, shifts);
    }
    return (a);

}

__global__ void testKernel(U_BN *A, U_BN *B, U_BN *C, int N){
    int i= blockIdx.x * blockDim.x + threadIdx.x;
    U_BN *TMP;
    //cu_dev_bn_usub(&A[i], &B[i], &C[i]);
    //cu_dev_bn_lshift(&C[i], 4);
    
    TMP = cu_dev_euclid(&A[i], &B[i]);
    //CUPRINTF("testKernel entrance by the global threadIdx= %d value: %u\n", i , TMP->d[0]);
    //CUPRINTF("testKernel entrance by the global threadIdx= %d value: %u\n", i , TMP->d[1]);
    C[i] = *TMP;
}

int main(void){
    int L = 5, //.Data length
        N = 100;

    unit_test(); //check all host bn functions
    print_mod_from_pem_file("keys_and_messages/1.pem");

    U_BN   *A;
    U_BN   *device_U_BN_A;
    U_BN   *B;
    U_BN   *device_U_BN_B;
    U_BN   *C;
    U_BN   *device_U_BN_C;

    cudaError_t cudaStatus;

    A =   (U_BN*)malloc(N*sizeof(U_BN));
    B =   (U_BN*)malloc(N*sizeof(U_BN));
    C =   (U_BN*)malloc(N*sizeof(U_BN));

    for(int i=0; i<N; i++){

        U_BN a;
        U_BN b;
        U_BN c;
        a.d = (unsigned*)malloc(L*sizeof(unsigned));
        b.d = (unsigned*)malloc(L*sizeof(unsigned));
        c.d = (unsigned*)malloc(L*sizeof(unsigned));
        a.top =   L;
        b.top =   L;
        c.top =   L;

        for(int j=0; j<L; j++)
            a.d[j]=0;

        for(int j=0; j<L; j++)
            b.d[j]=0;

        for(int j=0; j<L; j++)
            c.d[j]=0;

        A[i] = a;
        B[i] = b;
        C[i] = c;

        cu_BN_dec2bn(&A[i], "858238501677248042531768818944");
        cu_BN_dec2bn(&B[i], "8353015802438879251643065122143616");
    }

    L=A[0].top;
    L=B[0].top;
    L=C[0].top;


    cudaDeviceReset();
    cudaStatus = cudaMalloc((void**)&device_U_BN_A, N*sizeof(U_BN));    
    cudaStatus = cudaMalloc((void**)&device_U_BN_B, N*sizeof(U_BN));
    cudaStatus = cudaMalloc((void**)&device_U_BN_C, N*sizeof(U_BN));
    cudaStatus = cudaMemcpy(device_U_BN_A, A, N*sizeof(U_BN), cudaMemcpyHostToDevice);
    cudaStatus = cudaMemcpy(device_U_BN_B, B, N*sizeof(U_BN), cudaMemcpyHostToDevice);
    cudaStatus = cudaMemcpy(device_U_BN_C, C, N*sizeof(U_BN), cudaMemcpyHostToDevice);

    for(int i = 0; i < N; ++i) {
        unsigned long *out;
        cudaMalloc(&out, L*sizeof(unsigned));
        cudaMemcpy(out, A[i].d, L*sizeof(unsigned), cudaMemcpyHostToDevice);
        cudaMemcpy(&device_U_BN_A[i].d, &out, sizeof(void*), cudaMemcpyHostToDevice);

        cudaMalloc(&out, L*sizeof(unsigned));
        cudaMemcpy(out, B[i].d, L*sizeof(unsigned), cudaMemcpyHostToDevice);
        cudaMemcpy(&device_U_BN_B[i].d, &out, sizeof(void*), cudaMemcpyHostToDevice);

        cudaMalloc(&out, L*sizeof(unsigned));
        cudaMemcpy(out, C[i].d, L*sizeof(unsigned), cudaMemcpyHostToDevice);
        cudaMemcpy(&device_U_BN_C[i].d, &out, sizeof(void*), cudaMemcpyHostToDevice);

        // will re-allocate later, for simplicity sake
        free(A[i].d);
        free(B[i].d);
        free(C[i].d);
    }

    cudaPrintfInit();
    testKernel<<<1,N>>>(device_U_BN_A, device_U_BN_B, device_U_BN_C, N);//to test and see on a sigle thread
    cudaPrintfDisplay(stdout, true);
    cudaPrintfEnd();

    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "\n testKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
        return 1;
    }

    cudaStatus = cudaMemcpy(A, device_U_BN_A, N*sizeof(U_BN), cudaMemcpyDeviceToHost);
    cudaStatus = cudaMemcpy(B, device_U_BN_B, N*sizeof(U_BN), cudaMemcpyDeviceToHost);
    cudaStatus = cudaMemcpy(C, device_U_BN_C, N*sizeof(U_BN), cudaMemcpyDeviceToHost);


    for(int i = 0; i < N; ++i) {
        unsigned *array = (unsigned*)malloc(L*sizeof(unsigned));
        cudaMemcpy(array, A[i].d, L*sizeof(unsigned), cudaMemcpyDeviceToHost);
        A[i].d = array;

        cudaMemcpy(array, B[i].d, L*sizeof(unsigned), cudaMemcpyDeviceToHost);
        B[i].d = array;

        cudaMemcpy(array, C[i].d, L*sizeof(unsigned), cudaMemcpyDeviceToHost);
        C[i].d = array;
    }

    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "\n testKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
        return 1;
    }

    for(int i=0; i<N; i++){
        printf("cuda_kernel result c[%d]=%s\n", i, cu_bn_bn2hex(&C[i]));
    }

    cudaFree(device_U_BN_A);
    cudaFree(device_U_BN_B);
    cudaFree(device_U_BN_C);
    return 0;
}