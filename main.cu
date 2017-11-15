#include "cuPrintf.cu"
//#include "test.h"
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


__device__ int cu_BN_ucmp(const VQ_VECTOR *a, const VQ_VECTOR *b){

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

__device__ long cu_long_abs(long number){

    if(number<0)
        return -number;
    else
        return number;

}

__device__ int cu_bn_usub(const VQ_VECTOR *a, const VQ_VECTOR *b, VQ_VECTOR *r){

    unsigned max, min, dif;
    register unsigned t1, t2, *ap, *bp, *rp;
    int i, carry;

    if(NULL == a || NULL == b || NULL == r)
        return 0;

    if(NULL == a->d || NULL == b->d || NULL == r->d)
        return 0;

    max = a->top;
    min = b->top;
    dif = cu_long_abs(max - min);

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
            t1 = (t1 - t2 - 1) & BN_MASK2;
        } else {
            carry = (t1 < t2);
            t1 = (t1 - t2) & BN_MASK2;
        }
# if defined(IRIX_CC_BUG) && !defined(LINT)
        dummy = t1;
# endif
        *(rp++) = t1 & BN_MASK2;
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
            t2 = (t1 - 1) & BN_MASK2;
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


__device__ int cu_BN_rshift1(VQ_VECTOR *a){

    if(NULL == a)
        return 0;

    if(NULL == a->d)
        return 0;

    if (BN_is_zero(a))
        return 0;

    unsigned *ap, *rp , t, c;
    int i, j;

    i = a->top;
    DEBUG_PRINT("a->top: %u\n", a->top);
    ap = a->d;
    DEBUG_PRINT("a->d[0]: %u\n", a->d[0]);

    j = i - (ap[i - 1] == 1);
    DEBUG_PRINT("j: %d\n", j);

    rp = a->d;
    t = ap[--i];
    DEBUG_PRINT("t: %u\n", t);
    c = (t & 1) ? CU_BN_TBIT : 0;
    DEBUG_PRINT("c: %u\n", c);
    if (t >>= 1)
        rp[i] = t;
    while (i > 0) {
        t = ap[--i];
        DEBUG_PRINT("i: %d\n", i);
        DEBUG_PRINT("t: %u\n", t);
        rp[i] = ((t >> 1) & CU_BN_MASK2) | c;
        DEBUG_PRINT("rp[%d]: %u\n", i, rp[i]);
        c = (t & 1) ? CU_BN_TBIT : 0;
        DEBUG_PRINT("c: %u\n", c);
    }
    a->top = j;
    return (1);

}

__device__ int cu_BN_lshift(VQ_VECTOR *a, unsigned n){

    if(NULL == a)
        return 0;

    if(NULL == a->d)
        return 0;

    if (BN_is_zero(a))
        return 0;

    if (0 == n)
        return 0;

    unsigned nw = 0, lb, rb, l;
    int i;
    unsigned nwb = 0, c = 0;

    nw = (n / CU_BN_BITS2);
    DEBUG_PRINT("nw: %u\n", nw);
    lb = (n % CU_BN_BITS2);
    DEBUG_PRINT("lb: %u\n", lb);
    rb = (CU_BN_BITS2 - lb);
    DEBUG_PRINT("rb: %u\n", rb);

    for(i = 31; i>=0; i--){
        if((a->d[a->top-1]>>nwb)&1){ 
            nwb = ((rb>=i)?1:0);
            break;
        }
    }

    if(nw || nwb){
        return 0;
        //a->d = (unsigned*)cu_realloc(a->d, (a->top+ nw + nwb)*sizeof(unsigned)) ;
        //memset((a->d+a->top), 0, (nw + nwb)*sizeof(unsigned));
        //memset(a->d, 0, (nw + nwb)*sizeof(unsigned));
    }

    if (lb == 0 && nw != 0 ){
        for (i = a->top - 1; i >= 0; i--){
            a->d[nw + i] = a->d[i];
            DEBUG_PRINT("a->d[nw + i]: %u\n", a->d[nw + i]);
        }
    } else {
        for (i = 0; i < (a->top + nw + nwb); i++) {
            l = a->d[i];
            DEBUG_PRINT("l = a->d[%d]: %u\n", i, l);
            a->d[i] = (l << lb) | c;
            c = (l >> rb);
            DEBUG_PRINT("after lshift a->d[%d]: %u\n", i, a->d[i]);

        }

    }
    a->top += (nw + nwb);
    return (1);

}

__device__ VQ_VECTOR *cu_euclid(VQ_VECTOR *a, VQ_VECTOR *b){

    VQ_VECTOR *t = NULL;
    unsigned shifts = 0;
    while (!CU_BN_is_zero(b)) {
        if (cu_BN_is_odd(a)) {
            if (cu_BN_is_odd(b)) {
                DEBUG_PRINT("b id odd, a is equal: %s\n", cu_bn_bn2hex(a));
                cu_bn_usub(a, b, a);
                DEBUG_PRINT("b id odd, a-b is equal: %s\n", cu_bn_bn2hex(a));
                cu_BN_rshift1(a);
                DEBUG_PRINT("b id odd, a-b>>1 is equal: %s\n", cu_bn_bn2hex(a));
                if (cu_BN_ucmp(a, b) < 0) {
                    t = a;
                    a = b;
                    b = t;
                }
            } else {      
                DEBUG_PRINT("b id even, b is equal: %s\n", cu_bn_bn2hex(b));
                cu_BN_rshift1(b);
                DEBUG_PRINT("b id even, b>>1 is equal: %s\n", cu_bn_bn2hex(b));
                if (cu_BN_ucmp(a, b) < 0) {
                    t = a;
                    a = b;
                    b = t;
                }
            }
        } else {              
            if (cu_BN_is_odd(b)) {
                DEBUG_PRINT("a id even, b is odd, a is equal: %s\n", cu_bn_bn2hex(a));
                cu_BN_rshift1(a);
                DEBUG_PRINT("a id even, b is odd, a>>1 is equal: %s\n", cu_bn_bn2hex(a));
                if (cu_BN_ucmp(a, b) < 0) {
                    t = a;
                    a = b;
                    b = t;
                }
            } else {       
                DEBUG_PRINT("a id even, b is even, a is equal: %s\n", cu_bn_bn2hex(a));
                DEBUG_PRINT("a id even, b is even, b is equal: %s\n", cu_bn_bn2hex(b));
                cu_BN_rshift1(a);
                DEBUG_PRINT("a id even, b is even, a>>1 is equal: %s\n", cu_bn_bn2hex(a));
                cu_BN_rshift1(b);
                DEBUG_PRINT("a id even, b is even, b>>1 is equal: %s\n", cu_bn_bn2hex(b));
                shifts++;
            }
        }

    }

    if (shifts) {
        DEBUG_PRINT("a is equal: %s\n", cu_bn_bn2hex(a));
        DEBUG_PRINT("shifts is equal: %u\n", shifts);
        cu_BN_lshift(a, shifts);
        DEBUG_PRINT("a<<shifts is equal: %s\n", cu_bn_bn2hex(a));
    }
    //cu_BN_free(t);
    return (a);

}

__global__ void testKernel(VQ_VECTOR *A, VQ_VECTOR *B, VQ_VECTOR *C, int N){
    int i= blockIdx.x * blockDim.x + threadIdx.x;
    //int p;
    //for(int k=0; k<N; k++)
    VQ_VECTOR *TMP;

    CUPRINTF("testKernel entrance by the global threadIdx= %d value: %u\n", i , A[i].d[0]);
    CUPRINTF("testKernel entrance by the global threadIdx= %d value: %u\n", i , B[i].d[0]);
    cu_bn_usub(&A[i], &B[i], &C[i]);
    //CUPRINTF("testKernel entrance by the global threadIdx= %d value: %u\n", i , TMP->d[0]);
    CUPRINTF("testKernel entrance by the global threadIdx= %d value: %u\n", i , C[i].d[0]);
    //p = cu_bn_usub(dev_A[i], dev_B[i], dev_C[i]);
    //cuPrintf("testKernel: %d\n", p);
}

int main(void){
    int L = 128, //.Data length
        N = 1;

    VQ_VECTOR   *A;
    VQ_VECTOR   *device_VQ_VECTOR_A;
    VQ_VECTOR   *B;
    VQ_VECTOR   *device_VQ_VECTOR_B;
    VQ_VECTOR   *C;
    VQ_VECTOR   *device_VQ_VECTOR_C;

    cudaError_t cudaStatus;

    A =   (VQ_VECTOR*)malloc(N*sizeof(VQ_VECTOR));
    B =   (VQ_VECTOR*)malloc(N*sizeof(VQ_VECTOR));
    C =   (VQ_VECTOR*)malloc(N*sizeof(VQ_VECTOR));

    for(int i=0; i<N; i++){
        VQ_VECTOR a;
        VQ_VECTOR b;
        VQ_VECTOR c;
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
    }

    cu_BN_dec2bn(&A[0], "37358654342800865389770");
    cu_BN_dec2bn(&B[0], "34121356110372842586730");
    cu_BN_dec2bn(&C[0], "437768685634765");
    L=A[0].top;
    L=B[0].top;
    L=C[0].top;
    //Prinf of all the elements of A
    /*for(int i=0; i<N; i++){
        printf("\nA[%d]={", i);
        for(int j=0; j<L; j++)
            printf("%u ",A[i].d[j]);
        printf("}\n");
    }
    printf("\n\n");*/
    //I Allocate and Copy data from A to device_VQ_VECTORon the GPU memory

    cudaDeviceReset();
    cudaStatus = cudaMalloc((void**)&device_VQ_VECTOR_A, N*sizeof(VQ_VECTOR));    
    cudaStatus = cudaMalloc((void**)&device_VQ_VECTOR_B, N*sizeof(VQ_VECTOR));
    cudaStatus = cudaMalloc((void**)&device_VQ_VECTOR_C, N*sizeof(VQ_VECTOR));
    cudaStatus = cudaMemcpy(device_VQ_VECTOR_A, A, N*sizeof(VQ_VECTOR), cudaMemcpyHostToDevice);
    cudaStatus = cudaMemcpy(device_VQ_VECTOR_B, B, N*sizeof(VQ_VECTOR), cudaMemcpyHostToDevice);
    cudaStatus = cudaMemcpy(device_VQ_VECTOR_C, C, N*sizeof(VQ_VECTOR), cudaMemcpyHostToDevice);

    for(int i = 0; i != N; ++i) {
        unsigned long *out;
        cudaMalloc(&out, L*sizeof(unsigned));
        cudaMemcpy(out, A[i].d, L*sizeof(unsigned), cudaMemcpyHostToDevice);
        cudaMemcpy(&device_VQ_VECTOR_A[i].d, &out, sizeof(void*), cudaMemcpyHostToDevice);

        cudaMalloc(&out, L*sizeof(unsigned));
        cudaMemcpy(out, B[i].d, L*sizeof(unsigned), cudaMemcpyHostToDevice);
        cudaMemcpy(&device_VQ_VECTOR_B[i].d, &out, sizeof(void*), cudaMemcpyHostToDevice);

        cudaMalloc(&out, L*sizeof(unsigned));
        cudaMemcpy(out, C[i].d, L*sizeof(unsigned), cudaMemcpyHostToDevice);
        cudaMemcpy(&device_VQ_VECTOR_C[i].d, &out, sizeof(void*), cudaMemcpyHostToDevice);

        // will re-allocate later, for simplicity sake
        free(A[i].d);
        free(B[i].d);
        free(C[i].d);
    }

    cudaPrintfInit();
    testKernel<<<1,N>>>(device_VQ_VECTOR_A, device_VQ_VECTOR_B, device_VQ_VECTOR_C, N);//to test and see on a sigle thread
    cudaPrintfDisplay(stdout, true);
    cudaPrintfEnd();

    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "\n testKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
        return 1;
    }

    cudaStatus = cudaMemcpy(A, device_VQ_VECTOR_A, N*sizeof(VQ_VECTOR), cudaMemcpyDeviceToHost);
    cudaStatus = cudaMemcpy(B, device_VQ_VECTOR_B, N*sizeof(VQ_VECTOR), cudaMemcpyDeviceToHost);
    cudaStatus = cudaMemcpy(C, device_VQ_VECTOR_C, N*sizeof(VQ_VECTOR), cudaMemcpyDeviceToHost);


    for(int i = 0; i != N; ++i) {
        unsigned *array = (unsigned*)malloc(L*sizeof(unsigned));
        cudaMemcpy(array, A[i].d, L*sizeof(unsigned), cudaMemcpyDeviceToHost);
        A[i].d = array;

        cudaMemcpy(array, B[i].d, L*sizeof(unsigned),
        cudaMemcpyDeviceToHost);
        B[i].d = array;

        cudaMemcpy(array, C[i].d, L*sizeof(unsigned),
        cudaMemcpyDeviceToHost);
        C[i].d = array;
    }

    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "\n testKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
        return 1;
    }

    cudaFree(device_VQ_VECTOR_A);
    cudaFree(device_VQ_VECTOR_B);
    cudaFree(device_VQ_VECTOR_C);
    return 0;
}