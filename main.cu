#include "cuPrintf.cu"
#include "test.h"
#include "cuda_bignum.h"
#include "files_manager.h"
#include <time.h>
//#include <openssl/bn.h>


#if __CUDA_ARCH__ < 200     //Compute capability 1.x architectures
#define CUPRINTF cuPrintf
#else                       //Compute capability 2.x architectures
#define CUPRINTF(fmt, ...) printf("[%d, %d]:\t" fmt, \
                                  blockIdx.y*gridDim.x+blockIdx.x,\
                                  threadIdx.z*blockDim.x*blockDim.y+threadIdx.y*blockDim.x+threadIdx.x,\
                                  __VA_ARGS__)
#endif

__host__ __device__ int cu_dev_bn_ucmp(const U_BN *a, const U_BN *b){

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

__host__ __device__ long cu_dev_long_abs(long number){

    if(number<0)
        return -number;
    else
        return number;

}

__host__ __device__ int cu_dev_bn_usub(const U_BN *a, const U_BN *b, U_BN *r){

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


__host__ __device__ int cu_dev_bn_rshift1(U_BN *a){


    if(NULL == a)
        return 0;

    if(NULL == a->d)
        return 0;

    if (cu_bn_is_zero(a))
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

__host__ __device__ int cu_dev_bn_lshift(U_BN *a, unsigned n){

  if(NULL == a)
        return 0;

    if(NULL == a->d)
        return 0;

    if (cu_bn_is_zero(a))
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
    if( ((l << lb)%CU_BN_BITS2) > 0 ) nwb = 1;
    if (lb == 0)
        for (i = a->top - 1; i >= 0; i--){
            a->d[nw + i] = a->d[i];
        }
    else
        for (i = a->top - 1; i >= 0; i--) {
            l = a->d[i];
            a->d[nw + i + 1] |= (l >> rb) & CU_BN_MASK2;
            a->d[nw + i] = (l << lb) & CU_BN_MASK2;
        }


    a->top += (nw + nwb);
    char *buffer=NULL;
    buffer = (char *)malloc(sizeof(char) * (CU_BN_BYTES * a->top));
    free(buffer);
    return (1);

}


__host__ __device__ U_BN *cu_dev_binary_euclid(U_BN *a, U_BN *b){
    U_BN *t = NULL;
    unsigned shifts = 0;

    if (cu_dev_bn_ucmp(a, b) < 0) {
        t = a;
        a = b;
        b = t;
    }

    while (!cu_bn_is_zero(b)) {
        if (cu_bn_is_odd(a)) {
            if (cu_bn_is_odd(b)) {
                cu_dev_bn_usub(a, b, a);
                cu_dev_bn_rshift1(a);
                if (cu_dev_bn_ucmp(a, b) < 0) {
                    t = a;
                    a = b;
                    b = t;
                }
            } else {      
                cu_dev_bn_rshift1(b);
                if (cu_dev_bn_ucmp(a, b) < 0) {
                    t = a;
                    a = b;
                    b = t;
                }
            }
        } else {              
            if (cu_bn_is_odd(b)) {
                cu_dev_bn_rshift1(a);
                if (cu_dev_bn_ucmp(a, b) < 0) {
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


__host__ __device__ U_BN *cu_dev_fast_binary_euclid(U_BN *a, U_BN *b){
    U_BN *t = NULL;
    do {
        if (cu_dev_bn_ucmp(a, b) < 0) {
            t = a;
            a = b;
            b = t;
        }
        if(!cu_dev_bn_usub(a, b, a)) break;
        while(!(a->d[0]&1)) {
            if(!cu_dev_bn_rshift1(a)) break;
        }
    } while (!cu_bn_is_zero(b));
    return (a);
}

__host__ __device__ U_BN *cu_dev_classic_euclid(U_BN *a, U_BN *b){

    while (cu_dev_bn_ucmp(a, b) != 0) {
        if (cu_dev_bn_ucmp(a, b) > 0) {
            cu_dev_bn_usub(a, b, a); 
        }
        else {
            cu_dev_bn_usub(b, a, b);
        }
    }

    return (a);

}

void CPU_computation(unsigned number_of_keys, unsigned key_size, char *keys_directory){

    U_BN tmp;
    unsigned sum=0;
    BIGNUM *r;
    BN_CTX *ctx;
    EVP_PKEY* pPubKey  = NULL;
    FILE*     pemFile    = NULL;
    RSA* rsa;
    int i, j;
    int L = ((key_size+31) / (8*sizeof(unsigned)));
    unsigned k = 0;
    BIGNUM   *PEMs;
    char *tmp_path;

    PEMs = (BIGNUM*)malloc(number_of_keys*sizeof(BIGNUM));
    tmp.d = (unsigned*)malloc(L*sizeof(unsigned));
    tmp.top = L;

    for(j=0; j<L; j++){
        tmp.d[j]=0;
    }

    for(i=0; i<number_of_keys; i++){

        asprintf(&tmp_path, "%s/%d.pem", keys_directory, (i+1));
        if( !( (pemFile = fopen(tmp_path, "rt") ) && ( pPubKey = PEM_read_PUBKEY(pemFile,NULL,NULL,NULL) ) ) ) {
            fprintf(stderr,"Cannot read \"public key\".\n");
        }
        rsa = EVP_PKEY_get1_RSA(pPubKey);
        PEMs[i] = *(rsa->n);
        fclose(pemFile);

    }

    free(tmp.d);
    ctx = BN_CTX_new();
    r=BN_new();
    clock_t start = clock();
    for(i=0, k=0; i<number_of_keys; i++){
        for(j=(i+1); j<number_of_keys; j++, k++){
            BN_gcd(r, &PEMs[i], &PEMs[j], ctx);
            if(!BN_is_one(r)){
                //printf( "A[%d]: %s\nB[%d]: %s\neuclid: %s\n\n", i, BN_bn2hex(&PEMs[i]), j, BN_bn2hex(&PEMs[j]), BN_bn2hex(r));
                sum+=1;
            }
        }
    }
    clock_t stop = clock();
    double elapsed = (double)(stop - start) * 1000.0 / CLOCKS_PER_SEC;
    printf("[CPU] Time elapsed in ms: %f\n", elapsed);
    printf("[CPU] Weak keys: %u\n", sum);
    BN_CTX_free(ctx);
    BN_free(r);

}

__global__ void testKernel(U_BN *A, U_BN *B, U_BN *C, unsigned n) {
    int i= blockIdx.x * blockDim.x + threadIdx.x;
    U_BN *TMP=NULL;

    if(i<n){
        //cu_dev_bn_rshift1(&A[i]);
        //cu_dev_bn_lshift(&A[i], 1);
        //cu_dev_bn_usub(&A[i],&B[i],TMP);
        //TMP = cu_dev_binary_euclid(&A[i], &B[i]);
        //TMP = cu_dev_classic_euclid(&A[i], &B[i]);
        TMP = cu_dev_fast_binary_euclid(&A[i], &B[i]);
        if(TMP->d[0]!=1)
            cuPrintf("testKernel entrance by the global threadIdx= %d \n", i);
        C[i] = *TMP;
    }
}

int main(int argc, char* argv[]){
	
    unsigned number_of_keys;
    unsigned key_size;
    unsigned thread_per_block;
    unsigned number_of_comutations;
    char *keys_directory;
    int counter;
    if(argc==5) {
        for(counter=0;counter<argc;counter++){
            switch(counter){
                case 1:
                    printf("\nnumber_of_keys argv[%d]: %s\n",counter,argv[counter]);
                    number_of_keys=atoi(argv[counter]);
                    break;
                case 2:
                    printf("\nkey_size argv[%d]: %s\n",counter,argv[counter]);
                    key_size=atoi(argv[counter]);
                    break;
                case 3:
                    printf("\nthreads_per_block argv[%d]: %s\n",counter,argv[counter]);
                    thread_per_block=atoi(argv[counter]);
                    break;
                case 4:
                    printf("\nname of keys directory argv[%d]: %s\n",counter,argv[counter]);
                    keys_directory=argv[counter];
                    break;
                default:
                    break;
            }
        }
    } else {
        printf("\nFind weak keys\nUsage:\n ./GCD_RSA number_of_keys key_size threads_per_block keys_directory_name\n");
        return 0;
    }

    // simplified binomial coefficient
    number_of_comutations=((number_of_keys/2)*(number_of_keys-1));

    U_BN tmp;
    int L = ((key_size+31) / (8*sizeof(unsigned)));
    unsigned i, j;
    unsigned k = 0;
    U_BN   *A, *B, *C;
    U_BN   *device_U_BN_A, *device_U_BN_B, *device_U_BN_C;
    U_BN   *cu_PEMs;
    char *tmp_path;
    cudaError_t cudaStatus;

    //unit_test();
    //CPU_computation(number_of_keys, key_size, keys_directory);


    A    = (U_BN*)malloc(number_of_comutations*sizeof(U_BN));
    B    = (U_BN*)malloc(number_of_comutations*sizeof(U_BN));
    C    = (U_BN*)malloc(number_of_comutations*sizeof(U_BN));
    cu_PEMs = (U_BN*)malloc(number_of_keys*sizeof(U_BN));

    for(i=0; i<number_of_comutations; i++){

        U_BN a;
        U_BN b;
        U_BN c;
        a.d = (unsigned*)malloc(L*sizeof(unsigned));
        b.d = (unsigned*)malloc(L*sizeof(unsigned));
        c.d = (unsigned*)malloc(L*sizeof(unsigned));
        a.top =   L;
        b.top =   L;
        c.top =   L;

        for(j=0; j<L; j++)
            a.d[j]=0;

        for(j=0; j<L; j++)
            b.d[j]=0;

        for(j=0; j<L; j++)
            c.d[j]=0;

        A[i] = a;
        B[i] = b;
        C[i] = c;

    }

    for(i=0; i<number_of_keys; i++){
        cu_PEMs[i] = tmp;
        asprintf(&tmp_path, "%s/%d.pem", keys_directory, (i+1));
        get_u_bn_from_mod_PEM(tmp_path, &cu_PEMs[i]);
        //printf( "cu_PEM[%d]: %s\n", i, cu_bn_bn2hex(&cu_PEMs[i]));
    }

    for(i=0, k=0; i<number_of_keys; i++){
        for(j=(i+1); j<number_of_keys; j++, k++){
            A[k] = cu_PEMs[i];
            B[k] = cu_PEMs[j];
        }
    }



    /*cudaDeviceReset();
    cudaStatus = cudaMalloc((void**)&device_U_BN_A, number_of_comutations*sizeof(U_BN));    
    cudaStatus = cudaMalloc((void**)&device_U_BN_B, number_of_comutations*sizeof(U_BN));
    cudaStatus = cudaMalloc((void**)&device_U_BN_C, number_of_comutations*sizeof(U_BN));
    cudaStatus = cudaMemcpy(device_U_BN_A, A, number_of_comutations*sizeof(U_BN), cudaMemcpyHostToDevice);
    cudaStatus = cudaMemcpy(device_U_BN_B, B, number_of_comutations*sizeof(U_BN), cudaMemcpyHostToDevice);
    cudaStatus = cudaMemcpy(device_U_BN_C, C, number_of_comutations*sizeof(U_BN), cudaMemcpyHostToDevice);

    unsigned long *out;

    for(i = 0; i < number_of_comutations; i++) {

        cudaMalloc(&out, L*sizeof(unsigned));
        cudaMemcpy(out, A[i].d, L*sizeof(unsigned), cudaMemcpyHostToDevice);
        cudaMemcpy(&device_U_BN_A[i].d, &out, sizeof(void*), cudaMemcpyHostToDevice);
        cudaMalloc(&out, L*sizeof(unsigned));
        cudaMemcpy(out, B[i].d, L*sizeof(unsigned), cudaMemcpyHostToDevice);
        cudaMemcpy(&device_U_BN_B[i].d, &out, sizeof(void*), cudaMemcpyHostToDevice);
        cudaMalloc(&out, L*sizeof(unsigned));
        cudaMemcpy(out, C[i].d, L*sizeof(unsigned), cudaMemcpyHostToDevice);
        cudaMemcpy(&device_U_BN_C[i].d, &out, sizeof(void*), cudaMemcpyHostToDevice);
    }

    cudaPrintfInit();


    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start, 0);

    testKernel<<<((number_of_comutations + thread_per_block -1)/thread_per_block), thread_per_block>>>(device_U_BN_A, device_U_BN_B, device_U_BN_C, number_of_comutations);

    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    float time;
    cudaEventElapsedTime(&time, start, stop);
    printf("[GPU] Time elapsed in ms %fms\n", time);
    cudaPrintfDisplay(stdout, true);
    cudaPrintfEnd();

    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "\n testKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
        return 1;
    }

    cudaStatus = cudaMemcpy(A, device_U_BN_A, number_of_comutations*sizeof(U_BN), cudaMemcpyDeviceToHost);
    cudaStatus = cudaMemcpy(B, device_U_BN_B, number_of_comutations*sizeof(U_BN), cudaMemcpyDeviceToHost);
    cudaStatus = cudaMemcpy(C, device_U_BN_C, number_of_comutations*sizeof(U_BN), cudaMemcpyDeviceToHost);

    unsigned *array = (unsigned*)malloc(L*sizeof(unsigned));

    for(i = 0; i < number_of_comutations; ++i) {
        array = (unsigned*)malloc(L*sizeof(unsigned));
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
    }*/

    C = cu_classic_euclid(&A[i], &B[i]);

    /*char *buffer=NULL;
    buffer = (char *)malloc(sizeof(char) * (CU_BN_BYTES * C->top));
    for(i=0, k=0; i<number_of_keys; i++){
        for(j=(i+1); j<number_of_keys; j++, k++){
            if( !strcmp( "1", cu_bn_bn2hex(&C[k],buffer))){
                //printf("i: %d, j: %d, C[%d]: %s\n", j, i, k, cu_bn_bn2hex(&C[k]));
            	printf("1");
            }
        }
    }
    free(buffer);*/
    //cudaFree(out);
    //free(array);
    //cudaFree(device_U_BN_A);
    //cudaFree(device_U_BN_B);
    //cudaFree(device_U_BN_C);


    /*for(i=0; i<number_of_comutations; i++){
        free(A[i].d);
        free(B[i].d);
        free(C[i].d);
    }*/

    free(A);
    free(B);
    free(C);
    free(cu_PEMs);
    return (0);
}