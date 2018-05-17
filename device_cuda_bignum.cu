#include "device_cuda_bignum.h"


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

    if(!cu_dev_bn_ucmp(a, b)){  
        r->top = 1;
        r->d[0] = 0;
        return 1;
    }

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
    if( ((l << lb)/32) > 0 ) nwb = 1;

    if(nw || nwb){
        a->d = (unsigned*)realloc(a->d, (a->top + nw + nwb)*sizeof(unsigned)) ;
        a->d[a->top]=0;
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


__host__ __device__ U_BN *cu_dev_binary_gcd(U_BN *a, U_BN *b){
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



void OpenSSL_GCD(unsigned number_of_keys, unsigned key_size, char *keys_directory){

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


__global__ void orgEuclideanKernel(U_BN *A, U_BN *B, U_BN *C, unsigned n) {
    int i= blockIdx.x * blockDim.x + threadIdx.x;
    U_BN *TMP=NULL;

    if(i<n){
        TMP = cu_dev_classic_euclid(&A[i], &B[i]);
        C[i] = *TMP;
    }
}

__global__ void binEuclideanKernel(U_BN *A, U_BN *B, U_BN *C, unsigned n) {
    int i= blockIdx.x * blockDim.x + threadIdx.x;
    U_BN *TMP=NULL;

    if(i<n){
        TMP = cu_dev_binary_gcd(&A[i], &B[i]);
        C[i] = *TMP;
    }
}

__global__ void fastBinaryKernel(U_BN *A, U_BN *B, U_BN *C, unsigned n) {
    int i= blockIdx.x * blockDim.x + threadIdx.x;
    U_BN *TMP=NULL;

    if(i<n){
        TMP = cu_dev_fast_binary_euclid(&A[i], &B[i]);
        C[i] = *TMP;
    }
}