#include "cuPrintf.cu"
#include "test.h"
#include "cuda_bignum.h"
#include "files_manager.h"
//#include <openssl/bn.h>

#define N 4950
#define KEYS 100
#define THREADS_PER_BLOCK 512

#if __CUDA_ARCH__ < 200     //Compute capability 1.x architectures
#define CUPRINTF cuPrintf
#else                       //Compute capability 2.x architectures
#define CUPRINTF(fmt, ...) printf("[%d, %d]:\t" fmt, \
                                  blockIdx.y*gridDim.x+blockIdx.x,\
                                  threadIdx.z*blockDim.x*blockDim.y+threadIdx.y*blockDim.x+threadIdx.x,\
                                  __VA_ARGS__)
#endif

#define MAIN_COMPUTATIONS 1


__host__ __device__ int cu_dev_BN_ucmp(const U_BN *a, const U_BN *b){

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

__host__ __device__ int cu_dev_bn_lshift(U_BN *a, unsigned n){

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
        //a->d[a->top]=0;
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


__host__ __device__ U_BN *cu_dev_euclid(U_BN *a, U_BN *b){
    U_BN *t = NULL;
    unsigned shifts = 0;

    if (cu_dev_BN_ucmp(a, b) < 0) {
        t = a;
        a = b;
        b = t;
    }

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


__global__ void testKernel(U_BN *A, U_BN *B, U_BN *C, unsigned n) {
    int i= blockIdx.x * blockDim.x + threadIdx.x;
    //char Hex[] = "0123456789ABCDEF";

    U_BN *TMP=NULL;
    //int l;
    //cu_dev_bn_usub(&A[i], &B[i], &C[i]);
    //cu_dev_bn_lshift(&C[i], 4);
    
    //TMP = cu_dev_euclid(&A[i], &B[i]);
    //CUPRINTF("testKernel entrance by the global threadIdx= %d value: %u\n", i , TMP->d[0]);
    //CUPRINTF("testKernel entrance by the global threadIdx= %d value: %u\n", i , TMP->d[1]);
    //C[i] = *TMP;

    //int k, j, v, z = 0;




    if(i<n){
        //for(l=0; l<A->top; l++){
        //    cu_dev_bn_usub(&A[i], &B[i], &C[i]);
        //}
        TMP = cu_dev_euclid(&A[i], &B[i]);
        C[i] = *TMP;
        /*for (k = A[i].top - 1; k >= 0; k--) {
            for (j = CU_BN_BITS2 - 4; j >= 0; j -= 4) {
                v = ((unsigned)(A[i].d[k] >> j)) & 0x0f;
                if (z || (v != 0)) {
                    CUPRINTF("%c\n",(Hex[v & 0x0f]));
                    z = 1;
                }
            }
        }
        CUPRINTF("top: %u\n", A[i].top);*/
    }
}

int main(void){
    int L = 35;
    int i, j;
    unsigned k = 0;
    unsigned n = N;

    unit_test(); //check all host bn functions

    U_BN   *A;
    U_BN   *device_U_BN_A;
    U_BN   *B;
    U_BN   *device_U_BN_B;
    U_BN   *C;
    U_BN   *device_U_BN_C;
    U_BN   *cu_PEMs;
    BIGNUM   *PEMs;
    char *tmp_path;

    cudaError_t cudaStatus;

    A    = (U_BN*)malloc(N*sizeof(U_BN));
    B    = (U_BN*)malloc(N*sizeof(U_BN));
    C    = (U_BN*)malloc(N*sizeof(U_BN));
    cu_PEMs = (U_BN*)malloc(KEYS*sizeof(U_BN));
    PEMs = (BIGNUM*)malloc(KEYS*sizeof(BIGNUM));

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


        //cu_BN_dec2bn(&A[i], "132009813808533392577123110438741884286561400398429860761027919959189196549797215586297852825375342475728679074489933320371765026814849875692023263110656924146683347962741534495754902097930935910070831755220321417369411370818762253940133993629997648473607090782039210687337530507010114741418840031542303031081");
        //cu_BN_dec2bn(&B[i], "135472400918611757666622822789636901038207639581474006488496906937544113899819968264216470405393313301250508761651903965883874352772699986982247519612608840409853757718079903608120168842687889231898954817245684707914621259848016658887023606975529849256590875282759156328281549546230980205644358325222571914637");
        //cu_BN_dec2bn(&C[i], "136269636317215868658126726142543242028128679787201513621377420299644359247151157885793577216543689892988935986714087409150506883630386841292060595217129497897100280678153687017820663980404875865314501020301179267627899307057160787226214936662085381326053730017478234531591680965138499420169342895677786825700");
    }

    U_BN tmp;
    EVP_PKEY* pPubKey  = NULL;
    FILE*     pemFile    = NULL;
    RSA* rsa;

    tmp.d = (unsigned*)malloc(L*sizeof(unsigned));
    tmp.top = L;

    for(int j=0; j<L; j++){
        tmp.d[j]=0;
    }

   	for(i=0; i<KEYS; i++){
        cu_PEMs[i] = tmp;
		asprintf(&tmp_path, "keys_and_messages/%d.pem", (i+1));
    	get_u_bn_from_mod_PEM(tmp_path, &cu_PEMs[i]);

        //printf( "cu_PEM[%d]: %s\n", i, cu_bn_bn2hex(&cu_PEMs[i]));

        if( !( (pemFile = fopen(tmp_path, "rt") ) && ( pPubKey = PEM_read_PUBKEY(pemFile,NULL,NULL,NULL) ) ) ) {
            fprintf(stderr,"Cannot read \"public key\".\n");
        }

        rsa = EVP_PKEY_get1_RSA(pPubKey);
        //printf( "PEM[%d]: %s\n", i, BN_bn2dec(rsa->n));
        PEMs[i] = *(rsa->n);
        fclose(pemFile);
        	//BN_dec2bn(&cu_PEMs_ssl[i], tmp_path);
    }

    free(tmp.d);


	for(i=0, k=0; i<KEYS; i++){
		for(j=(i+1); j<KEYS; j++, k++){
			A[k] = cu_PEMs[i];
			B[k] = cu_PEMs[j];
		}
	}

    BN_CTX *ctx;
    ctx = BN_CTX_new();
    BIGNUM *r;
    r=BN_new();
    for(i=0, k=0; i<KEYS; i++){
		for(j=(i+1); j<KEYS; j++, k++){
            BN_gcd(r, &PEMs[i], &PEMs[j], ctx);
            if(!BN_is_one(r)){
    			//printf( "A[%d]: %s\nB[%d]: %s\neuclid: %s\n\neuclid_dev_cu: %s\n\n\n", i, BN_bn2hex(&PEMs[i]), j, BN_bn2hex(&PEMs[j]), BN_bn2hex(r), cu_bn_bn2hex(cu_dev_euclid(&A[k], &B[k])));
            }   
		}
	}
    
    BN_CTX_free(ctx);
    BN_free(r);
	printf("k: %d\n", k);
    //L=A[0].top;
    //L=B[0].top;
    //L=C[0].top;


    cudaDeviceReset();
    cudaStatus = cudaMalloc((void**)&device_U_BN_A, N*sizeof(U_BN));    
    cudaStatus = cudaMalloc((void**)&device_U_BN_B, N*sizeof(U_BN));
    cudaStatus = cudaMalloc((void**)&device_U_BN_C, N*sizeof(U_BN));
    cudaStatus = cudaMemcpy(device_U_BN_A, A, N*sizeof(U_BN), cudaMemcpyHostToDevice);
    cudaStatus = cudaMemcpy(device_U_BN_B, B, N*sizeof(U_BN), cudaMemcpyHostToDevice);
    cudaStatus = cudaMemcpy(device_U_BN_C, C, N*sizeof(U_BN), cudaMemcpyHostToDevice);

	unsigned long *out;

    for(int i = 0; i < N; ++i) {

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
        //free(A[i].d);
        //free(B[i].d);
        //free(C[i].d);

    }

    cudaPrintfInit();
    testKernel<<<((N + THREADS_PER_BLOCK -1)/THREADS_PER_BLOCK), THREADS_PER_BLOCK>>>(device_U_BN_A, device_U_BN_B, device_U_BN_C, n);//to test and see on a sigle thread
    cudaDeviceSynchronize();
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

    unsigned *array = (unsigned*)malloc(L*sizeof(unsigned));

    for(int i = 0; i < N; ++i) {
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
    }

	for(i=0, k=0; i<KEYS; i++){
		for(j=(i+1); j<KEYS; j++, k++){
			if( strcmp( "1", cu_bn_bn2hex(&C[k]))){
				printf("i: %d, j: %d, C[%d]: %s\n", j, i, k, cu_bn_bn2hex(&C[k]));
            }
		}
	}
    cudaFree(out);
	free(array);
    cudaFree(device_U_BN_A);
    cudaFree(device_U_BN_B);
    cudaFree(device_U_BN_C);
    return (0);
}