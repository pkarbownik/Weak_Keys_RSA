#include "device_cuda_bignum.h"
/*#include "cuPrintf.cu"

#if __CUDA_ARCH__ < 200     //Compute capability 1.x architectures
#define CUPRINTF cuPrintf
#else                       //Compute capability 2.x architectures
#define CUPRINTF(fmt, ...) printf("[%d, %d]:\t" fmt, \
                                  blockIdx.y*gridDim.x+blockIdx.x,\
                                  threadIdx.z*blockDim.x*blockDim.y+threadIdx.y*blockDim.x+threadIdx.x,\
                                  __VA_ARGS__)
#endif*/

typedef enum {
    EUCLIDEAN=0,
    BINARY_EUCLIDEAN,
    FAST_BINARY_EUCLIDEAN,
    UNKNOWN
} algorithms;

algorithms set_enum_algorithm(char * algorithm){
    if(!strcmp( "euclid", algorithm)){
        return EUCLIDEAN;
    } else if(!strcmp( "binary", algorithm)) {
        return BINARY_EUCLIDEAN;
    } else if(!strcmp( "fast", algorithm)) {
        return FAST_BINARY_EUCLIDEAN;
    } else {
        return UNKNOWN;
    }
}

int main(int argc, char* argv[]){
	
    unsigned number_of_keys;
    unsigned key_size;
    unsigned thread_per_block;
    unsigned number_of_comutations;
    char *keys_directory;
    int counter;
    algorithms gcd_kind;

    if(argc==6) {
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
                case 5:
                    printf("\nkind of algorithm argv[%d]: %s\n",counter,argv[counter]);
                    gcd_kind = set_enum_algorithm(argv[counter]);
                    break;
                default:
                    break;
            }
        }
    } else {
        printf("\nFind weak keys\n\rUsage:\n\r ./GCD_RSA number_of_keys key_size threads_per_block directory_name kind_of_algorithm\n\rAlgorithms:\n\r\t-\"euclid\"\n\r\t-\"binary\"\n\r\t-\"fast\"\n\r");
        return 0;
    }

    // simplified binomial coefficient
    number_of_comutations=((number_of_keys/2)*(number_of_keys-1));

    U_BN tmp;
    int L = ((key_size+31) / (8*sizeof(unsigned)));
    unsigned i, j;
    unsigned k = 0, l=0;
    U_BN   *A, *B, *C;
    U_BN   *device_U_BN_A, *device_U_BN_B, *device_U_BN_C;
    U_BN   *cu_PEMs;
    char *tmp_path;
    cudaError_t cudaStatus;

    //unit_test();
    //OpenSSL_GCD(number_of_keys, key_size, keys_directory);


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

        U_BN d;
        d.d = (unsigned*)malloc(L*sizeof(unsigned));
        d.top =   L;

        for(j=0; j<L; j++)
            d.d[j]=0;

        cu_PEMs[i] = d;

    }

    for(i=0; i<number_of_keys; i++){
        cu_PEMs[i] = tmp;
        asprintf(&tmp_path, "%s/%d.pem", keys_directory, (i+1));
        get_u_bn_from_mod_PEM(tmp_path, &cu_PEMs[i]);
    }

    for(i=0, k=0; i<number_of_keys; i++){
        for(j=(i+1); j<number_of_keys; j++, k++){
            A[k].top = cu_PEMs[i].top;
            B[k].top = cu_PEMs[j].top;
            for(l=0;l<L;l++){
                A[k].d[l] = cu_PEMs[i].d[l];
                B[k].d[l] = cu_PEMs[j].d[l];
            }
        }
    }


    int sum=0;
    clock_t start = clock();
    switch(gcd_kind){
        case EUCLIDEAN:
            for(i=0, k=0; i<number_of_keys; i++){
                for(j=(i+1); j<number_of_keys; j++, k++){
                    if( strcmp( "1", cu_bn_bn2hex(cu_dev_classic_euclid(&A[k], &B[k])))){
                        printf("[CPU] Euclidean Weak key: %s\n", cu_bn_bn2hex(cu_dev_classic_euclid(&A[k], &B[k])));
                        sum+=1;
                    }
                }
            }
            break;
        case BINARY_EUCLIDEAN:
            for(i=0, k=0; i<number_of_keys; i++){
                for(j=(i+1); j<number_of_keys; j++, k++){
                    if( strcmp( "1", cu_bn_bn2hex(cu_dev_binary_gcd(&A[k], &B[k])))){
                        printf("[CPU] Binary Weak key: %s\n", cu_bn_bn2hex(cu_dev_binary_gcd(&A[k], &B[k])));
                        sum+=1;
                    }
                }
            }
            break;
        case FAST_BINARY_EUCLIDEAN:
            for(i=0, k=0; i<number_of_keys; i++){
                for(j=(i+1); j<number_of_keys; j++, k++){
                    if( strcmp( "1", cu_bn_bn2hex(cu_dev_fast_binary_euclid(&A[k], &B[k])))){
                        printf("[CPU] Fast Weak key: %s\n", cu_bn_bn2hex(cu_dev_fast_binary_euclid(&A[k], &B[k])));
                        sum+=1;
                    }
                }
            }
            break;
        default:
            printf("[CPU] Unknown GCD algorithm");
            break;
    }
    clock_t stop = clock();
    double elapsed = (double)(stop - start) * 1000.0 / CLOCKS_PER_SEC;
    printf("[CPU] Time elapsed in ms: %f\n", elapsed);
    printf("[CPU] Weak keys: %d\n", sum);


    cudaDeviceReset();
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

    float time;
    cudaEvent_t start_cu, stop_cu;
    cudaEventCreate(&start_cu);
    cudaEventCreate(&stop_cu);
    cudaEventRecord(start_cu, 0);
    switch(gcd_kind){
        case EUCLIDEAN:
            printf("[GPU] Euclidean algorithm\n");
            orgEuclideanKernel<<<((number_of_comutations + thread_per_block -1)/thread_per_block), thread_per_block>>>(device_U_BN_A, device_U_BN_B, device_U_BN_C, number_of_comutations);
            break;
        case BINARY_EUCLIDEAN:
            printf("[GPU] Binary algorithm\n");
            binEuclideanKernel<<<((number_of_comutations + thread_per_block -1)/thread_per_block), thread_per_block>>>(device_U_BN_A, device_U_BN_B, device_U_BN_C, number_of_comutations);
            break;
        case FAST_BINARY_EUCLIDEAN:
            printf("[GPU] Fast Binary algorithm\n");
            fastBinaryKernel<<<((number_of_comutations + thread_per_block -1)/thread_per_block), thread_per_block>>>(device_U_BN_A, device_U_BN_B, device_U_BN_C, number_of_comutations);
            break;
        default:
            printf("[GPU] Unknown GCD algorithm\n");
            break;
    }

    cudaEventRecord(stop_cu, 0);
    cudaEventSynchronize(stop_cu);
    cudaEventElapsedTime(&time, start_cu, stop_cu);
    printf("[GPU] Time elapsed in ms %fms\n", time);

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
    }

    sum=0;
    for (k = 0; k < number_of_comutations; k++){
        if( strcmp( "1", cu_bn_bn2hex(&C[k]))){
            sum += 1;
        }
    }
    printf("[GPU] Weak keys: %d\n", sum);

    cudaFree(out);
    free(array);
    cudaFree(device_U_BN_A);
    cudaFree(device_U_BN_B);
    cudaFree(device_U_BN_C);
    free(A);
    free(B);
    free(C);
    free(cu_PEMs);
    return (0);
}