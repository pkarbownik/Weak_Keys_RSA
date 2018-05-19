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
    unsigned key_size, combinations;
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
    combinations=(number_of_keys*number_of_keys);

    U_BN tmp;
    int L = ((key_size+31) / (8*sizeof(unsigned)));
    unsigned i, j, l;
    unsigned k = 0;
    U_BN   *device_U_BN_A, *device_U_BN_B, *device_U_BN_R;
    U_BN   *KEYS, *R;
    char *tmp_path;
    U_BN   *A, *B;
    cudaError_t cudaStatus;

    //unit_test();
    OpenSSL_GCD(number_of_keys, key_size, keys_directory);

    R = (U_BN*)malloc(number_of_comutations*sizeof(U_BN));
    KEYS = (U_BN*)malloc(number_of_keys*sizeof(U_BN));
    A = (U_BN*)malloc(number_of_keys*sizeof(U_BN));
    B = (U_BN*)malloc(number_of_keys*sizeof(U_BN));

    for(i=0; i<number_of_keys; i++){

        U_BN a;
        U_BN b;
        U_BN keys;
        a.d = (unsigned*)malloc(L*sizeof(unsigned));
        b.d = (unsigned*)malloc(L*sizeof(unsigned));
        keys.d = (unsigned*)malloc(L*sizeof(unsigned));
        a.top =   L;
        b.top =   L;
        keys.top =   L;

        for(j=0; j<L; j++)
            a.d[j]=0;

        for(j=0; j<L; j++)
            b.d[j]=0;

        for(j=0; j<L; j++)
            keys.d[j]=0;

        A[i] = a;
        B[i] = b;
        KEYS[i] = keys;
    }


    for(i=0; i<number_of_comutations; i++){

        U_BN r;
        r.d = (unsigned*)malloc(L*sizeof(unsigned));
        r.top =   L;

        for(j=0; j<L; j++)
            r.d[j]=0;

        R[i] = r;
    }

    for(i=0; i<number_of_keys; i++){
        //KEYS[i] = tmp;
        asprintf(&tmp_path, "%s/%d.pem", keys_directory, (i+1));
        get_u_bn_from_mod_PEM(tmp_path, &A[i]);
        get_u_bn_from_mod_PEM(tmp_path, &B[i]);


    }


    /*int l;
    for(k=0; k<number_of_keys; k++){
        A[k].top = KEYS[k].top;
        B[k].top = KEYS[k].top;
        for(l=0;l<L;l++){
            A[k].d[l] = KEYS[k].d[l];
        }

        for(l=0;l<L;l++){
            B[k].d[l] = KEYS[k].d[l];
        }

    	printf("[CPU] A: %s\n", cu_bn_bn2hex(&A[k]));
    	printf("[CPU] B: %s\n", cu_bn_bn2hex(&B[k]));
    }*/

    U_BN *tmp_A = NULL;
    U_BN *tmp_B = NULL;
    tmp_A = (U_BN*)malloc(sizeof(U_BN));
    tmp_B = (U_BN*)malloc(sizeof(U_BN));

    tmp_A->d = (unsigned*)malloc(32*sizeof(unsigned));
    tmp_B->d = (unsigned*)malloc(32*sizeof(unsigned));

    int sum=0;
    clock_t start = clock();
    switch(gcd_kind){
        case EUCLIDEAN:
            for(i=0; i<number_of_keys; i++){
                for(j=(i+1); j<number_of_keys; j++){
				    tmp_A->top = A[i].top;
				    tmp_B->top = B[j].top;
				    for(l=0; l<L; l++){
				        tmp_A->d[l] = A[i].d[l];
				        tmp_B->d[l] = B[j].d[l];
    				}
                    if( strcmp( "1", cu_bn_bn2hex(cu_dev_classic_euclid(tmp_A, tmp_B)))){
                        printf("[CPU] Euclidean Weak key with selection: %s\n", cu_bn_bn2hex(cu_dev_classic_euclid(tmp_A, tmp_B)));
                        sum+=1;
                    }
                }
            }
            break;
        case BINARY_EUCLIDEAN:
            for(i=0; i<number_of_keys; i++){
                for(j=(i+1); j<number_of_keys; j++){
				    tmp_A->top = A[i].top;
				    tmp_B->top = B[j].top;
				    for(l=0; l<L; l++){
				        tmp_A->d[l] = A[i].d[l];
				        tmp_B->d[l] = B[j].d[l];
    				}
    				cu_dev_binary_gcd(tmp_A, tmp_B);
                    if( strcmp( "1", cu_bn_bn2hex(cu_dev_binary_gcd(tmp_A, tmp_B)))){
                        //printf("[CPU] Binary Weak key with selection: %s\n", cu_bn_bn2hex(cu_dev_binary_gcd(tmp_A, tmp_B)));
                        sum+=1;
                    }
                }
            }
            break;
        case FAST_BINARY_EUCLIDEAN:
            for(i=0, k=0; i<number_of_keys; i++){
                for(j=(i+1); j<number_of_keys; j++, k++){

				    tmp_A->top = A[i].top;
				    tmp_B->top = B[j].top;
				    for(l=0; l<L; l++){
				        tmp_A->d[l] = A[i].d[l];
				        tmp_B->d[l] = B[j].d[l];
    				}

                    if( strcmp( "1", cu_bn_bn2hex(cu_dev_fast_binary_euclid(tmp_A, tmp_B)))){
                        //printf("[CPU] Fast Weak key: %s\n", cu_bn_bn2hex(cu_dev_fast_binary_euclid(tmp_A, tmp_B)));
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
    cudaStatus = cudaMalloc((void**)&device_U_BN_A, number_of_keys*sizeof(U_BN));    
    cudaStatus = cudaMalloc((void**)&device_U_BN_B, number_of_keys*sizeof(U_BN));
    cudaStatus = cudaMalloc((void**)&device_U_BN_R, number_of_comutations*sizeof(U_BN));
    cudaStatus = cudaMemcpy(device_U_BN_A, A, number_of_keys*sizeof(U_BN), cudaMemcpyHostToDevice);
    cudaStatus = cudaMemcpy(device_U_BN_B, B, number_of_keys*sizeof(U_BN), cudaMemcpyHostToDevice);
    cudaStatus = cudaMemcpy(device_U_BN_R, R, number_of_comutations*sizeof(U_BN), cudaMemcpyHostToDevice);

    unsigned long *out;

    for(i = 0; i < number_of_keys; i++) {

        cudaMalloc(&out, L*sizeof(unsigned));
        cudaMemcpy(out, A[i].d, L*sizeof(unsigned), cudaMemcpyHostToDevice);
        cudaMemcpy(&device_U_BN_A[i].d, &out, sizeof(void*), cudaMemcpyHostToDevice);
        cudaMalloc(&out, L*sizeof(unsigned));
        cudaMemcpy(out, B[i].d, L*sizeof(unsigned), cudaMemcpyHostToDevice);
        cudaMemcpy(&device_U_BN_B[i].d, &out, sizeof(void*), cudaMemcpyHostToDevice);
    }

    for(i = 0; i < number_of_comutations; i++) {

        cudaMalloc(&out, L*sizeof(unsigned));
        cudaMemcpy(out, R[i].d, L*sizeof(unsigned), cudaMemcpyHostToDevice);
        cudaMemcpy(&device_U_BN_R[i].d, &out, sizeof(void*), cudaMemcpyHostToDevice);
    }

    float time;
    cudaEvent_t start_cu, stop_cu;
    cudaEventCreate(&start_cu);
    cudaEventCreate(&stop_cu);
    cudaEventRecord(start_cu, 0);
    /*switch(gcd_kind){
        case EUCLIDEAN:
            printf("[GPU] Euclidean algorithm with selection\n");
            orgEuclideanKernel_with_selection<<<((combinations + thread_per_block -1)/thread_per_block), thread_per_block>>>(device_U_BN_A, device_U_BN_B, device_U_BN_R, number_of_comutations, number_of_keys);
            break;
        case BINARY_EUCLIDEAN:
            printf("[GPU] Binary algorithm with selection\n");
            binEuclideanKernel_with_selection<<<((combinations + thread_per_block -1)/thread_per_block), thread_per_block>>>(device_U_BN_A, device_U_BN_B, device_U_BN_R, number_of_comutations, number_of_keys);
            break;
        case FAST_BINARY_EUCLIDEAN:
            printf("[GPU] Fast Binary algorithm with selection\n");
            fastBinaryKernel_with_selection<<<((combinations + thread_per_block -1)/thread_per_block), thread_per_block>>>(device_U_BN_A, device_U_BN_B, device_U_BN_R, number_of_comutations, number_of_keys);
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

    cudaStatus = cudaMemcpy(R, device_U_BN_R, number_of_comutations*sizeof(U_BN), cudaMemcpyDeviceToHost);

    unsigned *array = (unsigned*)malloc(L*sizeof(unsigned));

    for(i = 0; i < number_of_comutations; ++i) {
        array = (unsigned*)malloc(L*sizeof(unsigned));
        cudaMemcpy(array, R[i].d, L*sizeof(unsigned), cudaMemcpyDeviceToHost);
        R[i].d = array;

    }

    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "\n testKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
        return 1;
    }

    sum=0;
    for (k = 0; k < number_of_comutations; k++){
        if( strcmp( "1", cu_bn_bn2hex(&R[k]))){
            sum += 1;
        }
    }
    printf("[GPU] Weak keys: %d\n", sum);

    cudaFree(out);
    free(array);
    cudaFree(device_U_BN_A);
    cudaFree(device_U_BN_B);
    cudaFree(device_U_BN_R);*/
    free(A);
    free(B);
    free(KEYS);
    free(R);
    return (0);
}