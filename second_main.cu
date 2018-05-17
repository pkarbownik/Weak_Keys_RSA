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
    unsigned i, j, sum = 0;
    unsigned k = 0, l=0;
    U_BN   *A, *B, *C;
    U_BN   *device_U_BN_A, *device_U_BN_B, *device_U_BN_C;
    U_BN   *cu_PEMs;
    char *tmp_path;
    cudaError_t cudaStatus;

    unit_test();
    CPU_computation(number_of_keys, key_size, keys_directory);


    A    = (U_BN*)malloc(number_of_comutations*sizeof(U_BN));
    B    = (U_BN*)malloc(number_of_comutations*sizeof(U_BN));
    C    = (U_BN*)malloc(number_of_comutations*sizeof(U_BN));
    cu_PEMs = (U_BN*)malloc(number_of_keys*sizeof(U_BN));

    for(i=0; i<number_of_comutations; i++){

        U_BN a;
        U_BN b;
        U_BN c;
        U_BN d;
        a.d = (unsigned*)malloc(L*sizeof(unsigned));
        b.d = (unsigned*)malloc(L*sizeof(unsigned));
        c.d = (unsigned*)malloc(L*sizeof(unsigned));
        d.d = (unsigned*)malloc(L*sizeof(unsigned));
        a.top =   L;
        b.top =   L;
        c.top =   L;
        d.top =   L;

        for(j=0; j<L; j++)
            a.d[j]=0;

        for(j=0; j<L; j++)
            b.d[j]=0;

        for(j=0; j<L; j++)
            c.d[j]=0;

        for(j=0; j<L; j++)
            c.d[j]=0;

        A[i] = a;
        B[i] = b;
        C[i] = c;
        C[i] = d;

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

    cudaPrintfInit();

    float time;
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start, 0);

    testKernel<<<((number_of_comutations + thread_per_block -1)/thread_per_block), thread_per_block>>>(device_U_BN_A, device_U_BN_B, device_U_BN_C, number_of_comutations);

    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
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
    }

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
    //free(C);
    free(cu_PEMs);
    return (0);
}