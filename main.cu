#include "GCD.cuh"

int N  = 100;	
EVP_PKEY* pPubKey_first  = NULL;
FILE*     pemFile_first    = NULL;
EVP_PKEY* pPubKey_second  = NULL;
FILE*     pemFile_second    = NULL;
BIGNUM* first_modulus = NULL;
BIGNUM* second_modulus = NULL;
BIGNUM* GCD_result = NULL;
RSA* first_rsa = NULL;
RSA* second_rsa = NULL;
unsigned int i, j;
char *nr_first = NULL;
char *nr_second = NULL;

void init_variables(){
	pPubKey_first  = EVP_PKEY_new();
	pPubKey_second  = EVP_PKEY_new();
	first_modulus = BN_new();
	second_modulus = BN_new();
	GCD_result = BN_new();
}

void free_variables(){
	BN_free(first_modulus);
	BN_free(second_modulus);
	BN_free(GCD_result);
	EVP_PKEY_free(pPubKey_first);
    EVP_PKEY_free(pPubKey_second);
    fclose(pemFile_first);
    fclose(pemFile_second);
}

/*__global__ void GCD(int n, float a, float *x, float *y){
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  if (i < n) y[i] = a*x[i] + y[i];
}*/



__global__ void
GCD(const float *A, const float *B, float *C, int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < numElements)
    {
        euclid_modulo(GCD_result, first_modulus, second_modulus);
    }
}



int main(void)
{
	// Error code to check return values for CUDA calls
    cudaError_t err = cudaSuccess;
	int numElements = 50000;
	size_t size = numElements * sizeof(float);
	printf("[Vector addition of %d elements]\n", numElements);

    // Allocate the host input vector A
    float *h_A = (float *)malloc(size);

    // Allocate the host input vector B
    float *h_B = (float *)malloc(size);

    // Allocate the host output vector C
    float *h_C = (float *)malloc(size);

    // Verify that allocations succeeded
    if (h_A == NULL || h_B == NULL || h_C == NULL)
    {
        fprintf(stderr, "Failed to allocate host vectors!\n");
        exit(EXIT_FAILURE);
    }

    // Initialize the host input vectors
    for (i = 0; i < numElements; ++i)
    {
        h_A[i] = rand()/(float)RAND_MAX;
        h_B[i] = rand()/(float)RAND_MAX;
    }

    // Allocate the device input vector A
    float *d_A = NULL;
    err = cudaMalloc((void **)&d_A, size);

    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device vector A (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    // Allocate the device input vector B
    float *d_B = NULL;
    err = cudaMalloc((void **)&d_B, size);

    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device vector B (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    // Allocate the device output vector C
    float *d_C = NULL;
    err = cudaMalloc((void **)&d_C, size);

    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device vector C (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    // Copy the host input vectors A and B in host memory to the device input vectors in
    // device memory
    printf("Copy input data from the host memory to the CUDA device\n");
    err = cudaMemcpy(d_A, h_A, size, cudaMemcpyHostToDevice);

    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to copy vector A from host to device (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(d_B, h_B, size, cudaMemcpyHostToDevice);

    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to copy vector B from host to device (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    // Launch the Vector Add CUDA Kernel
    int threadsPerBlock = 512;
    int blocksPerGrid =(numElements + threadsPerBlock - 1) / threadsPerBlock;
    printf("CUDA kernel launch with %d blocks of %d threads\n", blocksPerGrid, threadsPerBlock);
    GCD<<<blocksPerGrid, threadsPerBlock>>>(d_A, d_B, d_C, numElements);
    err = cudaGetLastError();

    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to launch vectorAdd kernel (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    // Copy the device result vector in device memory to the host result vector
    // in host memory.
    printf("Copy output data from the CUDA device to the host memory\n");
    err = cudaMemcpy(h_C, d_C, size, cudaMemcpyDeviceToHost);

    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to copy vector C from device to host (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    // Verify that the result vector is correct
    for (i = 0; i < numElements; ++i)
    {
        if (fabs(h_A[i] + h_B[i] - h_C[i]) > 1e-5)
        {
            fprintf(stderr, "Result verification failed at element %d!\n", i);
            exit(EXIT_FAILURE);
        }
    }

    printf("Test PASSED\n");

    // Free device global memory
    err = cudaFree(d_A);

    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to free device vector A (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    err = cudaFree(d_B);

    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to free device vector B (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    err = cudaFree(d_C);

    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to free device vector C (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    // Free host memory
    free(h_A);
    free(h_B);
    free(h_C);

    // Reset the device and exit
    // cudaDeviceReset causes the driver to clean up all state. While
    // not mandatory in normal operation, it is good practice.  It is also
    // needed to ensure correct operation when the application is being
    // profiled. Calling cudaDeviceReset causes all profile data to be
    // flushed before the application exits
    err = cudaDeviceReset();

    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to deinitialize the device! error=%s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    printf("Done\n");


	clock_t start = clock();
	for(i=1;i<=N;i++){
		for(j=(i+1);j<=N;j++){
			init_variables();
			asprintf(&nr_first, "keys_and_messages/%d.pem", i);
			asprintf(&nr_second, "keys_and_messages/%d.pem", j);


			if((pemFile_first = fopen(nr_first,"rt")) && (pPubKey_first = PEM_read_PUBKEY(pemFile_first,NULL,NULL,NULL))){
		        //fprintf(stderr,"Public key read.\n");
		    }
		    else
		    {
		        fprintf(stderr,"Cannot read \"public key\".\n");

		    }

		    if((pemFile_second = fopen(nr_second,"rt")) && (pPubKey_second = PEM_read_PUBKEY(pemFile_second,NULL,NULL,NULL))){
		        //fprintf(stderr,"Public key read.\n");
		    } else {
		        fprintf(stderr,"Cannot read \"public key\".\n");

		    }

			first_rsa = EVP_PKEY_get1_RSA(pPubKey_first);
			second_rsa = EVP_PKEY_get1_RSA(pPubKey_second);
			first_modulus = first_rsa->n;
			second_modulus = second_rsa->n;
			euclid_modulo(GCD_result, first_modulus, second_modulus);
			if(!BN_is_one(GCD_result)){
				printf("%s and %s:\n", nr_first, nr_second);
				fprintf(stdout, "GCD result:\n%s\n", BN_bn2dec(GCD_result));
			}

			free_variables();
		}
	}
	clock_t stop = clock();
	double elapsed = (double)(stop - start) * 1000.0 / CLOCKS_PER_SEC;
	printf("Time elapsed in ms: %f\n", elapsed);

	return 0;
}

