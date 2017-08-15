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

__global__
void saxpy(int n, float a, float *x, float *y)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  if (i < n) y[i] = a*x[i] + y[i];
}

int main(void)
{

	int M = 1<<20;
	float *x, *y, *d_x, *d_y;
	x = (float*)malloc(M*sizeof(float));
	y = (float*)malloc(M*sizeof(float));

	cudaMalloc(&d_x, M*sizeof(float)); 
	cudaMalloc(&d_y, M*sizeof(float));

	for (i = 0; i < M; i++) {
	x[i] = 1.0f;
	y[i] = 2.0f;
	}

	cudaMemcpy(d_x, x, M*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_y, y, M*sizeof(float), cudaMemcpyHostToDevice);

	// Perform SAXPY on 1M elements
	saxpy<<<(M+255)/256, 256>>>(M, 2.0f, d_x, d_y);

	cudaMemcpy(y, d_y, M*sizeof(float), cudaMemcpyDeviceToHost);

	float maxError = 0.0f;
	for (i = 0; i < M; i++)
	maxError = max(maxError, abs(y[i]-4.0f));
	printf("Max error: %f\n", maxError);

	cudaFree(d_x);
	cudaFree(d_y);
	free(x);
	free(y);


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

