#include "GCD.h"


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


int main(){

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
				fprintf(stdout, "GCD result:\n %s\n", BN_bn2dec(GCD_result));
			}
			
			free_variables();
		}
	}
}

