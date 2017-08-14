#include "GCD.h"


#include <openssl/rsa.h>
#include <openssl/pem.h>
#include <openssl/bn.h>


const int kBits = 1024;
const int kExp = 3;


int main(){

    EVP_PKEY* pPubKey  = NULL;
    FILE*     pemFile    = NULL;


    if((pemFile = fopen("keys_and_messages/1.pem","rt")) && 
       (pPubKey = PEM_read_PUBKEY(pemFile,NULL,NULL,NULL)))
    {
        fprintf(stderr,"Public key read.\n");
    }
    else
    {
        fprintf(stderr,"Cannot read \"public key\".\n");

    }

	RSA* rsa = EVP_PKEY_get1_RSA(pPubKey);
	//BN_print_fp(stdout, rsa->n);
	printf("Public modulus:\n");
	fprintf(stdout, "%s", BN_bn2dec(rsa->n));
	printf("\n");
	//function();

}


