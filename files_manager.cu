#include "files_manager.h"

void print_mod_from_pem_file(char * filePath){

	EVP_PKEY* pPubKey  = NULL;
    FILE*     pemFile    = NULL;

	if( !( (pemFile = fopen(filePath, "rt") ) && ( pPubKey = PEM_read_PUBKEY(pemFile,NULL,NULL,NULL) ) ) ) {
        fprintf(stderr,"Cannot read \"public key\".\n");
	}

	RSA* rsa = EVP_PKEY_get1_RSA(pPubKey);
	//BN_print_fp(stdout, rsa->n);
	printf("Public modulus:\n");
	fprintf(stdout, "%s", BN_bn2dec(rsa->n));
	printf("\n");
}