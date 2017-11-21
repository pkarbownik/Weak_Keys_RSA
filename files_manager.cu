#include "files_manager.h"

void print_mod_from_pem_file(char * filePath){

	EVP_PKEY* pPubKey  = NULL;
    FILE*     pemFile    = NULL;

	if( !( (pemFile = fopen(filePath, "rt") ) && ( pPubKey = PEM_read_PUBKEY(pemFile,NULL,NULL,NULL) ) ) ) {
        fprintf(stderr,"Cannot read \"public key\".\n");
        return;
	}

	RSA* rsa = EVP_PKEY_get1_RSA(pPubKey);
	BN_print_fp(stdout, rsa->n);
	printf("\n");
}

int get_u_bn_from_mod_PEM(char * filePath, U_BN* bignum){

    if(NULL == bignum)
        return 0;

    if(NULL == bignum->d)
        return 0;

	EVP_PKEY* pPubKey  = NULL;
    FILE*     pemFile    = NULL;

	if( !( (pemFile = fopen(filePath, "rt") ) && ( pPubKey = PEM_read_PUBKEY(pemFile,NULL,NULL,NULL) ) ) ) {
        fprintf(stderr,"Cannot read \"public key\".\n");
        return 0;
	}

	RSA* rsa = EVP_PKEY_get1_RSA(pPubKey);
	bignum2u_bn(rsa->n, bignum);
	fclose(pemFile);
	//INFO("%s\n", cu_bn_bn2hex(bignum));
	return (1);
}