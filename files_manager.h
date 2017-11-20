#ifndef FILES_MANAGER_H
#define FILES_MANAGER_H

#include <openssl/rsa.h>
#include <openssl/pem.h>
#include <openssl/bn.h>
#include "cuda_bignum.h"

void print_mod_from_pem_file(char * filePath);
int get_u_bn_from_mod_PEM(char * filePath, U_BN* bignum);

#endif /* CUDA_BIGNUM_H */