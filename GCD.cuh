#ifndef GCD_H
#define GCD_H 

#include <stdio.h>
#include <openssl/rsa.h>
#include <openssl/pem.h>
#include <openssl/bn.h>
#include <string.h>
#include <time.h>


void euclid_modulo(BIGNUM *r, const BIGNUM *a, const BIGNUM *b);


#endif /* GCD_H */
