#ifndef GCD_H
#define GCD_H 

#define _GNU_SOURCE 
#include <stdio.h>
#include <openssl/rsa.h>
#include <openssl/pem.h>
#include <openssl/bn.h>
#include <string.h>
 void euclid_modulo(BIGNUM *r, const BIGNUM *a, const BIGNUM *b);

#endif /* GCD_H */