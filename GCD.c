#include "GCD.cuh"

void euclid_modulo(BIGNUM *r, const BIGNUM *a, const BIGNUM *b){
	BN_CTX *ctx;
	ctx = BN_CTX_new();
	if(!BN_gcd(r, a, b, ctx)){
		printf("BN_gcd_error");
	}
	BN_CTX_free(ctx);
}

