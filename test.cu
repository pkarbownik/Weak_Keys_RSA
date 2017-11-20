#include "test.h"

void unit_test(void){
	INFO("tests start...\n");
	cu_BN_new_test();
	Hw_test();
	Lw_test();
	cu_BN_mul_words_test();
	cu_BN_mul_word_test();
	cu_BN_set_word_test();
	cu_BN_add_word_test();
	cu_BN_dec2bn_test();
	cu_BN_bn2hex_test();
	cu_BN_ucmp_test();
	cu_long_abs_test();
	cu_bn_usub_test();
	cu_bn_num_bits_word_test();
	cu_bn_num_bits_test();
	string_num_add_test();
	number_of_digits_test();
	long2string_test();
	string_num_add_long_test();
	cu_bn_rshift1_test();
	cu_BN_lshift_test();
	cu_euclid_test();
	bignum2u_bn_test();
	get_u_bn_from_mod_PEM_test();
	INFO("tests completed\n");
}

void cu_BN_new_test(void){
	U_BN *bn = NULL;
	unsigned i;
	for (i=0; i<100; i++){
		bn = cu_BN_new();
		cu_BN_free(bn);
	}
}

void Hw_test(void){
	assert(10 == Hw(42949672988));
	INFO("Test passed\n");
}

void Lw_test(void){
	assert(28 == Lw(42949672988));
	INFO("Test passed\n");
}

void cu_BN_mul_words_test(void){
	unsigned w = 4294967295;
	unsigned *ap = ((unsigned*)malloc(100*sizeof(unsigned)));
	unsigned *rp = ((unsigned*)malloc(101*sizeof(unsigned)));
	unsigned i;
	for(i=0; i<100; i++){
		ap[i]=4294967295;
	}
	rp[100]=cu_BN_mul_words(rp, ap, 100, w);
	assert(4294967294 == rp[100]);
	assert(1 == rp[0]);
	for(i=1; i<100; i++){
		assert(4294967295 == rp[i]);
	}
	free(ap);
	free(rp);	
	INFO("Test passed\n");
}

void cu_BN_mul_word_test(void){
	unsigned i;
	U_BN   *A = NULL;
	A = cu_BN_new();
	unsigned w = 4294967295;
	A->top=100;
	A->d = ((unsigned*)malloc(A->top * sizeof(unsigned)));
	for(i=0; i<A->top; i++){
		A->d[i] = 4294967295;
	}
	assert(1==cu_BN_mul_word(A, w));
	assert(1 == A->d[0]);
	assert(101 == A->top);
	assert(4294967294 == A->d[100]);
	cu_BN_free(A);
	INFO("Test passed\n");
}

void cu_BN_set_word_test(void){
	U_BN   *A = NULL;
	A = cu_BN_new();
	unsigned w = 0;
	A->top=2;
	A->d = ((unsigned*)malloc(A->top*sizeof(unsigned)));
	A->d[0] = 1;
	A->d[1] = 1;
	assert(1==cu_BN_set_word(A, w));
	assert(w == A->d[0]);
	assert(1 == A->top);
	cu_BN_free(A);
	INFO("Test passed\n");
}
void cu_BN_add_word_test(void){
	unsigned i;
	U_BN   *A = NULL;
	A = cu_BN_new();
	A->top=100;
	A->d = ((unsigned*)malloc(A->top*sizeof(unsigned)));
	unsigned w = 4294967295;
	for(i=0; i<A->top; i++){
		A->d[i] = 4294967295;
	}
	assert(1==cu_BN_add_word(A, w));
	assert(4294967294 == A->d[0]);
	for(i=1; i<100; i++){
		assert(0 == A->d[i]);
	}
	assert(1 == A->d[100]);
	cu_BN_free(A);
	INFO("Test passed\n");
}


void cu_BN_dec2bn_test(void){
	U_BN   *A = NULL;
	unsigned i;
	A = cu_BN_new();
	assert(1==cu_BN_dec2bn(A, "32317006071311007300714876688669951960"\
		"4441026697154840321303454275246551388678908931972014115229134"\
		"6368871796092189801949411955915049092109508815238644828312063"\
		"0877367300996091750197750389652106796057638384067568276792218"\
		"6426197561618380943384761704705816458520363050428875758915410"\
		"6580860755239912393038552191433338966834242068497478656456949"\
		"4856176035326322058077805659331026192708460314150258592864177"\
		"1167259436037184618573575983511523016459044036976132332872312"\
		"2712568471082020972515710172693132346967854258065669793504599"\
		"7268352998638215525166389437335543602135433229604645318478604"\
		"952148193555853611059596230656")); //2^2048
	for(i=0; i<64; i++){
		assert(0 ==A->d[i]);
	}
	assert(1 ==A->d[64]);
	cu_BN_free(A);
	INFO("Test passed\n");
}

void cu_BN_bn2hex_test(void){
	U_BN   *A = NULL;
	A = cu_BN_new();
	assert(1==cu_BN_dec2bn(A, "1848764767645778789"));
	assert(!strcmp("19A821C2D0C8C365", cu_bn_bn2hex(A)));
	cu_BN_free(A);
	INFO("Test passed\n");
}

void cu_BN_ucmp_test(void){
	U_BN   *A = NULL, *B = NULL;
	A = cu_BN_new();
	B = cu_BN_new();
	assert(1 == cu_BN_dec2bn(A, "1848764767645778788"));
	assert(1 == cu_BN_dec2bn(B, "1848764767645778788"));
	assert(0 == cu_BN_ucmp(A, B));
	cu_BN_free(A);
	cu_BN_free(B);
	INFO("Test passed\n");
}

void cu_long_abs_test(void){
	assert(2147483649999L == cu_long_abs(-2147483649999L));
	assert(2147483649999L == cu_long_abs(2147483649999L));
	INFO("Test passed\n");
}



void cu_bn_usub_test(void){
	U_BN   *A = NULL, *B = NULL, *C = NULL;
	A = cu_BN_new();
	B = cu_BN_new();
	C = cu_BN_new();
	cu_BN_dec2bn(A, "184876476346363755645778788");
	cu_BN_dec2bn(B, "184876476346363755644778788");
	assert(1 == cu_bn_usub(A, B, C));
	assert(!strcmp("F4240", cu_bn_bn2hex(C)));
	cu_BN_free(A);
	cu_BN_free(B);
	cu_BN_free(C);
	INFO("Test passed\n");
} 

void cu_bn_num_bits_word_test(void){
	assert(16 == cu_bn_num_bits_word(0b1111000011110000));
	INFO("Test passed\n");
} 

void cu_bn_num_bits_test(void){
	U_BN   *A = NULL;
	A = cu_BN_new();
	assert(1 == cu_BN_dec2bn(A, "1848764763497967886363755645778788"));
	assert(111 == cu_bn_num_bits(A));
	cu_BN_free(A);
	INFO("Test passed\n");
} 

void string_num_add_test(void){
	assert(!strcmp("37558289099012189180223099189280269177576", \
	string_num_add("37556443534534534534577453543634634638788", "1845564477654645645645645645634538788")));
	INFO("Test passed\n");
}

void number_of_digits_test(void){
	assert(9 == number_of_digits(123456789));
	INFO("Test passed\n");
} 

void long2string_test(void){
	assert(!strcmp("123456789", long2string(123456789)));
	INFO("Test passed\n");
}

void string_num_add_long_test(void)
{
	assert(!strcmp("37556443534534534534577455389199112293433", \
	string_num_add_long("37556443534534534534577453543634634638788", 1845564477654645)));
	INFO("Test passed\n");
}

void cu_bn_rshift1_test(void){
	U_BN   *W = NULL;
	W = cu_BN_new();
	assert(1==cu_BN_dec2bn(W, "27928727520532098560054510086934803266769027328779773633"));//2^2048
	assert(1 == cu_BN_rshift1(W));
	assert(!strcmp("91CB756A91160F4177795203E8ECFB5C5E6C6A03223360", cu_bn_bn2hex(W)));
	cu_BN_free(W);
	INFO("Test passed\n");
}

void cu_BN_lshift_test(void){
	U_BN   *A = NULL;
	A = cu_BN_new();
	cu_BN_dec2bn(A, "231622341");
	assert(1 == cu_BN_lshift(A, 7));
	assert(!strcmp("6E7236280", cu_bn_bn2hex(A)));
	cu_BN_free(A);
	INFO("Test passed\n");
}

void cu_euclid_test(void){
	U_BN   *A = NULL, *B = NULL;
	unsigned L=5, N=1; 
	/*A = cu_BN_new();
	B = cu_BN_new();
	cu_BN_dec2bn(A, "8353015802438879251643065122143616");
	cu_BN_dec2bn(B, "858238501677248042531768818944");
	A = cu_euclid(A, B);
	INFO("6E7236280=%s\n", cu_bn_bn2hex(A));
	cu_BN_free(A);*/


    A =   (U_BN*)malloc(N*sizeof(U_BN));
    B =   (U_BN*)malloc(N*sizeof(U_BN));

    for(int i=0; i<N; i++){
        U_BN a;
        U_BN b;
        a.d = (unsigned*)malloc(L*sizeof(unsigned));
        b.d = (unsigned*)malloc(L*sizeof(unsigned));
        a.top =   L;
        b.top =   L;

        for(int j=0; j<L; j++)
            a.d[j]=0;

        for(int j=0; j<L; j++)
            b.d[j]=0;

        A[i] = a;
        B[i] = b;
    }

    cu_BN_dec2bn(&A[0], "858238501677248042531768818944");
    cu_BN_dec2bn(&B[0], "8353015802438879251643065122143616");

    A = cu_euclid(&A[0], &B[0]);
	INFO("6E7236280=%s\n", cu_bn_bn2hex(A));

    INFO("Test passed\n");
}

void bignum2u_bn_test(void){
	BIGNUM *bn;
	U_BN *u_bn;
	bn = BN_new();
	u_bn = cu_BN_new();
	assert( 0 < BN_dec2bn(&bn, "54635484657846634") );
	assert( 1 == bignum2u_bn(bn, u_bn) );
	assert(!strcmp("C21AAF0F29896A", cu_bn_bn2hex(u_bn)));
	cu_BN_free(u_bn);
	BN_free(bn);
	INFO("Test passed\n");

}

void get_u_bn_from_mod_PEM_test(void){
	U_BN *u_bn;
	u_bn = cu_BN_new();
	//assert( 1 == get_u_bn_from_mod_PEM("keys_and_messages/1.pem", u_bn));
	//printf("from PEM file: %s\n", cu_bn_bn2hex(u_bn));
	INFO("Test passed\n");
}