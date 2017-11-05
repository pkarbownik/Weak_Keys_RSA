#include "test.h"

void unit_test(void){
	INFO("test start...\n");
	cu_BN_new_test();
	Hw_test();
	Lw_test();
	mul_test();
	cu_BN_mul_words_test();
	cu_BN_mul_word_test();
	cu_BN_set_word_test();
	//cu_BN_add_word_test();
	//cu_BN_dec2bn_test();
	//BN_bn2hex_test();
	//cu_BN_ucmp_test();
	//cu_long_abs_test();
	//cu_bn_usub_test();
	//cu_bn_num_bits_word_test();
	//cu_bn_num_bits_test();
	//string_num_add_test();
	//number_of_digits_test();
	//long2string_test();
	//string_num_add_long_test();
	//cu_BN_rshift1_test();
	//cu_BN_lshift_test();
	//cu_euclid_test();
	INFO("test stop\n");
}

void cu_BN_new_test(void){
	VQ_VECTOR *bn = NULL;
	bn = cu_BN_new();
	bn -> top = 1;
	INFO("bn -> top: %u\n", bn -> top);
	cu_BN_free(bn);
}

void mul_test(void){
	unsigned r,a,w,c;
	r=0;
	a=4;
	w=CU_BN_DEC_CONV;
	c=294967299;
	mul(r,a,w,c);
	INFO("a: %u\n", a);
	INFO("r: %u\n", r);
	INFO("w: %u\n", w);
	INFO("c: %u\n", c);
}

void Hw_test(void){
	unsigned t;
	t=Hw(42949672988);
	INFO("t: %u\n", t);

}

void Lw_test(void){
	unsigned t;
	t=Lw(42949672988);
	INFO("t: %u\n", t);
}

void cu_BN_mul_words_test(void){
	unsigned w = 429496725, i;
	unsigned *ap = ((unsigned*)malloc(2*sizeof(unsigned)));
	unsigned *rp = ((unsigned*)malloc(3*sizeof(unsigned)));
	ap[0]=429496725;
	ap[1]=429496725;
	INFO("cu_BN_mul_words returns: %u\n", cu_BN_mul_words(rp, ap, 3, w));
	INFO("w: %u\n", w);
	for(i=0; i<2; i++){
		INFO("ap[%u]: %u\n", i, ap[i]);
	}
	for(i=0; i<3; i++){
		INFO("rp[%u]: %u\n", i, rp[i]);
	}
}

void cu_BN_mul_word_test(void){
	unsigned i;
	VQ_VECTOR   *A = NULL;
	A = cu_BN_new();
	unsigned size =2;
	unsigned w = 429496725;
	A->top=size;
	A->d = ((unsigned*)malloc(size*sizeof(unsigned)));
	A->d[0] = 429496725;
	A->d[1] = 0;
	assert(1==cu_BN_mul_word(A, w));
	assert(Lw(429496725L*429496725) == A->d[0]);
	assert(Hw(429496725L*429496725) == A->d[1]);
	INFO("Test passed\n");
	cu_BN_free(A);
}

void cu_BN_set_word_test(void){
	
	unsigned i;
	VQ_VECTOR   *A = NULL;
	A = cu_BN_new();
	unsigned size =2;
	unsigned w = 429496725;
	A->top=size;
	A->d = ((unsigned*)malloc(size*sizeof(unsigned)));
	A->d[0] = 1;
	A->d[1] = 1;
	assert(1==cu_BN_set_word(A, w));
	assert(w == A->d[0]);
	assert(1 == A->d[1]);
	assert(1 == A->top);
	INFO("Test passed\n");
	cu_BN_free(A);
}
void cu_BN_add_word_test(void){
	unsigned i;
	VQ_VECTOR   *A = NULL;
	A = cu_BN_new();
	unsigned size =2;
	unsigned w = 4294967295;
	A->top=size;
	A->d = ((unsigned*)malloc(size*sizeof(unsigned)));
	A->d[0] = 4294967295;
	A->d[1] = 4;
	INFO("cu_BN_add_word returns: %u\n", cu_BN_add_word(A, w));
	cu_BN_free(A);
}


void cu_BN_dec2bn_test(void){
	unsigned i;
	VQ_VECTOR   *A = NULL;
	char *a = "34023467455236674558236"; //2^128
	A = cu_BN_dec2bn(a);
	INFO("a: %s\n", a);
	for(i=0; i<A->top; i++){
		INFO("A->d[%u]: %u\n", i, A->d[i]);
	}
}

void BN_bn2hex_test(void)
{
	unsigned i;
	VQ_VECTOR   *A = NULL;
	char *a = "1848764767645778789"; //2^32
	char *b;
	A = cu_BN_dec2bn(a);
	INFO("a: %s\n", a);
	INFO("A->top: %u\n", A->top);
	for(i=0; i<A->top; i++){
		INFO("A->d[%u]: %u\n", i, A->d[i]);
	}
	INFO("b: %s\n", cu_bn_bn2hex(A));
}

void cu_BN_ucmp_test(void)
{
	unsigned i;
	VQ_VECTOR   *A = NULL, *B = NULL;
	char *a = "1848764767645778788"; //2^32
	char *b = "1848764767645778788";
	A = cu_BN_dec2bn(a);
	B = cu_BN_dec2bn(b);
	INFO("a: %s\n", a);
	INFO("b: %s\n", b);
	INFO("A->top: %u\n", A->top);
	INFO("B->top: %u\n", B->top);
	INFO("cu_BN_ucmp: %d\n", cu_BN_ucmp(A, B));
}

void cu_long_abs_test(void)
{
	long negative;
	unsigned long positive;
	negative = -2147483649999;
	positive = cu_long_abs(negative);
	INFO("negative: %li\n", negative);
	positive = cu_long_abs(positive);
	INFO("positive %li\n", positive);
}



void cu_bn_usub_test(void)
{
	unsigned size =2;
	VQ_VECTOR   *A = NULL, *B = NULL, *C = NULL;
	char *a = "184876476346363755645778788"; //2^32
	char *b = "184876476346363755644778788";
	C =   (VQ_VECTOR*)malloc(sizeof(VQ_VECTOR));
	C->top=size;
	C->d = ((unsigned*)malloc(size*sizeof(unsigned)));
	A = cu_BN_dec2bn(a);
	B = cu_BN_dec2bn(b);
	C = cu_bn_usub(A, B);
	INFO("size: %u\n", size);
	INFO("a: %s\n", a);
	INFO("b: %s\n", b);
	INFO("A->top: %u\n", A->top);
	INFO("B->top: %u\n", B->top);
	INFO("c in hex: %s\n", cu_bn_bn2hex(C));
} 

void cu_bn_num_bits_word_test(void)
{
	long number =23457635347456;
	INFO("number: %lu\n", number);
	INFO("cu_bn_num_bits_word of word : %d\n", cu_bn_num_bits_word(number));
} 

void cu_bn_num_bits_test(void)
{
	unsigned size =2;
	VQ_VECTOR   *A = NULL;
	char *a = "1848764763497967886363755645778788"; //2^32
	A = cu_BN_dec2bn(a);
	INFO("size: %u\n", size);
	INFO("a: %s\n", a);
	INFO("A->top: %u\n", A->top);
	INFO("cu_bn_num_bits: %d\n", cu_bn_num_bits(A));
} 

void string_num_add_test(void)
{
	char *a = "37556443534534534534577453543634634638788"; //2^32
	char *b = "1845564477654645645645645645634538788";
	INFO("a: %s\n", a);
	INFO("b: %s\n", b);
	INFO("a+b: %s\n", string_num_add(a, b));
} 

void number_of_digits_test(void)
{
	long number = 123456789;
	INFO("number: %lu\n", number);
	INFO("number of digits: %u\n", number_of_digits(number));
} 

void long2string_test(void){
	long number = 123456789;
	INFO("number: %lu\n", number);
	INFO("long2string: %s\n", long2string(number));
}

void string_num_add_long_test(void)
{
	char *a = "37556443534534534534577453543634634638788"; //2^32
	long b = 1845564477654645;
	INFO("a (string): %s\n", a);
	INFO("b (long): %lu\n", b);
	INFO("a+b: %s\n", string_num_add_long(a, b));
}

void cu_BN_rshift1_test(void){
	unsigned i;
	VQ_VECTOR   *A = NULL;
	char *a = "664683265346895634856568568568678867688687565685686785858585857857857857576676345865647588346856346";
	A = cu_BN_dec2bn(a);
	A = cu_BN_rshift1(A);
	INFO("a: %s\n", a);
	//for(i=0; i<A->top; i++){
	//	INFO("A->d[%u]: %u\n", i, A->d[i]);
	//}
	//INFO("result a=%s rshift1: %s\n", a, cu_bn_bn2hex(A));
}

void cu_BN_lshift_test(void){
	unsigned i;
	VQ_VECTOR   *A = NULL;
	char *a = "6646832653468945893465834685634657834685683467856347563465862347623974563475723562340956348568346856346";
	A = cu_BN_dec2bn(a);
	A = cu_BN_lshift(A, 4);
	INFO("a: %s\n", a);
	//for(i=0; i<A->top; i++){
	//	INFO("A->d[%u]: %u\n", i, A->d[i]);
	//}
	//INFO("result a=%s lshift (A,4): %s\n", a, cu_bn_bn2hex(A));
}

void cu_euclid_test(void){
	VQ_VECTOR   *A = NULL, *B = NULL;
	char *a = "211319228244187486513184688950596901432020552950859782";
	char *b = "37724566494969212902300091545866760828606124226";
	A = cu_BN_dec2bn(a);
	INFO("a: %s\n", a);
	B = cu_BN_dec2bn(b);
	INFO("b: %s\n", b);
	A = cu_euclid(A, B);
	INFO("15417B103640DD7A2752917A=%s\n", cu_bn_bn2hex(A));
}