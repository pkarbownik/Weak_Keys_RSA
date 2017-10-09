#include "test.h"

void unit_test(void){
	INFO("test start...\n");
	Hw_test();
	Lw_test();
	mul_test();
	cu_BN_mul_words_test();
	cu_BN_mul_word_test();
	cu_BN_add_word_test();
	INFO("test stop\n");
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
	unsigned long t;
	t=Hw(42949672988);
	INFO("t: %u\n", t);

}

void Lw_test(void){
	unsigned t;
	t=Lw(42949672988);
	INFO("t: %u\n", t);
}

void cu_BN_mul_words_test(void){
	unsigned w = 429496725;
	unsigned i, num = 2;
	unsigned *ap = ((unsigned*)malloc(num*sizeof(unsigned)));
	ap[0]=429496725;
	ap[1]=415346545;
	unsigned *rp = ((unsigned*)malloc(num*sizeof(unsigned)));
	INFO("cu_BN_mul_words returns: %u\n", cu_BN_mul_words(rp, ap, num, w));
	INFO("num: %u\n", num);
	INFO("w: %u\n", w);
	for(i=0; i<num; i++){
		INFO("ap[%u]: %u\n", i, ap[i]);
		INFO("rp[%u]: %u\n", i, rp[i]);
	}
}

void cu_BN_mul_word_test(void){
	unsigned i;
	VQ_VECTOR   *A = NULL;
	A =   (VQ_VECTOR*)malloc(sizeof(VQ_VECTOR));
	unsigned size =2;
	unsigned w = 429496725;
	A->top=size;
	A->d = ((unsigned*)malloc(size*sizeof(unsigned)));
	A->d[0] = 429496725;
	A->d[1] = 415346545;
	INFO("cu_BN_mul_word returns: %u\n", cu_BN_mul_word(A, w));
	INFO("size: %u\n", size);
	INFO("w: %u\n", w);
	INFO("A->top: %u\n", A->top);
	for(i=0; i<A->top; i++){
		INFO("A->d[%u]: %u\n", i, A->d[i]);
	}
}


void cu_BN_add_word_test(void){
	unsigned i;
	VQ_VECTOR   *A = NULL;
	A =   (VQ_VECTOR*)malloc(sizeof(VQ_VECTOR));
	unsigned size =2;
	unsigned w = 4294967295;
	A->top=size;
	A->d = ((unsigned*)malloc(size*sizeof(unsigned)));
	A->d[0] = 4294967295;
	A->d[1] = 429496725;
	INFO("cu_BN_add_word returns: %u\n", cu_BN_add_word(A, w));
	INFO("size: %u\n", size);
	INFO("w: %u\n", w);
	INFO("A->top: %u\n", A->top);
	for(i=0; i<A->top; i++){
		INFO("A->d[%u]: %u\n", i, A->d[i]);
	}
}