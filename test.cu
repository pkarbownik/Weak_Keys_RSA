#include "test.h"

void unit_test(void){
	INFO("test start...\n");
	Hw_test();
	Lw_test();
	mul_test();
	cu_bn_mul_words_test();
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

void cu_bn_mul_words_test(void){
	unsigned w = 429496725;
	unsigned i, num = 2;
	unsigned *ap = ((unsigned*)malloc(num*sizeof(unsigned)));
	ap[0]=429496725;
	ap[1]=415346545;
	unsigned *rp = ((unsigned*)malloc(num*sizeof(unsigned)));
	INFO("cu_bn_mul_words returns: %u\n", cu_bn_mul_words(rp, ap, num, w));
	INFO("num: %u\n", num);
	INFO("w: %u\n", w);
	for(i=0; i<num; i++){
		INFO("ap[%u]: %u\n", i, ap[i]);
		INFO("rp[%u]: %u\n", i, rp[i]);
	}
}