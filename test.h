#ifndef TEST_H
#define TEST_H

#include "cuda_bignum.h"
#include <assert.h>

#define TRACE 1

#if defined(TRACE) && TRACE > 0
 #define INFO(fmt, args...) fprintf(stderr, "INFO: %s:%d:%s(): " fmt, \
    __FILE__, __LINE__, __func__, ##args)
#else
 #define INFO(fmt, args...) /* Don't do anything in release builds */
#endif

void unit_test(void);
void mul_test(void);
void Hw_test();
void Lw_test();
void cu_bn_mul_words_test(void);

#endif /* TEST_H */