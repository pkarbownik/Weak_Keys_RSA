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
void cu_BN_new_test(void);
void mul_test(void);
void Hw_test();
void Lw_test();
void cu_BN_mul_words_test(void);
void cu_BN_mul_word_test(void);
void cu_BN_add_word_test(void);
void cu_BN_dec2bn_test(void);
void BN_bn2hex_test(void);
void cu_BN_ucmp_test(void);
void cu_long_abs_test(void);
void cu_bn_usub_test(void);
void cu_bn_bn2dec_test(void);
void cu_bn_num_bits_word_test(void);
void cu_bn_num_bits_test(void);
void string_num_add_test(void);
void number_of_digits_test(void);
void long2string_test(void);
void string_num_add_long_test(void);
void cu_bn_copy_test(void);
#endif /* TEST_H */