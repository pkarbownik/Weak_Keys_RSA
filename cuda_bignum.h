#ifndef CUDA_BIGNUM_H
#define CUDA_BIGNUM_H

#include <stdio.h>
#include <openssl/rsa.h>
#include <openssl/pem.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <ctype.h>
#include "cuda_runtime.h"

#define DEBUG 0

#if defined(DEBUG) && DEBUG > 0
 #define DEBUG_PRINT(fmt, args...) fprintf(stderr, "DEBUG: %s:%d:%s(): " fmt, \
    __FILE__, __LINE__, __func__, ##args)
#else
 #define DEBUG_PRINT(fmt, args...) /* Don't do anything in release builds */
#endif

struct   __Q_VECTOR__{
    unsigned* d;       
    int     top;    
};

typedef struct __Q_VECTOR__     VQ_VECTOR;

#define debug(fmt, ...) printf("%s:%d: " fmt, __FILE__, __LINE__, __VA_ARGS__);


# define cu_bn_correct_top(a) \
        { \
        unsigned *ftl; \
        int tmp_top = (a)->top; \
        if (tmp_top > 0) \
                { \
                for (ftl= &((a)->d[tmp_top-1]); tmp_top > 0; tmp_top--) \
                        if (*(ftl--)) break; \
                (a)->top = tmp_top; \
                } \
        }

/************************64bit version*********************/
//# define CU_BN_MASK2        (0xffffffffffffffffL)
//# define CU_BN_DEC_NUM      19
//# define CU_BN_DEC_CONV     (10000000000000000000UL)
//# define CU_BN_MASK2l       (0xffffffffL)

/************************32bit version*********************/
#define cu_BN_zero(a)      (cu_BN_set_word((a),0))
#define cu_BN_is_odd(a)        (((a)->top > 0) && ((a)->d[0] & 1))
#define CU_BN_is_zero(a)       ((a)->top == 0)
#define cu_BN_is_initialized() 
#define CU_BN_BITS2        32
#define CU_BN_BITS4        16
#define CU_BN_BYTES        8
#define CU_BN_MASK2        (0xffffffffL)
#define CU_BN_MASK2l       (0xffff)
#define CU_BN_MASK2h1      (0xffff8000L)
#define CU_BN_MASK2h       (0xffff0000L)
#define CU_BN_DEC_CONV     (1000000000L)
#define CU_BN_DEC_NUM      9
#define CU_BN_DEC_FMT1     "%u"
#define CU_BN_DEC_FMT2		"%09u"
#define CU_BN_TBIT         (0x80000000L)

#define Lw(t)    (((unsigned)t))
#define Hw(t)    ((unsigned)((t)>>CU_BN_BITS2))

static const char Hex[] = "0123456789ABCDEF";

//extern __global__ void testKernel(VQ_VECTOR *X, int N);

//unsigned bn_mul_add_words(unsigned *rp, const unsigned *ap, int num, unsigned w);
char *strrev(char *str);
VQ_VECTOR *cu_BN_new();
void cu_BN_free(VQ_VECTOR *a);
int cu_BN_set_word(VQ_VECTOR *a, unsigned w);
int cu_BN_mul_word(VQ_VECTOR *a, unsigned w);
int cu_BN_add_word(VQ_VECTOR *a, unsigned w);
int cu_BN_dec2bn(VQ_VECTOR * ret, const char *a);
unsigned  cu_BN_mul_words(unsigned  *rp, const unsigned  *ap, int num, unsigned  w);
char *cu_bn_bn2hex(const VQ_VECTOR *a);
//__device__ int cu_BN_ucmp(const VQ_VECTOR *a, const VQ_VECTOR *b);
//__device__ long cu_long_abs(long number);
//__device__ int cu_bn_usub(const VQ_VECTOR *a, const VQ_VECTOR *b, VQ_VECTOR *c);
int cu_bn_num_bits_word(long l);
int cu_bn_num_bits(const VQ_VECTOR *a);
unsigned number_of_digits(long number);
char *long2string(long number);
char *string_num_add(const char *a, const char *b);
char *string_num_add_long(const char *a, long b);
int cu_bn_copy(VQ_VECTOR *a, const  VQ_VECTOR *b);
//__device__ int cu_BN_rshift1(VQ_VECTOR *a);
//__device__ int cu_BN_lshift(VQ_VECTOR *a, unsigned n);
//__device__ VQ_VECTOR *cu_euclid(VQ_VECTOR *a, VQ_VECTOR *b);
#endif /* CUDA_BIGNUM_H */