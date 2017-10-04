#include <stdio.h>
#include <openssl/rsa.h>
#include <openssl/pem.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <ctype.h>
#include "cuda_runtime.h"

struct   __Q_VECTOR__{
    unsigned long* d;       
    int     top;    
};

typedef struct __Q_VECTOR__     VQ_VECTOR;



# define cu_bn_correct_top(a) \
        { \
        unsigned long *ftl; \
        int tmp_top = (a)->top; \
        if (tmp_top > 0) \
                { \
                for (ftl= &((a)->d[tmp_top-1]); tmp_top > 0; tmp_top--) \
                        if (*(ftl--)) break; \
                (a)->top = tmp_top; \
                } \
        }


# define CU_BN_MASK2        (0xffffffffffffffffL)
# define CU_BN_is_zero(a)       ((a)->top == 0)
# define CU_BN_DEC_NUM      19
# define CU_BN_DEC_CONV     (10000000000000000000UL)
# define CU_BN_BITS2        64
# define CU_BN_BITS4        32
# define CU_BN_MASK2l       (0xffffffffL)
# define cu_BN_zero(a)      (cu_BN_set_word((a),0))


# define BN_UMULT_HIGH(a,b)   ({      \
        register unsigned long  ret,discard;  \
        asm ("mulq      %3"             \
             : "=a"(discard),"=d"(ret)  \
             : "a"(a), "g"(b)           \
             : "cc");                   \
        ret;                  })            


# define mul(r,a,w,c)    {               \
        unsigned long high,low,ret,ta=(a);   \
        low =  (w) * ta;                \
        high=  BN_UMULT_HIGH(w,ta);     \
        ret =  low + (c);               \
        (c) =  high;                    \
        (c) += (ret<low)?1:0;           \
        (r) =  ret;                     \
        }

//extern __global__ void testKernel(VQ_VECTOR *X, int N);

//unsigned long bn_mul_add_words(unsigned long *rp, const unsigned long *ap, int num, unsigned long w);
VQ_VECTOR *cu_BN_new(void);
void cu_BN_free(VQ_VECTOR *a);
int cu_BN_set_word(VQ_VECTOR *a, unsigned long w);
int cu_BN_mul_word(VQ_VECTOR *a, unsigned long w);
int cu_BN_add_word(VQ_VECTOR *a, unsigned long w);
int cu_BN_dec2bn(VQ_VECTOR *bn, const char *a);
