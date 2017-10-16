#include "cuda_bignum.h"

#ifndef MAX
#define MAX(a,b) (a > b ? a : b)
#endif


VQ_VECTOR *cu_BN_new(void)
{
    VQ_VECTOR *ret=NULL;
    ret->top = 0;
    ret->d = NULL;
    return (ret);
}

void cu_BN_free(VQ_VECTOR *a)
{
    if (a == NULL)
        return;
    a->d = NULL;
}

int cu_BN_set_word(VQ_VECTOR *a, unsigned w)
{
    a->d[0] = w;
    a->top = (w ? 1 : 0);
    return (1);
}

unsigned cu_BN_mul_words(unsigned  *rp, const unsigned  *ap, int num, unsigned  w)
{
    unsigned  c1 = 0;
    assert(num >= 0);
    if (num <= 0)
        return (c1);
    while (num) {
        mul(rp[0], ap[0], w, c1);
        ap++;
        rp++;
        num--;
    }
    return (c1);
}

int cu_BN_mul_word(VQ_VECTOR *a, unsigned w)
{
    unsigned ll;
    w &= CU_BN_MASK2;
    if (a->top) {
        if (w == 0)
            cu_BN_zero(a);
        else {
            ll = cu_BN_mul_words(a->d, a->d, a->top, w);
            if (ll) {
                a->d[a->top++] = ll;
            }
        }
    }
    return (1);
}

int cu_BN_add_word(VQ_VECTOR *a, unsigned w)
{
    unsigned l;
    int i;

    w &= CU_BN_MASK2;

    /* degenerate case: w is zero */
    if (!w)
        return 1;
    /* degenerate case: a is zero */
    if (CU_BN_is_zero(a))
        return cu_BN_set_word(a, w);

    for (i = 0; w != 0 && i < a->top; i++) {
        a->d[i] = l = (a->d[i] + w) & CU_BN_MASK2;
        w = (w > l) ? 1 : 0;
    }
    if (w && i == a->top) {
        a->top++;
        a->d[i] = w;
    }
    return (1);
}

int cu_BN_dec2bn(VQ_VECTOR *bn, const char *a)
{
    VQ_VECTOR *ret = NULL;
    unsigned l = 0;
    int i, j;
    int num;


    if ((a == NULL) || (*a == '\0'))
        return (0);

    for (i = 0; i <= (INT_MAX/4) && isdigit((unsigned char)a[i]); i++)
        continue;

    if (i > INT_MAX/4){
        if (bn == NULL)
            cu_BN_free(ret);
        return (0);
    }
    num = i;
    if (bn == NULL)
        return (num);

    if (bn == NULL) {
        if ((ret = cu_BN_new()) == NULL)
            return (0);
    } else {
        ret = bn;
        cu_BN_zero(ret);
    }
    j = CU_BN_DEC_NUM - (i % CU_BN_DEC_NUM);
    if (j == CU_BN_DEC_NUM)
        j = 0;
    l = 0;
    while (--i >= 0) {
        l *= 10;
        l += *a - '0';
        a++;
        if (++j == CU_BN_DEC_NUM) {
            cu_BN_mul_word(ret, CU_BN_DEC_CONV);
            cu_BN_add_word(ret, l);
            l = 0;
            j = 0;
        }
    }

    cu_bn_correct_top(ret);
    bn = ret;
    return (num);
}

/* Must 'OPENSSL_free' the returned data */
char *BN_bn2hex(const VQ_VECTOR *a)
{
    int i, j, v, z = 0;
    char *buf;
    char *p;
    if (BN_is_zero(a))
        return "0";
    buf = (char *)malloc(a->top * CU_BN_BYTES * 2 + 2);
    assert(buf != NULL);
    p = buf;
    for (i = a->top - 1; i >= 0; i--) {
        for (j = CU_BN_BITS2 - 8; j >= 0; j -= 8) {
            /* strip leading zeros */
            v = ((unsigned)(a->d[i] >> j)) & 0xff;
            DEBUG_PRINT("when i is: %u and j is: %u v is equal: %u\n", i, j, v);
            if (z || (v != 0)) {
                *(p++) = Hex[v >> 4];
                DEBUG_PRINT("when i is: %u and j is: %u p high is equal: %u\n", i, j, p);
                *(p++) = Hex[v & 0x0f];
                DEBUG_PRINT("when i is: %u and j is: %u p low is equal: %u\n", i, j, p);
                z = 1;
            }
        }
    }
    *p = '\0';
    return (buf);
}

int cu_BN_ucmp(const VQ_VECTOR *a, const VQ_VECTOR *b)
{
    int i;
    unsigned t1, t2, *ap, *bp;

    i = a->top - b->top;
    if (i != 0)
        return (i);
    ap = a->d;
    bp = b->d;
    for (i = a->top - 1; i >= 0; i--) {
        t1 = ap[i];
        t2 = bp[i];
        if (t1 != t2)
            return ((t1 > t2) ? 1 : -1);
    }
    return (0);
}
