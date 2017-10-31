#include "cuda_bignum.h"

#ifndef MAX
#define MAX(a,b) (a > b ? a : b)
#endif

char *strrev(char *str)
{
      char *p1, *p2;

      if (! str || ! *str)
            return str;
      for (p1 = str, p2 = str + strlen(str) - 1; p2 > p1; ++p1, --p2)
      {
            *p1 ^= *p2;
            *p2 ^= *p1;
            *p1 ^= *p2;
      }
      return str;
}


VQ_VECTOR *cu_BN_new(void)
{
    VQ_VECTOR *ret = NULL;
    ret =  (VQ_VECTOR*)malloc(sizeof(VQ_VECTOR));
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
    DEBUG_PRINT("num: %u\n", num);
    while (num) {
        mul(rp[0], ap[0], w, c1);
        DEBUG_PRINT("rp[0]: %u\n", rp[0]);
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
    DEBUG_PRINT("w: %u\n", w);
    DEBUG_PRINT("a->top: %u\n", a->top);
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

VQ_VECTOR *cu_BN_dec2bn(const char *a)
{
    VQ_VECTOR *ret = NULL;
    unsigned l = 0;
    int j, i;
    ret = cu_BN_new();
    ret->top=0;
    ret->d = ((unsigned*)malloc(sizeof(unsigned)));
    if ((a == NULL) || (*a == '\0'))
        return (0);

    for (i = 0; i <= (INT_MAX/4) && isdigit((unsigned char)a[i]); i++)
        continue;

    if (i > INT_MAX/4)
        return 0;

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
            DEBUG_PRINT("ret->d[0]: %u\n", ret->d[0]);
            cu_BN_add_word(ret, l);
            DEBUG_PRINT("ret->d[0]: %u\n", ret->d[0]);
            l = 0;
            j = 0;
        }
    }
    cu_bn_correct_top(ret);
    DEBUG_PRINT("ret->top: %u\n", ret->top);
    return (ret);
}

/* Must 'OPENSSL_free' the returned data */
char *cu_bn_bn2hex(const VQ_VECTOR *a)
{
    int i, j, v, z = 0;
    char *buf;
    char *p;
    if (BN_is_zero(a))
        return 0;
    buf = (char *)malloc(a->top * CU_BN_BYTES * 2 + 2);
    assert(buf != NULL);
    p = buf;
    for (i = a->top - 1; i >= 0; i--) {
        for (j = CU_BN_BITS2 - 8; j >= 0; j -= 8) {
            /* strip leading zeros */
            v = ((unsigned)(a->d[i] >> j)) & 0xff;
            if (z || (v != 0)) {
                *(p++) = Hex[v >> 4];
                *(p++) = Hex[v & 0x0f];
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

long cu_long_abs(long number){
    if(number<0)
        return -number;
    else
        return number;
}

/* unsigned subtraction of b from a, a must be larger than b. */
int cu_bn_usub(VQ_VECTOR *r, const VQ_VECTOR *a, const VQ_VECTOR *b)
{
    unsigned max, min, dif;
    register unsigned t1, t2, *ap, *bp, *rp;
    int i, carry;

    max = a->top;
    min = b->top;
    dif = cu_long_abs(max - min);

    ap = a->d;
    bp = b->d;
    rp = r->d;

#if 1
    carry = 0;
    for (i = min; i != 0; i--) {
        t1 = *(ap++);
        t2 = *(bp++);
        if (carry) {
            carry = (t1 <= t2);
            t1 = (t1 - t2 - 1) & BN_MASK2;
        } else {
            carry = (t1 < t2);
            t1 = (t1 - t2) & BN_MASK2;
        }
# if defined(IRIX_CC_BUG) && !defined(LINT)
        dummy = t1;
# endif
        *(rp++) = t1 & BN_MASK2;
    }
#else
    carry = bn_sub_words(rp, ap, bp, min);
    ap += min;
    bp += min;
    rp += min;
#endif
    if (carry) {                /* subtracted */
        if (!dif)
            /* error: a < b */
            return 0;
        while (dif) {
            dif--;
            t1 = *(ap++);
            t2 = (t1 - 1) & BN_MASK2;
            *(rp++) = t2;
            if (t1)
                break;
        }
    }
#if 0
    memcpy(rp, ap, sizeof(*rp) * (max - i));
#else
    if (rp != ap) {
        for (;;) {
            if (!dif--)
                break;
            rp[0] = ap[0];
            if (!dif--)
                break;
            rp[1] = ap[1];
            if (!dif--)
                break;
            rp[2] = ap[2];
            if (!dif--)
                break;
            rp[3] = ap[3];
            rp += 4;
            ap += 4;
        }
    }
#endif

    r->top = max;
    cu_bn_correct_top(r);
    return (1);
}

int cu_bn_num_bits_word(long l)
{
     static const unsigned char bits[256] = {
        0, 1, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4,
        5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
        6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
        6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
        7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
        8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
        8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
        8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
        8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
        8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
        8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
        8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    };

#if defined(SIXTY_FOUR_BIT_LONG)
    if (l & 0xffffffff00000000L) {
        if (l & 0xffff000000000000L) {
            if (l & 0xff00000000000000L) {
                return (bits[(int)(l >> 56)] + 56);
            } else
                return (bits[(int)(l >> 48)] + 48);
        } else {
            if (l & 0x0000ff0000000000L) {
                return (bits[(int)(l >> 40)] + 40);
            } else
                return (bits[(int)(l >> 32)] + 32);
        }
    } else
#else
# ifdef SIXTY_FOUR_BIT
    if (l & 0xffffffff00000000LL) {
        if (l & 0xffff000000000000LL) {
            if (l & 0xff00000000000000LL) {
                return (bits[(int)(l >> 56)] + 56);
            } else
                return (bits[(int)(l >> 48)] + 48);
        } else {
            if (l & 0x0000ff0000000000LL) {
                return (bits[(int)(l >> 40)] + 40);
            } else
                return (bits[(int)(l >> 32)] + 32);
        }
    } else
# endif
#endif
    {
#if defined(THIRTY_TWO_BIT) || defined(SIXTY_FOUR_BIT) || defined(SIXTY_FOUR_BIT_LONG)
        if (l & 0xffff0000L) {
            if (l & 0xff000000L)
                return (bits[(int)(l >> 24L)] + 24);
            else
                return (bits[(int)(l >> 16L)] + 16);
        } else
#endif
        {
#if defined(THIRTY_TWO_BIT) || defined(SIXTY_FOUR_BIT) || defined(SIXTY_FOUR_BIT_LONG)
            if (l & 0xff00L)
                return (bits[(int)(l >> 8)] + 8);
            else
#endif
                return (bits[(int)(l)]);
        }
    }
}

int cu_bn_num_bits(const VQ_VECTOR *a)
{
    int i = a->top - 1;

    if (CU_BN_is_zero(a))
        return 0;
    return ((i * CU_BN_BITS2) + cu_bn_num_bits_word(a->d[i]));
}

#define c2d(c) (c-'0')
#define d2c(c) (c+'0')

char *string_num_add(const char *a, const char *b){
    int alen, blen;
    int i, carry=0;
    char *wk = NULL;
    char *awk=strdup(a);
    char *bwk=strdup(b);

    alen=strlen(strrev(awk));
    blen=strlen(strrev(bwk));
    if(alen<blen){
        alen ^= blen;blen ^= alen;alen ^= blen;//swap
        wk = awk ; awk = bwk ; bwk = wk;
    }
    char *ans = ((char *) malloc (alen * sizeof(char)));
    ans[alen+1]=ans[alen]='\0';
    for(i=0;i<alen;++i){
        int sum = c2d(awk[i])+(i<blen ? c2d(bwk[i]): 0)+carry;
        ans[i] = d2c(sum % 10);
        carry = sum / 10;
    }
    if(carry){
        ans[i++]='1';
    }
    free(awk);
    free(bwk);
    return strrev(ans);
}

unsigned number_of_digits(long number){
    unsigned count = 0;
    while(number != 0){
        number /= 10;
        ++count;
    }
    return count;
}

char *long2string(long number){
    char * result = ((char *) malloc (number_of_digits(number)*sizeof(char)));
    sprintf(result, "%lu", number);
    return result;
}

char *string_num_add_long(const char *a, long word){
    return string_num_add(a, long2string(word));
}

int cu_bn_copy(VQ_VECTOR *a, const  VQ_VECTOR *b){
    a = (VQ_VECTOR*)malloc(sizeof(VQ_VECTOR));
    a->top = b->top;
    int i;
    a->d = (unsigned *)malloc(sizeof(unsigned)*a->top);
    for(i=0; i<b->top; i++){
        a->d[i]=b->d[i];
    }
    return 1;
}