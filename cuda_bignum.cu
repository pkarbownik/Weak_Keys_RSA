#include "cuda_bignum.h"

#ifndef MAX
#define MAX(a,b) (a > b ? a : b)
#endif

char *strrev(char *str){
      char *p1, *p2;

      if (! str || ! *str)
            return str;

      for (p1 = str, p2 = str + strlen(str) - 1; p2 > p1; ++p1, --p2){
            *p1 ^= *p2;
            *p2 ^= *p1;
            *p1 ^= *p2;
      }
      return str;
}

VQ_VECTOR *cu_BN_new(void)
{
    VQ_VECTOR *ret;

    if ((ret = (VQ_VECTOR *)malloc(sizeof(*ret))) == NULL) {
        return (NULL);
    }
    ret->top = 0;
    ret->d = (unsigned *) malloc(sizeof(unsigned));
    return (ret);
}

void cu_BN_free(VQ_VECTOR *a)
{
    if (a == NULL)
        return;
    if (a->d != NULL)
        free(a->d);
    if (a == NULL)
        free(a);
}

unsigned cu_BN_mul_words(unsigned  *rp, const unsigned  *ap, int num, unsigned  w){

    unsigned i=0;
    unsigned tmp;
    unsigned  c1 = 0;
    if (num <= 0)
        return 0;
    while (num) {
        tmp= Lw(((unsigned long)(ap[i]) * (w) + (c1)));
        c1= Hw(((unsigned long)(ap[i]) * (w) + (c1)));
        DEBUG_PRINT("ap[%u]: %u\n", i, ap[i]);
        DEBUG_PRINT("rp[%u]: %u\n", i, tmp);
        rp[i] = tmp;
        i++;
        num--;
    }
    DEBUG_PRINT("w: %u\n", w);
    return (c1);

}

int cu_BN_mul_word(VQ_VECTOR *a, unsigned w){

    unsigned ll;
    if(NULL == a)
        return 0;

    if(NULL == a->d)
        return 0;

    w &= CU_BN_MASK2;

    if (a->top) {
        if (w == 0)
            cu_BN_zero(a);
        else {
            ll = cu_BN_mul_words(a->d, a->d, a->top, w);
            if (ll) {
                a->d = ((unsigned*)realloc(a->d, (++a->top)*sizeof(unsigned))); 
                a->d[(a->top-1)] = ll;
                DEBUG_PRINT("a->d[(a->top-1)]: %u\n", a->d[(a->top-1)]);
                DEBUG_PRINT("a->top: %u\n", a->top);
            }
        }
    }
    return (1);

}

int cu_BN_set_word(VQ_VECTOR *a, unsigned w){

    if(NULL == a)
        return 0;

    if(NULL == a->d)
        return 0;

    a->d = ((unsigned *)malloc(sizeof(unsigned)));
    a->d[0] = w;
    a->top = 1;
    return (1);

}

int cu_BN_add_word(VQ_VECTOR *a, unsigned w){

    unsigned l;
    int i;

    if(NULL == a)
        return 0;

    if(NULL == a->d)
        return 0;

    w &= CU_BN_MASK2;

    /* degenerate case: w is zero */
    if (!w)
        return (1);

    /* degenerate case: a is zero */
    if (CU_BN_is_zero(a))
        return cu_BN_set_word(a, w);

    for (i = 0; i < a->top; i++) {
        a->d[i] = l = (a->d[i] + w) & CU_BN_MASK2;
        w = (w > l) ? 1 : 0;
    }

    if ((w) && (i == a->top)) {
        a->d = ((unsigned*)realloc(a->d, (++a->top)*sizeof(unsigned)));
        a->d[i] = w;
    }
    return (1);

}

int cu_BN_dec2bn(VQ_VECTOR *ret, const char *a){

    unsigned l;
    int j, i;

    if(NULL == ret)
        return 0;

    if(NULL == ret->d)
        return 0;

    if ((a == NULL) || (*a == '\0')){
        return (0);
    }

    if ('0'==*a)
        return cu_BN_set_word(ret, 0);

    for (i = 0; i <= (INT_MAX/4) && isdigit((unsigned char)a[i]); i++)
        continue;


    DEBUG_PRINT("number of digits: %u\n", i);

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
            DEBUG_PRINT("l: %u\n", l);
            cu_BN_mul_word(ret, CU_BN_DEC_CONV);
            cu_BN_add_word(ret, l);
            l = 0;
            j = 0;
        }
    }
    cu_bn_correct_top(ret);
    return (1);

}

/* Must 'OPENSSL_free' the returned data */
char *cu_bn_bn2hex(const VQ_VECTOR *a){

    int i, j, v, z = 0;
    char *buf = NULL, *p = NULL;

    if(NULL == a)
        return 0;

    if(NULL == a->d)
        return 0;

    if (BN_is_zero(a))
        return 0;

    buf = (char *)malloc(sizeof(char) * (CU_BN_BYTES * a->top + 1));
    p=buf;
    for (i = a->top - 1; i >= 0; i--) {
        for (j = CU_BN_BITS2 - 4; j >= 0; j -= 4) {
            v = ((unsigned)(a->d[i] >> j)) & 0x0f;
            if (z || (v != 0)) {
                //DEBUG_PRINT("v: %d\n", v);
                *(buf++) = Hex[v & 0x0f];
                //DEBUG_PRINT("*buf: %c\n", *(buf-1));
                z = 1;
            }
        }
    }
    *buf = '\0';
    //DEBUG_PRINT("p: %s\n", p);
    return (p);

}

int cu_BN_ucmp(const VQ_VECTOR *a, const VQ_VECTOR *b){

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
int cu_bn_usub(const VQ_VECTOR *a, const VQ_VECTOR *b, VQ_VECTOR *r){

    unsigned max, min, dif;
    register unsigned t1, t2, *ap, *bp, *rp;
    int i, carry;

    if(NULL == a || NULL == b || NULL == r)
        return 0;

    if(NULL == a->d || NULL == b->d || NULL == r->d)
        return 0;

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

int cu_bn_num_bits_word(long l){

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

int cu_bn_num_bits(const VQ_VECTOR *a){

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


int cu_BN_rshift1(VQ_VECTOR *a){

    if(NULL == a)
        return 0;

    if(NULL == a->d)
        return 0;

    if (BN_is_zero(a))
        return 0;

    unsigned *ap, *rp , t, c;
    int i, j;

    i = a->top;
    DEBUG_PRINT("a->top: %u\n", a->top);
    ap = a->d;
    DEBUG_PRINT("a->d[0]: %u\n", a->d[0]);

    j = i - (ap[i - 1] == 1);
    DEBUG_PRINT("j: %d\n", j);

    rp = a->d;
    t = ap[--i];
    DEBUG_PRINT("t: %u\n", t);
    c = (t & 1) ? CU_BN_TBIT : 0;
    DEBUG_PRINT("c: %u\n", c);
    if (t >>= 1)
        rp[i] = t;
    while (i > 0) {
        t = ap[--i];
        DEBUG_PRINT("i: %d\n", i);
        DEBUG_PRINT("t: %u\n", t);
        rp[i] = ((t >> 1) & CU_BN_MASK2) | c;
        DEBUG_PRINT("rp[%d]: %u\n", i, rp[i]);
        c = (t & 1) ? CU_BN_TBIT : 0;
        DEBUG_PRINT("c: %u\n", c);
    }
    a->top = j;
    return (1);

}

int cu_BN_lshift(VQ_VECTOR *a, unsigned n){

    if(NULL == a)
        return 0;

    if(NULL == a->d)
        return 0;

    if (BN_is_zero(a))
        return 0;

    if (0 == n)
        return 0;

    unsigned nw = 0, lb, rb, l;
    int i;
    unsigned nwb = 0, c = 0;

    nw = (n / CU_BN_BITS2);
    DEBUG_PRINT("nw: %u\n", nw);
    lb = (n % CU_BN_BITS2);
    DEBUG_PRINT("lb: %u\n", lb);
    rb = (CU_BN_BITS2 - lb);
    DEBUG_PRINT("rb: %u\n", rb);

    for(i = 31; i>=0; i--){
        if((a->d[a->top-1]>>nwb)&1){ 
            nwb = ((rb>=i)?1:0);
            break;
        }
    }

    if(nw || nwb){
        a->d = (unsigned*)realloc(a->d, (a->top+ nw + nwb)*sizeof(unsigned)) ;
        memset((a->d+a->top), 0, (nw + nwb)*sizeof(unsigned));
        memset(a->d, 0, (nw + nwb)*sizeof(unsigned));
    }

    if (lb == 0 && nw != 0 ){
        for (i = a->top - 1; i >= 0; i--){
            a->d[nw + i] = a->d[i];
            DEBUG_PRINT("a->d[nw + i]: %u\n", a->d[nw + i]);
        }
    } else {
        for (i = 0; i < (a->top + nw + nwb); i++) {
            l = a->d[i];
            DEBUG_PRINT("l = a->d[%d]: %u\n", i, l);
            a->d[i] = (l << lb) | c;
            c = (l >> rb);
            DEBUG_PRINT("after lshift a->d[%d]: %u\n", i, a->d[i]);

        }

    }
    a->top += (nw + nwb);
    return (1);

}

VQ_VECTOR *cu_euclid(VQ_VECTOR *a, VQ_VECTOR *b){

    VQ_VECTOR *t = NULL;
    unsigned shifts = 0;
    while (!CU_BN_is_zero(b)) {
        if (cu_BN_is_odd(a)) {
            if (cu_BN_is_odd(b)) {
                DEBUG_PRINT("b id odd, a is equal: %s\n", cu_bn_bn2hex(a));
                cu_bn_usub(a, b, a);
                DEBUG_PRINT("b id odd, a-b is equal: %s\n", cu_bn_bn2hex(a));
                cu_BN_rshift1(a);
                DEBUG_PRINT("b id odd, a-b>>1 is equal: %s\n", cu_bn_bn2hex(a));
                if (cu_BN_ucmp(a, b) < 0) {
                    t = a;
                    a = b;
                    b = t;
                }
            } else {            /* a odd - b even */
                DEBUG_PRINT("b id even, b is equal: %s\n", cu_bn_bn2hex(b));
                cu_BN_rshift1(b);
                DEBUG_PRINT("b id even, b>>1 is equal: %s\n", cu_bn_bn2hex(b));
                if (cu_BN_ucmp(a, b) < 0) {
                    t = a;
                    a = b;
                    b = t;
                }
            }
        } else {                /* a is even */
            if (cu_BN_is_odd(b)) {
                DEBUG_PRINT("a id even, b is odd, a is equal: %s\n", cu_bn_bn2hex(a));
                cu_BN_rshift1(a);
                DEBUG_PRINT("a id even, b is odd, a>>1 is equal: %s\n", cu_bn_bn2hex(a));
                if (cu_BN_ucmp(a, b) < 0) {
                    t = a;
                    a = b;
                    b = t;
                }
            } else {            /* a even - b even */
                DEBUG_PRINT("a id even, b is even, a is equal: %s\n", cu_bn_bn2hex(a));
                DEBUG_PRINT("a id even, b is even, b is equal: %s\n", cu_bn_bn2hex(b));
                cu_BN_rshift1(a);
                DEBUG_PRINT("a id even, b is even, a>>1 is equal: %s\n", cu_bn_bn2hex(a));
                cu_BN_rshift1(b);
                DEBUG_PRINT("a id even, b is even, b>>1 is equal: %s\n", cu_bn_bn2hex(b));
                shifts++;
            }
        }
        /* 0 <= b <= a */
    }

    if (shifts) {
        DEBUG_PRINT("a is equal: %s\n", cu_bn_bn2hex(a));
        DEBUG_PRINT("shifts is equal: %u\n", shifts);
        cu_BN_lshift(a, shifts);
        DEBUG_PRINT("a<<shifts is equal: %s\n", cu_bn_bn2hex(a));
    }
    cu_BN_free(t);
    return (a);

}