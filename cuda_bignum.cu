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

U_BN *cu_bn_new(void)
{
    U_BN *ret;

    if ((ret = (U_BN *)malloc(sizeof(*ret))) == NULL) {
        return (NULL);
    }
    ret->top = 0;
    ret->d = (unsigned *) malloc(sizeof(unsigned));
    return (ret);
}

void cu_bn_free(U_BN *a)
{
    if (a == NULL)
        return;
    if (a->d != NULL)
        free(a->d);
    if (a == NULL)
        free(a);
}

unsigned cu_bn_mul_words(unsigned  *rp, const unsigned  *ap, int num, unsigned  w){

    unsigned i=0;
    unsigned tmp;
    unsigned  c1 = 0;
    if (num <= 0)
        return 0;
    while (num) {
        tmp= Lw(((unsigned long)(ap[i]) * (w) + (c1)));
        c1= Hw(((unsigned long)(ap[i]) * (w) + (c1)));
        //DEBUG_PRINT("ap[%u]: %u\n", i, ap[i]);
        //DEBUG_PRINT("rp[%u]: %u\n", i, tmp);
        rp[i] = tmp;
        i++;
        num--;
    }
    //DEBUG_PRINT("w: %u\n", w);
    return (c1);

}

/*__device__ void * cu_realloc(void * old_p, unsigned new_s)
{
    void* new_p = (void*) malloc (new_s);
    unsigned old_s = sizeof(old_p);
    if(old_s <= new_s){
        memcpy(old_p, new_p, old_s);
    } else {
        memcpy(old_p, new_p, new_s);
    }
    free(old_p);
    return (new_p);
}*/

int cu_bn_mul_word(U_BN *a, unsigned w){

    unsigned ll;
    if(NULL == a)
        return 0;

    if(NULL == a->d)
        return 0;

    w &= CU_BN_MASK2;

    if (a->top) {
        if (w == 0)
            cu_bn_zero(a);
        else {
            ll = cu_bn_mul_words(a->d, a->d, a->top, w);
            if (ll) {
                a->d = ((unsigned*)realloc(a->d, (++a->top)*sizeof(unsigned))); 
                a->d[(a->top-1)] = ll;
                //DEBUG_PRINT("a->d[(a->top-1)]: %u\n", a->d[(a->top-1)]);
                //DEBUG_PRINT("a->top: %u\n", a->top);
            }
        }
    }
    return (1);

}

int cu_bn_set_word(U_BN *a, unsigned w){

    if(NULL == a)
        return 0;

    if(NULL == a->d)
        return 0;

    a->d = ((unsigned *)malloc(sizeof(unsigned)));
    a->d[0] = w;
    a->top = 1;
    return (1);

}

int cu_bn_add_word(U_BN *a, unsigned w){

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
    if (cu_bn_is_zero(a))
        return cu_bn_set_word(a, w);

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

int cu_bn_dec2bn(U_BN *ret, const char *a){

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
        return cu_bn_set_word(ret, 0);

    for (i = 0; i <= (INT_MAX/4) && isdigit((unsigned char)a[i]); i++)
        continue;


    //DEBUG_PRINT("number of digits: %u\n", i);

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
            //DEBUG_PRINT("l: %u\n", l);
            cu_bn_mul_word(ret, CU_BN_DEC_CONV);
            cu_bn_add_word(ret, l);
            l = 0;
            j = 0;
        }
    }
    cu_bn_correct_top(ret);
    return (1);

}

/* Must 'OPENSSL_free' the returned data */
char *cu_bn_bn2hex(const U_BN *a){

    int i, j, v, z = 0;
    char *buf = NULL, *p = NULL;
    char Hex[] = "0123456789ABCDEF";

    if(NULL == a)
        return 0;

    if(NULL == a->d)
        return 0;

    if (cu_bn_is_zero(a))
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

int cu_bn_ucmp(const U_BN *a, const U_BN *b){

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
int cu_bn_usub(const U_BN *a, const U_BN *b, U_BN *r){

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
            t1 = (t1 - t2 - 1);// & BN_MASK2;
        } else {
            carry = (t1 < t2);
            t1 = (t1 - t2);// & BN_MASK2;
        }
# if defined(IRIX_CC_BUG) && !defined(LINT)
        dummy = t1;
# endif
        *(rp++) = t1;// & BN_MASK2;
    }
#else
    carry = bn_sub_words(rp, ap, bp, min);
    ap += min;
    bp += min;
    rp += min;
#endif
    if (carry) {              /* subtracted */
        if (!dif)
            /* error: a < b */
            return 0;
        while (dif) {
            dif--;
            t1 = *(ap++);
            t2 = (t1 - 1);// & BN_MASK2;
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

/* unsigned subtraction of b from a, a must be larger than b. */
int cu_bn_usub_optimized(const U_BN *a, const U_BN *b, U_BN *r){

   unsigned max, min, dif;
    register unsigned t1, t2, *ap, *bp, *rp;
    int i, carry;

    max = a->top;
    min = b->top;
    dif = max - min;

    ap = a->d;
    bp = b->d;
    rp = r->d;
    carry = 0;
    for (i = 0; i < min; i++) {
        t1 = *(ap++);
        t2 = *(bp++);
        *(rp++) = (t1 - t2 - carry);
        carry = (t1 < (t2 + carry));
    }
    while (dif) {
        t1 = *(ap++);
        t2 = (t1 - carry);
        carry = (carry > t1);
        *(rp++) = t2;
        dif--;
    }
    r->top = max;
    cu_bn_correct_top(r);
    return (1);

}

int cu_bn_num_bits_word(long l){

     unsigned char bits[256] = {
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

int cu_bn_num_bits(const U_BN *a){

    int i = a->top - 1;

    if (cu_bn_is_zero(a))
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

int cu_bn_rshift1(U_BN *a){

    if(NULL == a)
        return 0;

    if(NULL == a->d)
        return 0;

    if (cu_bn_is_zero(a))
        return 0;

    unsigned *ap, *rp , t, c;
    int i, j;

    i = a->top;
    //DEBUG_PRINT("a->top: %u\n", a->top);
    ap = a->d;
    //DEBUG_PRINT("a->d[0]: %u\n", a->d[0]);

    j = i - (ap[i - 1] == 1);
    //DEBUG_PRINT("j: %d\n", j);

    rp = a->d;
    t = ap[--i];
    //DEBUG_PRINT("t: %u\n", t);
    c = (t & 1) ? CU_BN_TBIT : 0;
    //DEBUG_PRINT("c: %u\n", c);
    if (t >>= 1)
        rp[i] = t;
    while (i > 0) {
        t = ap[--i];
        //DEBUG_PRINT("i: %d\n", i);
        //DEBUG_PRINT("t: %u\n", t);
        rp[i] = ((t >> 1) & CU_BN_MASK2) | c;
        //DEBUG_PRINT("rp[%d]: %u\n", i, rp[i]);
        c = (t & 1) ? CU_BN_TBIT : 0;
        //DEBUG_PRINT("c: %u\n", c);
    }
    a->top = j;
    return (1);

}

int cu_bn_lshift(U_BN *a, unsigned n){

    if(NULL == a)
        return 0;

    if(NULL == a->d)
        return 0;

    if (cu_bn_is_zero(a))
        return 0;

    if (0 == n)
        return 0;

    unsigned nw = 0, lb, rb, l;
    int i;
    unsigned nwb = 0, c = 0;

    nw = (n / CU_BN_BITS2);
    //DEBUG_PRINT("nw: %u\n", nw);
    lb = (n % CU_BN_BITS2);
    //DEBUG_PRINT("lb: %u\n", lb);
    rb = (CU_BN_BITS2 - lb);
    //DEBUG_PRINT("rb: %u\n", rb);

    l=a->d[a->top-1];
    if( (l >> rb) > 0 ) nwb = 1;

    if(nw || nwb){
        //a->d = (unsigned*)realloc(a->d, (a->top + nw + nwb)*sizeof(unsigned)) ;
        //a->d[a->top]=0;
        //memset((a->d+a->top-1), 0, (nw + nwb));
        //memset(a->d, 0, (nw + nwb)*sizeof(unsigned));
    }
    //DEBUG_PRINT("nwb: %u\n", nwb);

    if (lb == 0 && nw != 0 ){
        for (i = a->top - 1; i >= 0; i--){
            a->d[nw + i] = a->d[i];
            //DEBUG_PRINT("a->d[nw + i]: %u\n", a->d[nw + i]);
        }
    } else {
        for (i = 0; i < (a->top + nw + nwb); i++) {
            l = a->d[i];
            //DEBUG_PRINT("l = a->d[%d]: %u\n", i, l);
            a->d[i] = (l << lb) | c;
            c = (l >> rb);
            //DEBUG_PRINT("after lshift a->d[%d]: %u\n", i, a->d[i]);

        }

    }
    a->top += (nw + nwb);
    //DEBUG_PRINT("a is equal: %s\n", cu_bn_bn2hex(a));
    return (1);

}

U_BN *cu_euclid(U_BN *a, U_BN *b){

    U_BN *t = NULL;

    if (cu_bn_ucmp(a, b) < 0) {
        t = a;
        a = b;
        b = t;
    }

    unsigned shifts = 0;
    while (!cu_bn_is_zero(b)) {
        if (cu_bn_is_odd(a)) {
            if (cu_bn_is_odd(b)) {
                DEBUG_PRINT("b id odd, a is odd, a is equal: %s\n", cu_bn_bn2hex(a));
                cu_bn_usub(a, b, a);
                DEBUG_PRINT("b id odd, a is odd, a-b is equal: %s\n", cu_bn_bn2hex(a));
                cu_bn_rshift1(a); // beacuse subtraction is even
                DEBUG_PRINT("b id odd, a is odd, a-b>>1 is equal: %s\n", cu_bn_bn2hex(a));
                if (cu_bn_ucmp(a, b) < 0) {
                    t = a;
                    a = b;
                    b = t;
                }
            } else {   
                DEBUG_PRINT("b id even, b is equal: %s\n", cu_bn_bn2hex(b));
                cu_bn_rshift1(b);
                DEBUG_PRINT("b id even, b>>1 is equal: %s\n", cu_bn_bn2hex(b));
                if (cu_bn_ucmp(a, b) < 0) {
                    t = a;
                    a = b;
                    b = t;
                }
            }
        } else {   
            if (cu_bn_is_odd(b)) {
                DEBUG_PRINT("a id even, b is odd, a is equal: %s\n", cu_bn_bn2hex(a));
                cu_bn_rshift1(a);
                DEBUG_PRINT("a id even, b is odd, a>>1 is equal: %s\n", cu_bn_bn2hex(a));
                if (cu_bn_ucmp(a, b) < 0) {
                    t = a;
                    a = b;
                    b = t;
                }
            } else {     
                DEBUG_PRINT("a id even, b is even, a is equal: %s\n", cu_bn_bn2hex(a));
                DEBUG_PRINT("a id even, b is even, b is equal: %s\n", cu_bn_bn2hex(b));
                cu_bn_rshift1(a);
                DEBUG_PRINT("a id even, b is even, a>>1 is equal: %s\n", cu_bn_bn2hex(a));
                cu_bn_rshift1(b);
                DEBUG_PRINT("a id even, b is even, b>>1 is equal: %s\n", cu_bn_bn2hex(b));
                shifts++;
            }
        }
    }

    if (shifts) {
        DEBUG_PRINT("a is equal: %s\n", cu_bn_bn2hex(a));
        DEBUG_PRINT("shifts is equal: %u\n", shifts);
        cu_bn_lshift(a, shifts);
        DEBUG_PRINT("a<<shifts is equal: %s\n", cu_bn_bn2hex(a));
    }
    //cu_bn_free(t);
    return (a);

}

int bignum2u_bn(BIGNUM* bignum, U_BN *u_bn){

    if(NULL == bignum)
        return 0;

    if(NULL == bignum->d)
        return 0;

    if(NULL == u_bn)
        return 0;

    if(NULL == u_bn->d)
        return 0;

    u_bn->top = ( (2 * bignum->top) + 1);
    u_bn->d = (unsigned *) malloc ( sizeof(unsigned) * u_bn->top );
    memcpy(u_bn->d, bignum->d, ( sizeof(unsigned) * u_bn->top ));

    return (1);
}

U_BN *cu_fast_binary_euclid(U_BN *a, U_BN *b){
    U_BN *t = NULL;
    do {
        if (cu_bn_ucmp(a, b) < 0) {
            t = a;
            a = b;
            b = t;
        }
        DEBUG_PRINT("a is equal: %s\n", cu_bn_bn2hex(a));
        DEBUG_PRINT("b is equal: %s\n", cu_bn_bn2hex(b));
        if(!cu_bn_usub(a, b, a)) break;
        DEBUG_PRINT("a is equal: %s\n", cu_bn_bn2hex(a));
        DEBUG_PRINT("b is equal: %s\n", cu_bn_bn2hex(b));
        while(!(a->d[0]&1)) {
            if(!cu_bn_rshift1(a)) break;
        }
    } while (!cu_bn_is_zero(b));
    return (a);
}

U_BN *cu_classic_euclid(U_BN *a, U_BN *b){

    while (cu_bn_ucmp(a, b) != 0) {
        if (cu_bn_ucmp(a, b) > 0) {
            cu_bn_usub(a, b, a); 
        }
        else {
            cu_bn_usub(b, a, b);
        }
        DEBUG_PRINT("a is equal: %s\n", cu_bn_bn2hex(a));
    }

    return (a);

}


int cu_ubn_copy(U_BN *a, const U_BN *b)
{
    int i;
    unsigned a0;

    if (a == b)
        return (0);

    for (i = 0; i < b->top; i++) {
        a0 = b->d[i];
        a->d[i] = a0;
    }
    //memcpy(a->d, b->d, sizeof(b->d[0]) * b->top);
    a->top = b->top;
    return (1);
}


unsigned cu_ubn_add_words(unsigned *r, const unsigned *a, const unsigned *b, int n)
{
    long ll = 0;
    int i;
    if (n <= 0)
        return (0);

    for (i = 0; i < n; i++) {
        ll += (BN_ULLONG) a[i] + b[i];
        r[i] = (long)ll & CU_BN_MASK2;
        ll >>= CU_BN_BITS2;
    }
    return (ll);
}

int cu_ubn_uadd(const U_BN *a, const U_BN *b, U_BN *r)
{
    int max, min, dif;
    unsigned *ap, *bp, *rp, carry, t1, t2;
    const U_BN *tmp;

    if (a->top < b->top) {
        tmp = a;
        a = b;
        b = tmp;
    }
    max = a->top;
    min = b->top;
    dif = max - min;

    r->top = max;

    ap = a->d;
    bp = b->d;
    rp = r->d;

    carry = cu_ubn_add_words(rp, ap, bp, min);
    //rp += min;
    //ap += min;
    //bp += min;

    if (carry) {
        while (dif) {
            dif--;
            t1 = *(ap++);
            t2 = (t1 + 1) & CU_BN_MASK2;
            *(rp++) = t2;
            if (t2) {
                carry = 0;
                break;
            }
        }
        if (carry) {
            /* carry != 0 => dif == 0 */
            *rp = 1;
            r->top++;
        }
    }
    if (dif && rp != ap)
        while (dif--)
            /* copy remaining words if ap != rp */
            *(rp++) = *(ap++);

    return (1);
}

U_BN *q_algorithm_PM(U_BN *a, U_BN *b){
    int q=0;
    U_BN *t = NULL;
    while (!cu_bn_is_zero(b)) {
        DEBUG_PRINT("q: %d b: %s\n", q, cu_bn_bn2hex(b));
        while(!cu_bn_is_odd(b)){
            cu_bn_rshift1(b);
            q++;
        }
        if(q >= 0) {
            t = a;
            a = b;
            b = t;
            q=-q;
        }
        if((a->d[0]+b->d[0])&3==0) {
            cu_ubn_uadd(a, b, b);
            cu_bn_rshift1(b);
        } else {
            cu_bn_usub(a, b, b);
            cu_bn_rshift1(b);
        }
    }
    return (a);

}

U_BN *algorithm_PM(U_BN *a, U_BN *b, unsigned keysize){
    int beta=keysize, alfa=keysize, tmp;
    U_BN *t = NULL;
    while (!cu_bn_is_zero(b)) {
        //DEBUG_PRINT("b: %s\n", cu_bn_bn2hex(b));
        while(!cu_bn_is_odd(b)){
            cu_bn_rshift1(b);
            beta--;
        }
        if(alfa >= beta) {
            t = a;
            a = b;
            b = t;
            tmp = alfa;
            alfa = beta;
            beta = tmp;
            
        }
        if((a->d[0]+b->d[0])&3==0) {
            cu_ubn_uadd(a, b, b);
            cu_bn_rshift1(b);
        } else {
            cu_bn_usub(a, b, b);
            cu_bn_rshift1(b);
        }
    }
    return (a);

}

U_BN *cu_approximate_euclid(U_BN *a, U_BN *b, unsigned keysize){
    int beta=keysize, alfa=keysize, tmp;
    U_BN *t = NULL;
    do{
        if(beta==0){
            if(!(alfa&1)) alfa --;
            cu_bn_mul_word(b, alfa);
            cu_bn_usub(a, b, a);
        } else {
            
        }
    } while(!cu_bn_is_odd(b));


    while (cu_bn_ucmp(a, b) != 0) {
        if (cu_bn_ucmp(a, b) > 0) {
            cu_bn_usub(a, b, a); 
        }
        else {
            cu_bn_usub(b, a, b);
        }
        DEBUG_PRINT("a is equal: %s\n", cu_bn_bn2hex(a));
    }

    return (a);

}