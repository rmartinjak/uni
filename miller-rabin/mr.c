#include "mr.h"

#include <stdlib.h>
#include <limits.h>

#define EVEN(n) (!(n % 2))

#define WITNESSES 200

#define ULLONG_BIT (sizeof(unsigned long long) * CHAR_BIT)


unsigned long long shift_mod(unsigned long long a, unsigned int sh, unsigned long long mod)
{
    while (sh--) {
        a <<= 1;
        a %= mod;
    }
    return a;
}

unsigned long long mult_mod(unsigned long long a, unsigned long long b, unsigned long long mod)
{
    unsigned long long ret = 0;
    unsigned int sh = 0;

    while (b) {
        if (b & 1) {
            ret += shift_mod(a, sh, mod);
            ret %= mod;
        }
        b >>= 1;
        sh++;
    }
    return ret;
}

unsigned long long pow_mod(unsigned long long base, unsigned long long exp, unsigned long long mod)
{
    unsigned long long ret = 1;
    unsigned long long b;

    b = base;

    while (exp)
    {
        if (exp & 1) {
            ret = mult_mod(ret, b, mod);
        }
        exp >>= 1;
        b = mult_mod(b, b, mod);
    }
    return (unsigned long long)ret;
}

unsigned long long randull(unsigned long long min, unsigned long long max)
{
    static size_t rand_bits = 0;
    size_t bits;
    unsigned long long rnd = 0;

    if (!rand_bits)
    {
        int r;
        for (r = RAND_MAX; r > 0; r >>= 1)
            rand_bits++;
    }

    for (bits = 0; bits < ULLONG_BIT; bits += rand_bits)
    {
        rnd <<= rand_bits;
        rnd |= rand();
    }

    rnd = min + max * ((double)rnd / ULLONG_MAX);
    return rnd;
}

int millerrabin(unsigned long long n)
{
    unsigned long long d, s, i, a;
    int k;
    int r1, r2;

    if (EVEN(n) || n <= 2)
        return 0;


    d = n-1;
    s = 0;
    while (EVEN(d)) {
        s++;
        d >>= 1;
    }

    for (k = 0; k < WITNESSES; k++) {
        a = randull(2, n-2);
        r1 = 0;
        r2 = 0;

        /* first test */
        if (pow_mod(a, d, n) == 1) {
            r1 = 1;
        }

        /* second test */
        for (i = 0; i < s; i++) {
            if (pow_mod(a, (1 << i) * d, n) == n-1) {
                r2 = 1;
                break;
            }
        }
        if (!r1 && !r2) {
            return 0;
        }
    }

    return 1;
}
