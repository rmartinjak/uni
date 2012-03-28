#include "seq.h"

#include <stdlib.h>

void seqfree(seq_t *s)
{
    if (!s)
        return;

    seqfree(s->next);

    free(s->id);
    free(s->data);
    free(s);
}

void fprintaa(FILE *f, const aa_t *a, size_t offset, size_t len)
{
    a += offset;
    while (len--)
        fputc(aatoc(*a++), f);
}

aa_t ctoaa(char c)
{
#define CASE(aa) if (c == SQ(aa)) return aa
    CASE(C);
    CASE(S);
    CASE(T);
    CASE(P);
    CASE(A);
    CASE(G);
    CASE(N);
    CASE(D);
    CASE(E);
    CASE(Q);
    CASE(H);
    CASE(R);
    CASE(K);
    CASE(M);
    CASE(I);
    CASE(L);
    CASE(V);
    CASE(F);
    CASE(Y);
    CASE(W);
    return -1;
#undef CASE
}

char aatoc(aa_t a)
{
#define CASE(aa) if (a == aa) return SQ(aa)
    if (a == GAP)
        return GAP_CHAR;
    CASE(C);
    CASE(S);
    CASE(T);
    CASE(P);
    CASE(A);
    CASE(G);
    CASE(N);
    CASE(D);
    CASE(E);
    CASE(Q);
    CASE(H);
    CASE(R);
    CASE(K);
    CASE(M);
    CASE(I);
    CASE(L);
    CASE(V);
    CASE(F);
    CASE(Y);
    CASE(W);
    return -1;
#undef CASE
}
