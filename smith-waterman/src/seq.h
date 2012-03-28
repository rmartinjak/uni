#ifndef SEQ_H
#define SEQ_H

#include <stdio.h>
#include <sys/types.h>

#define SQ(x) ((#x)[0])

typedef int aa_t;
typedef struct seq_t
{
    char *id;
    aa_t *data;
    size_t len;
    struct seq_t *next;
}
seq_t;

#define GAP_CHAR '-'
#define AA_COUNT 20
enum aminoacids { C, S, T, P, A, G, N, D, E, Q, H, R, K, M, I, L, V, F, Y, W, GAP };

#define SEQINIT(p) { (p)->id = NULL; (p)->data = NULL; (p)->len = 0; (p)->next = NULL; }
void seqfree(seq_t *s);

aa_t ctoaa(char c);
char aatoc(aa_t a);

void fprintaa(FILE *f, const aa_t *a, size_t offset, size_t len);

#endif
