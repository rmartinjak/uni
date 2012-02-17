#ifndef ALIGN_H
#define ALIGN_H

#include "seq.h"

#ifndef PRINTLEN
#define PRINTLEN 70
#endif

typedef short score_t;
typedef short op_t;

typedef struct align_t {
	aa_t *s1;
	aa_t *s2;
	size_t len;
	score_t score;
} align_t;

enum ops { ALIGN, GAPX, GAPY };

align_t *align(score_t *submat, int gap, seq_t *s1, seq_t *s2);
void alignfree(align_t *p);
void fprintalign(FILE *f, const align_t *a);
#define printalign(a) fprintalign(stdout, a)

#endif
