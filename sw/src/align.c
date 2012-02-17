#include "align.h"

#include "seq.h"
#include "matrix.h"

#include <sys/types.h>
#include <errno.h>
#include <stdlib.h>

void fprintalign(FILE *f, const align_t *a) {
	size_t off=0;
	size_t len = a->len;

	fprintf(f, "    Score: %d, Length: %d\n\n", a->score, a->len);

	while (off < a->len) {
		len = (a->len - off > PRINTLEN) ? PRINTLEN : a->len - off;

		fprintf(f, "%5d  ", off+1);
		fprintaa(f, a->s1, off, len);
		fprintf(f, "\n");

		fprintf(f, "       ");
		fprintaa(f, a->s2, off, len);
		fprintf(f, "\n");

		fprintf(f, "\n");
		off += PRINTLEN;
	}
}

void alignfree(align_t *p) {
	if (!p)
		return;
	free(p->s1);
	free(p->s2);
	free(p);
}

#define SUBMAT(x, y) submat[(x << 5) + y]
align_t *align(score_t *submat, int gap, seq_t *s1, seq_t *s2) {
	align_t *a;

	aa_t *a1, *a2;
	int align_len;

	score_t **dp_mat;

	int i, k;
	int imax=0, kmax=0, scmax=0;
	score_t sc, tmp;

	if ((a = malloc(sizeof(align_t))) == NULL) {
		errno = ENOMEM;
		return NULL;
	}

	if ((dp_mat = (score_t **)matrixalloc(sizeof(score_t), s1->len+1, s2->len+1)) == NULL) {
		errno = ENOMEM;
		return NULL;
	}

	/* erste zeile/spalte auf 0 */
	for (i=0; i <= s1->len; ++i)
		dp_mat[i][0] = 0;
	for (i=0; i <= s2->len; ++i)
		dp_mat[0][i] = 0;

	/* restliche felder anhand der rekursionsvorschrift berechnen */
	for (i=1; i <= s1->len; ++i) {
		for (k=1; k <= s2->len; ++k) {
			/* aligned? */
			sc = dp_mat[i-1][k-1] + SUBMAT(s1->data[i-1], s2->data[k-1]);

			/* gap in s1 */
			if ((tmp = dp_mat[i][k-1] - gap) > sc)
				sc = tmp;

			/* gap in s2 */
			if ((tmp = dp_mat[i-1][k] - gap) > sc)
				sc = tmp;

			if (sc < 0)
				sc = 0;

			dp_mat[i][k] = sc;

			if (sc > scmax) {
				imax = i;
				kmax = k;
				scmax = sc;

			}
		}
	}


	/* BACKTRACE */

	i = imax;
	k = kmax;
	
	/* laenge (= benoetigter speicher) feststellen */
	align_len = 0;
	while ((sc = dp_mat[i][k]) > 0) {
		align_len++;

		if (sc == dp_mat[i-1][k-1] + SUBMAT(s1->data[i-1], s2->data[k-1])) {
			i--;
			k--;
		}
		else if (sc == dp_mat[i][k-1] - gap)
			k--;
		else
			i--;
	}

	a->s1 = calloc(align_len, sizeof(aa_t));
	a->s2 = calloc(align_len, sizeof(aa_t));

	if (!a->s1 || !a->s2) {
		free(a->s1);
		free(a->s2);
		errno = ENOMEM;
		return NULL;
	}

	/* laenge und score stehen schon fest */
	a->len = align_len;
	a->score = scmax;

	/* sequenzdaten werden von hinten eingefuegt */
	a1 = a->s1 + align_len-1;
	a2 = a->s2 + align_len-1;

	i = imax;
	k = kmax;

	/* pfad nochmal durchgehen und alignmentdaten schreiben */
	while ((sc = dp_mat[i][k]) > 0) {

		if (sc == dp_mat[i-1][k-1] + SUBMAT(s1->data[i-1], s2->data[k-1])) {
			*a1 = s1->data[i-1];
			*a2 = s2->data[k-1];
			i--;
			k--;
		}
		else if (sc == dp_mat[i][k-1] - gap) {
			*a1 = GAP;
			*a2 = s2->data[k-1];
			k--;
		}
		else {
			*a1 = s1->data[i-1];
			*a2 = GAP;
			i--;
		}

		a1--;
		a2--;
	}

	matrixfree((void **)dp_mat, s1->len+1);

	return a;
}
