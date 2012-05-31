#include "align.h"

#include "seq.h"
#include "matrix.h"

#include <sys/types.h>
#include <errno.h>
#include <stdlib.h>

#define SUBMAT(x, y) submat[(x << 5) + y]

typedef struct dpmat_t
{
    score_t m;
    score_t ix;
    score_t iy;
    score_t max;
}
dpmat_t;

void fprintalign(FILE *f, const align_t *a)
{
    size_t off=0;
    size_t len = a->len;

    fprintf(f, "    Score: %d, Length: %zu\n\n", a->score, a->len);

    while (off < a->len) {
        len = (a->len - off > PRINTLEN) ? PRINTLEN : a->len - off;

        fprintf(f, "%5zu  ", off+1);
        fprintaa(f, a->s1, off, len);
        fprintf(f, "\n");

        fprintf(f, "       ");
        fprintaa(f, a->s2, off, len);
        fprintf(f, "\n");

        fprintf(f, "\n");
        off += PRINTLEN;
    }
}

void alignfree(align_t *p)
{
    if (!p)
        return;
    free(p->s1);
    free(p->s2);
    free(p);
}

align_t *align(score_t *submat, score_t gap_start, score_t gap_cont, seq_t *s1, seq_t *s2)
{
    align_t *a;

    aa_t *a1, *a2;
    size_t align_len;

    dpmat_t **dp_mat;

    int i, k;
    int imax=0, kmax=0, scmax=0;
    score_t sc, tmp_max, tmp;

    if ((a = malloc(sizeof(align_t))) == NULL) {
        errno = ENOMEM;
        return NULL;
    }

    if ((dp_mat = (dpmat_t **)matrixalloc(sizeof(dpmat_t), s1->len+1, s2->len+1)) == NULL) {
        errno = ENOMEM;
        return NULL;
    }

    /* erste zeile/spalte auf 0 */
    for (i=0; i <= s1->len; ++i) {
        dp_mat[i][0].m = 0;
        dp_mat[i][0].ix = 0;
        dp_mat[i][0].iy = 0;
        dp_mat[i][0].max = 0;
    }
    for (i=0; i <= s2->len; ++i) {
        dp_mat[0][i].m = 0;
        dp_mat[0][i].ix = 0;
        dp_mat[0][i].iy = 0;
        dp_mat[0][i].max = 0;
    }

    /* restliche felder anhand der rekursionsvorschrift berechnen */
    for (i=1; i <= s1->len; ++i) {
        for (k=1; k <= s2->len; ++k) {

            /*---*/
            /* m */
            /*---*/
            sc = dp_mat[i-1][k-1].m + SUBMAT(s1->data[i-1], s2->data[k-1]);

            tmp = dp_mat[i-1][k-1].ix + SUBMAT(s1->data[i-1], s2->data[k-1]);
            if (tmp > sc) sc = tmp;

            tmp = dp_mat[i-1][k-1].iy + SUBMAT(s1->data[i-1], s2->data[k-1]);
            if (tmp > sc) sc = tmp;

            if (sc < 0) sc = 0;
            dp_mat[i][k].m = sc;
            tmp_max = sc;

            /*----------------*/
            /* ix (gap in s2) */
            /*----------------*/
            sc = dp_mat[i-1][k].m - gap_start;

            tmp = dp_mat[i-1][k].ix - gap_cont;
            if (tmp > sc) sc = tmp;

            /* score >= 0 */
            if (sc < 0) sc = 0;

            dp_mat[i][k].ix = sc;
            if (sc > tmp_max) tmp_max = sc;

            /*----------------*/
            /* iy (gap in s1) */
            /*----------------*/
            sc = dp_mat[i][k-1].m - gap_start;

            tmp = dp_mat[i][k-1].iy - gap_cont;
            if (tmp > sc) sc = tmp;

            /* score >= 0 */
            if (sc < 0) sc = 0;
            dp_mat[i][k].iy = sc;
            if (sc > tmp_max) tmp_max = sc;

            dp_mat[i][k].max = tmp_max;

            if (tmp_max > scmax) {
                imax = i;
                kmax = k;
                scmax = tmp_max;
            }
        }
    }


    /* BACKTRACE */

    i = imax;
    k = kmax;

    /* laenge (= benoetigter speicher) feststellen */
    align_len = 0;
    while ((sc = dp_mat[i][k].max) > 0) {
        align_len++;

        if (sc == dp_mat[i][k].m) {
            i--;
            k--;
        }
        else if (sc == dp_mat[i][k].ix)
            i--;
        else
            k--;
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
    while ((sc = dp_mat[i][k].max) > 0) {
        if (sc == dp_mat[i][k].m) {
            *a1 = s1->data[i-1];
            *a2 = s2->data[k-1];
            i--;
            k--;
        }
        else if (sc == dp_mat[i][k].ix) {
            *a1 = s1->data[i-1];
            *a2 = GAP;
            i--;
        }
        else {
            *a1 = GAP;
            *a2 = s2->data[k-1];
            k--;
        }

        a1--;
        a2--;
    }

    matrixfree((void **)dp_mat, s1->len+1);

    return a;
}
