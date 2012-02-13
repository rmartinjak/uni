#include "matrix.h"
#include "seq.h"
#include "fasta.h"
#include "align.h"

#include <stdlib.h>
#include <stdio.h>

void doalign(score_t *submat, int gap, seq_t *s1, seq_t *s2);
score_t *readsubmat(const char *path);

void doalign(score_t *submat, int gap, seq_t *s1, seq_t *s2) {
	align_t *a;

	printf("Alignment of \"%s\" (Length %d) and \"%s\" (Length %d):\n", s1->id, s1->len, s2->id, s2->len);
	a = align(submat, gap, s1, s2);

	if (!a) {
		perror("aligning sequences");
		return;
	}

	printalign(a);
	alignfree(a);
}

score_t *readsubmat(const char *path) {
	FILE *f;
	char line[8];
	int aa1;
	int aa2;
	score_t sc;
	score_t *mat;
	
	if ((mat = malloc((AA_COUNT << 5) * sizeof(score_t))) == NULL)
		return NULL;

	
	f = fopen(path, "r");
	if (!f)
		return NULL;

	while (fgets(line, 8, f)) {
		aa1 = (int)ctoaa(line[0]);
		aa2 = (int)ctoaa(line[2]);
		sc = (score_t)strtol((line+4), NULL, 10);

		if (aa1 < 0 || aa2 < 0) {
			fclose(f);
			return NULL;
		}

		mat[(aa1 << 5) + aa2] = sc;
	}

	fclose(f);

	return mat;
}

int main(int argc, char **argv) {
	int gap;
	char *endptr;

	score_t *submat;

	seq_t *seq;
	seq_t *p1;
	seq_t *p2;

	/* argumente: 
		argv[1] ist gap penalty
		argv[2] sequenzdatei
	*/

	if (argc != 3) {
		printf("usage: %s gap sequence_file\n", *argv);
		return EXIT_FAILURE;
	}

	/* argv[1] in zahl umwandeln */
	gap = (int)strtol(argv[1], &endptr, 10);
	if (*endptr != '\0') {
		printf("\"%s\" is not a number\n", argv[1]);
		return EXIT_FAILURE;
	}

	/* sequenzen lesen */
	if ((seq = fasta_read(argv[2])) == NULL) {
		perror("reading sequences");
		return EXIT_FAILURE;
	}
	/* nur eine sequenz */
	if (!seq->next) {
		printf("at least 2 sequences are needed\n");
		return EXIT_FAILURE;
	}

	/* substitutionsmatrix einlesen */
	if ((submat = readsubmat("BLOSUM62_NICE")) == NULL) {
		perror("reading substitution matrix");
		return EXIT_FAILURE;
	}

	/* 2 sequenzen */
	if (!seq->next->next) {
		doalign(submat, gap, seq, seq->next);
		printf("\n");
	}

	/* mehr als 2 -> jede mit jeder */
	else {
		p1 = seq;
		while (p1->next) {
			p2 = p1->next;

			for (p2 = p1->next; p2; p2 = p2->next) {
				doalign(submat, gap, p1, p2);
				printf("\n");
			}
			p1 = p1->next;
		}
	}

	seqfree(seq);
	free(submat);

	return EXIT_SUCCESS;
}
