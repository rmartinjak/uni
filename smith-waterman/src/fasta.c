#include "fasta.h"
#include "seq.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>

static int fasta_read_id(FILE *f, seq_t *s);
static int fasta_read_data(FILE *f, seq_t *s);

/* id lesen und kommentare ueberspringen */
static int fasta_read_id(FILE *f, seq_t *s)
{
    int c;
    char line[1024];

    /* id muss mit > beginnen */
    if ((c = fgetc(f)) != '>') {
        errno = EINVAL;
        return -1;
    }

    /* zeile einlesen und zeilenumbruch abschneiden */
    fgets(line, sizeof(line), f);
    if (line[strlen(line)-1] == '\n')
        line[strlen(line)-1] = '\0';

    /* "leere" ID ist nicht erlaubt */
    if (*line == '\0') {
        errno = EINVAL;
        return -1;
    }
    /* id aus zeilenpuffer nach s kopieren */
    if ((s->id = strdup(line)) == NULL) {
        errno = ENOMEM;
        return -1;
    }

    /* eventuelle kommentarzeilen ueberspringen */
    while ((c = fgetc(f)) == ';') {
        fgets(line, sizeof(line), f);
    }

    /* dateiende */
    if (c == EOF)
        return EOF;

    /* das letzte gelesene zeichen enthaelt daten! wieder auf den stream schieben */
    ungetc(c, f);
    return 0;
}

/* sequenzdaten lesen */
static int fasta_read_data(FILE *f, seq_t *s)
{
    int c;
    fpos_t fpos;
    int seq_len;
    aa_t *a;

    /* anfangsposition merken */ 
    fgetpos(f, &fpos);
    seq_len = 0;

    /* laenge der sequenz ermitteln */
    while ((c = fgetc(f)) != '>' && c != EOF) {
        if (c != '\n')
            seq_len++;
    }

    a = calloc((seq_len+1), sizeof(aa_t));
    if (!a) {
        errno = ENOMEM;
        return -1;
    }

    s->len = seq_len;
    s->data = a;

    fsetpos(f, &fpos);

    /* sequenzdaten einlesen */
    while (seq_len) {
        if ((c = fgetc(f)) != '\n') {
            if ((*a++ = ctoaa(c)) == -1) {
                free(a);
                s->len = 0;
                errno = EINVAL;
                return -1;
            }
            seq_len--;
        }
    }
    *a = 0;

    /* zeilenumbruch der letzten datenzeile */
    fgetc(f);

    /* anfang der neuen zeile */
    c = fgetc(f);

    if (c == '>') {
        ungetc(c, f);
        return 0;
    }
    else
        return EOF;
}

/* sequenzen aus datei im fasta-format lesen (gibt verkettete liste zurueck) */
seq_t *fasta_read(const char *path)
{
    FILE *f;
    seq_t *ret;
    seq_t *s;

    ret = malloc(sizeof(seq_t));
    if (!ret) {
        errno = ENOMEM;
        return NULL;
    }
    SEQINIT(ret);
    s = ret;

    f = fopen(path, "r");
    if (!f)
        return NULL;

    /* abwechselnd id und daten lesen, bei fehler oder EOF abbruch */
    while (1) {
        /* fehler oder EOF beim lesen von id fuehrt zum abbruch (nach einer id muessen immer daten kommen!) */
        if (fasta_read_id(f, s) != 0) {
            seqfree(ret);
            fclose(f);
            return NULL;
        }

        errno = 0;
        if (fasta_read_data(f, s) != 0) {
            /* fehler */
            if (errno) {
                seqfree(ret);
                fclose(f);
                return NULL;
            }
            /* EOF, einlesen korrekt beendet */
            else {
                break;
            }
        }
        /* speicher fuer naechste sequenz bereitstellen */
        else {
            if ((s->next = malloc(sizeof(seq_t))) == NULL) {
                seqfree(ret);
                fclose(f);
                errno = ENOMEM;
                return NULL;
            }
            s = s->next;
            SEQINIT(s);
        }
    }

    fclose(f);
    return ret;
}
