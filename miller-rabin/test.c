#include "mr.c"

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>

int main(int argc, char **argv) {
    char *endptr;
    unsigned long long min, max, n;

    if (argc < 2) {
        printf("usage: %s <min> <max>\n", argv[0]);
        return EXIT_FAILURE;
    }

    min = strtoull(argv[1], &endptr, 10);
    if (*endptr || *argv[1] == '\0') {
        printf("min has wrong format\n");
        return EXIT_FAILURE;
    }

    max = strtoull(argv[2], &endptr, 10);
    if (*endptr || *argv[2] == '\0') {
        printf("max has wrong format\n");
        return EXIT_FAILURE; }

    if (min > max) {
        printf("min must be <= max\n");
        return EXIT_FAILURE;
    }

    for (n = min; n <= max; n++) {
        if (millerrabin(n)) printf("%llu\n", n);
    }

    return EXIT_SUCCESS;
}
