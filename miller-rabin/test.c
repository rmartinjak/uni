#include "mr.c"

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <time.h>

int main(int argc, char **argv) {
    char *endptr;
    unsigned long long min, max, n;
    unsigned long long x;

    srand(time(NULL));

    /*
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

    for (n = min | 1; n <= max; n += 2) {
        if (millerrabin(n)) printf("%llu\n", n);
    }
    */
    min = ULLONG_MAX;
    max = 0;

    n = 10000000;
    while (n--) {
        x = randull(0, ULLONG_MAX);
        if (x < min)
            min = x;
        if (x > max)
            max = x;
    }

    printf("min: %llu\n", min);
    printf("max: %llu\n", max);
    printf("max: %llu (diff)\n", ULLONG_MAX - max);
    printf("ULLONG_MAX: %llu\n", ULLONG_MAX);

    return EXIT_SUCCESS;
}
