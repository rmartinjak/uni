#include "mr.c"

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <time.h>
#include <inttypes.h>

int main(int argc, char **argv)
{
    uintmax_t min, max, n, x;

    srand(time(NULL));

    min = UINTMAX_MAX;
    max = 0;

    n = 10000000;
    while (n--) {
        x = randumax(0, UINTMAX_MAX);
        if (x < min)
            min = x;
        if (x > max)
            max = x;
    }

    printf("min: %" PRIuMAX "\n", min);
    printf("max: %" PRIuMAX "\n", max);
    printf("max: %" PRIuMAX " (diff)\n", UINTMAX_MAX - max);
    printf("UINTMAX_MAX: %" PRIuMAX "\n", UINTMAX_MAX);

    return EXIT_SUCCESS;
}
