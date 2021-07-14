#include <stdio.h>
#include <stdlib.h>
#include "solvers.h"
#include "physics.h"
#include <time.h>

void orbit(double tSpan[], double y0[], char output[], double ATOL) {
    /*
    Propagates orbit using RK810 solver.
    tSpan: t0, tf of times.
    y0: initial x y z vx vy vz state vector
    output: filename of output csv data file
    ATOL: tolerance
    */
    clock_t start, end;

    start = clock();
    solution* result = RK810vec(combined_perturbations, tSpan, y0, ATOL);
    end = clock();
    double cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    printf("Integration using RK810 solver completed in %.8f seconds. Steps taken: %d\n", cpu_time_used, result->n);

    FILE *fptr = fopen(output, "w");

    for (int i = 0; i < result->n; i++) {
        fprintf(fptr, "%.15f, %.15f, %.15f, %.15f, %.15f, %.15f, %.15f\n", result->t[i], result->y[i][0], result->y[i][1], result->y[i][2], result->y[i][3], result->y[i][4], result->y[i][5]);
        free(result->y[i]);
    }
    fclose(fptr);

    free(result->t);
    free(result);
}