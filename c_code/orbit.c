#include <stdio.h>
#include <stdlib.h>
#include "solvers.h"
#include "physics.h"
#include <time.h>

void orbit(double tSpan[], double y0[], char integrator[], int ord, char output[]) {
    /*
    tSpan: t0, tf of times.
    y0: initial x y z vx vy vz state vector
    integrator: path to a butcher tableau eg. ./tables/RK8.txt
    ord: order of integrator eg. 5
    output: filename of output csv data file
    */
    clock_t start, end;

    start = clock();
    solution* result = RKvec_adaptive(combined_perturbations, tSpan, y0, 6, ord, integrator);
    end = clock();
    double cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    printf("Integration completed in %f seconds. Steps taken: %d", cpu_time_used, result->n);

    FILE *fptr = fopen(output, "w");

    for (int i = 0; i < result->n; i++) {
        fprintf(fptr, "%.15f, %.15f, %.15f, %.15f, %.15f, %.15f, %.15f\n", result->t[i], result->y[i][0], result->y[i][1], result->y[i][2], result->y[i][3], result->y[i][4], result->y[i][5]);
        free(result->y[i]);
    }
    fclose(fptr);

    free(result->t);
    free(result);
}
