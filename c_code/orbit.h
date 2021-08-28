#ifndef ORBIT_H
#define ORBIT_H
#include "solvers.h"
#include "physics.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>

void orbit(char solver[], double tSpan[], double y0[], char output[], double ATOL) {
    /*
    Propagates orbit using solver of choice.
    solver: name, either "RK810" or "RK1012"
    tSpan: t0, tf of times.
    y0: initial x y z vx vy vz state vector
    output: filename of output csv data file
    ATOL: tolerance
    */
    clock_t start, end;
    solution* result;

    start = clock();
    if (!strcmp(solver, "RK810")) result = RK810vec(combined_perturbations, tSpan, y0, ATOL);
    else if (!strcmp(solver, "RK1012")) result = RK1012vec(combined_perturbations, tSpan, y0, ATOL);
    else {
        printf("Invalid solver chosen.\n");
        return;
    }
    end = clock();
    double cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    printf("Integration using %s solver completed in %.8f seconds. Steps taken: %d\n", solver, cpu_time_used, result->n);

    FILE *fptr = fopen(output, "w");

    for (int i = 0; i < result->n; i++) {
        fprintf(fptr, "%.15f, %.15f, %.15f, %.15f, %.15f, %.15f, %.15f\n", result->t[i], result->y[i][0], result->y[i][1], result->y[i][2], result->y[i][3], result->y[i][4], result->y[i][5]);
    }
    fclose(fptr);

    free(result->t);
    free(result);
}

solution* orbitPTR(char solver[], double tSpan[], double y0[], double ATOL) {
    /*
    Propagates orbit using solver of choice and returns pointer to solution
    solver: name, either "RK810" or "RK1012"
    tSpan: t0, tf of times.
    y0: initial x y z vx vy vz state vector
    ATOL: tolerance
    */
    if (!strcmp(solver, "RK810")) return RK810vec(combined_perturbations, tSpan, y0, ATOL);
    else if (!strcmp(solver, "RK1012")) return RK1012vec(combined_perturbations, tSpan, y0, ATOL);
    else {
        printf("Invalid solver chosen.\n");
        return NULL;
    }
}

#endif