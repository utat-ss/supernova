#include <stdio.h>
#include <stdlib.h>
#include "solvers.h"
#include "physics.h"

void orbit(int n, double y0[], char integrator[], char output[]) {
    /*
    n: number of time steps (10 seconds each)
    y0: initial x y z vx vy vz state vector
    integrator: path to a butcher tableau eg. ./tables/RK8.txt
    output: filename of output csv data file
    */
    double* t = malloc(n * sizeof(double));

    for (int i = 0; i < n; i++) t[i] = i*10;

    double** result = RKvec(combined_perturbations, t, y0, 6, n, integrator);

    FILE *fptr = fopen(output, "w");

    for (int i = 0; i < n; i++) {
        fprintf(fptr, "%f, %f, %f, %f, %f, %f, %f\n", t[i], result[i][0], result[i][1], result[i][2], result[i][3], result[i][4], result[i][5]);
        free(result[i]);
    }
    fclose(fptr);

    free(t);
    free(result);
}
