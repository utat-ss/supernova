// High level interface functions

#ifndef ORBIT_H
#define ORBIT_H

#include "solvers.h"
#include "forces.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>

void orbit(char solver[], char model[], double tSpan[], double y0[], char output[], double ATOL) {
    /*
    Propagates orbit using solver of choice.
    solver: name, either "RK810" or "RK1012"
    model: perturbation model, either "simplified" or "advanced"
    tSpan: t0, tf of times in seconds.
    y0: initial x y z vx vy vz state vector
    output: filename of output csv data file
    ATOL: tolerance
    */
    clock_t start, end;
    solution* result;

    start = clock();

    void (*f)(double, double *, double *);

    if (!strcmp(model, "simplified")) f = simplified_perturbations;
    else if (!strcmp(model, "advanced")) f = advanced_perturbations;
    else {
        printf("Invalid model chosen. Pick either 'simplified' or 'advanced'\n");
        return;
    }

    if (!strcmp(solver, "RK810")) result = RK810vec(f, tSpan, y0, ATOL);
    else if (!strcmp(solver, "RK1012")) result = RK1012vec(f, tSpan, y0, ATOL);
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

solution* orbitPTR(char solver[], char model[], double tSpan[], double y0[], double ATOL) {
    /*
    Propagates orbit using solver of choice and returns pointer to solution
    solver: name, either "RK810" or "RK1012"
    model: perturbation model, either "simplified" or "advanced"
    tSpan: t0, tf of times.
    y0: initial x y z vx vy vz state vector
    ATOL: tolerance
    */
    void (*f)(double, double *, double *);

    if (!strcmp(model, "simplified")) f = simplified_perturbations;
    else if (!strcmp(model, "advanced")) f = advanced_perturbations;
    else {
        printf("Invalid model chosen. Pick either 'simplified' or 'advanced'\n");
        return NULL;
    }

    if (!strcmp(solver, "RK810")) return RK810vec(f, tSpan, y0, ATOL);
    else if (!strcmp(solver, "RK1012")) return RK1012vec(f, tSpan, y0, ATOL);
    else {
        printf("Invalid solver chosen. Pick either 'RK810' or 'RK1012'\n");
        return NULL;
    }
}

double* nextStatefromCurrent(double x, double y, double z, double vx, double vy, double vz, double h) {
    // Advance forward in your orbit by h seconds.
    double* state =  malloc(6 * sizeof(double));
    double y0[6] = {x, y, z, vx, vy, vz};
    RK1012Step(simplified_perturbations, h, y0, state);
    return state;
}

double* stateFromKeplerian(double sma, double ecc, double inc, double raan, double aop, double m) {
    // Initial ECI position/velocity using Keplerian elements

    double (*q)[3] = QpX(aop, inc, raan); // Perifocal to ECI matrix
    double* Perifocal_pos_vel = PFE(sma, ecc, inc, aop, m);
    // Perifocal coordinates

    double* y0 = malloc(6 * sizeof(double));

    double ECI_position[3];
    double ECI_vel[3];
    matXvec(q, Perifocal_pos_vel, ECI_position);
    matXvec(q, Perifocal_pos_vel+3, ECI_vel);

    for (int i = 0; i < 3; i++) {
        y0[i] = ECI_position[i];
        y0[i+3] = ECI_vel[i];
    }
    return y0;
}

void freemem(void* ptr) {
    free(ptr);
}

#endif