#ifndef SOLVERS_H
#define SOLVERS_H

typedef struct AdaptiveSolution {
    double* t; // timesteps
    double** y; // function value
    int n; // number of steps taken
} solution;

solution* RK810vec(void (*f)(double, double[], int, double*), double tSpan[], double y0[], int m, double ATOL);

#endif