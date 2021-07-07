#ifndef SOLVERS_H
#define SOLVERS_H

typedef struct AdaptiveSolution {
    double* t; // timesteps
    double** y; // function value
    int n; // number of steps taken
} solution;

double* euler(double (*f)(double, double), double* t, double y0, int n);
double* RK4(double (*f)(double, double), double* t, double y0, int n);
solution* RKvec_adaptive(void (*f)(double, double[], int, double*), double tSpan[], double y0[], int m, int ord, char fname[], double ATOL);
double** RKvec(void (*f)(double, double[], int, double*), double* t, double y0[], int m, int n, char fname[]);

#endif