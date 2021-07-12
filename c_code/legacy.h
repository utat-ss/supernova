#ifndef LEGACY_H
#define LEGACY_H

typedef struct AdaptiveSolution {
    double* t; // timesteps
    double** y; // function value
    int n; // number of steps taken
} solution;
void butcher(double*** a, double** b, double** c, char filename[], int* n);
void butcher_adaptive(double*** a, double** b, double** c1, double** c2, char filename[], int* n);
double* euler(double (*f)(double, double), double* t, double y0, int n);
double* RK4(double (*f)(double, double), double* t, double y0, int n);
solution* RKvec_adaptive(void (*f)(double, double[], int, double*), double tSpan[], double y0[], int m, int ord, char fname[], double ATOL);
double** RKvec(void (*f)(double, double[], int, double*), double* t, double y0[], int m, int n, char fname[]);
void orbitRK(double tSpan[], double y0[], char integrator[], int ord, char output[], double ATOL);

#endif