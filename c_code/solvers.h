#ifndef SOLVERS_H
#define SOLVERS_H

double* euler(double (*f)(double, double), double* t, double y0, int n);
double* RK4(double (*f)(double, double), double* t, double y0, int n);
double** RK4vec(void (*f)(double, double[], int, double*), double* t, double y0[], int m, int n);
double** RKvec(void (*f)(double, double[], int, double*), double* t, double y0[], int m, int n, char fname[]);

#endif