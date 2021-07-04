#ifndef SOLVERS_H
#define SOLVERS_H


double* euler(double (*f)(double, double), double* t, double y0, int n);
double* RK4(double (*f)(double, double), double* t, double y0, int n);

#endif