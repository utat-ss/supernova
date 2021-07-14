#ifndef ORBIT_H
#define ORBIT_H
#include "solvers.h"

void orbit(double tSpan[], double y0[], char output[], double ATOL);
void orbitgeneral(solution* (*solver)(void (*)(double, double[], double*), double[], double[], double), double tSpan[], double y0[], char output[], double ATOL);

#endif