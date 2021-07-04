#include <stdio.h>
#include "solvers.h"

double standard(double t, double y) {
    // Standard testing a DE solver
    return y;
}

int main(void) {
    double t[11] = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};

    double* result = RK4(standard, t, 1, 11);

    for (int i = 0; i < 11; i++) {
        printf("%f ", result[i]);
    }
	return 0;
}