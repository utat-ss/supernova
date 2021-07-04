#include <stdio.h>
#include <stdlib.h>
#include "solvers.h"
#include "physics.h"
#include "butcher.h"

double standard(double t, double y) {
    // Standard testing a DE solver
    return y;
}

void orbit(int n, char integrator[], char output[]) {
    // orbit propagator test
    double* t = calloc(n, sizeof(double));

    for (int i = 0; i < n; i++) {
        t[i] = i*10;
    }

    double y0[6] = {(6378 + 500)*1e3, 0, 0, 0, 8.5e3, 0};

    double** result = RKvec(earth_gravity, t, y0, 6, n, integrator);

    FILE *fptr = fopen(output, "w");

    for (int i = 0; i < n; i++) {
        fprintf(fptr, "%f, %f, %f, %f, %f, %f, %f\n", t[i], result[i][0], result[i][1], result[i][2], result[i][3], result[i][4], result[i][5]);
    }
}

void showbutcher(char integrator[]) {
    double** a;
    double *b, *c;
    int n;

    butcher(&a, &b, &c, integrator, &n);
    for (int i = 0; i <= n; i++) {
        printf("%f ", c[i]);
        for (int j = 0; j < i; j++) {
            printf("%f ", a[i][j]);
        }
        printf("\n");
    }

    for (int i = 0; i <= n; i++) {
        printf("%f ", b[i]);

        free(a[i]);
    }

    free(a);
    free(b);
    free(c);
}

int main(void) {
    char RK4[] = "./tables/RK4.txt";
    char RK8[] = "./tables/RK8.txt";
    char euler[] = "./tables/euler.txt";
    orbit(1200, euler, "euler.csv");
    orbit(1200, RK4, "RK4.csv");
    orbit(1200, RK8, "RK8.csv");

	return 0;
}