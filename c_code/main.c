#include <stdio.h>
#include <stdlib.h>
#include "orbit.h"
#include "butcher.h"

double standard(double t, double y) {
    // Standard testing a DE solver
    return y;
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
    double y0[6] = {(6378 + 500)*1e3, 0, 0, 0, 8.5e3, 0};

    char RK4[] = "./tables/RK4.txt";
    char RK8[] = "./tables/RK8.txt";
    char euler[] = "./tables/euler.txt";
    orbit(1200, y0, euler, "euler.csv");
    orbit(1200, y0, RK4, "RK4.csv");
    orbit(1200, y0, RK8, "RK8.csv");

	return 0;
}