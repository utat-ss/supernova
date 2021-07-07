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

    printf("\n");

    for (int i = 0; i <= n; i++) {
        printf("%f ", b[i]);

        free(a[i]);
    }

    free(a);
    free(b);
    free(c);
}

void showbutcher_adaptive(char integrator[]) {
    double** a;
    double *b1, *b2, *c;
    int n;

    butcher_adaptive(&a, &b1, &b2, &c, integrator, &n);
    for (int i = 0; i < n; i++) {
        printf("%f: ", c[i]);
        for (int j = 0; j < i; j++) {
            printf("%f ", a[i][j]);
        }
        printf("\n");
    }

    printf("\n");

    for (int i = 0; i < n; i++) {
        printf("%f ", b1[i]);
        free(a[i]);
    }

    printf("\n");

    for (int i = 0; i < n; i++) {
        printf("%f ", b2[i]);
    }

    free(a);
    free(b1);
    free(c);
}

int main(void) {
    double y0[6] = {(6378 + 500)*1e3, 0, 0, 0, 8.5e3, 0};
    double tSpan[2] = {0, 50};

    showbutcher_adaptive("./tables/RKF78.txt");

    // orbit(tSpan, y0, "./tables/RKF56.txt", 6, "RKF56.csv");
    orbit(tSpan, y0, "./tables/RKF78.txt", 7, "RKF78.csv", 1e-2);

	return 0;
}