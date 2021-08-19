#ifndef ECEFsupergravity_H
#define ECEFsupergravity_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define GM 3.986004405e14
#define R_e 6.378137e6

// for propagation in ECEF coordinates
double*** VW(int size, double u[]) {
    // Zonal terms

    double r2 = u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
    double x = u[0];
    double y = u[1];
    double z = u[2];

    double** v = (double**) malloc(size * sizeof(double*));
    double** w = (double**) malloc(size * sizeof(double*));
    
    for (int i = 0; i < size; i++) {
        v[i] = (double*) malloc(size * sizeof(double));
        w[i] = (double*) malloc(size * sizeof(double));
    }
    w[0][0] = 0;
    v[0][0] = R_e/sqrt(r2);

    for (int m = 1; m < size; m++) {
        // eq. 3.29
        v[m][m] = (2.0*m - 1) * (x * R_e / r2 * v[m-1][m-1] - y * R_e / r2 * w[m-1][m-1]);
        w[m][m] = (2.0*m - 1) * (x * R_e / r2 * w[m-1][m-1] + y * R_e / r2 * v[m-1][m-1]);
    }
    
    for (int m = 1; m < size; m++) {

        v[m+1][m] = (2.0 *(m+1)-1) * z * R_e / r2 * v[m][m];
        w[m+1][m] = (2.0 *(m+1)-1) * z * R_e / r2 * w[m][m];

        for (int n = m+2; n < size; n++) {
            v[n][m] = (2.0 *n-1)/((double)(n-m)) * z * R_e / r2 * v[n-1][m] - ((double)(n+m-1))/((double)(n-m) * R_e * R_e / r2 * v[n-2][m]);
            w[n][m] = (2.0 *n-1)/((double)(n-m)) * z * R_e / r2 * w[n-1][m] - ((double)(n+m-1))/((double)(n-m) * R_e * R_e / r2 * w[n-2][m]);
        }
    }
    
}
#endif