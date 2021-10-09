// Gravity model used by Supernova

#ifndef GRAVITY_H
#define GRAVITY_H

// number of gravity terms to take (2-21)
#define JGM3SIZE 21

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "constants.h"
#include "vecmath.h"

typedef struct HarmonicMatrices {
    double (*V)[JGM3SIZE+1];
    double (*W)[JGM3SIZE+1];
} harmonicmatrices;

// for propagation in ECEF coordinates
harmonicmatrices* VW(double r[3]) {
    harmonicmatrices* mats = (harmonicmatrices*) malloc(sizeof(harmonicmatrices));

    // Variables to make it easier
    double x = r[0];
    double y = r[1];
    double z = r[2];
    double r2 = x*x + y*y + z*z; // square of magnitude of position

    // Initialize arrays
    double (*V)[JGM3SIZE+1] = (double (*)[JGM3SIZE+1]) malloc((JGM3SIZE+1) * sizeof(double[JGM3SIZE+1]));
    double (*W)[JGM3SIZE+1] = (double (*)[JGM3SIZE+1]) malloc((JGM3SIZE+1) * sizeof(double[JGM3SIZE+1]));

    W[0][0] = 0;
    V[0][0] = R_e/sqrt(r2);

    // Step 1: fill in Vertical column 0
    for (int n = 1; n < JGM3SIZE+1; n++) {
        // eq. 3.30
        V[n][0] = (2.0*n-1)/((double)(n)) * z * R_e/r2 * V[n-1][0] - ((double)(n-1))/((double)(n)) * R_e * R_e/r2 * V[n-2][0];
        W[n][0] = 0;
        //printf("Processed index: %d 0 %f\n", n, V[n][0]);
    }

    // Step 2: Fill in rest of the diagonals
    for (int m = 1; m < JGM3SIZE+1; m++) {
        // eq. 3.29
        V[m][m] = (2.0*m - 1) * (x * R_e / r2 * V[m-1][m-1] - y * R_e / r2 * W[m-1][m-1]);
        W[m][m] = (2.0*m - 1) * (x * R_e / r2 * W[m-1][m-1] + y * R_e / r2 * V[m-1][m-1]);
        //printf("Processed index: %d %d %f\n", m, m, V[m][m]);

        // Step 3: Fill in columns under diagonal
        for (int n = m+1; n < JGM3SIZE+1; n++) {
            // eq. 3.30
            V[n][m] = (2.0 *n-1)/((double)(n-m)) * z * R_e / r2 * V[n-1][m] - ((double)(n+m-1))/((double)(n-m)) * R_e * R_e / r2 * V[n-2][m];
            W[n][m] = (2.0 *n-1)/((double)(n-m)) * z * R_e / r2 * W[n-1][m] - ((double)(n+m-1))/((double)(n-m)) * R_e * R_e / r2 * W[n-2][m];
            //printf("Processed index: %d %d %f\n", n, m, V[n][m]);
        }
    }

    mats->V = V;
    mats->W = W;
    return mats;
    
}

void JGM_gravity(double t, double r[3], double accel[3]) {
    // Gravitational acceleration vector for spacecraft in ECEF coordinates
    // accelerations vector will be modified directly

    harmonicmatrices* h = VW(r);

    double (*V)[JGM3SIZE+1] = h->V;
    double (*W)[JGM3SIZE+1] = h->W;

    // C_n,m = CS(n,m);
    // S_n,m = CS(m-1,n);

    // m = 0 accelerations
    for (int n = 0; n < JGM3SIZE; n++) {
        accel[0] += GMdivREsq * (-CS[n][0] * V[n+1][1]);
        accel[1] += GMdivREsq * (-CS[n][0] * W[n+1][1]);
        accel[2] += GMdivREsq * ((n+1) * (-CS[n][0] * V[n+1][0]));
        // other term with W[n+1][0] removed since these are equally zero
    }
    // printf("Base: ");
    // printMagVec(accel);


    // m > 0 accelerations
    for (int m = 1; m < JGM3SIZE; m++) {
        for (int n = m; n < JGM3SIZE; n++) {
            accel[0] += GMdivREsq * 0.5 * ((-CS[n][m] * V[n+1][m+1] - CS[m-1][n] * W[n+1][m+1]) + factorials[n][m] * (CS[n][m] * V[n+1][m-1] + CS[m-1][n] * W[n+1][m-1]));
            accel[1] += GMdivREsq * 0.5 * ((-CS[n][m] * W[n+1][m+1] + CS[m-1][n] * V[n+1][m+1]) + factorials[n][m] * (-CS[n][m] * W[n+1][m-1] + CS[m-1][n] * V[n+1][m-1]));

            accel[2] += GMdivREsq * ((n-m+1) * (-CS[n][m] * V[n+1][m] - CS[m-1][n] * W[n+1][m]));
        }
    }
    // printf("Extra: ");
    // printMagVec(accel);

    free(V);
    free(W);
    free(h);
}

void J2_gravity(double t, double r[6], double accel[3]) {
    // accelerations vector will be modified directly

    double mag_r = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);  // distance of spacecraft from centre of Earth
    double factor = (3.0/2.0) * GM * R_e * R_e * J2 / pow(mag_r, 5); // J2 Factor

    // 1. Earth Gravity
    for (int i = 0; i < 3; i++) accel[i] = -GM * r[i] / (mag_r * mag_r * mag_r); // v' = -GM r / |r|^3

    // 2. J2 effect
    // J2 Perturbation, from Curtis 12.30
    accel[0] += factor * r[0] * (5 * r[2] * r[2] / (mag_r * mag_r) - 1);
    accel[1] += factor * r[1] * (5 * r[2] * r[2] / (mag_r * mag_r) - 1);
    accel[2] += factor * r[2] * (5 * r[2] * r[2] / (mag_r * mag_r) - 3);

    printVec(accel);

}
#endif
