#ifndef ECEFsupergraVity_H
#define ECEFsupergraVity_H

// size of graVity model
#define JGM3SIZE 21

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "constants.h"
#include "vecmath.h"
#include "JGM3Constants.h"

typedef struct HarmonicMatrices {
    double (*V)[JGM3SIZE+1];
    double (*W)[JGM3SIZE+1];
} harmonicmatrices;

// for propagation in ECEF coordinates
harmonicmatrices* VW(double u[]) {

    harmonicmatrices* mats = (harmonicmatrices*) malloc(sizeof(harmonicmatrices));

    // Variables to make it easier
    double x = u[0];
    double y = u[1];
    double z = u[2];
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

void JGM_gravity(double t, double u[6], double output[6]) {
    harmonicmatrices* h = VW(u);

    double (*V)[JGM3SIZE+1] = h->V;
    double (*W)[JGM3SIZE+1] = h->W;

    double accel[3]; // accelerations

    // C_n,m = CS(n,m);
    // S_n,m = CS(m-1,n);

    // m = 0 accelerations
    for (int n = 0; n < JGM3SIZE; n++) {
        accel[0] += GM / (R_e * R_e) * (-CS[n][0] * V[n+1][1]);
        accel[1] += GM / (R_e * R_e) * (-CS[n][0] * W[n+1][1]);

        accel[2] += GM / (R_e * R_e) * ((n+1) * (-CS[n][0] * V[n+1][0]));
        // other term with W[n+1][0] removed since these are equally zero
    }

    // m > 0 accelerations
    for (int m = 1; m < JGM3SIZE; m++) {
        for (int n = m; n < JGM3SIZE; n++) {
            accel[0] += GM / (R_e * R_e) * 0.5 * ((-CS[n][m] * V[n+1][m+1] - CS[m-1][n] * W[n+1][m+1]) + factorials[n][m] * (CS[n][m] * V[n+1][m-1] + CS[m-1][n] * W[n+1][m-1]));
            accel[1] += GM / (R_e * R_e) * 0.5 * ((-CS[n][m] * W[n+1][m+1] + CS[m-1][n] * V[n+1][m+1]) + factorials[n][m] * (-CS[n][m] * W[n+1][m-1] + CS[m-1][n] * V[n+1][m-1]));

            accel[2] += GM / (R_e * R_e) * ((n-m+1) * (-CS[n][m] * V[n+1][m] - CS[m-1][n] * W[n+1][m]));
        }
    }

    printVec(accel);

}
#endif
