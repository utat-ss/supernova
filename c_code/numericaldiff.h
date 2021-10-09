#ifndef NUMERICALDIFF_H
#define NUMERICALDIFF_H

#include <math.h>
#include <stdio.h>
#include "constants.h"
#include "vecmath.h"

// Numerical differentiation testing project
void GMgradient(double position[3]) {
    // central difference formula
    double h = 1;
    double output[3];

    for (int i = 0; i < 3; i++) {
        // i: direction
        double r1[3];
        double r2[3];

        for (int j = 0; j < 3; j++) {
            r1[j] = position[j];
            r2[j] = position[j];
        }
        r1[i] -= h;
        r2[i] += h;

        double magr1 = sqrt(r1[0] * r1[0] +r1[1] * r1[1] + r1[2] * r1[2]);
        double magr2 = sqrt(r2[0] * r2[0] +r2[1] * r2[1] + r2[2] * r2[2]);

        output[i] = GM * (1/magr2 - 1/magr1) / (2*h);
    }
    printVec(output);
}
#endif