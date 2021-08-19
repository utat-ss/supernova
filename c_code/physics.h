#ifndef PHYSICS_H
#define PHYSICS_H

#define GM 3.986004405e14
#define R_e 6.378137e6
#define J2 1.08262668e-3
#include <math.h>
#include "vecmath.h"

// Physics models
// operates on 6D state vector
const double earthomega[3] = {0, 0, 72.9211e-6};

void combined_perturbations(double t, double u[6], double output[6]) {
    /*
    force model on spacecraft
    t: time, not used
    u: vector of size 6 containing x y z vx vy vz
    output: vector of size m which stores x' y' z' vx' vy' vz'
    */
    double mag_r = sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
    double factor = (3.0/2.0) * GM * R_e * R_e * J2 / pow(mag_r, 5);
    //double mag_v = sqrt(u[3]*u[3] + u[4]*u[4] + u[5]*u[5]);

    //double atmvel[3];
    //cross(u, earthomega, atmvel);

    // 0. Velocity update
    for (int i = 0; i < 3; i++) output[i] = u[i+3];

    // 1. Earth Gravity
    for (int i = 0; i < 3; i++) output[i+3] = -GM * u[i] / (mag_r * mag_r * mag_r); // v' = -GM r / |r|^3

    // 2. J2 effect
    // J2 Perturbation, from Curtis 12.30
    output[3] += factor * u[0] * (5 * u[2] * u[2] / (mag_r * mag_r) - 1);
    output[4] += factor * u[1] * (5 * u[2] * u[2] / (mag_r * mag_r) - 1);
    output[5] += factor * u[2] * (5 * u[2] * u[2] / (mag_r * mag_r) - 3);

    // 3. Atmospheric Drag
    // double rho = 1;
    // for (int i = 3; i < 6; i++) output[i] -= rho * CD * AM * u[i] * mag_v/2.0;
    return;
}

#endif