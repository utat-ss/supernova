#ifndef PHYSICS_H
#define PHYSICS_H

// Aero configuration (static in code for now)
#define CD 2.2
#define Am 0.01

// Model configurations
#define J2_ON 1
#define AERO_ON 1

#include <math.h>
#include "constants.h"
#include "vecmath.h"

// Physics models
// operates on 6D state vector

double atmosphere(double z) {
    // Atmospheric density based on Curtis D.41
    double h[] = {0, 25e3, 30e3, 40e3, 50e3, 60e3, 70e3, 80e3, 90e3, 100e3, 110e3, 120e3, 130e3, 140e3, 150e3, 180e3, 200e3, 250e3, 300e3, 350e3, 400e3, 450e3, 500e3, 600e3, 700e3, 800e3, 900e3, 1000};
    // geometric heights in m
    
    double r[] = {1.225, 4.008e-2, 1.841e-2, 3.996e-3, 1.027e-3, 3.097e-4, 8.283e-5, 1.846e-5, 3.416e-6, 5.606e-7, 9.708e-8, 2.222e-8, 8.152e-9, 3.831e-9, 2.076e-9, 5.194e-10, 2.541e-10, 6.073e-11, 1.916e-11, 7.014e-12, 2.803e-12, 1.184e-12, 5.215e-13, 1.137e-13, 3.070e-14, 1.136e-14, 5.759e-15, 3.561e-15};
    // densities in kg/m^3
    
    double H[] = {7.310e3, 6.427e3, 6.546e3, 7.360e3, 8.342e3, 7.583e3, 6.661e3, 5.927e3, 5.533e3, 5.703e3, 6.782e3, 9.973e3, 13.243e3, 16.322e3, 21.652e3, 27.974e3, 34.934e3, 43.342e3, 49.755e3, 54.513e3, 58.019e3, 60.980e3, 65.654e3, 76.377e3, 100.587e3, 147.203e3, 208.020e3, 0};
    // scale heights in m

    if (z < 0) z = 0;
    else if (z > 1000e3) z = 1000e3;

    int i;
    for (i = 0; i < 27; i++) {
        if (z >= h[i] && z < h[i+1]) break; // we found the correct altitude
    }

    return r[i] * exp(-(z-h[i])/H[i]);
}

void combined_perturbations(double t, double u[6], double output[6]) {
    /*
    force model on spacecraft
    t: time, not used
    u: vector of size 6 containing x y z vx vy vz
    output: vector of size m which stores x' y' z' vx' vy' vz'
    */
    double mag_r = sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);  // distance of spacecraft from centre of Earth
    double factor = (3.0/2.0) * GM * R_e * R_e * J2 / pow(mag_r, 5); // J2 Factor

    // 0. Velocity update
    for (int i = 0; i < 3; i++) output[i] = u[i+3];

    // 1. Earth Gravity
    for (int i = 0; i < 3; i++) output[i+3] = -GM * u[i] / (mag_r * mag_r * mag_r); // v' = -GM r / |r|^3

    #if J2_ON == 1
    // 2. J2 effect
    // J2 Perturbation, from Curtis 12.30
    output[3] += factor * u[0] * (5 * u[2] * u[2] / (mag_r * mag_r) - 1);
    output[4] += factor * u[1] * (5 * u[2] * u[2] / (mag_r * mag_r) - 1);
    output[5] += factor * u[2] * (5 * u[2] * u[2] / (mag_r * mag_r) - 3);
    #endif

    
    #if AERO_ON == 1
    // 3. Atmospheric Drag
    // from Curtis 12.12

    double v_rel[3] = {u[3]+u[1]*W_e, u[4]-W_e*u[0], u[5]};
    // determine atmospheric velocity and subsequent relative velocity by subtracting atmospheric velocity

    //printVec(v_rel);

    double mag_v_rel = sqrt(v_rel[0]*v_rel[0] + v_rel[1]*v_rel[1] + v_rel[2]*v_rel[2]);  // magnitude of relative velocity
    
    double rho = atmosphere(mag_r-R_e); // atmospheric density at altitude

    for (int i = 3; i < 6; i++) {
        output[i] -= 0.5 * rho * CD * Am * mag_v_rel * v_rel[i-3];
        //printf("Acceleration from drag in direction %d: %.10f m/s^2\t", i-3, -0.5 * rho * CD * Am * mag_v_rel * v_rel[i-3]);
    }
    #endif
    
    return;
}

#endif