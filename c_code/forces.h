// Force model for Supernova

#ifndef FORCES_H
#define FORCES_H

// Aero configuration (static in code for now)
#define CD 2.2
#define Am 0.01

// Model configurations
#define J2_ON 1
#define AERO_ON 1

#include <math.h>
#include <stdio.h>
#include "constants.h"
#include "vecmath.h"
#include "gravity.h"

// Physics models
// operates on 6D state vector

double atmosphere(double z) {
    // Atmospheric density based on Curtis D.41
    // Model: US Standard 76

    if (z < 0) z = 0;
    else if (z > 1000e3) z = 1000e3;

    int i;
    for (i = 0; i < 27; i++) {
        if (z >= ATMh[i] && z < ATMh[i+1]) break; // we found the correct altitude
    }

    return ATMr[i] * exp(-(z-ATMh[i])/ATMH[i]);
}

void simplified_perturbations(double t, double u[6], double output[6]) {
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

void advanced_perturbations(double t, double u[6], double output[6]) {
    /*
    advanced force model on spacecraft
    t: mission elapsed time in seconds
    u: vector of size 6 containing x y z vx vy vz
    output: vector of size m which stores x' y' z' vx' vy' vz'
    */

    double mag_r = sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);  // distance of spacecraft from centre of Earth
    double curr_jd = t/86400;
    //printf("%f\n", mag_r - 6378e3);

    // 0. Velocity update
    for (int i = 0; i < 3; i++) output[i] = u[i+3];

    // 1. JGM3 Gravity
    double ECEF[3];
    double Q[3][3];
    double accel[3] = {0, 0, 0}; // store temporary accel in ECEF
    ECI2ECEF(curr_jd, Q); // get conversion matrix ECI -> ECEF => Q
    matXvec(Q, u, ECEF);
    JGM_gravity(curr_jd, ECEF, accel); // determine gravitational acceleration
    transpose(Q); // get conversion matrix ECI <- ECEF => Q^T
    matXvec(Q, accel, output+3); // Accelerations in ECI
    
    
    #if AERO_ON == 1
    // 2. Atmospheric Drag
    // from Curtis 12.12
    double v_rel[3] = {u[3]+u[1]*W_e, u[4]-W_e*u[0], u[5]};
    // determine atmospheric velocity and subsequent relative velocity by subtracting atmospheric velocity

    double mag_v_rel = sqrt(v_rel[0]*v_rel[0] + v_rel[1]*v_rel[1] + v_rel[2]*v_rel[2]);  // magnitude of relative velocity
    double rho = atmosphere(mag_r-R_e); // atmospheric density at altitude
    for (int i = 3; i < 6; i++) output[i] -= 0.5 * rho * CD * Am * mag_v_rel * v_rel[i-3];
    #endif
    
    return;
}

#endif