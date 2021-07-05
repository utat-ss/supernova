#include <math.h>

const double GM = 3.986004405e14;
const double R_e = 6.378137e6;
const double J2 = 1.08262668e-3;
const double flat = 1/298.257;

void earth_gravity(double t, double u[], int m, double* output) {
    /*
    Uses keplers law to do physics
    t: time, not used
    u: vector of size m containing x y z vx vy vz
    m: size
    output: vector of size m which stores x' y' z' vx' vy' vz'
    */

    double mag_r3 = pow(u[0]*u[0] + u[1]*u[1] + u[2]*u[2], 1.5); // |r|^(3);

    for (int i = 0; i < 3; i++) {
        // r' = v
        output[i] = u[i+3];

        // v' = -GM r / |r|^3
        output[i+3] = -GM * u[i] / mag_r3;
    }
}

void J2_accel(double t, double u[], int m, double* output) {
    // J2 Perturbation, from Curtis 12.30
    double mag_r = sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
    double factor = (3.0/2.0) * J2 * GM * R_e * R_e / pow(mag_r, 5);

    output[3] = factor * u[0] * (5 * u[2] * u[2] / (mag_r * mag_r) - 1);
    output[4] = factor * u[1] * (5 * u[2] * u[2] / (mag_r * mag_r) - 1);
    output[5] = factor * u[2] * (5 * u[2] * u[2] / (mag_r * mag_r) - 3);
}

void combined_perturbations(double t, double u[], int m, double* output) {
    double temp[6]; // temporary variable to store J2 perturbation effects

    // 1. Earth Gravity
    earth_gravity(t, u, m, output);

    // 2. J2 effect
    J2_accel(t, u, m, temp);

    for (int i = 3; i < 6; i++) output[i] += temp[i];
}