#include <math.h>

const double G = 6.67e-11;
const double M_e = 5.98e24;

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
        output[i+3] = -G*M_e * u[i] / mag_r3;
    }
}