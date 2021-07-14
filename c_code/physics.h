#ifndef PHYSICS_H
#define PHYSICS_H

void earth_gravity(double t, double u[], double* output);
void J2_accel(double t, double u[], double* output);
void combined_perturbations(double t, double u[], double* output);

#endif