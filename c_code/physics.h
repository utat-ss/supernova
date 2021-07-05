#ifndef PHYSICS_H
#define PHYSICS_H

void earth_gravity(double t, double u[], int m, double* output);
void J2_accel(double t, double u[], int m, double* output);
void combined_perturbations(double t, double u[], int m, double* output);

#endif