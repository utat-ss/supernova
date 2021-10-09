// Vector Math Utilities for Supernova

#ifndef VECMATH_H
#define VECMATH_H

#include <math.h>
#include <stdio.h>
#include "constants.h"

void cross(double A[3], double B[3], double U[3]) {
    // Crosses vectors A and B then stores result in vector U
    U[0] = A[1]*B[2] - B[1]*A[2];
    U[1] = A[2]*B[0] - B[2]*A[0];
    U[2] = A[0]*B[1] - B[0]*A[1];
}

double dot(double A[3], double B[3]) {
    // Returns dot product of Vectors A and B as scalar
    return A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
}

void matXvec(double (*A)[3], double V[3], double U[3]) {
    // Multiplies matrix A by vector V, storing result in vector U
    for (int i=0; i<3; i++) U[i] = A[i][0] * V[0] + A[i][1] * V[1] + A[i][2] * V[2];
}

void matXmat(double (*A)[3], double (*B)[3], double (*C)[3]) {
    // Multiplies matrix A by matrix B, storing result in matrix C
    for (int i=0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            C[i][j] = A[i][0] * B[0][j] + A[i][1] * B[1][j] + A[i][2] * B[2][j];
        }
    }
}

void transpose(double (*A)[3]) {
    // in place matrix transpose of A
    double tmp;
    for (int i = 0; i < 2; i++) {
        for (int j = i+1; j < 3; j++) {
            tmp = A[i][j];
            A[i][j] = A[j][i];
            A[j][i] = tmp;
        }
    }
}

/*
Specialty Matrices
*/

double (*QXp(double aop, double inc, double raan))[3] {
    // ECI to Perifocal matrix
    // aop: AOP
    // raan: RAAN
    // inc: inclination
    // all angles in radians please
    double (*matrix)[3] = malloc(3 * sizeof(double[3]));

    matrix[0][0] = -sin(raan) * cos(inc) * sin(aop) + cos(raan) * cos(aop);
    matrix[0][1] = cos(raan) * cos(inc) * sin(aop) + sin(raan) * cos(aop);
    matrix[0][2] = sin(inc) * sin(aop);

    matrix[1][0] = -sin(raan) * cos(inc) * cos(aop) - cos(raan) * sin(aop);
    matrix[1][1] = cos(raan) * cos(inc) * cos(aop) - sin(raan) * sin(aop);
    matrix[1][2] = sin(inc) * cos(aop);

    matrix[2][0] = sin(raan) * sin(inc);
    matrix[2][1] = -cos(raan) * sin(inc);
    matrix[2][2] = cos(inc);
    return matrix;
}

double (*QpX(double aop, double inc, double raan))[3] {
    // Perifocal to ECI matrix    
    // aop: AOP
    // raan: RAAN
    // inc: inclination
    // all angles in radians please
    double (*matrix)[3] = malloc(3 * sizeof(double[3]));

    matrix[0][0] = -sin(raan) * cos(inc) * sin(aop) + cos(raan) * cos(aop);
    matrix[1][0] = cos(raan) * cos(inc) * sin(aop) + sin(raan) * cos(aop);
    matrix[2][0] = sin(inc) * sin(aop);

    matrix[0][1] = -sin(raan) * cos(inc) * cos(aop) - cos(raan) * sin(aop);
    matrix[1][1] = cos(raan) * cos(inc) * cos(aop) - sin(raan) * sin(aop);
    matrix[2][1] = sin(inc) * cos(aop);

    matrix[0][2] = sin(raan) * sin(inc);
    matrix[1][2] = -cos(raan) * sin(inc);
    matrix[2][2] = cos(inc);
    return matrix;
}

void ECI2ECEF(double jd, double (*matrix)[3]) {
    // ECI to ECEF matrix, given time in Julian

    // Get Earth rotation angle given time using Celest algorithm.
    // https://github.com/JaiWillems/Celest/blob/develop-v0.2.0/celest/satellite/coordinate.py
    double ERA = -fmod(6.30038748702467 * (jd - 2451545) + 4.895, 2*PI);
    //printf("ERA: %f\n", ERA * 180 / PI);

    matrix[0][0] = cos(ERA);
    matrix[0][1] = -sin(ERA);
    matrix[0][2] = 0;

    matrix[1][0] = sin(ERA);
    matrix[1][1] = cos(ERA);
    matrix[1][2] = 0;

    matrix[2][0] = 0;
    matrix[2][1] = 0;
    matrix[2][2] = 1;

}

void ECEF2ECI(double jd, double (*matrix)[3]) {
    // ECEF to ECI matrix, given time in Julian

    double ERA = fmod(6.30038748702467 * (jd - 2451545) + 4.895, 2*PI);

    matrix[0][0] = cos(ERA);
    matrix[0][1] = -sin(ERA);
    matrix[0][2] = 0;

    matrix[1][0] = sin(ERA);
    matrix[1][1] = cos(ERA);
    matrix[1][2] = 0;

    matrix[2][0] = 0;
    matrix[2][1] = 0;
    matrix[2][2] = 1;

}

double* PFE(double a, double e, double i, double aop, double m) {
    // Returns perifocal position and velocity as 6-vector

    double* perifocals = malloc(6 * sizeof(double));

    double E = m + ((e*sin(m))/(1-sin(m + e) + sin(m))); // Eccentric anomaly 
    double theta = 2 * atan(sqrt((1+e)/(1-e))*tan(E/2)); // total anomaly
    double r = a * (1-e*e) / (1+ e*cos(theta)); // semi latus rectum
    double v = sqrt(GM/(a*(1 - e*e))); // velocity magnitude

    // PERIFOCAL R/V
    perifocals[0] = r * cos(theta);
    perifocals[1] = r * sin(theta);
    perifocals[2] = 0;

    perifocals[3] = v * -sin(theta);
    perifocals[4] = v * (e + cos(theta));
    perifocals[5] = 0;

    return perifocals;
}

/*
FORMATTERS
*/
void printVec(double vec[3]) {
    printf("[%f, %f, %f]\n", vec[0], vec[1], vec[2]);
}

void printMatrix(double (*matrix)[3]) {
    for (int i = 0; i < 3; i++) {
        printf("[%f, %f, %f]\n", matrix[0][i], matrix[1][i], matrix[2][i]);
    }
}
void printMagVec(double vec[3]) {
    printf("%f\n", sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]));
}

#endif