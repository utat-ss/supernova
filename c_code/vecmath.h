#ifndef VECMATH_H
#define VECMATH_H

// vector math library
#include <math.h>

void cross(double a[3], double b[3], double u[3]) {
    // Crosses a, b and stores result in u
    // 231 312
    // 120 201
    u[0] = a[1]*b[2] - b[1]*a[2];
    u[1] = a[2]*b[0] - b[2]*a[0];
    u[2] = a[0]*b[1] - b[0]*a[1];
}

double dot(double a[3], double b[3]) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

double** QMatrix(double omega, double inc, double OMEGA) {
    // Perifocal to ECI matrix
    double** matrix = malloc(3 * sizeof(double*));
    for (int i=0; i<3; i++) matrix[i] = malloc(3 * sizeof(double));

    matrix[0][0] = -sin(OMEGA) * cos(inc) * sin(omega) + cos(OMEGA) * cos(omega);
    matrix[0][1] = cos(OMEGA) * cos(inc) * sin(omega) + sin(OMEGA) * cos(omega);
    matrix[0][2] = sin(inc) * sin(omega);

    matrix[1][0] = -sin(OMEGA) * cos(inc) * cos(omega) - cos(OMEGA) * sin(omega);
    matrix[1][1] = cos(OMEGA) * cos(inc) * cos(omega) - sin(OMEGA) * sin(omega);
    matrix[1][2] = sin(inc) * cos(omega);

    matrix[2][0] = sin(OMEGA) * sin(inc);
    matrix[2][1] = -cos(OMEGA) * sin(inc);
    matrix[2][2] = cos(inc);
    return matrix;
}

void vecmatmul(double** mat, double vec[3], double u[3]) {
    // vector multiply 3d matrix to 3d vector
    for (int i=0; i<3; i++) u[i] = mat[i][0] * vec[0] + mat[i][1] * vec[1] + mat[i][2] * vec[2];
}

void matmatmul(double** A, double** B, double** U) {
    // vector multiply 3d matrix to 3d matrix
    
}

void matT(double** A) {
    // in place matrix transpose
    double tmp;
    for (int i = 0; i < 2; i++) {
        for (int j = i+1; j < 3; j++) {
            tmp = A[i][j];
            A[i][j] = A[j][i];
            A[j][i] = tmp;
        }
    }
}


#endif