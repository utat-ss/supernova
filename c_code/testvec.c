#include <stdio.h>
#include <stdlib.h>
#include "vecmath.h"
#include <math.h>

int main() {
    // Validation from Mathworks (MATLAB documentation)

    // 1. CROSS PRODUCT
    double A1[3] = {4, -2, 1};
    double B1[3] = {1, -1, 3};
    double C1[3];
    cross(A1, B1, C1);
    printVec(C1);
    // Result: [-5, -11, -2]

    // 2. DOT PRODUCT
    printf("%f\n", dot(A1, C1));
    // Result: 0
    printf("%f\n", dot(B1, C1));
    // Result: 0

    double A2[3] = {4, -1, 2};
    double B2[3] = {2, -2, -1};
    printf("%f\n", dot(A2, B2));
    // Result: 8

    // 3. Matrices
    double A3[3][3] = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    double B3[3][3] = {{0, 1, 0}, {1, 0, 0}, {0, 0, 1}};
    double (*C3)[3] = malloc(3 * sizeof(double[3]));

    matXmat(A3, B3, C3);
    printMatrix(C3);

    double (*q)[3] = QpX(0, 97.49588373 * PI/180.0, 227.05566156 * PI/180.0);
    printMatrix(q);

    // 4. Matrix times Vec
    double A4[3][3] = {{-0.681288, -0.095495, -0.725760}, {-0.732016, 0.088877, 0.675465}, {0.000000, 0.991454, -0.130455}};
    double B4[3] = {6877999.998558, 0.000000, 0.000000};
    double C4[3];
    matXvec(A4, B4, C4);
    printVec(C4);

    // 5. Orbit alteration
    double* Perifocal_pos_vel = PFE(6903e3, 0.003621614, 97.49588373 * PI/180.0, 0, 0);
    // Perifocal coordinates

    double ECI_position[3];
    double ECI_vel[3];
    matXvec(q, Perifocal_pos_vel, ECI_position);
    matXvec(q, Perifocal_pos_vel+3, ECI_vel);

    printVec(ECI_position);
    printVec(ECI_vel);

}