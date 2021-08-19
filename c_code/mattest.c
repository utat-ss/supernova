#include <stdio.h>
#include <stdlib.h>
#include "vecmath.h"
#include <math.h>
#define PI 	3.14159265358979323846

int main() {
    double** q = QMatrix(60.0 * PI/180.0, 30.0 * PI/180.0, 40.0 * PI/180.0);

    matT(q);

    for (int i = 0; i < 3; i++) {
        for (int j = 0 ; j < 3; j++) {
            printf(" %f ", q[i][j]);
        }
        printf("\n");
    }
    double pos[3] = {6285, 3628.6, 0};
    double reci[3];
    vecmatmul(q, pos, reci);

    for (int i = 0; i < 3; i++) printf(" %f ", reci[i]);

}