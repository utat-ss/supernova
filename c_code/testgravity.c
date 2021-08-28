#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gravity.h"

int main(void) {
    
    double u[3] = {5803949.216150, 3690671.686340, 0.000000};
    double accel[3];
    JGM_gravity(1.0, u, accel);
    J2_gravity(1.0, u, accel);

}