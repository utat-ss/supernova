#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ECEFsupergravity.h"

int main(void) {
    
    double u[6] = {5000e3, 0, 1000e3, 0, 0, 0};
    JGM_gravity(1.0, u, NULL);

}