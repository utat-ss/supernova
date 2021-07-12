#include <stdio.h>
#include <stdlib.h>
#include "orbit.h"
#include "legacy.h"

int main(void) {
    double y0[6] = {(6378 + 500)*1e3, 0, 0, 0, 7.5e3, 0};
    double tSpan[2] = {0, 1000};

    orbitRK(tSpan, y0, "./tables/RK45.txt", 5, "RK45.csv", 1e-5);
    orbitRK(tSpan, y0, "./tables/RKF56.txt", 6, "RKF56.csv", 1e-5);
    orbit(tSpan, y0, "RK810.csv", 1e-12);

	return 0;
}