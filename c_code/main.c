#include <stdio.h>
#include <stdlib.h>
#include "orbit.h"

int main(void) {
    double y0[6] = {(6378 + 500)*1e3, 0, 0, 0, 7.5e3, 0};
    double tSpan[2] = {0, 100000};
    double tol = 1e-12;
    orbit(tSpan, y0, "RK810.csv", tol);

	return 0;
}