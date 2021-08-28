#include <stdio.h>
#include <stdlib.h>
#include "orbit.h"

int main(void) {
    double y0[6] = {(6378 + 500)*1e3, 0, 0, 0, 7.6e3, 0};
    double tSpan[2] = {0, 86400*10};
    double tol = 1e-6;
    orbit("RK810", tSpan, y0, "RK810.csv", tol);
    orbit("RK1012", tSpan, y0, "RK1012.csv", tol);
    
	return 0;
}