#include <stdio.h>
#include <stdlib.h>
#include "orbit.h"

int main(void) {
    double y0[6] = {-4749231.102296294, -4975106.82469687, 0.0, -719.6538503589323, 686.980716081442, 7561.282735263496};
    double tSpan[2] = {0, 86400*10};
    double tol = 1e-5;
    orbit("RK810", "simplified", tSpan, y0, "simplified.csv", tol);
    orbit("RK810", "advanced", tSpan, y0, "advanced.csv", tol);
    
	return 0;
}