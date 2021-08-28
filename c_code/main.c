#include <stdio.h>
#include <stdlib.h>
#include "orbit.h"

int main(void) {
    double y0[6] = {-4749231.102296294, -4975106.82469687, 0.0, -719.6538503589323, 686.980716081442, 7561.282735263496};
    double tSpan[2] = {0, 86400*10};
    double tol = 1e-12;
    orbit("RK810", "advanced", tSpan, y0, "RK810.csv", tol);
    orbit("RK1012", "advanced", tSpan, y0, "RK1012.csv", tol);
    
	return 0;
}