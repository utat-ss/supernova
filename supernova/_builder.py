from cffi import FFI

ffibuilder = FFI()

# two public functions: orbit and stateFromKeplerian
ffibuilder.cdef(
    """
typedef struct AdaptiveSolution {
    double* t; // timesteps
    double (*y)[6]; // function value
    int n; // number of steps taken
} solution;

solution* orbitPTR(char solver[], char model[], double tSpan[], double y0[], double ATOL);
double* stateFromKeplerian(double sma, double ecc, double inc, double raan, double aop, double m);
void free(void* ptr); // free memory allocated by C code
"""
)


ffibuilder.set_source(
    "supernova.backend",
    """
                      #include "constants.h"
                      #include "forces.h"
                      #include "gravity.h"
                      #include "numericaldiff.h"
                      #include "orbit.h"
                      #include "solvers.h"
                      #include "vecmath.h"
                      """,
    include_dirs=["./c_code"],
    extra_compile_args=["-std=c99", "-O3", "-march=native", "-ffast-math", "-lm"],
)

if __name__ == "__main__":
    ffibuilder.compile(verbose=False)
