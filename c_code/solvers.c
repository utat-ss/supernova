#include <stdlib.h>

double* euler(double (*f)(double, double), double* t, double y0, int n) {
    // Integrates from one side of the interval to another using Euler's method
    // f: function for which y' = f(t, y)
    // t: evaluation points
    // y0: initial condition
    // n: number of steps

    double* result = calloc(n, sizeof(double));

    if (n <= 1) return NULL; // error case

    double dt;
    double deriv;

    result[0] = y0;

    for (int step = 1; step < n; step++) {
        // Euler Step
        // // y_n = y_n-1 + f(t-1, y_n-1) * dt

        dt = t[step] - t[step-1];        
        deriv = f(t[step-1], result[step-1]);

        result[step] = result[step-1] + deriv * dt;
    }

    return result;
}

double* RK4(double (*f)(double, double), double* t, double y0, int n) {
    // Integrates from one side of the interval to another using RK4 method
    // f: function for which y' = f(t, y)
    // t: evaluation points
    // y0: initial condition
    // n: number of steps

    double* result = calloc(n, sizeof(double));

    if (n <= 1) return NULL; // error case

    double dt, k1, k2, k3, k4;

    result[0] = y0;

    for (int step = 0; step < n-1; step++) {
        // RK Step

        dt = t[step+1] - t[step];
        k1 = f(t[step], result[step]);
        k2 = f(t[step] + dt/2, result[step] + dt*k1/2);
        k3 = f(t[step] + dt/2, result[step] + dt*k2/2);
        k4 = f(t[step] + dt, result[step] + dt*k3);

        result[step+1] = result[step] + (k1 + 2*k2 + 2*k3 + k4)*dt/6;
    }

    return result;

}