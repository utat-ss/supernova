#include <stdlib.h>
#include "butcher.h"

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

double** RK4vec(void (*f)(double, double[], int, double*), double* t, double y0[], int m, int n) {
    /*
    Integrates between two time points using RK4 method using n timesteps on a vector of size m

    t: evaluation points
    y0: array of length m containing initial conditions
    m, n: dimensions of vector, number of timesteps
    f: function for which y' = f(t, y, m, ARR), where m is size, and ARR stores the resulting vector
    */

    if (n <= 1) return NULL; // error case

    // Allocate result array
    double** result = calloc(n, sizeof(double*));
    for (int i = 0; i < n; i++) {
        result[i] = calloc(m, sizeof(double));
    }

    // Step variables
    double h; // stepsize
    double* k1 = calloc(m, sizeof(double));
    double* k2 = calloc(m, sizeof(double));
    double* k2calc = calloc(m, sizeof(double));
    double* k3 = calloc(m, sizeof(double));
    double* k3calc = calloc(m, sizeof(double));
    double* k4 = calloc(m, sizeof(double));
    double* k4calc = calloc(m, sizeof(double));
    // Temporary step variables for vector values
    
    // Init y0
    for (int j = 0; j < m; j++) {
        result[0][j] = y0[j];
    }

    for (int i = 0; i < n-1; i++) {
        // RK Step
        h = t[i+1] - t[i];

        // K1 Step
        f(t[i], result[i], m, k1);

        // K2 Step
        for (int j = 0; j < m; j++) k2calc[j] = result[i][j] + h * k1[j] / 2;
        f(t[i] + h/2, k2calc, m, k2);

        // K3 Step
        for (int j = 0; j < m; j++) k3calc[j] = result[i][j] + h * k2[j] / 2;
        f(t[i] + h/2, k3calc, m, k3);

        // K4 Step
        for (int j = 0; j < m; j++) k4calc[j] = result[i][j] + h * k3[j];
        f(t[i] + h, k4calc, m, k4);

        // Append to result
        for (int j = 0; j < m; j++) {
            result[i+1][j] = result[i][j] + h * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j]) / 6;
        }
    }

    // Free memory
    free(k1);
    free(k2);
    free(k3);
    free(k4);
    free(k2calc);
    free(k3calc);
    free(k4calc);

    return result;

}

double** RKvec(void (*f)(double, double[], int, double*), double* t, double y0[], int m, int n, char fname[]) {
    /*
    Integrates along time points using a butcher formatted RK method using n timesteps on a vector of dimension m

    t: evaluation points
    y0: array of length m containing initial conditions
    m, n: dimensions of vector, number of timesteps
    f: function for which y' = f(t, y, m, ARR), where m is size, and ARR stores the resulting vector
    */

    if (n <= 1) return NULL; // error case

    ////// Allocate result array
    double** result = calloc(n, sizeof(double*));
    for (int i = 0; i < n; i++) result[i] = calloc(m, sizeof(double));

    ////// Prepare Step variables
    double h; // stepsize

    // Read RK Step Parameters
    double** a;
    double *b, *c;
    int s; // number of RK steps
    butcher(&a, &b, &c, fname, &s);

    // Create Step Variables
    double** k = calloc(s+1, sizeof(double*)); // k values
    double** y_i = calloc(s+1, sizeof(double*)); // function inputs at each step

    for (int i = 0; i <= s; i++) {
        k[i] = calloc(m, sizeof(double));
        y_i[i] = calloc(m, sizeof(double));
    }
    
    // Init y0
    for (int j = 0; j < m; j++) result[0][j] = y0[j];

    /////// INTEGRATION STEPS
    for (int i = 0; i < n-1; i++) {
        h = t[i+1] - t[i]; // stepsize

        // RK Step 0
        f(t[i], result[i], m, k[0]);

        // Perform all RK Steps [1, s]
        for (int r = 1; r <= s; r++) {
            // Prepare input vector
            for (int j = 0; j < m; j++) {
                y_i[r][j] = result[i][j];

                // Add previous steps
                for (int w = 0; w < r; w++) y_i[r][j] += h * k[w][j] * a[r][w];
            }
            // evaluate next k
            f(t[i] + h * c[r], y_i[r], m, k[r]);
        }      

        // Append to result
        for (int j = 0; j < m; j++) {
            result[i+1][j] = result[i][j];     
            // Add all weights
            for (int r = 0; r <= s; r++) result[i+1][j] += h * b[r] * k[r][j];
        }
    }

    ////// Free memory
    for (int i = 0; i < s; i++) {
        free(k[i]);
        free(y_i[i]);
    }
    free(k);
    free(y_i);

    return result;

}