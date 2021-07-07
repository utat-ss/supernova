#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "butcher.h"
#include "solvers.h"

const double H0 = 20; // initial step

double* euler(double (*f)(double, double), double* t, double y0, int n) {
    // Integrates from one side of the interval to another using Euler's method
    // f: function for which y' = f(t, y)
    // t: evaluation points
    // y0: initial condition
    // n: number of steps

    double* result = malloc(n * sizeof(double));

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

    double* result = malloc(n * sizeof(double));

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
    double** result = malloc(n * sizeof(double*));
    for (int i = 0; i < n; i++) result[i] = malloc(m * sizeof(double));

    ////// Prepare Step variables
    double h; // stepsize

    // Read RK Step Parameters
    double** a;
    double *b, *c;
    int s; // number of RK steps
    butcher(&a, &b, &c, fname, &s);

    // Create Step Variables
    double** k = malloc((s+1) * sizeof(double*)); // k values
    double* y_i = malloc(m * sizeof(double)); // function inputs at each step

    for (int i = 0; i <= s; i++) k[i] = malloc(m * sizeof(double)); // init k values
    
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
                y_i[j] = result[i][j];

                // Add previous steps
                for (int w = 0; w < r; w++) y_i[j] += h * k[w][j] * a[r][w];
            }
            // evaluate next k
            f(t[i] + h * c[r], y_i, m, k[r]);
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
    }
    free(k);
    free(y_i);

    return result;

}

solution* RKvec_adaptive(void (*f)(double, double[], int, double*), double tSpan[], double y0[], int m, int ord, char fname[], double ATOL) {
    /*
    Integrates along time points using a butcher formatted RK method using n timesteps on a vector of dimension m

    f: function for which y' = f(t, y, m, ARR), where m is dimension of the vector, and ARR stores the resulting vector
    tSpan: stores (t0, tf)
    y0: array of length m containing initial conditions
    m: dimension of vector
    ord: order of method
    fname: integrator filename
    ATOL: tolerance
    */

    if (tSpan[1] <= tSpan[0]) return NULL; // error case

    ////// Initial estimate for n
    int n = (int)((tSpan[1] - tSpan[0])/H0); // 1 step every 10 seconds
    int step = 0;

    ////// Allocate result struct
    solution* result = (solution*) malloc(sizeof(solution));
    result->y = malloc(n * sizeof(double*));
    result->t = malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) result->y[i] = malloc(m * sizeof(double));

    ////// Prepare Step variables
    double h = H0; // initial stepsize
    double h_old;
    double t = tSpan[0]; // intial time

    // Init y0, t0 (step 0)
    for (int j = 0; j < m; j++) result->y[0][j] = y0[j];
    result->t[0] = t;

    // Read RK Step Parameters
    double** a;
    double *b1, *b2, *c;
    int s; // number of RK stages in the method
    butcher_adaptive(&a, &b1, &b2, &c, fname, &s);

    // Create Step Variables
    double** k = malloc((s+1) * sizeof(double*)); // k values
    double* y_i = malloc(m * sizeof(double)); // function inputs at each step
    for (int i = 0; i <= s; i++) k[i] = malloc(m * sizeof(double)); // init k values

    double* sol1 = malloc(m * sizeof(double)); // temporary solution for error term
    double* sol2 = malloc(m * sizeof(double)); // temporary solution for error term
    double err; // magnitude of error

    /////// INTEGRATION STEPS
    while (t-h < tSpan[1]) {

        // RK Stage 0
        f(t, result->y[step], m, k[0]);

        // Perform all RK Stages [1, s)
        for (int r = 1; r < s; r++) {
            // Prepare input vector
            for (int j = 0; j < m; j++) {
                y_i[j] = result->y[step][j];    
                for (int w = 0; w < r; w++) {
                    y_i[j] += h * k[w][j] * a[r][w]; // Add previous steps

                    //if (j == 0) printf("RK Weight for stage %d is %f\n", w+1, a[r][w]);
                } 
            }
            f(t + h * c[r], y_i, m, k[r]); // evaluate next k
        } 

        // Calculate error
        for (int j = 0; j < m; j++) {
            sol1[j] = 0; // zero out vector from last time
            sol2[j] = 0;
        }
        for (int j = 0; j < m; j++) {
            // Add all weights
            for (int r = 0; r < s; r++) {
                sol1[j] += h * b1[r] * k[r][j];
                sol2[j] += h * b2[r] * k[r][j];
            }
        }
        //printf("Sol1: %f %f   Sol2: %f %f \n", sol1[0], sol1[1], sol2[0], sol2[1]);

        // Determine RMS error of vector
        err = 0;
        for (int j = 0; j < m; j++) err += pow(sol1[j] - sol2[j], 2);
        err = pow(err, 0.5);

        //printf("Step %d, Error: %f Stepsize: %f\n", step, err, h);

        //// Step size adjustment
        // Determine new step size
        h_old = h;
        h = 0.9 * h * pow(ATOL/err, 1./((float) ord)); // next step size

        if (err < ATOL) {
            // step within tolerance

            // check array size
            if ((step + 1) >= n) {
                // reallocate array if needed
                n *= 2;
                result->y = realloc(result->y, n * sizeof(double*));
                result->t = realloc(result->t, n * sizeof(double));
                for (int i = step + 1; i < n; i++) result->y[i] = malloc(m * sizeof(double));
            }

            // Append to result
            for (int j = 0; j < m; j++) {
                result->y[step+1][j] = result->y[step][j];     
                // Add all weights
                for (int r = 0; r < s; r++) result->y[step+1][j] += h * b1[r] * k[r][j];
            }
            step++; // advance step
            t += h_old; // advance time
            result->t[step] = t; // record time
        }
        // Otherwise, retry step
    }

    // record # of steps
    result->n = step;

    ////// Free memory
    for (int i = 0; i < s; i++) free(k[i]);
    free(k);
    free(y_i);

    return result;
}