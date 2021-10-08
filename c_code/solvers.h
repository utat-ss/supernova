// Numerical Runge Kutta Solvers for Supernova

#ifndef SOLVERS_H
#define SOLVERS_H

#define VEC_SIZE 6 // dimension of vector (kept static at 6 for orbit prop)
#define H0810 60 // initial step
#define H01012 200 // initial step
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "constants.h"

typedef struct AdaptiveSolution {
    double* t; // timesteps
    double (*y)[VEC_SIZE]; // function value
    int n; // number of steps taken
} solution;

solution* RK810vec(void (*f)(double, double[], double*), double tSpan[2], double y0[VEC_SIZE], double ATOL) {
    /*
    Integrates along time points using an RK810 method using n timesteps on a vector of dimension m
    Source: https://sce.uhcl.edu/feagin/courses/rk10.pdf

    f: function for which y' = f(t, y, ARR), where ARR stores the resulting vector
    tSpan: stores (t0, tf)
    y0: array of length [VEC_SIZE] containing initial conditions
    ATOL: tolerance
    */

    if (tSpan[1] <= tSpan[0]) return NULL; // error case

    ////// Initial estimate for n = number of rows in solution
    int n = (int)((tSpan[1] - tSpan[0])/H0810); // 1 step every 10 seconds
    int step = 0;

    ////// Allocate result struct
    solution* result = (solution*) malloc(sizeof(solution));
    result->y = (double(*)[6]) malloc(n * sizeof(double[VEC_SIZE]));
    result->t = (double(*)) malloc(n * sizeof(double));

    ////// Prepare Step variables
    double h = H0810; // initial stepsize
    double h_old; // previous step variable
    double t = tSpan[0]; // intial time

    // Init y0, t0 (step 0)
    for (int j = 0; j < VEC_SIZE; j++) result->y[0][j] = y0[j];
    result->t[0] = t;

    // Create Step Variables
    double k[17][VEC_SIZE]; // k values
    double y_i[VEC_SIZE]; // function inputs at each step
    double y_curr[VEC_SIZE]; // current solution array, and also used to update itself for next step

    double err; // magnitude of error

    /////// INTEGRATION STEPS
    while (t-h < tSpan[1]) {
        for (int j = 0; j < VEC_SIZE; j++) y_curr[j] = result->y[step][j]; // get current y (for vectorization purposes)

        f(t, y_curr, k[0]); // RK Stage 0

        // Perform all RK Stages [1, s)
        for (int r = 1; r < 17; r++) {
            //// Prepare input vector
            for (int j = 0; j < VEC_SIZE; j++) {             
                y_i[j] = y_curr[j]; // take current sol      
                for (int w = 0; w < r; w++) y_i[j] += h * k[w][j] * RK810a[r][w]; // Add previous steps
            }
            f(t + h * RK810c[r], y_i, k[r]); // evaluate next k
        }

        // Calculate error using error estimate formula from paper
        // take abs value of all errors as convervative estimate
        err = 0;
        for (int j = 0; j < VEC_SIZE; j++) err += fabs(1.0/360.0 * h * (k[1][j] - k[15][j]));

        //// Step size adjustment
        // Determine new step size
        h_old = h;
        h = 0.9 * h * pow(ATOL/err, 1.0/9.0); // next step size (based on ninth order local error)

        if (err < ATOL) {
            /// step within tolerance, append solution
            // check array size and increase size if needed
            if ((step + 1) >= n) {
                n *= 2; // double size of array and reallocate memory
                result->y = (double(*)[6]) realloc(result->y, n * sizeof(double[VEC_SIZE]));
                result->t = (double(*)) realloc(result->t, n * sizeof(double));
            }

            // Append to result
            for (int j = 0; j < VEC_SIZE; j++) {
                for (int r = 0; r < 17; r++) y_curr[j] += h_old * RK810b[r] * k[r][j]; // Add all weights
            }
            for (int j = 0; j < VEC_SIZE; j++) result->y[step+1][j] = y_curr[j]; // write (separated for vectorization)

            step++; // advance step
            t += h_old; // advance time
            result->t[step] = t; // record time
        }
        // Otherwise, retry step
    }
    // record # of steps after finish
    result->n = step;

    return result;
}

solution* RK1012vec(void (*f)(double, double[], double*), double tSpan[2], double y0[VEC_SIZE], double ATOL) {
    /*
    Integrates along time points using an RK1012 method using n timesteps on a vector of dimension m
    Source: https://sce.uhcl.edu/rungekutta/rk1210.txt

    f: function for which y' = f(t, y, ARR), where ARR stores the resulting vector
    tSpan: stores (t0, tf)
    y0: array of length [VEC_SIZE] containing initial conditions
    ATOL: tolerance
    */

    if (tSpan[1] <= tSpan[0]) return NULL; // error case

    ////// Initial estimate for n
    int n = (int)((tSpan[1] - tSpan[0])/H01012); // 1 step every 10 seconds
    int step = 0;

    ////// Allocate result struct
    solution* result = (solution*) malloc(sizeof(solution));
    result->y = (double(*)[6]) malloc(n * sizeof(double[VEC_SIZE]));
    result->t = (double(*)) malloc(n * sizeof(double));


    ////// Prepare Step variables
    double h = H01012; // initial stepsize
    double h_old; // previous step variable
    double t = tSpan[0]; // intial time

    // Init y0, t0 (step 0)
    for (int j = 0; j < VEC_SIZE; j++) result->y[0][j] = y0[j];
    result->t[0] = t;

    // Create Step Variables
    double k[25][VEC_SIZE]; // k values
    double y_i[VEC_SIZE]; // function inputs at each step
    double y_curr[VEC_SIZE]; // current solution array, and also used to update itself for next step

    double err; // magnitude of error

    /////// INTEGRATION STEPS
    while (t-h < tSpan[1]) {
        for (int j = 0; j < VEC_SIZE; j++) y_curr[j] = result->y[step][j]; // get current y (for vectorization purposes)

        f(t, y_curr, k[0]); // RK Stage 0

        // Perform all RK Stages [1, 25)
        for (int r = 1; r < 25; r++) {
            //// Prepare input vector
            for (int j = 0; j < VEC_SIZE; j++) {             
                y_i[j] = y_curr[j]; // take current sol      
                for (int w = 0; w < r; w++) y_i[j] += h * k[w][j] * RK1012a[r][w]; // Add previous steps
            }
            f(t + h * RK1012c[r], y_i, k[r]); // evaluate next k
        }

        // Calculate error using error estimate formula from paper
        // take abs value of all errors as convervative estimate
        err = 0;
        for (int j = 0; j < VEC_SIZE; j++) err += fabs(49.0/640.0 * h * (k[1][j] - k[23][j]));

        //// Step size adjustment
        // Determine new step size
        h_old = h;
        if (err == 0) h *= 1.1; // if error is too small to be detected, increase step size
        else h = 0.9 * h * pow(ATOL/err, 1.0/11.0); // next step size (based on eleventh order local error)

        if (err < ATOL) {
            /// step within tolerance, append solution
            // check array size and increase size if needed
            if ((step + 1) >= n) {
                n *= 2; // double size of array and reallocate memory
                result->y = (double(*)[6]) realloc(result->y, n * sizeof(double[VEC_SIZE]));
                result->t = (double(*)) realloc(result->t, n * sizeof(double));
            }

            // Append to result
            for (int j = 0; j < VEC_SIZE; j++) {
                for (int r = 0; r < 25; r++) y_curr[j] += h_old * RK1012b[r] * k[r][j]; // Add all weights
            }
            for (int j = 0; j < VEC_SIZE; j++) result->y[step+1][j] = y_curr[j]; // write (separated for vectorization)

            step++; // advance step
            t += h_old; // advance time
            result->t[step] = t; // record time
        }
        // Otherwise, retry step
    }
    // record # of steps after finish
    result->n = step;

    return result;
}

void RK1012Step(void (*f)(double, double[], double*), double h, double y_curr[VEC_SIZE], double y_out[VEC_SIZE]) {
    /*
    Integrates a single step of RK1012, agnostic to time

    f: function for which y' = f(t, y, ARR), where ARR stores the resulting vector
    h: stepsize (make sure this is small enough)
    y_curr: initial state
    y_out: where final state is stored
    */

    double t = 0; // filler time

    // Create Step Variables
    double k[25][VEC_SIZE]; // k values
    double y_i[VEC_SIZE]; // function inputs at each step

    for (int j = 0; j < VEC_SIZE; j++) y_out[j] = y_curr[j]; // Null out vector for now

    f(t, y_curr, k[0]); // RK Stage 0

    // Perform all RK Stages [1, 25)
    for (int r = 1; r < 25; r++) {
        //// Prepare input vector
        for (int j = 0; j < VEC_SIZE; j++) {             
            y_i[j] = y_curr[j]; // take current sol      
            for (int w = 0; w < r; w++) y_i[j] += h * k[w][j] * RK1012a[r][w]; // Add previous steps
        }
        f(t + h * RK1012c[r], y_i, k[r]); // evaluate next k
    }

    // Append to result
    for (int j = 0; j < VEC_SIZE; j++) {
        for (int r = 0; r < 25; r++) y_out[j] += h * RK1012b[r] * k[r][j]; // Add all weights
    }
}

#endif