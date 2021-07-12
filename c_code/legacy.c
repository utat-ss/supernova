#include <stdlib.h>
#include <stdio.h>
#include "solvers.h"
#include "physics.h"
#include <time.h>
#include <math.h>

/*
Legacy code that is no longer used by supernova, but is here for sake of learning
*/

// Parser for butcher table
void butcher(double*** a, double** b, double** c, char filename[], int* num) {
    /*
    a (pointer to 2D array) = steps: lower triangular array of size n, should be sent UNALLOCATED
    b (pointer to 1D array) = sums: single array of size n, should be sent UNALLOCATED
    c (pointer to 1D array) = weights: single array of size n, should be sent UNALLOCATED
    
    filename: filename of butcher tableau, space delimited with all numbers written as rationals
    n: length of tableau
    */
    FILE *fp = fopen(filename, "r");

    // PREAMBLE
    int lines = 0;
    
    // scan number of lines
    while(!feof(fp)) if(fgetc(fp) == '\n') lines++;
    rewind(fp); // rewind

    // PARSE
    // parse steps
    int n = lines-1;
    *num = n; // assign length

    *a = malloc(n * sizeof(double*));
    *b = malloc(n * sizeof(double));
    *c = malloc(n * sizeof(double));

    // BUTCHER TABLE CONSTRUCTION
    // first row of butcher table
    (*c)[0] = 0;
    (*a)[0] = NULL; // not needed, so make it null

    // dump first line of buffer
    float f1, f2; // numerator and denominator
    // MUST BE FLOAT TO PARSE PROPERLY

    fscanf(fp, "%f/%f\n", &f1, &f2);

    // start from the second matrix row, since the first is just zero.
    for (int i = 1; i <= n; i++) {        
        
        fscanf(fp, "%f/%f", &f1, &f2);
        (*c)[i] = f1/f2; // parse c coefficient

        // prepare list for a coefficients in this row
        (*a)[i] = malloc(i * sizeof(double));

        for (int j = 0; j < i; j++) {
            fscanf(fp, " %f/%f", &f1, &f2);
            (*a)[i][j] = f1/f2;
        }
        
        fscanf(fp, "\n"); // Skip line break
    }

    // SUMS
    fscanf(fp, "%f/%f", &f1, &f2);
    (*b)[0] = f1/f2;
    for (int i = 1; i <= n; i++) {
        fscanf(fp, " %f/%f", &f1, &f2);
        (*b)[i] = f1/f2;
    }
}

// Parser for adaptive butcher table
void butcher_adaptive(double*** a, double** b1, double** b2, double** c, char filename[], int* num) {
    /*
    a (pointer to 2D array) = steps: lower triangular array of size n, should be sent UNALLOCATED
    b1, b2 (pointer to 1D array) = sums: single array of size n, should be sent UNALLOCATED
    b1 should contain the higher order method terms
    c (pointer to 1D array) = weights: single array of size n, should be sent UNALLOCATED
    
    filename: filename of butcher tableau, space delimited with all numbers written as rationals
    n: length of tableau
    */
    FILE *fp = fopen(filename, "r");

    // PREAMBLE
    int lines = 0;
    
    // scan number of lines
    while(!feof(fp)) if(fgetc(fp) == '\n') lines++;
    lines++;
    rewind(fp); // rewind

    // PARSE
    // parse steps
    int n = lines-2; // since last two rows are sums
    *num = n; // assign length

    printf("Integrator Loaded. Name: %s Stages: %d\n", filename, n);

    *a = malloc(n * sizeof(double*));
    *b1 = malloc(n * sizeof(double));
    *b2 = malloc(n * sizeof(double));
    *c = malloc(n * sizeof(double));

    // BUTCHER TABLE CONSTRUCTION
    // first row of butcher table
    (*c)[0] = 0;
    (*a)[0] = NULL; // not needed, so make it null

    // dump first line of buffer
    float f1, f2; // numerator and denominator
    // MUST BE FLOAT TO PARSE PROPERLY

    fscanf(fp, "%f/%f\n", &f1, &f2);

    // start from the second matrix row, since the first is just zero.
    for (int i = 1; i < n; i++) {        
        
        fscanf(fp, "%f/%f", &f1, &f2);
        (*c)[i] = f1/f2; // parse c coefficient

        // prepare list for a coefficients in this row
        (*a)[i] = malloc(i * sizeof(double));

        for (int j = 0; j < i; j++) {
            fscanf(fp, " %f/%f", &f1, &f2);
            (*a)[i][j] = f1/f2;
        }
        
        fscanf(fp, "\n"); // Skip line break
    }

    // SUMS
    fscanf(fp, "%f/%f", &f1, &f2);
    (*b1)[0] = f1/f2;
    for (int i = 1; i < n; i++) {
        fscanf(fp, " %f/%f", &f1, &f2);
        (*b1)[i] = f1/f2;
    }
    fscanf(fp, "\n"); // Skip line break

    fscanf(fp, "%f/%f", &f1, &f2);
    (*b2)[0] = f1/f2;
    for (int i = 1; i < n; i++) {
        fscanf(fp, " %f/%f", &f1, &f2);
        (*b2)[i] = f1/f2;
    }
}

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

    // INITIAL STEPSIZE OF 20

    ////// Initial estimate for n
    int n = (int)((tSpan[1] - tSpan[0])/20); // 1 step every 10 seconds
    int step = 0;

    ////// Allocate result struct
    solution* result = (solution*) malloc(sizeof(solution));
    result->y = malloc(n * sizeof(double*));
    result->t = malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) result->y[i] = malloc(m * sizeof(double));

    ////// Prepare Step variables
    double h = 20; // initial stepsize
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
        //if (step == 0) printf("Sol1: %f %f   Sol2: %f %f \n", sol1[0], sol1[1], sol2[0], sol2[1]);

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
                for (int r = 0; r < s; r++) result->y[step+1][j] += h_old * b1[r] * k[r][j];
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

void orbitRK(double tSpan[], double y0[], char integrator[], int ord, char output[], double ATOL) {
    /*
    USES CUSTOM RK METHODS FOR ORBIT PROP
    tSpan: t0, tf of times.
    y0: initial x y z vx vy vz state vector
    integrator: path to a butcher tableau eg. ./tables/RK8.txt
    ord: order of integrator eg. 5
    output: filename of output csv data file
    */
    clock_t start, end;

    start = clock();
    solution* result = RKvec_adaptive(combined_perturbations, tSpan, y0, 6, ord, integrator, ATOL);
    end = clock();
    double cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    printf("Integration completed in %f seconds. Steps taken: %d\n", cpu_time_used, result->n);

    FILE *fptr = fopen(output, "w");

    for (int i = 0; i < result->n; i++) {
        fprintf(fptr, "%.15f, %.15f, %.15f, %.15f, %.15f, %.15f, %.15f\n", result->t[i], result->y[i][0], result->y[i][1], result->y[i][2], result->y[i][3], result->y[i][4], result->y[i][5]);
        free(result->y[i]);
    }
    fclose(fptr);

    free(result->t);
    free(result);
}
