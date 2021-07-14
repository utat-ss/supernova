#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "solvers.h"
#define VEC_SIZE 6
#define H0 60 // initial step


solution* RK810vec(void (*f)(double, double[], double*), double tSpan[], double y0[], double ATOL) {
    /*
    Integrates along time points using an RK810 method using n timesteps on a vector of dimension m
    Source: https://sce.uhcl.edu/feagin/courses/rk10.pdf

    f: function for which y' = f(t, y, ARR), where ARR stores the resulting vector
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
    for (int i = 0; i < n; i++) result->y[i] = malloc(VEC_SIZE * sizeof(double));

    ////// Prepare Step variables
    double h = H0; // initial stepsize
    double h_old; // previous step variable
    double t = tSpan[0]; // intial time

    // Init y0, t0 (step 0)
    for (int j = 0; j < VEC_SIZE; j++) result->y[0][j] = y0[j];
    result->t[0] = t;

    // RK Step Parameters from https://sce.uhcl.edu/rungekutta/rk108.txt
    // a coefficients (written at beta in the paper) [intermediate solution weights]
    double a[17][16] ={
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0.1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {-0.9151765613752915, 1.4545344021782731, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0.20225919030111816, 0.0, 0.6067775709033545, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0.18402471470864357, 0.0, 0.19796683122719236, -0.07295478473136326, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0.08790073402066813, 0.0, 0.0, 0.41045970252026065, 0.4827137536788665, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0.08597005049024603, 0.0, 0.0, 0.3308859630407222, 0.4896629573094502, -0.07318563750708508, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0.12093044912533372, 0.0, 0.0, 0.0, 0.2601246757582956, 0.032540262154909134, -0.0595780211817361, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0.11085437958039149, 0.0, 0.0, 0.0, 0.0, -0.06057614882550056, 0.3217637056017784, 0.510485725608063, 0, 0, 0, 0, 0, 0, 0, 0},
    {0.112054414752879, 0.0, 0.0, 0.0, 0.0, -0.14494277590286592, -0.3332697190962567, 0.4992692295568801, 0.5095046089296861, 0, 0, 0, 0, 0, 0, 0},
    {0.11397678396418598, 0.0, 0.0, 0.0, 0.0, -0.07688133642033569, 0.23952736032439065, 0.3977746623680946, 0.010755895687360746, -0.3277691241640189, 0, 0, 0, 0, 0, 0},
    {0.07983145282801961, 0.0, 0.0, 0.0, 0.0, -0.052032968680060306, -0.05769541461685489, 0.19478191571210415, 0.14538492318832508, -0.07829427103516708, -0.11450329936109892, 0, 0, 0, 0, 0},
    {0.9851156101648573, 0.0, 0.0, 0.3308859630407222, 0.4896629573094502, -1.3789648657484357, -0.8611641950276356, 5.784288136375372, 3.2880776198510357, -2.386339050931364, -3.254793424836439, -2.16343541686423, 0, 0, 0, 0},
    {0.8950802957716328, 0.0, 0.19796683122719236, -0.07295478473136326, 0.0, -0.8512362396620076, 0.3983201123185333, 3.639372631810356, 1.5482287703983033, -2.122217147040537, -1.5835039854532618, -1.7156160828593627, -0.024403640575012746, 0, 0, 0},
    {-0.9151765613752915, 1.4545344021782731, 0.0, 0.0, -0.7773336436449683, 0.0, -0.0910895662155176, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0910895662155176, 0.7773336436449683, 0, 0},
    {0.1, 0.0, -0.15717866579977116, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.15717866579977116, 0},{0.1817813007000953, 0.675, 0.3427581598471898, 0.0, 0.25911121454832275, -0.35827896671795206, -1.0459489594088331, 0.930327845415627, 1.7795095943170811, 0.1, -0.2825475695390441, -0.15932735011997254, -0.14551589464700151, -0.25911121454832275, -0.3427581598471898, -0.675}};
    
    // c coefficients (written as a in the paper) [intermediate timestep values]
    double c[17] = { 0.000000000000000000000000000000000000000000000000000000000000,0.100000000000000000000000000000000000000000000000000000000000,0.539357840802981787532485197881302436857273449701009015505500,0.809036761204472681298727796821953655285910174551513523258250,0.309036761204472681298727796821953655285910174551513523258250,0.981074190219795268254879548310562080489056746118724882027805,0.833333333333333333333333333333333333333333333333333333333333,0.354017365856802376329264185948796742115824053807373968324184,0.882527661964732346425501486979669075182867844268052119663791,0.642615758240322548157075497020439535959501736363212695909875,0.357384241759677451842924502979560464040498263636787304090125,0.117472338035267653574498513020330924817132155731947880336209,0.833333333333333333333333333333333333333333333333333333333333,0.309036761204472681298727796821953655285910174551513523258250,0.539357840802981787532485197881302436857273449701009015505500,0.100000000000000000000000000000000000000000000000000000000000,1.00000000000000000000000000000000000000000000000000000000000};

    // b coefficeints (written as c in the paper) [solution weights]
    double b[17] = {0.0333333333333333333333333333333333333333333333333333333333333,0.0250000000000000000000000000000000000000000000000000000000000,0.0333333333333333333333333333333333333333333333333333333333333,0.000000000000000000000000000000000000000000000000000000000000,0.0500000000000000000000000000000000000000000000000000000000000,0.000000000000000000000000000000000000000000000000000000000000,0.0400000000000000000000000000000000000000000000000000000000000,0.000000000000000000000000000000000000000000000000000000000000,0.189237478148923490158306404106012326238162346948625830327194,0.277429188517743176508360262560654340428504319718040836339472,0.277429188517743176508360262560654340428504319718040836339472,0.189237478148923490158306404106012326238162346948625830327194,-0.0400000000000000000000000000000000000000000000000000000000000,-0.0500000000000000000000000000000000000000000000000000000000000,-0.0333333333333333333333333333333333333333333333333333333333333,-0.0250000000000000000000000000000000000000000000000000000000000,0.0333333333333333333333333333333333333333333333333333333333333};

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
                for (int w = 0; w < r; w++) y_i[j] += h * k[w][j] * a[r][w]; // Add previous steps
            }
            f(t + h * c[r], y_i, k[r]); // evaluate next k
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
                n *= 2;
                result->y = realloc(result->y, n * sizeof(double*));
                result->t = realloc(result->t, n * sizeof(double));
                for (int i = step + 1; i < n; i++) result->y[i] = malloc(VEC_SIZE * sizeof(double));
            }

            // Append to result
            for (int j = 0; j < VEC_SIZE; j++) {
                for (int r = 0; r < 17; r++) y_curr[j] += h_old * b[r] * k[r][j]; // Add all weights
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