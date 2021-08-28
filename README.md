# Supernova 🌟
Supernova is an orbit propagator coded in C for the purpose of accelerating calculations, both on the ground and aboard the FINCH spacecraft.

Supernova propagates orbits using Cowell's method, integrating Newton's laws.

# Components
## Solvers
Supernova uses an [adaptive 10th order method](https://sce.uhcl.edu/feagin/courses/rk10.pdf) with an truncation error estimate of eigth order called RK810. It outperforms an equivalent DOP853 based Python solver by a factor of 100x in terms of runtime when tested using only J2 perturbations.

## Physics
Supernova defines functions for gravitational attraction of the Earth, with advanced perturbation effects coming soon. Currently, only J2 and atmospheric drag are supported.

Aero characteristics are currently hard-coded for FINCH (drag coefficient = 2.2, Area to mass ratio = 0.01 m^2 / kg)

### Orbit
Accepts input for initial orbit state vector and serves as a wrapper for the RK810 solver.

## Python interface
Supernova can be compiled as a `.so` shared library, and imported into Python as `ctypes` for easy interfacing between the programs.

### Plotter
A simple Python script is provided to visualize the results of the integration.

# Building
Run the makefile. `make program` for C, and `make python` for the Python interface. Enjoy!



