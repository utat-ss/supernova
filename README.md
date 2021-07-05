# Supernova 🌟
Supernova is an orbit propagator coded in C for the purpose of accelerating calculations, both on the ground and abord the FINCH spacecraft.

Supernova propagates orbits using Cowell's method, integrating Newton's laws.

# Components
## Solvers
Supernova uses a custom arbitrary explicit Runge-Kutta method to solve vector valued ordinary differential equations. The solver accepts [Butcher Tableaus](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods#Explicit_Runge%E2%80%93Kutta_methods) as definitions of different solver types.

Supernova implements an 8th order method as outlined in [this book, p.288](https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19760017203.pdf), which is on par with what is available with `scipy.solve_ivp`.

### Butcher
The Butcher module interprets Butcher tables and feeds them into the RK solver to specify the steps being taken by the solver.

## Physics
Supernova defines functions for gravitational attraction of the Earth, with perturbation effects coming soon.

### Orbit
Accepts input for initial orbit state vector.

## Python interface
Supernova can be compiled as a `.so` shared library, and imported into Python as `ctypes` for easy interfacing between the programs.

### Plotter
A simple Python script is provided to visualize the results of the integration.

# Building
Run the makefile. `make program` for C, and `make python` for the Python interface. Enjoy!



