# JGM3 Geopotential Gravity Model

def JGM_gravity(t, r):
    '''
    Calculates gravitational forces on spacecraft
    using 20th order spherical harmonic gravity model
    '''
    # generate VW coefficients
    
    