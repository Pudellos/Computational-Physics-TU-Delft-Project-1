import numpy as np
import sympy as sym

def k_B(P,V,T,N):
    ''' Calculation of Boltzmann constant for the simulation
    
    parameters:
    P = pressure [Pa]
    V= volume [m^3]
    T = temperature [K]
    N = number of gas molecules in the simulation '''
    return((P*V)/(T*N))

def LJ_potential(r, N, P=10*5, V=1, T=273.15, symbol_mode=False, differenciate=False):
    '''Calculating Lennar Jones potential between particle 1 and 2
    
    parameters:
    r : np.array = positions of molecules 1 and 2 [m]
    P : float = pressure [Pa]
    V : float = volume [m^3]
    T : float = temperature [K]
    N : int = number of gas molecules in the simulation
    symbol_mode : bool 
    differenciate : bool

    how to use:
    function returns Lennard Jones potential value
    if symbol_mode == True : symbolic equation returned
    if differenciate == True : derivative of Lennard Jones potential is returned 
    
    parameters used by symbol_mode:
    sigma: float = 3.405 [A] = 3.405e-10 [m]
    eps : float = 119.8/k_B
    r_abs : float = distance between molecules A and B
    '''
    def equation(sigma,eps,r_abs):
        return((4*eps*((sigma/r_abs)**12-(sigma/r_abs)**6)))
    
    sigma = sym.Symbol('sigma')
    eps = sym.Symbol('eps')
    r_abs=sym.Symbol('r_abs')
    
    sigma_value=3.405e-10
    eps_value=119.8/k_B(P,V,T,N) 
    try:
        if len(r)==2:
            r_abs_value=abs(r[0]-r[1])
    except TypeError: 
        if len(r[0])==2:
            r_abs_value_x=abs(r[0][0]-r[1][0])
            r_abs_value_y=abs(r[0][1]-r[1][1])
            r_abs_value=(r_abs_value_x**2+r_abs_value_y**2)**(1/2)
    
    if not differenciate:
        eqn=equation(sigma,eps,r_abs)
        if not symbol_mode:
            eqn=eqn.evalf(subs={sigma: sigma_value, eps: eps_value, r_abs:r_abs_value})
            
    if differenciate:
        eqn=sym.diff(equation(sigma,eps,r_abs), r_abs)
        if not symbol_mode:
            eqn=eqn.evalf(subs={sigma: sigma_value, eps: eps_value, r_abs:r_abs_value})
            
    return(eqn) 