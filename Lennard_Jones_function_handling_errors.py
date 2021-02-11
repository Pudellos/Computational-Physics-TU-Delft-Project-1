import numpy as np
import sympy as sym

quiet=True #if True all prited warnings are supressed


def k_B(N,P=10*5, V=1, T=273.15):
    ''' Calculation of Boltzmann constant for the simulation
    
    parameters:
    P = pressure [Pa]
    V= volume [m^3]
    T = temperature [K]
    N = number of gas molecules in the simulation '''
    try:
        return((P*V)/(T*N))
    except ZeroDivisionError:
        if not quiet:
            print('T, N cannot be equal to zero, I return/use k_B at stp for 1 gas molecule')
        return((P*V)/T) 

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
    sigma: float = 3.405
    eps : float = 119.8/k_B
    r_abs : float = distance between molecules A and B
    '''
    def equation(sigma,eps,r_abs):
            return((4*eps*((sigma/r_abs)**12-(sigma/r_abs)**6)))
    
    sigma = sym.Symbol('sigma')
    eps = sym.Symbol('eps')
    r_abs=sym.Symbol('r_abs')
    
    sigma_value=3.405
    try:
        eps_value=119.8/k_B(N,P,V,T)
    except ZeroDivisionError:
        if not quiet:
            print('P, V cannot be equal to zero, I use k_B at stp')
        eps_value=119.8/k_B(N)
    except TypeError:
        if not quiet:
            print('something went from, check your inputs!')
    
    try:
        r_abs_value=abs(r[0]-r[1])
    except:
        if not quiet:
            print('something went wrong, check your inputs')
    if r_abs_value==0:
        if not quiet:
            print('the molecules have identical position, impossible, I return zero')
    
    if not differenciate:
        eqn=equation(sigma,eps,r_abs)
        if not symbol_mode:
            eqn=eqn.evalf(subs={sigma: sigma_value, eps: eps_value, r_abs:r_abs_value})
            
    if differenciate:
        eqn=sym.diff(equation(sigma,eps,r_abs), r_abs)
        if not symbol_mode:
            eqn=eqn.evalf(subs={sigma: sigma_value, eps: eps_value, r_abs:r_abs_value})
            
    try:
        return(eqn)  
    except:
        if not quiet:
            print('something went wrong, check your inputs')
        return(0)