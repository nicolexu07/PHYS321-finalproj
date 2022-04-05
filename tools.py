import pandas as pd
import numpy as np
from io import StringIO
import re
from scipy.optimize import fsolve
from astropy import constants as const

def get_data(filename, path='data/'):
    """ (str, str) -> (pd.DataFrame)
    Takes a file name as input 
    Returns a pandas dataframe with 3 cols corresponding
    to days, radial velocity, and uncertainty
    Converts Julian days to seconds setting the earliest time as 0
    """
    lines = "0,1,2\n" + "".join([re.sub(r"\s+", ',', line)[1::]+'\n' for line in open(f'./{path}{filename}') 
                     if not (line.startswith('\ '[0]) or line.startswith('|'))])
    #print(lines)
    df = pd.read_csv(StringIO(lines), sep=',', index_col=False)
    df = df.rename(columns={'0':0, '1':1, '2':2})
    df[0] = (df[0]-min(df[0]))*86400
    return df
        
    
 
def get_obs_info(filename, path='data/'):
    """ (str, str) -> (np.array)
    Takes file name as input and returns numpy array with
    observation information as entries 
    """
    lines = '0,1,2\n' + "".join([re.sub(r"\s+", ',', line)[1::]+'\n' for line in open(f'./{path}{filename}') 
                     if line.startswith('\ '[0])])
    #print(lines)
    df = pd.read_csv(StringIO(lines), sep=',', index_col=False)
    df = df.drop('1', axis=1)
    df =df.rename(columns={'0':0, '2':1})
    
    return df



# function to use in scipy.optimize.fsolve
def funct(u, t, tau, T, e):
    return u-e*np.sin(u)-((2*np.pi/T)*(t-tau))



def solve_for_u(t, tau, T, e):
    # since esinu < u, using RHS of eqn as guess 
    u_guess = (2*np.pi/T)*(t-tau)
    root = fsolve(funct, u_guess, args=(t, tau, T, e))
    
    return root 



def radial_velocity(t, m, M, T, I, e, v_0, omega, tau):
    kappa = ((2*np.pi*const.G.value)**(1/3)*m*np.sin(I))/(T**(1/3)*(M+m)**(2/3)*np.sqrt(1-e**2))
    u = solve_for_u(t, tau, T, e)
    f = 2*np.arctan(np.sqrt((1+e)/(1-e))*np.tan(u/2))
    rad_vel = kappa*(np.cos(f+omega)+e*np.cos(omega))+v_0
    
    return rad_vel



if __name__ == "__main__":
    import doctest
    doctest.testmod()