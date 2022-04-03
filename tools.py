from importlib.resources import path
import pandas as pd
import numpy as np
from io import StringIO
import re
from scipy.optimize import fsolve
from astropy import constants as const

def get_data(filename, path='data/'):
    """ (str, str) -> (pd.DataFrame)
    Takes a file name as input. 
    Returns a pandas dataframe with 3 cols corresponding to days, radial velocity, and uncertainty.
    Converts Julian days to seconds setting the earliest time as 0.
    
    >>> temp = get_data('test_data.tbl', '')
    >>> np.all(temp[0]==pd.Series([0, 86400]))
    True
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
    Takes file name as input and returns numpy array with observation information as entries. 
    """
    lines = '0 1\n' + "".join([re.sub('=', ' ', re.sub('\'', '', re.sub(r"\s+", '', line)))+'\n' 
                    for line in open(f'./{path}{filename}') if line.startswith('\ '[0])])
    #print(lines)

    df = pd.read_csv(StringIO(lines), delimiter=' ', index_col=False)
    #df = df.drop('1', axis=1)
    df =df.rename(columns={'0':0, '1':1})

    return df



# function to use in scipy.optimize.fsolve
def func_u(u, t, tau, T, e):
  
    return u-e*np.sin(u)-((2*np.pi/T)*(t-tau))



def solve_for_u(t, tau, T, e):
    """ (np.array, num, num, num) -> (np.array)
    Numerically solves u-e*np.sin(u)-((2*np.pi/T)*(t-tau))=0 for u.


    When t=tau we expect u=0
    Numerical solver gives approximately 0
    >>> e = np.random.uniform(0, 1)
    >>> T = np.random.uniform(3282.3503, 3.46896e13)
    >>> tau = np.random.uniform(3282.3503, 3.46896e13)    
    >>> temp = solve_for_u(tau, tau, T, e)
    >>> abs(temp) < 1e-7
    array([ True])

    When e=0 the function is linear
    >>> e = 0
    >>> t = np.random.uniform(0, T, 10)
    >>> temp = solve_for_u(t, tau, T, e)
    >>> np.all(abs(2*np.pi*(t-tau)/T - temp) < 1e-7)
    True
    """
    # since esinu < u, using RHS of eqn as guess 
    u_guess = (2*np.pi/T)*(t-tau)
    root = fsolve(func_u, u_guess, args=(t, tau, T, e))
    
    return root 



def radial_velocity(t, m, M, T, I, e, v_0, omega, tau):
    kappa = ((2*np.pi*const.G.value)**(1/3)*m*np.sin(I))/(T**(1/3)*(M+m)**(2/3)*np.sqrt(1-e**2))
    u = solve_for_u(t, tau, T, e)
    f = 2*np.arctan(np.sqrt((1+e)/(1-e))*np.tan(u/2))
    rad_vel = kappa*(np.cos(f+omega)+e*np.cos(omega))+v_0
    
    return rad_vel



def list_files(dir):
    """ (str) -> (list)
    Lists the names of all files in the dir directory found in the current working directory.
    """
    return os.listdir("./"+dir)


def get_star_id(filename, path='data/'):
    """ (str, str) -> (str)
    Returns the star id of the star in filename.

    >>> get_star_id('UID_0302504_RVC_002.tbl')
    'HATS-2'

    >>> get_star_id('test_data.tbl', '')
    'TestStar'

    >>> get_star_id('UID_0250999_RVC_001.tbl')
    'GJ3634'
    """
    return get_obs_info(filename, path)[1][0]


def find_files_for_star(star_id):
    """ (str) -> (list)
    Finds all the files in /data/ with with data pertaining to star_id.
    
    >>> find_files_for_star('CrazyFrog')
    Traceback (most recent call last):
        ...
    ValueError: No files with this star id found.

    >>> find_files_for_star('HD 142')
    ['UID_0000522_RVC_001.tbl', 'UID_0000522_RVC_002.tbl']
    """
    star_id = re.sub(r'\s+', '', star_id)
    files = list_files('data')
    ans = []
    for file in files:
        if get_star_id(file) == star_id:
            ans.append(file)
    if len(ans) == 0:
        raise ValueError('No files with this star id found.')
    return ans

    
if __name__ == "__main__":
    import doctest
    doctest.testmod()