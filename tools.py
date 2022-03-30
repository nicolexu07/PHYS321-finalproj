import wget
import os
import requests
import pandas as pd
import numpy as np
from io import StringIO
import re
from scipy.optimize import fsolve

cwd = os.getcwd()



def get_links(site):
    """ (str) -> (list)
    Returns a list of strings, each entry in the list is a http:// link found on the website site.
    Ignores links beginning with http://irsa as these do not link to .tbl files.
    """
    site_str = str(requests.get(site).content) #contents of site as string
    links = []

    i = 0
    while i < len(site_str)-8:
        if site_str[i:i+7]=='http://' and site_str[i+7:i+7+4]!='irsa':
            link = 'http://'
            j = 0
            while site_str[i+7+j] != ' ':
                link += site_str[i+7+j]
                j += 1
            links.append(link)
            i += j
        i += 1
    return links



def download_from_link(link):
    """ (str) -> ()
    Downloads .tbl files found through the links at:
    https://exoplanetarchive.ipac.caltech.edu/bulk_data_download/wget_RADIAL.bat
    Saves them to a folder in the current working directory called data.
    """
    wget.download(link, cwd+'\data')



def download_radial():
    """ () -> ()
    Downloads all the .tbl files linked to in the exoplanetarchive.
    """
    links = get_links('https://exoplanetarchive.ipac.caltech.edu/bulk_data_download/wget_RADIAL.bat')
    for link in links:
        download_from_link(link)



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
def func(u, tau, T, e):
    return u-e*np.sin(u)-((2*np.pi/T)*(t-tau))



def solve_for_u(t, tau, T, e):
    # since esinu < u, using RHS of eqn as guess 
    u_guess = (2*np.pi/T)*(t-tau)
    root = fsolve(func, u_guess, args=(tau, T, e))
    
    return root 



def radial_velocity(t, m, M, T, I, e, v_0, omega, tau):
    kappa = ((2*np.pi*G)**(1/3)*m*np.sin(I))/(T**(1/3)*(M+m)**(2/3)*np.sqrt(1-e**2))
    u = solve_for_u(t, tau, T, e)
    f = 2*arctan(np.sqrt((1+e)/(1-e))*tan(u/2))
    rad_vel = kappa*(np.cos(f+omega)+e*np.cos(omega))+v_0
    
    return rad_vel



if __name__ == "__main__":
    import doctest
    doctest.testmod()