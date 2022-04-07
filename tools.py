import pandas as pd
import numpy as np
from io import StringIO
import re
from scipy.optimize import fsolve
from astropy import constants as const
import emcee
import corner
import matplotlib.pyplot as plt
import os

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


def radial_velocity(t, mu, T, I, e, v_0, omega, tau):
    kappa = mu*np.sin(I)/(T**(1/3)*np.sqrt(1-e**2))
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



def get_telescope(filename, path='data/'):
    """ (str, str) -> (str)
    Returns the telescope of the star in filename.
    """
    return get_obs_info(filename, path)[1][12]


def find_files_for_telescope(telescope):
    """ (str) -> (list)
    Finds all the files in /data/ with with data pertaining to telescope.
    """
    telescope = re.sub(r'\s+', '', telescope)
    files = list_files('data')
    ans = []
    for file in files:
        if get_telescope(file) == telescope:
            ans.append(file)
    if len(ans) == 0:
        raise ValueError('No files with this telescope found.')
    
    return ans



class BinarySystem:
    """
    Represents a Binary System
    """
    def __init__(self, data=None, parameters=None, num_points=None):
        """ (self, pd.DataFrame, np.array(), int)


        If no inputs are given, then raise a ValueError. 
        If data is given, then use that data (assumes proper format).
        If parameters and num_points are given, then generates num_points radial velocity data (adds Gaussian noise).
        If only num_points is given, then generates random parameters and num_points radial velocity data (adds Gaussian noise).
        """
        if (data is None) and (parameters is None) and (num_points is None):
            raise ValueError('At least one initializing argument must be specified.')
        elif (data is None) and (parameters is None):
            # generating values for parameters
            self.mu = np.random.uniform(0, 1.246059e6) # in kg
            self.e = np.random.uniform(0, 1)
            self.I = np.random.uniform(-np.pi, np.pi)
            self.omega = np.random.uniform(0, np.pi/2)
            self.T = np.random.uniform(3282.3503, 3.46896e13) # in seconds
            self.tau = np.random.uniform(3282.3503, 3.46896e13) # in seconds
            self.v_0 = np.random.uniform(-10000, 10000) # in m/s
            
            # generating radial velocity data from 
            t = np.linspace(0, 3e8, num_points)
            self.time = t
            radial_velocities = radial_velocity(t, self.mu, self.T, 
                                                      self.I, self.e, self.v_0, self.omega, self.tau)
            # adding random Gaussian noise
            radial_velocities += np.random.normal(20, 10, len(t))
            self.radial_velocity = radial_velocities
            self.uncertainty = np.array([]) #how do we determine uncertainty when generating data?

            self.sampler = None
            self.samples = None
        elif (data is None):
            #generate random values from given parameters
            self.mu = param[0]
            self.e = param[1]
            self.I = param[2]
            self.omega = param[3]
            self.T = param[4]
            self.tau = param[5]
            self.v_0 = param[6]

            t = np.linspace(0, 3e8, num_points)
            self.time = t
            radial_velocities = radial_velocity(t, self.m, self.M, self.T, 
                                                      self.I, self.e, self.v_0, self.omega, self.tau)
            # adding random Gaussian noise
            radial_velocities += np.random.normal(20, 10, len(t))
            self.radial_velocity = radial_velocities
            self.uncertainty = np.array([]) #how do we determine uncertainty when generating data?
        
            self.sampler = None
            self.samples = None
        elif (num_points is None) and (parameters is None):
            #no known parameters
            self.mu = None
            self.e = None
            self.I = None
            self.omega = None
            self.T = None
            self.tau = None
            self.v_0 = None
            
            #data is given 
            self.data = data #dataframe with all data
            self.time = self.data[0].to_numpy()
            self.radial_velocity = self.data[1].to_numpy()
            self.uncertainty = self.data[2].to_numpy()

            self.sampler = None
            self.samples = None
        else:
            raise ValueError('Only certain combinations of inputs are accepted when defining a BinarySystem.')
    
    labels = [r"$\mu$", r"$e$", r"$I$", r"$\omega$", r"$T$", r"$\tau$", r"$v_0$"]

    def truth(self):
        """ (self) -> (np.array)
        If the parameters of the system are known, then return these. 
        """
        ans = [self.mu,
                self.e,
                self.I,
                self.omega,
                self.T,
                self.tau,
                self.v_0]
        return np.array(ans)
    
    @staticmethod
    def log_likelihood():
        pass

    @staticmethod
    def log_prior():
        pass

    @staticmethod
    def log_post():
        pass
    
    def initialize_mcmc(self, nwalkers, ndim=7):
        """ (self, int, int) -> ()
        Sets up the MCMC in self.sampler with nwalkers walkers and ndim dimension.
        """
        self.sampler = emcee.EnsembleSampler(nwalkers, ndim, BinarySystem.log_post, args=(self.time, self.radial_velocity, self.uncertainty))
    
    def run_mcmc(self, num_iter, nwalkers=None, ndim=7):
        """ (self, int, int, int) -> ()
        Runs the MCMC with num_iter iterations.
        If the MCMC is not initialized, sets it up with nwalkers walkers and ndim dimension.
        """
        if (self.sampler is None) and (nwalkers is None):
            raise ValueError('EnsembleSampler is not initialized, pass nwalkers as an argument or run self.initialize_mcmc.')
        elif self.sampler is None: #need to initialize the mcmc first
            self.initialize_mcmc(nwalkers, ndim)
        
        initial_pos = np.array(np.random.uniform(0, 1.246059e6), #mu
                                np.random.uniform(0, 1), #e
                                np.random.uniform(-np.pi, np.pi), #I
                                np.random.uniform(0, np.pi/2), #omega
                                np.random.uniform(3282.3503, 3.46896e13), #T
                                np.random.uniform(3282.3503, 3.46896e13), #tau
                                np.random.uniform(-10000, 10000),) #v_0
        
        self.sampler.run_mcmc(initial_pos, num_iter, progress=True)

    def get_samples(self, flat=False, thin=1, discard=0):
        """ (self, boolean, int, int) -> (np.array)
        Returns the chain of Monte Carlo steps.
        flat : determines whether or not to flatten the array.
        thin : take only every thin steps from the chain.
        discard : discards the first given number of steps. 
        """
        if self.sampler is None:
            raise ValueError('Must run the MCMC before accessing the samples.')
        return self.sampler.get_chain(flat=flat, thin=thin, discard=discard)
    
    def save_samples(self, flat=False, thin=1, discard=0):
        """ (self, boolean, int, int) -> ()
        Saves the chain of Monte Carlo steps to an instance attribute.
        flat : determines whether or not to flatten the array.
        thin : take only every thin steps from the chain.
        discard : discards the first given number of steps. 
        """
        if self.sampler is None:
            raise ValueError('Must run the MCMC before accessing the samples.')
        self.samples = self.sampler.get_chain(flat=flat, thin=thin, discard=discard)

    def plot_data(self, ls='', marker='.', color='blue'):
        plt.figure()

        plt.errorbar(self.time, self.data, self.uncertainty, ls=ls, marker=marker, color=color)

        plt.xlabel('Time (s)')
        plt.ylabel('Radial Velocity (m/s)')
        plt.show()

    def trace(self, ndim=7, flat=False, thin=1, discard=0):
        """ (self, int, boolean, int, int) -> ()
        Displays trace plot of the MCMC.
        """
        f, axes = plt.subplots(ndim, figsize=(10, 7), sharex=True)
        samples = self.get_samples(flat=flat, thin=thin, discard=discard)
        for i in range(ndim):
            ax = axes[i]
            ax.plot(samples[:, :, i], alpha=0.3)
            ax.set_xlim(0, len(samples))
            ax.set_ylabel(BinarySystem.labels[i])
            #ax.yaxis.set_label_coords(-0.1, 0.5)
        axes[-1].set_xlabel("Step number");

    def corner(self, thin=15, discard=100, quantiles=None):
        """ (self) -> ()
        Displays corner plot of the MCMC results
        """
        flat_samples = self.get_samples(flat=True, thin=thin, discard=discard)
        if quantiles is None:
            fig = corner.corner(flat_samples, labels=BinarySystem.labels);
        else:
            fig = corner.corner(flat_samples, labels=BinarySystem.labels, quantiles=[0.16, 0.5, 0.84]);

if __name__ == "__main__":
    import doctest
    doctest.testmod()
