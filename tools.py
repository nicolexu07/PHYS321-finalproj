import pandas as pd
import numpy as np
from io import StringIO
import re
import scipy
from scipy.optimize import fsolve, curve_fit
from astropy import constants as const
import emcee
import corner
import matplotlib.pyplot as plt
import os

from sympy import li


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
    """
    Function to use in scipy.optimize.fsolve
    """
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


def radial_velocity(t, mu, T, e, v_0, omega, tau):
    kappa = mu/(T**(1/3)*np.sqrt(1-e**2))
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


def get_uncertainties_telescope(telescope):
    """ (str) -> (list)
    Returns list of uncertainty values for a given telescope across all files in data
    """
    files = find_files_for_telescope(telescope)
    uncertainties = []
    for file in files:
        df = get_data(f'{file}')
        uncertainties += list(df.iloc[:, 2])
        
    return uncertainties


def plot_telescope_uncertainty_hist(telescope, bins=100):
    """(str, int, tup) -> ()
    Plots histogram of uncertainty values for a given telescope
    """
    unc = get_uncertainties_telescope(telescope)
    plt.hist(unc, bins=bins)
    plt.xlabel('Uncertainty Range')
    plt.ylabel('Number of Samples')
    plt.title(f'Histogram from Uncertainties of {telescope} Telescope')

    
    
    
def get_instrument(filename, path='data/'):
    """ (str, str) -> (str)
    Returns the instrument of the star in filename.
    """
    return get_obs_info(filename, path)[1][13]


def find_files_for_instrument(instrument):
    """ (str) -> (list)
    Finds all the files in /data/ with with data pertaining to instrument.
    """
    instrument = re.sub(r'\s+', '', instrument)
    files = list_files('data')
    ans = []
    for file in files:
        if get_instrument(file) == instrument:
            ans.append(file)
    if len(ans) == 0:
        raise ValueError('No files with this instrument found.')
    
    return ans


def get_uncertainties_instrument(instrument):
    """ (str) -> (list)
    Returns list of uncertainty values for a given instrument across all files in data
    """
    files = find_files_for_instrument(instrument)
    uncertainties = []
    for file in files:
        df = get_data(f'{file}')
        uncertainties += list(df.iloc[:, 2])
        
    return uncertainties


def plot_instrumental_uncertainty_hist(instrument, bins=100):
    """(str, int, tup) -> ()
    Plots histogram of uncertainty values
    """
    unc = get_uncertainties_instrument(instrument)
    plt.hist(unc, bins=bins)
    plt.xlabel('Uncertainty Range')
    plt.ylabel('Number of Samples')
    plt.title(f'Histogram from Uncertainties of {instrument} Instrument')    
    
    
    
def plot_uncertainty_for_file(filename):
    """ (str) -> ()
    Plots radial velocity vs uncertainty for a dataset 'filename'
    """
    df = get_data(filename)
    radial_velocities = list(df.iloc[:, 1])
    uncertainties = list(df.iloc[:, 2])
    
    star_id = get_star_id(filename)
    
    plt.scatter(radial_velocities, uncertainties)
    plt.xlabel('Uncertainty (m/s)')
    plt.ylabel('Radial Velocity (m/s)')
    plt.title(f'Relationship between Radial Velocity and {star_id}')
    
    
def plot_uncertainty_correlation():
    """ () -> ()
    Plots radial velocity vs associated uncertainty across all systems in data
    """

    files = list_files('data')
    radial_velocities = []
    uncertainties = []
    for file in files[:200]:
        df = get_data(file)
        radial_velocities += list(df.iloc[:, 1])
        uncertainties += list(df.iloc[:, 2])

    plt.scatter(radial_velocities, uncertainties, s=0.1)
    plt.xlim(-100, 100)
    plt.ylim(0, 50)  
    
    
def plot_uncertainty_correlation_telescope(telescope):
    """ (str) -> ()
    Plots radial velocity vs uncertainty for all data points of a given telescope
    """
    files = find_files_for_telescope(telescope)
    radial_velocities = []
    uncertainties = []
    for file in files[:200]:
        df = get_data(file)
        radial_velocities += list(df.iloc[:, 1])
        uncertainties += list(df.iloc[:, 2])

    plt.scatter(radial_velocities, uncertainties)
    
    
    
def gen_uncertainty(radial_velocities, instrument):
    """ (np.array, str) -> (np.array)
    Returns an array of uncertainty values associated with each 
    radial velocity value in radial_velocities array
    """
    # the uncertainty distribution we will draw from 
    uncertainty_dist = get_uncertainties_instrument(instrument)
    
    uncertainties = []
    #assigning uncertainties by drawing randomly from distribution 
    for v in radial_velocities:
        index = np.random.randint(0, len(uncertainty_dist)-1)
        uncertainties.append(uncertainty_dist[index])

    return np.array(uncertainties)
    
    
def plot_ellipse(e, a=1, points=1000, title=None):
    """ (num, num, int, str) -> ()
    Plots an ellipse with eccentricity e and semimajor axis a.
    """
    b = a*np.sqrt(1-e**2)
    theta = np.linspace(0, 2*np.pi, points)
    r = a*b/np.sqrt((b*np.cos(theta))**2+(a*np.sin(theta))**2)
    x = r*np.cos(theta)
    y = r*np.sin(theta)

    plt.figure(figsize=(8,8))
    plt.plot(x, y)

    plt.xlabel('x')
    plt.ylabel('y')
    if title is not None:
        plt.title(title)
    plt.ylim(-a, a)
    plt.xlim(-a, a)
    plt.show()


class BinarySystem:
    """
    Represents a Binary System
    """
    def __init__(self, data=None, parameters=None, num_points=None, tmax=None):
        """ (self, pd.DataFrame, np.array(), int, num)

        If no inputs are given, raise a ValueError. 

        If data is given, use that data (assumes proper format).

        If parameters and num_points are given, generates num_points radial velocity data between times 0 and tmax (adds noise). 
        The parameters ought to have a certain order:
        [mu, e, omega, T, tau, v_0]
        T and tau are expected to be in units of seconds, this method handles taking their logarithm.

        If only num_points is given, generates random parameters and num_points radial velocity data between times 0 and tmax (adds noise).

        If tmax is not specified, defaults to 3e8.
        """
        if (data is None) and (parameters is None) and (num_points is None):
            raise ValueError('At least one initializing argument must be specified.')
        elif (data is None) and (parameters is None):
            # assigning a random instrument 
            # selecting a random file from directory
            all_files = list_files('data')
            index = np.random.randint(0, len(all_files)-1)
            # finding associated instrument 
            self.instrument = get_instrument(all_files[index])
            
            # generating values for parameters
            self.mu = np.random.choice(np.linspace(-1.246059e6, 1.246059e6, 20000)) # in kg
            self.e = np.random.uniform(0, 1)
            self.omega = np.random.uniform(0, np.pi*2)
            self.log_T = np.random.uniform(3.516184928, 13.54019929) # in seconds
            self.log_tau = np.random.uniform(3.516184928, self.log_T) # in seconds
            self.v_0 = np.random.choice(np.linspace(-9000, 9000, 10000)) # in m/s
            
            # generating radial velocity data from parameters
            if tmax is None:
                t = np.linspace(0, 3e8, num_points)
            else:
                t = np.linspace(0, tmax, num_points)
            self.time = t
            radial_velocities = radial_velocity(t, self.mu, 10**self.log_T, 
                                                      self.e, self.v_0, self.omega, 10**self.log_tau)

            self.uncertainty = np.minimum(50, gen_uncertainty(radial_velocities, self.instrument))
            
            # we should add noise to model's radial velocity (based on generated parameters)
            # to get our final simulated data for radial velocity 
            
            # add Gaussian noise with stdev of 1 which we scale with our uncertainty values
            noise = np.random.normal(0, 1, size=len(radial_velocities)) * self.uncertainty
            self.radial_velocity = radial_velocities + noise
                
            self.sampler = None
            self.samples = None
        elif (data is None):
            # assigning a random instrument 
            # selecting a random file from directory
            all_files = list_files('data')
            index = np.random.randint(0, len(all_files)-1)
            # finding associated instrument 
            self.instrument = get_instrument(all_files[index])
            
            #generate random values from given parameters
            self.mu = parameters[0]
            self.e = parameters[1]
            self.omega = parameters[2]
            self.log_T = np.log10(parameters[3])
            self.log_tau = np.log10(parameters[4])
            self.v_0 = parameters[5]

            if tmax is None:
                t = np.linspace(0, 3e8, num_points)
            else:
                t = np.linspace(0, tmax, num_points)
            self.time = t
            radial_velocities = radial_velocity(t, self.mu, 10**self.log_T, 
                                                      self.e, self.v_0, self.omega, 10**self.log_tau)

            self.uncertainty = gen_uncertainty(radial_velocities, self.instrument)
            
            # we should add noise to model's radial velocity (based on generated parameters)
            # to get our final simulated data for radial velocity 
            
            # add Gaussian noise with stdev of 1 which we scale with our uncertainty values
            noise = np.random.normal(0, 1, size=len(radial_velocities)) * self.uncertainty
            self.radial_velocity = radial_velocities + noise
        
            self.sampler = None
            self.samples = None
        elif (num_points is None) and (parameters is None):
            #no known parameters
            self.mu = None
            self.e = None
            self.omega = None
            self.log_T = None
            self.log_tau = None
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
    
    labels = [r"$\mu$", r"$e$", r"$\omega$", r"$log_{10}T$", r"$log_{10}\tau$", r"$v_0$"]

    def truth(self):
        """ (self) -> (np.array)
        If the parameters of the system are known, then return these. 
        """
        ans = [self.mu,
                self.e,
                self.omega,
                self.log_T,
                self.log_tau,
                self.v_0]
        return np.array(ans)
    
    
    
    def log_likelihood(self, theta):
        """ (self, np.array) -> (num)
        The log probability of measuring self.data given a system with parameters theta. The probability is determined
        by a Gaussian distribution where each data point is assumed to be independent of the others.         
        """
        # note theta = [mu, e, I, omega, T, tau, v_0]
        mu, e, omega, log_T, log_tau, v_0 = theta 
        model = radial_velocity(self.time, mu, 10**log_T, e, v_0, omega, 10**log_tau)
        
        return -0.5*np.sum((self.radial_velocity - model)**2 / self.uncertainty**2 + np.log(2*np.pi*self.uncertainty**2))
        

    def log_prior(self, theta):
        """ (self, np.array) -> (np.inf/float)
        tau must be less than T, this inequality follows to the logarithms of these values.
        """
        mu, e, omega, log_T, log_tau, v_0 = theta 
        
        # based on prior research on allowed values 
        if -1.246059e6> mu or 1.246059e6 <= mu:
            return -np.inf
        elif 3.516184928 > log_T or 13.54019929 < log_T:
            return -np.inf
        elif 0 > e or 1 <= e:
            return -np.inf
        elif -100000 > v_0 or 100000 < v_0:
            return -np.inf
        elif 0 > omega or np.pi*2 < omega:
            return -np.inf
        elif 3.516184928 > log_tau or log_T < log_tau:
            return -np.inf
        return 0


    def log_post(self, theta):
        """ (self, np.array) -> (np.inf/num)
        log(posterior) = log(likelihood) + log(prior) + const.
        """
        lp = self.log_prior(theta)
        if not np.isfinite(lp):
            return -np.inf
        return lp + self.log_likelihood(theta)
    
    
    
    def initialize_mcmc(self, nwalkers, ndim=6):
        """ (self, int, int) -> ()
        Sets up the MCMC in self.sampler with nwalkers walkers and ndim dimension.
        """
        self.sampler = emcee.EnsembleSampler(nwalkers, ndim, self.log_post)
        self.nwalkers = nwalkers

    def run_mcmc(self, num_iter, nwalkers=None, ndim=6):
        """ (self, int, int, int) -> ()
        Runs the MCMC with num_iter iterations and random initial conditions.
        If the MCMC is not initialized, sets it up with nwalkers walkers and ndim dimension.
        """
        if (self.sampler is None) and (nwalkers is None):
            raise ValueError('EnsembleSampler is not initialized, pass nwalkers as an argument or run self.initialize_mcmc.')
        elif self.sampler is None: #need to initialize the mcmc first
            self.initialize_mcmc(nwalkers, ndim)
        
        # note [mu, e, I, omega, T, tau, v_0] = theta
        initial_pos = [] 
        for i in range(self.nwalkers):
            temp = [np.random.uniform(-1.246059e6, 1.246059e6), #mu
                    np.random.uniform(0, 1), #e
                    np.random.uniform(0, np.pi*2), #omega
                    np.random.uniform(3.516184928, 13.54019929)] #log_T
            temp.append(np.random.uniform(3.516184928, temp[-1])) #log_tau
            temp.append(np.random.uniform(-9999, 9999))#v_0
            initial_pos.append(temp)
        initial_pos = np.array(initial_pos)

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

    def plot_data(self, ls='', marker='.', color='blue', title=None, real=False):
        """ (self, str, str, str, str) -> ()
        Shows a plot of the BinarySystem's radial velocity data as a function of time.
        If real is set to True, will display the actual curve generated by the true parameters if possible.
        """
        plt.figure()

        plt.errorbar(self.time, self.radial_velocity, self.uncertainty, 
                        ls=ls, marker=marker, color=color, elinewidth=0.5, label='Data')

        if real:
            t = np.linspace(0, max(self.time))
            y = radial_velocity(t, self.mu, 10**self.log_T, self.e, self.v_0, self.omega, 10**self.log_tau)
            plt.plot(t, y, label='Actual', color='red')
            plt.legend()
        if title is not None:
            plt.title(title)
        plt.xlabel('Time (s)')
        plt.ylabel('Radial Velocity (m/s)')
        plt.show()

    def trace(self, ndim=6, flat=False, thin=1, discard=0, title=None):
        """ (self, int, boolean, int, int) -> ()
        Displays trace plot of the MCMC.
        """
        f, axes = plt.subplots(ndim, figsize=(10, 15), sharex=True)
        samples = self.get_samples(flat=flat, thin=thin, discard=discard)
        if title is not None:
            plt.title(title)
        for i in range(ndim):
            ax = axes[i]
            ax.plot(samples[:, :, i], alpha=0.3)
            ax.set_xlim(0, len(samples))
            ax.set_ylabel(BinarySystem.labels[i])
            #ax.yaxis.set_label_coords(-0.1, 0.5)
        axes[-1].set_xlabel("Step number");


    def corner(self, thin=1, discard=250, quantiles=None):
        """ (self, int, int, list) -> ()
        Displays corner plot of the MCMC results
        """
        flat_samples = self.get_samples(flat=True, thin=thin, discard=discard)
        if quantiles is None:
            fig = corner.corner(flat_samples, labels=BinarySystem.labels);
        else:
            fig = corner.corner(flat_samples, labels=BinarySystem.labels, quantiles=quantiles);

    def param_hist(self, param, bins, limits=None, flat=True, thin=1, discard=250):
        """ (self, int, int, tuple, boolean, int, int) -> (np.array, np.array)
        Returns the bin edges and corresponding counts of a histogram of the parameter of the system 
        index by the param input.
        """
        if limits is None:
            val = self.get_samples(flat=flat, thin=thin, discard=discard)
            ans = []
            for entry in val:
                ans.append(entry[param])
            return np.histogram(ans, bins=bins)
        else:
            val = self.get_samples(flat=flat, thin=thin, discard=discard)
            ans = []
            for entry in val:
                ans.append(entry[param])
            return np.histogram(ans, bins=bins, range=limits)

    def param_hist_plot(self, param, bins, limits=None, flat=True, thin=1, discard=250, actual_value=None, title=None):
        """ (self, int, int, tuple, boolean, int, int, num) -> ()
        Plots the histogram of the parameter indexed by param.
        """
        val = self.get_samples(flat=flat, thin=thin, discard=discard)
        ans = []
        for entry in val:
            ans.append(entry[param])
        plt.figure()
        plt.hist(ans, bins=bins)
        if actual_value is not None:
            count, temp = np.histogram(ans, bins=bins)
            plt.vlines(actual_value, ymin=0, ymax=1.05*max(count), color='black')
        if title is not None:
            plt.title(title)
        if limits is not None:
            plt.xlim(limits)
        plt.xlabel(BinarySystem.labels[param])
        plt.ylabel("Count")
        plt.show()
        """
        else:
            val = self.get_samples(flat=flat, thin=thin, discard=discard)
            ans = []
            for entry in val:
                ans.append(entry[param])
            plt.figure()
            plt.hist(ans, bins=bins)
            plt.vlines(actual_value, ymin=0, ymax=1.25*max(ans), color='black')
            if title is not None:
                plt.title(title)
            plt.xlabel(BinarySystem.labels[param])
            plt.ylabel("Count")
            plt.show()
        """

    def plot_samples(self, num_samples, thin=1, discard=250, alpha=0.2, title=None, least_squares=False, real=False):
        """ (self, int, int, int, int, str, boolean, boolean) -> ()
        Plots num_samples sample curves from the results of the MCMC if least_squares is False. If least_squares is  
        True then only plots the result from the nonlinear least-squares fitting of radial_velocity to the data.
        """
        plt.figure()

        t = np.linspace(min(self.time), max(self.time), 1000)

        if least_squares:
            param, param_cov = self.least_squares_fit()
            mu, e, omega, T, tau, v_0 = param
            y = radial_velocity(t, mu, T, e, v_0, omega, tau)
            plt.plot(t, y, color='yellow', label='Fit')
       
       #plots the samples
        samples = self.get_samples(flat=True, thin=thin, discard=discard)
        indices = np.random.randint(0, len(samples), size=num_samples)
        for i in range(len(indices)-1):
            mu, e, omega, log_T, log_tau, v_0 = samples[indices[i]]
            y = radial_velocity(t, mu, 10**log_T, e, v_0, omega, 10**log_tau)
            plt.plot(t, y, color='red', alpha=alpha)
        #last one separate to make the legend proper
        mu, e, omega, log_T, log_tau, v_0 = samples[indices[-1]]
        y = radial_velocity(t, mu, 10**log_T, e, v_0, omega, 10**log_tau)
        plt.plot(t, y, color='red', alpha=alpha, label='Samples')

        if real:
            t = np.linspace(0, max(self.time))
            y = radial_velocity(t, self.mu, 10**self.log_T, self.e, self.v_0, self.omega, 10**self.log_tau)
            plt.plot(t, y, label='Actual', color='green')
            
        plt.errorbar(self.time, self.radial_velocity, self.uncertainty, 
                        ls='', marker='.', color='blue', elinewidth=0.5, label='Data')

        if title is not None:
            plt.title(title)
        plt.xlabel('Time (s)')
        plt.legend()
        plt.ylabel('Radial Velocity (m/s)')
        plt.ylim(np.amin(self.radial_velocity-self.uncertainty)-5, np.amax(self.radial_velocity+self.uncertainty)+5)
        plt.show()


    def least_squares_fit(self):
        """ (self) -> (np.array, np.array)
        Uses scipy.optimize.curve_fit to fit the radial_velocity function to the data.
        """
        return curve_fit(radial_velocity, self.time, self.radial_velocity, sigma=self.uncertainty)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
