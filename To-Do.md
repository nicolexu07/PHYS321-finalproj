March 25th to 31st
- [ ] Research range on priors (William)
  - [ ] M
  - [ ] m
  - [ ] &mu;
  - [ ] T
  - [X] _e_
  - [X] I
  - [X] &omega;
  - [ ] &tau;
  - [ ] v<sub>0</sub> 
- [ ] Code simulations to be used to test MCMC
  - [ ] t
- [ ] Implement radial velocity function v(t)=&kappa;[cos(_f_+&omega;)+_e_*cos(&omega;)]+v<sub>0</sub>
  - [ ] Make function for &kappa;
  - [ ] Solve u-_e_*sin(u)=2&pi;(t-&tau;)/_T_ for u using scipy.optimize.fsolve
  - [ ] Make function for _f_=2*arctan[ sqrt((1+_e_)/(1-_e_)) * tan(u/2)]
- [ ] Download and load data (William & Nicole)
  - [X] Function(s) to download all .tbl files (William)
  - [X] Function to load radial velocity data from .tbl to pandas DataFrame or numpy array (Nicole)
  - [X] Function to load info from .tbl files (each file has lines about the observation such as star name, telescope used, etc. that come before the actual radial velocity data) (Nicole)
  - [ ] Improve load data function to take file path as parameter and format data
- [ ] Format data
  - [ ] Dates are in Julian days, find out how to convert to seconds
  - [ ] Double-check proper data types in each column  
- [ ] Test functions https://docs.python.org/3/library/unittest.html
  - [ ] Test loading data
  - [ ] Test loading info
  - [ ] Test formatting data

Apr. 1st to 8th
- [ ] Code MCMC
  * _The way I was thinking about this was to create a class that would hold all the data and methods we could need. The class would take radial velocity data as an input and store it as an instance attribute. It would also have an emcee EnsembleSampler object as an attribute, and the sampler's chain once it has been run as another attribute. It would also have methods to run the sampler, get the chain, plot trace plots, plot corner plots, plot data, plot example fits to the data, git MAP parameters, etc. Having all this contained in a class would make it easier to run MCMCs for many different files and store the results in a way that is easily accesible. -William_ 
  - [ ] Make new class called BinarySystem
  - [ ] \_\_init\_\_ method should take radial velocity data as input  
