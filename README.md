# Parameter Inference for Binary Systems
Generally, exoplanets can't be directly observed due to the high brightness contrast between the two objects, and limited telescope resolution. However, it's relatively easy to obtain radial velocity measurements by observing Doppler Shift. These radial velocity measurements can be used to recover information about the binary system, including orbital period, eccentricity, and mean velocity of the system's center of mass. We modeled the posterior (parameters given data) with Bayes' Theorem and sampled this posterior distribution with a Markov Chain Monte Carlo (MCMC). We found that the orbital parameters can be accurately recovered for both simulated data (including noise) as well as for real radial velocity data from existing binary systems. 

## Relevant Files
__tools.py__ contains code for functions to format data, create simulations, run the MCMC, and plot our results.

__Exoplanet_Parameter_Detection.ipynb__ calls functions from tools.py to run our simulations and MCMC and visualize results.

__Radial_Velocity_MCMC_Paper.pdf__ our final report submission. 

## Other Files (that you don't need to look at):

__data__ is a directory containing radial velocity data from the CalTech Exoplanet Archive. 

__import_tar_files.py__ Code used to import files in the 'data' directory from an external source. 

__test_data.tbl__ Test data used in one of our unit tests. 

__To-Do.md__ was used to track our progress and assign responsibilities. Can be a more detailed reference of the distribution of work.

__research.md__ was just for us to keep all our research in one place! All the relevant information in this file is also included in our paper.

## Work Distribution

Nicole Xu: Implementation of radial velocity function, posterior, likelihood, priors. Data simulation (exploring uncertainty, generating parameter values, adding noise). Report writing (abstract, methods, conclusion).

William Pugsley: Retrieving, formating, and searching through data files. Graphing and emcee methods for BinarySystem class. Report writing (introduction, methods, results, conclusion).
