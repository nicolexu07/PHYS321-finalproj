# Parameter Inference for Binary Systems

__tools.py__ contains code for functions to format data, create simulations, run the MCMC, and plot our results.

__Exoplanet_Parameter_Detection.ipynb__ calls functions from tools.py to run our simulations and MCMC and visualize results.

## Work Distribution Summary

Nicole Xu: Implementation of radial velocity function, posterior, likelihood, priors. Data simulation (exploring uncertainty, generating parameter values, adding noise). Report writing (abstract, methods, conclusion).

William Pugsley: Retrieving, formating, and searching through data files. Graphing and emcee methods for BinarySystem class. Report writing (introduction, methods, results, conclusion).

## Other Files (that you don't need to look at):

__data__ is a directory containing radial velocity data from the CalTech Exoplanet Archive. 

__import_tar_files.py__ Code used to import files in the 'data' directory from an external source. 

__test_data.tbl__ Test data used in one of our unit tests. 

__To-Do.md__ was used to track our progress and assign responsibilities. Can be a more detailed reference of the distribution of work.

__research.md__ was just for us to keep all our research in one place! All the relevant information in this file is also included in our paper.

