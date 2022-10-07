# Bayesian Statistics.

Repository for the Assignment of Bayesian Statistics (2022). Utrecht University - Methodology and Statistics. 

It contains 3 files:    
1- manusccript.pdf: report to be graded.     
2- df.csv: the data used for the analysis.  
3- Villalobos_bayes.R: script.  

About Sofía_Villalobos.R:  
Lines 1-15 need to be run to load data, packages and define variables.  
Lines 16-213 specify informative and informative priors, number of iterations, initial values,  
and runs and storages the results for the Gibbs sampler with the Metropolis Hasting step for B1.  
Lines 214-264 creates trace plots to check convergence of chains for informative priors.   
First 900000 iterations of each chain are deleted and check new plots.  
Lines 265-276 Generates autocorrelation plots to evaluate convergence  
Lines 277-307 compute Monte Carlo error to check convergence   
Lines 308-375 computes predictive posterior distribution and check normality assumption  
with bayesian p-value  
Lines 376-385 create tables with Mean, SD and Credible Intervals for both models.  
This tables can be observed in line 385 "output_inf" and in line 385 "output_unin".  
Lines 397-448 computes Deviation Information Criteria  
Lines 450-463 consists in hypotheses testing with Bayes Factor using the package  
BAIN. Results were obtain with set.seed(43838), and set.seed(78874) to test stability of results.  
