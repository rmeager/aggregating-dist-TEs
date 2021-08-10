README for the Meager Distributional Effects Paper Repository

This repository contains all scripts needed to generate the main results of the paper entitled “Aggregating Distributional Treatment Effects: A Bayesian Hierarchical Analysis of the Microcredit Literature.” It does not contain all scripts for all Appendices, but will be updated with the relevant Appendices during the revision process. 

The repo contains: 
1.  A data subfolder, with the combined microcredit data file I used throughout the analysis.
2. A stan-models subfolder, with all the model scripts written in the stan language.
3. An R code subfolder with all R scripts for the analysis and graphics/tables. 
4. An output folder with the manuscript tex file. 

The project, which began when I was in graduate school, has a somewhat bespoke R script structure. To reproduce the paper one need only interact with the R code subfolder and then compile the tex file in the output folder. However, the R code subfolder contains two distinct types of files that must be used as follows. 

1. ANALYSIS AND RESULTS GENERATION SCRIPTS 

All R files beginning with “tailored-hierarchical…” run the hierarchical models that underlie the main results of the papers and produce the output processed by the graphics scripts. 

All R files beginning with “tailored-no-pooling…” and “tailored-full-pooling…” run these versions of the lognormal model for comparison purposes, informing table 2 for profit and results comparison more broadly.  

There are 2 important features of the analysis scripts above, worth noting: 

(A) The above files call the stan models in the stan-models subfolder into R. RStan performs Hamiltonian Monte Carlo simulations to characterise posterior distributions. It requires the machine running the operation to have a C++ compiler that R can talk to installed on it. This is common now on modern machines but may not be common on remote servers, and the code will not run otherwise. Please see “Getting started in RStan” for information on this: https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started

(B) These files fit hierarchical distributional models for data sets variously comprised of 10,000-40,000 households (depending on the outcome variable), and as such, each model takes several hours to run even on a powerful server. The paper is composed of the results of more than 20 of these model runs. Running this entire repository will take several days at minimum. I have therefore additionally provided the model output from each file as RData files in the Dropbox output folder of this repository here: 

https://www.dropbox.com/sh/9n1qgjeokxg9tun/AADe_Mcg7v95ldadD3Tdk0Vva?dl=0

2. GRAPHICS AND TABLE PROCESSING SCRIPTS

The file “graphics and computation of quantiles lognormal models business variables.R” generates the business variable results shown in figures 1, 2 and prepares .rds objects for input into table 2.  

The file “graphics and computation of quantiles lognormal consumption variables.R” generates the consumption variable results shown in figures 1, 2 and prepares .rds objects for input into table 1. 

The file “graphics for no pooling models and site specific output lognormal.R” generates the no pooling table results for profit, and some deprecated figures for all business variables which you may enjoy should you wish to see the no pooling results with the partial pooling results side by side. 

The file “tables-full-results-profit-consumption.R” generates the results displayed in tables 1 and 2. It can only be run AFTER the first 3 graphics files listed above, because it exploits the subtables created and stored as .RDS objects by those scripts. 

The file “Bayesian tests of equality.R” generates table 3. 

The file “posterior-predictive-graphic-comparison-lognormal-pareto-full-sim.R” generates figure 3. 

The file “graphics and computation of quantiles composite tail model.R” generates panel 1 of figure 4. The file “graphics and computation of quantiles from PLN model.R” generates panel 2 of figure 4. Panel 3 of figure 4 is the profit graphic from figure 1. 

The file “graphics and computation of quantiles from tailored hierarchical pdf lognormal models consumption types pb split.R” generates figure 5. 

The file “graphics and computation of quantiles from tailored hierarchical pdf lognormal models pb split.R” generates figure 6.

