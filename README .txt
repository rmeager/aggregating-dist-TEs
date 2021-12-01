README for Meager (2021) "Aggregating Distributional Treatment Effects"

1. OVERVIEW 


This repository contains all R scripts needed to generate the main results of the paper entitled “Aggregating Distributional Treatment Effects: A Bayesian Hierarchical Analysis of the Microcredit Literature.” In its default state the masterfile.R generates the paper's tables and figures from saved MCMC output which takes about 15 minutes. If rStan is installed, then the masterfile can be toggled to run the MCMC scripts from the raw data, which takes about 72 hours on a high performance computing server. 

2. DATA AVAILABILITY AND PROVENANCE STATEMENTS

This paper uses external data from 7 sources:

(a) Angelucci, Manuela, Dean Karlan, and Jonathan Zinman. 2015. "Microcredit Impacts: Evidence from a Randomized Microcredit Program Placement Experiment by Compartamos Banco." American Economic Journal: Applied Economics, 7 (1): 151-82.
DOI: 10.1257/app.20130537 
URL: https://www.aeaweb.org/articles?id=10.1257/app.20130537


(b) Attanasio, Orazio, Britta Augsburg, Ralph De Haas, Emla Fitzsimons, and Heike Harmgart. 2015. "The Impacts of Microfinance: Evidence from Joint-Liability Lending in Mongolia." American Economic Journal: Applied Economics, 7 (1): 90-122.
DOI: 10.1257/app.20130489
URL: https://www.aeaweb.org/articles?id=10.1257/app.20130489

(c) Augsburg, Britta, Ralph De Haas, Heike Harmgart, and Costas Meghir. 2015. "The Impacts of Microcredit: Evidence from Bosnia and Herzegovina." American Economic Journal: Applied Economics, 7 (1): 183-203.
DOI: 10.1257/app.20130272
URL: https://www.aeaweb.org/articles?id=10.1257/app.20130272


(d) Banerjee, Abhijit, Esther Duflo, Rachel Glennerster, and Cynthia Kinnan. 2015. "The Miracle of Microfinance? Evidence from a Randomized Evaluation." American Economic Journal: Applied Economics, 7 (1): 22-53.
DOI: 10.1257/app.20130533
URL: https://www.aeaweb.org/articles?id=10.1257/app.20130533

(e) Crépon, Bruno, Florencia Devoto, Esther Duflo, and William Parienté. 2015. "Estimating the Impact of Microcredit on Those Who Take It Up: Evidence from a Randomized Experiment in Morocco." American Economic Journal: Applied Economics, 7 (1): 123-50.
DOI: 10.1257/app.20130535
URL: https://www.aeaweb.org/articles?id=10.1257/app.20130535

(f) Karlan D, Zinman J. Microcredit in theory and practice: using randomized credit scoring for impact evaluation. Science. 2011 Jun 10;332(6035):1278-84. 
DOI: 10.1126/science.1200138.
URL: https://www.science.org/doi/10.1126/science.1200138

(g) Tarozzi, Alessandro, Jaikishan Desai, and Kristin Johnson. 2015. "The Impacts of Microcredit: Evidence from Ethiopia." American Economic Journal: Applied Economics, 7 (1): 54-89.
DOI: 10.1257/app.20130475
URL: https://www.aeaweb.org/articles?id=10.1257/app.20130475


All data are publicly available. 

3. DATASET LIST 

Wrangling the 7 sources above is formidable; you should not spend months of your life on it like I did. Moreover, it uses Stata, software to which not all users can afford access. Therefore, I have provided a single dataset file of the relevant components of the data used in my analysis, cleaned and processed. 

The data file is /data/microcredit_project_data.RData

As a codebook for interpreting the data file, please consult /data/import_organise_data_tracking_script.R 

This is the file I actually used to generate the dataset, but I do not intend for you to run this script; for it to run, you have to first run all of the authors' do files from the original sources, because in many cases they stitched their headline variables together within those do files themselves. This takes many days for a person to replicate because it has not been automated by any of the original papers' replication packages. Please use the RData file I have provided.

4. COMPUTATIONAL REQUIREMENTS

In the current repository I supply you with the output of the MCMC analysis run in rStan so that you do not need to use that package yourself. Therefore, the script /R code/masterfile.R can be run in R on any machine. This takes about 15-25 minutes on my current machine which runs MacOS Catalina 10.15.7: 
Processor Name:	Quad-Core Intel Core i5
  Processor Speed:	2 GHz
  Number of Processors:	1
  Total Number of Cores:	4
  L2 Cache (per Core):	512 KB
  L3 Cache:	6 MB
  Hyper-Threading Technology:	Enabled
  Memory:	32 GB

Re-running all the MCMC would take about 4-5 days on a high performance server, and likely 2-3 weeks on the machine above. 

5. DESCRIPTION OF PROGRAMS / CODE 

This enumeration can be found in /R code/masterfile.R 

6. INSTRUCTIONS TO REPLICATORS 

There are two ways to interact with the paper. I recommend the first. 

6.A. - Fully Automated (Taking the MCMC output stored in the subfolder "output" as given) 

1. Open the file aggregation-list-TEs.Rproj in R. 
2. Run source("/R code/masterfile.R")

6.B - May require manual adjustments (to configure rStan)

1. Install and check rstan's operation using the instructions https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
2. Open the file aggregation-list-TEs.Rproj in R. 
3. Open the file /R code/masterfile.R in R.
4. Go to line 11 of the masterfile and set rerun_MCMC <- TRUE 
5. Run source("/R code/masterfile.R"). 

Note: The reason this may not be possible to implement fully automated is that step 1 can sometimes require a small number of tailored interventions to ensure configuration is correct depending on the particular machine and OS. The reason for this is that rStan requires the machine to have on it a C++ compiler that can be called from R. Most machines do by default, but not all. 


7. TABLES AND FIGURES 

This enumeration can be found in /R code/masterfile.R 

8. REFERENCES 

Contained in section 2. 




