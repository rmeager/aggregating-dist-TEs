# Load data for Microcredit Meta-Analysis Project 
# Rachael Meager
# First Version: August 28th, 2014
# This Version: Feb 2015

### Notes ###

# This code loads variables from .dta files from 7 RCTs of microcredit
# All data was publicly available from the authors and journals, without whom this project would not be possible. Hooray for open science!
# I use "i" as the general index for household-level variables and "k" as the index for site-level variables


# A note about NAs: I have tried to keep NAs as NAs, except in TWO CASES. 
# 1. Profit/Business expenses/Revenues. An NA in here can indicate not operating a business and we want to capture that in a zero.
# 2. Variables constructed by summing other variables should not be NA unless *ALL* of the constituent variables are NA.


### Preliminaries ###

rm(list = ls())

chooseCRANmirror(ind=90)

installer <- FALSE
if(installer == TRUE){
  install.packages('foreign')
  install.packages('Hmisc')
  install.packages('xtable')
  install.packages('coda')
  install.packages('stargazer')
  install.packages('dummies')
  install.packages('zoo')
  install.packages('lmtest')
  install.packages('plm')
  install.packages('sandwich')
  install.packages('Matrix')
  install.packages("ggplot2")
  install.packages("gridExtra")
  install.packages("taRifx")
install.packages("rstan")}


library(dummies)
library(stargazer)
library(foreign)
library(Hmisc)
library(xtable)
library(coda)
library(zoo)
library(plm)
library(lmtest)
library(sandwich)
library(ggplot2)
library(gridExtra)
library(Matrix)
library(taRifx)
library(rstan)
# stan2codafit function declaration

# function to convert stan to coda
stan2coda <- function(stanfit) {
  mcmc.list(lapply(1:ncol(stanfit), function(x) mcmc(as.array(stanfit)[,x,])))
}


### Load and Organise the Data ###

# Load Angelucci et al 2015 #

# If unit is monetary, it is pesos unless otherwise specified

raw_data_angelucci_2015 <- read.dta("//bbkinghome/rmeager/winprofile/mydocs/Research work/bayesian meta analysis/microcredit analysis/microcredit rct data/angelucci_et_al_2015.dta")

angelucci_treatment <- raw_data_angelucci_2015$Treatment # not sure how much compliance/defiance this was...
angelucci_treatmentassigned <- raw_data_angelucci_2015$BTreatment # I have asked the authors to confirm interpretation of this variable but did not receive a reply.
angelucci_profit_non_imputed <- raw_data_angelucci_2015$Q10_9_toprof # fortnight
angelucci_profit <- raw_data_angelucci_2015$Q10_9_toprof #fortnight
angelucci_profit[is.na(angelucci_profit)] <- 0
angelucci_revenues <- raw_data_angelucci_2015$Q10_9_totrev #fortnight
angelucci_revenues_non_imputed <- raw_data_angelucci_2015$Q10_9_totrev #fortnight 
angelucci_revenues[is.na(angelucci_revenues)] <- 0
angelucci_expenditures <- raw_data_angelucci_2015$Q10_8_totexp #fortnight 
angelucci_expenditures_non_imputed <- raw_data_angelucci_2015$Q10_8_totexp
angelucci_expenditures[is.na(angelucci_expenditures)] <- 0
angelucci_assets <- raw_data_angelucci_2015$Q9_totvalue
angelucci_income <- raw_data_angelucci_2015$Q13_3_amount + raw_data_angelucci_2015$Q13_jobincome + raw_data_angelucci_2015$remitandtrans + raw_data_angelucci_2015$Q13_5_amount_mo
angelucci_businessincome <- raw_data_angelucci_2015$Q13_3_amount
angelucci_businessincome_non_imputed <- raw_data_angelucci_2015$Q13_3_amount #I should check how this differs from profit though
angelucci_businessincome[is.na(angelucci_businessincome)] <- 0
angelucci_existingbusiness <- as.numeric(raw_data_angelucci_2015$business_before)
angelucci_temptation <- (7/30)*raw_data_angelucci_2015$Q7_tempt # now this is weekly too

# consumption is all WEEKLY 
angelucci_consumption_source_dataframe <- data.frame(angelucci_temptation, raw_data_angelucci_2015$foodconsump, raw_data_angelucci_2015$Q18_2_1_amount_wk, raw_data_angelucci_2015$Q18_2_2_amount_wk, raw_data_angelucci_2015$Q18_2_3_amount_wk, raw_data_angelucci_2015$Q7_nonfood_nondur)
angelucci_consumption <- rowSums( angelucci_consumption_source_dataframe, na.rm = TRUE)  

angelucci_district <- raw_data_angelucci_2015$supercluster_xi
angelucci_womenempower1 <- raw_data_angelucci_2015$Q2_3_sing_any
angelucci_womenempower2 <- raw_data_angelucci_2015$Q16_tot_somepower_cond
angelucci_eduspend <- raw_data_angelucci_2015$Q18_2_1_amount_wk
angelucci_healthspend <- raw_data_angelucci_2015$Q18_2_2_amount_wk
angelucci_foodspend <- raw_data_angelucci_2015$foodconsump

angelucci_urban <- raw_data_angelucci_2015$urban
angelucci_indicator <- rep(1, length(angelucci_treatment))
angelucci_cluster_id <- raw_data_angelucci_2015$cluster 

angelucci_takeup <- raw_data_angelucci_2015$Q21_3_comp
angelucci_loan_amount <- raw_data_angelucci_2015$Q21_5_comp

angelucci_informal_loan_indicator <- raw_data_angelucci_2015$Q21_3_informal2
angelucci_informal_loan_amount <- raw_data_angelucci_2015$Q21_5_informal2

angelucci_other_formal_loan_indicator <- raw_data_angelucci_2015$Q21_3_othformal
angelucci_other_formal_loan_amount <- raw_data_angelucci_2015$Q21_5_othformal
  
angelucci_bank_loan_indicator <- raw_data_angelucci_2015$Q21_3_bank
angelucci_bank_loan_amount <- raw_data_angelucci_2015$Q21_5_bank
  
angelucci_other_MFI_loan_indicator <- raw_data_angelucci_2015$Q21_3_othmfi
angelucci_other_MFI_loan_amount <- raw_data_angelucci_2015$Q21_5_othmfi

angelucci_unspecified_other_loan_indicator <- raw_data_angelucci_2015$Q21_3_oth #this may include family loans
angelucci_unspecified_other_loan_amount <- raw_data_angelucci_2015$Q21_5_oth

angelucci_targetwomen <- 1
angelucci_promotion <- 1
angelucci_grouploan <- rep(1, length(angelucci_treatment))
angelucci_collateralised <- 0 
angelucci_averageloansize <- 538 #in USD (but not PPP)
angelucci_APR <- 100 
angelucci_timegap <- 27 #months
angelucci_currentmarket <- 2 #scale is integers 0-1-2-3 marking no/little/some/full existing formal credit market activity 
angelucci_village_rand <- 1
angelucci_loanterm <- 14 #months
angelucci_loansize_percentincome <- 6 #percent of income in country per year


# Load Attanasio et al 2015 # 

# have to make a choice here about whether to use the _r  variable (just respondent) OR to use _j (joint with partner) ASK ESTHER
# I have chosen _j for now but perhaps we should change it 
raw_data_attanasio_2015 <- read.dta( "//bbkinghome/rmeager/winprofile/mydocs/Research work/bayesian meta analysis/microcredit analysis/microcredit rct data/attanasio et al 2015/Analysis files/data/attanasio_processed_for_rm_analysis.dta") 
# NOW DROP EVERYONE WHO IS BASELINE HERE 
raw_data_attanasio_2015$followup <- as.numeric(raw_data_attanasio_2015$followup)
raw_data_attanasio_2015 <- raw_data_attanasio_2015[raw_data_attanasio_2015$followup==2,]

additional_loans_data_attanasio_2015 <- read.dta("//bbkinghome/rmeager/winprofile/mydocs/Research work/bayesian meta analysis/microcredit analysis/microcredit rct data/attanasio et al 2015/Analysis files/data/Followup/f_Section K -Loans.dta")
if (identical(additional_loans_data_attanasio_2015$hhid, raw_data_attanasio_2015$hhid) == FALSE){stop("DANGER! IDs in Attanasio data files are NOT aligned!")}



attanasio_treatment <- raw_data_attanasio_2015$treatment 
attanasio_treatment <- as.numeric(raw_data_attanasio_2015$treatment)
# we need to split out assignment to individual versus group loans
attanasio_grouploan <- rep(0, length(attanasio_treatment))
attanasio_grouploan[attanasio_treatment == 3] <- 1
attanasio_individualloan <- rep(0, length(attanasio_treatment))
attanasio_individualloan[attanasio_treatment == 2] <- 1
# I now transform it into a binary indicator
attanasio_treatment[attanasio_treatment == 1] <- 0
attanasio_treatment[attanasio_treatment != 0] <- 1
attanasio_profit <- raw_data_attanasio_2015$profit_j
attanasio_profit_non_imputed <- raw_data_attanasio_2015$profit_j
attanasio_profit[is.na(attanasio_profit)] <- 0
attanasio_revenues <- raw_data_attanasio_2015$rev_j 
attanasio_revenues_non_imputed <- raw_data_attanasio_2015$rev_j 
attanasio_revenues[is.na(attanasio_revenues)] <- 0
attanasio_expenditures <- raw_data_attanasio_2015$exp_j
attanasio_expenditures_non_imputed <- raw_data_attanasio_2015$exp_j
attanasio_expenditures[is.na(attanasio_expenditures)] <- 0
attanasio_assets <- raw_data_attanasio_2015$assets_all
attanasio_businessassets <- raw_data_attanasio_2015$BAI
attanasio_businessassets_non_imputed <- raw_data_attanasio_2015$BAI
attanasio_businessassets[is.na(attanasio_businessassets)] <- 0
attanasio_assets <- raw_data_attanasio_2015$assets_all
attanasio_income <- raw_data_attanasio_2015$hhincome
attanasio_consumption <- raw_data_attanasio_2015$totalc #monthly
attanasio_saving <- raw_data_attanasio_2015$totalsav_hh
attanasio_foodspend <- raw_data_attanasio_2015$foodexp_m
attanasio_eduspend <- raw_data_attanasio_2015$schoo12
attanasio_existingbusiness <- as.numeric(raw_data_attanasio_2015$BLenterprise)-1
attanasio_indicator <- rep(2, length(attanasio_treatment))
attanasio_cluster_id <- raw_data_attanasio_2015$soum
attanasio_temptation <- raw_data_attanasio_2015$tempt # monthly
attanasio_consumerdurables <- raw_data_attanasio_2015$durconm #monthly


attanasio_takeup <- raw_data_attanasio_2015$nloans_x
attanasio_takeup[attanasio_takeup !=0] <- 1 
attanasio_takeup_number_of_loans <- additional_loans_data_attanasio_2015$f_nloans_x
attanasio_loan_amount_source_dataframe <- data.frame(additional_loans_data_attanasio_2015$f_amount1_x,additional_loans_data_attanasio_2015$f_amount2_x, additional_loans_data_attanasio_2015$f_amount3_x)
attanasio_loan_amount <- rowSums(attanasio_loan_amount_source_dataframe, na.rm=TRUE) 

attanasio_bank_loans_amount_1_only <- additional_loans_data_attanasio_2015$f_amount1_e
attanasio_bank_loans_amount_1_only[additional_loans_data_attanasio_2015$f_source1_e != "Bank"] <- NA
attanasio_bank_loans_amount_2_only <- additional_loans_data_attanasio_2015$f_amount2_e
attanasio_bank_loans_amount_2_only[additional_loans_data_attanasio_2015$f_source2_e != "Bank"] <- NA
attanasio_bank_loans_amount_3_only <- additional_loans_data_attanasio_2015$f_amount3_e
attanasio_bank_loans_amount_3_only[additional_loans_data_attanasio_2015$f_source3_e != "Bank"] <- NA

attanasio_bank_loan_amount_source_dataframe <- data.frame(attanasio_bank_loans_amount_1_only,
                                                          attanasio_bank_loans_amount_2_only,
                                                          attanasio_bank_loans_amount_3_only)
                                                                                                          
attanasio_bank_loan_amount <- rowSums(attanasio_bank_loan_amount_source_dataframe, na.rm=TRUE) 

attanasio_bank_loan_indicator <- attanasio_bank_loan_amount
attanasio_bank_loan_indicator[attanasio_bank_loan_indicator > 0 ] <- 1

attanasio_other_MFI_loan_indicator <- rep(0, length=dim(additional_loans_data_attanasio_2015)[1])
attanasio_other_MFI_loan_amount <- rep(NA, length=dim(additional_loans_data_attanasio_2015)[1])

attanasio_family_loans_amount_1_only <- additional_loans_data_attanasio_2015$f_amount1_e
attanasio_family_loans_amount_1_only[additional_loans_data_attanasio_2015$f_source1_e != "Relative"] <- NA
attanasio_family_loans_amount_2_only <- additional_loans_data_attanasio_2015$f_amount2_e
attanasio_family_loans_amount_2_only[additional_loans_data_attanasio_2015$f_source2_e != "Relative"] <- NA
attanasio_family_loans_amount_3_only <- additional_loans_data_attanasio_2015$f_amount3_e
attanasio_family_loans_amount_3_only[additional_loans_data_attanasio_2015$f_source3_e != "Relative"] <- NA

attanasio_family_loan_amount_source_dataframe <- data.frame(attanasio_family_loans_amount_1_only,
                                                          attanasio_family_loans_amount_2_only,
                                                          attanasio_family_loans_amount_3_only)

attanasio_family_loan_amount <- rowSums(attanasio_family_loan_amount_source_dataframe, na.rm=TRUE) 

attanasio_family_loan_indicator <- attanasio_family_loan_amount
attanasio_family_loan_indicator[attanasio_family_loan_indicator > 0 ] <- 1

attanasio_informal_loans_amount_1_only <- additional_loans_data_attanasio_2015$f_amount1_e
attanasio_informal_loans_amount_1_only[additional_loans_data_attanasio_2015$f_source1_e != "Local Shop"] <- NA
attanasio_informal_loans_amount_2_only <- additional_loans_data_attanasio_2015$f_amount2_e
attanasio_informal_loans_amount_2_only[additional_loans_data_attanasio_2015$f_source2_e != "Local Shop"] <- NA
attanasio_informal_loans_amount_3_only <- additional_loans_data_attanasio_2015$f_amount3_e
attanasio_informal_loans_amount_3_only[additional_loans_data_attanasio_2015$f_source3_e != "Local Shop"] <- NA

attanasio_informal_loan_amount_source_dataframe <- data.frame(attanasio_informal_loans_amount_1_only,
                                                          attanasio_informal_loans_amount_2_only,
                                                          attanasio_informal_loans_amount_3_only)

attanasio_informal_loan_amount <- rowSums(attanasio_informal_loan_amount_source_dataframe, na.rm=TRUE) 

attanasio_informal_loan_indicator <- attanasio_informal_loan_amount
attanasio_informal_loan_indicator[attanasio_informal_loan_indicator > 0 ] <- 1

# now I have to bundle the informal loan and the family loan since most sites did NOT separate them

attanasio_informal_loan_amount_source_dataframe <- data.frame(attanasio_informal_loan_amount,
                                                              attanasio_family_loan_amount)

attanasio_informal_loan_amount <- rowSums(attanasio_informal_loan_amount_source_dataframe, na.rm=TRUE) 

attanasio_informal_loan_indicator <- attanasio_informal_loan_amount
attanasio_informal_loan_indicator[attanasio_informal_loan_indicator > 0 ] <- 1


attanasio_other_formal_loans_amount_1_only <- additional_loans_data_attanasio_2015$f_amount1_e
attanasio_other_formal_loans_amount_1_only[additional_loans_data_attanasio_2015$f_source1_e != "Savings/Credit Co-op"] <- NA
attanasio_other_formal_loans_amount_2_only <- additional_loans_data_attanasio_2015$f_amount2_e
attanasio_other_formal_loans_amount_2_only[additional_loans_data_attanasio_2015$f_source2_e != "Savings/Credit Co-op"] <- NA
attanasio_other_formal_loans_amount_3_only <- additional_loans_data_attanasio_2015$f_amount3_e
attanasio_other_formal_loans_amount_3_only[additional_loans_data_attanasio_2015$f_source3_e != "Savings/Credit Co-op"] <- NA

attanasio_other_formal_loan_amount_source_dataframe <- data.frame(attanasio_other_formal_loans_amount_1_only,
                                                          attanasio_other_formal_loans_amount_2_only,
                                                          attanasio_other_formal_loans_amount_3_only)

attanasio_other_formal_loan_amount <- rowSums(attanasio_other_formal_loan_amount_source_dataframe, na.rm=TRUE) 

attanasio_other_formal_loan_indicator <- attanasio_other_formal_loan_amount
attanasio_other_formal_loan_indicator[attanasio_other_formal_loan_indicator > 0 ] <- 1


attanasio_rosca_loan_indicator <- attanasio_other_formal_loan_indicator
attanasio_rosca_loan_amount <- attanasio_other_formal_loan_amount 


# no healthspend, no women empowerment

attanasio_urban <- rep(0, length(attanasio_treatment))
attanasio_targetwomen <- 1
attanasio_promotion <- 0
attanasio_collateralised <- 1 
attanasio_averageloansize <- 435  #in USD (but not PPP)
attanasio_APR <-  120
attanasio_timegap <- 19 #months
attanasio_currentmarket <-  1 #scale is integers 0-1-2-3 marking no/little/some/full existing credit market activity 
attanasio_village_rand <- 1
attanasio_loanterm <- 7 # CAREFUL it's different for group versus individual loans
attanasio_loansize_percentincome <- 36 #yearly income in the country as a whole




# Load Augsberg et al 2015 #
# has to be put together from multiple dta files, so we will check that the hhids line up
# NOTE I am having an issue with imputed zeroes here from adding things and dropping NAs... 
# I should change this to reflect the NAs as NA...

raw_data_augsberg_2015 <- read.dta("//bbkinghome/rmeager/winprofile/mydocs/Research work/bayesian meta analysis/microcredit analysis/microcredit rct data/augsberg et al 2015/Analysis files_AEJApp-2013-0272/data/Followup/SECTION 6 - Business cl.dta")
augsberg_treatment <- raw_data_augsberg_2015$treatment
augsberg_expenditures_source_dataframe <- data.frame(raw_data_augsberg_2015$bm_expenses, raw_data_augsberg_2015$bs_expenses)
augsberg_expenditures <- rowSums(augsberg_expenditures_source_dataframe, na.rm = TRUE)
augsberg_expenditures_is_missing <- which(rowSums(is.na(augsberg_expenditures_source_dataframe))==ncol(augsberg_expenditures_source_dataframe))
augsberg_expenditures_non_imputed <- rowSums(augsberg_expenditures_source_dataframe, na.rm = TRUE)
augsberg_expenditures_non_imputed[augsberg_expenditures_is_missing] <- NA

augsberg_revenues_source_dataframe <- data.frame(raw_data_augsberg_2015$bm_revenue, raw_data_augsberg_2015$bs_revenue)
augsberg_revenues <- rowSums(augsberg_revenues_source_dataframe, na.rm = TRUE)
augsberg_revenues_is_missing <- which(rowSums(is.na(augsberg_revenues_source_dataframe))==ncol(augsberg_revenues_source_dataframe))
augsberg_revenues_non_imputed <- rowSums(augsberg_revenues_source_dataframe, na.rm = TRUE)
augsberg_revenues_non_imputed[augsberg_revenues_is_missing] <- NA

augsberg_profit_source_dataframe <- data.frame(raw_data_augsberg_2015$bm_profit, raw_data_augsberg_2015$bs_profit)
augsberg_profit <- rowSums(augsberg_profit_source_dataframe, na.rm = TRUE)
augsberg_profit_is_missing <- which(rowSums(is.na(augsberg_profit_source_dataframe))==ncol(augsberg_profit_source_dataframe))
augsberg_profit_non_imputed <- rowSums(augsberg_profit_source_dataframe, na.rm = TRUE)
augsberg_profit_non_imputed[augsberg_profit_is_missing] <- NA

augsberg_profit_bm_only_non_imputed <- raw_data_augsberg_2015$bm_profit
augsberg_profit_bm_only_imputed <- raw_data_augsberg_2015$bm_profit
augsberg_profit_bm_only_imputed[is.na(augsberg_profit_bm_only_imputed)] <- 0

# now let us use the bm only here
augsberg_profit <- augsberg_profit_bm_only_imputed

augsberg_profit_non_imputed <- augsberg_profit_bm_only_non_imputed

augsberg_id_1 <- raw_data_augsberg_2015$intervid

raw_data_augsberg_2015 <- read.dta("//bbkinghome/rmeager/winprofile/mydocs/Research work/bayesian meta analysis/microcredit analysis/microcredit rct data/augsberg et al 2015/Analysis files_AEJApp-2013-0272/data/Followup/SECTION 4 - Assets cl.dta")
augsberg_assets <- raw_data_augsberg_2015$assetvalue
augsberg_id_2 <- raw_data_augsberg_2015$intervid

if (identical(augsberg_id_1, augsberg_id_2) == FALSE){stop("DANGER! IDs in Augsberg data files are NOT aligned!")}

raw_data_augsberg_2015 <- read.dta("//bbkinghome/rmeager/winprofile/mydocs/Research work/bayesian meta analysis/microcredit analysis/microcredit rct data/augsberg et al 2015/Analysis files_AEJApp-2013-0272/data/Followup/SECTION 5 - Income cl.dta")
augsberg_income <- raw_data_augsberg_2015$y_tot #this is yearly
augsberg_businessincome <- raw_data_augsberg_2015$yval_selfempl
augsberg_businessincome[is.na(augsberg_businessincome)] <- 0
augsberg_businessincome_non_imputed <- raw_data_augsberg_2015$yval_selfempl
augsberg_id_3 <- raw_data_augsberg_2015$intervid

if (identical(augsberg_id_1, augsberg_id_3) == FALSE){stop("DANGER! IDs in Augsberg data files are NOT aligned!")}

raw_data_augsberg_2015 <- read.dta("//bbkinghome/rmeager/winprofile/mydocs/Research work/bayesian meta analysis/microcredit analysis/microcredit rct data/augsberg et al 2015/Analysis files_AEJApp-2013-0272/data/Followup/SECTION 3 - Consumption cl.dta")
raw_data_augsberg_2015 <- raw_data_augsberg_2015[order(raw_data_augsberg_2015$intervid),]

augsberg_consumption_source_dataframe <- data.frame(raw_data_augsberg_2015$durablec, 12*(raw_data_augsberg_2015$nondc)) #yearly
designates_missing_augsberg_consumption <- max(augsberg_consumption_source_dataframe, na.rm=TRUE)+100
augsberg_consumption_source_dataframe[!rowSums(!is.na(augsberg_consumption_source_dataframe)),] <- rep(designates_missing_augsberg_consumption,dim(augsberg_consumption_source_dataframe)[2])
augsberg_consumption <- rowSums(augsberg_consumption_source_dataframe, na.rm = TRUE)  
augsberg_consumption[augsberg_consumption== (dim(augsberg_consumption_source_dataframe)[2])*designates_missing_augsberg_consumption ] <- NA
# consumption is yearly here

augsberg_consumerdurables <- raw_data_augsberg_2015$durablec #per year 

# temptation goods are here per WEEK though which is different so convert to year
augsberg_temptation <- 52*(raw_data_augsberg_2015$cigalcc) # now it's per year
augsberg_id_4 <- raw_data_augsberg_2015$intervid

if (identical(augsberg_id_1, augsberg_id_4) == FALSE){stop("DANGER! IDs in Augsberg data files are NOT aligned!")}

raw_data_augsberg_2015 <- read.dta("//bbkinghome/rmeager/winprofile/mydocs/Research work/bayesian meta analysis/microcredit analysis/microcredit rct data/augsberg et al 2015/Analysis files_AEJApp-2013-0272/data/Baseline/BL - SECTION 6 - Business cl merged with followup by RM.dta")
raw_data_augsberg_2015 <- raw_data_augsberg_2015[order(raw_data_augsberg_2015$intervid),]
augsberg_existingbusiness <- raw_data_augsberg_2015$bm_own
augsberg_id_5 <- raw_data_augsberg_2015$intervid

if (identical(augsberg_id_1, augsberg_id_5) == FALSE){stop("DANGER! IDs in Augsberg data files are NOT aligned!")}

raw_data_augsberg_2015 <- read.dta("//bbkinghome/rmeager/winprofile/mydocs/Research work/bayesian meta analysis/microcredit analysis/microcredit rct data/augsberg et al 2015/Analysis files_AEJApp-2013-0272/data/Followup/SECTION 2 - Loans cl.dta")
raw_data_augsberg_2015 <- raw_data_augsberg_2015[order(raw_data_augsberg_2015$intervid),]
augsberg_takeup <- raw_data_augsberg_2015$l1_eki
augsberg_id_6 <- raw_data_augsberg_2015$intervid

if (identical(augsberg_id_1, augsberg_id_6) == FALSE){stop("DANGER! IDs in Augsberg data files are NOT aligned!")}

raw_data_augsberg_2015 <- read.dta("//bbkinghome/rmeager/winprofile/mydocs/Research work/bayesian meta analysis/microcredit analysis/microcredit rct data/augsberg et al 2015/Analysis files_AEJApp-2013-0272/data/Followup/SECTION 7 - Savingscl.dta")
raw_data_augsberg_2015 <- raw_data_augsberg_2015[order(raw_data_augsberg_2015$intervid),]
augsberg_savings <- raw_data_augsberg_2015$savings_avg 
augsberg_id_7 <- raw_data_augsberg_2015$intervid

if (identical(augsberg_id_1, augsberg_id_7) == FALSE){stop("DANGER! IDs in Augsberg data files are NOT aligned!")}

raw_data_augsberg_2015 <- read.dta("//bbkinghome/rmeager/winprofile/mydocs/Research work/bayesian meta analysis/microcredit analysis/microcredit rct data/augsberg et al 2015/Analysis files_AEJApp-2013-0272/data/locations RM generated.dta")
raw_data_augsberg_2015 <- raw_data_augsberg_2015[order(raw_data_augsberg_2015$intervid),]
augsberg_id_8 <- raw_data_augsberg_2015$intervid

if (identical(augsberg_id_1, augsberg_id_8) == FALSE){stop("DANGER! IDs in Augsberg data files are NOT aligned!")}

augsberg_urban <- raw_data_augsberg_2015$urban
augsberg_indicator <- rep(3, length(augsberg_treatment))
augsberg_district <- raw_data_augsberg_2015$branch # yeah this is just the branch not an administrative thing...but that's a great subdivision tactic!

# we don't have villages here... so yeah..

augsberg_targetwomen <- 0
augsberg_promotion <- 0
augsberg_grouploan <- rep(0, length(augsberg_treatment))
augsberg_collateralised <- 1 
augsberg_averageloansize <- 1012 #in USD (but not PPP)
augsberg_APR <- 22
augsberg_timegap <- 14 #months
augsberg_currentmarket <-  2 #scale is integers 0-1-2-3 marking no/little/some/full existing credit market activity 
augsberg_village_rand <- 0 
augsberg_loansize_percentincome <- 9 #percent of income in country per year
augsberg_loanterm <- 14 #months




# Load Banerjee et al 2015 #
# wait but which endline should I use? I say first endline, as second has control hh with access to MFI.

raw_data_banerjee_2015 <- read.dta("//bbkinghome/rmeager/winprofile/mydocs/Research work/bayesian meta analysis/microcredit analysis/microcredit rct data/banerjee et al 2015/2013-0533_data (TO SUBMIT)/2013-0533_data_endlines1and2_stata12.dta")
banerjee_treatment <- as.numeric(raw_data_banerjee_2015$treatment)-1
banerjee_profit <- raw_data_banerjee_2015$bizprofit_1
banerjee_profit[is.na(banerjee_profit)] <- 0
banerjee_profit_non_imputed <- raw_data_banerjee_2015$bizprofit_1
banerjee_revenues <- raw_data_banerjee_2015$bizrev_1
banerjee_revenues[is.na(banerjee_revenues)] <- 0
banerjee_revenues_non_imputed <- raw_data_banerjee_2015$bizrev_1
banerjee_expenditures <- raw_data_banerjee_2015$bizexpense_1
banerjee_expenditures[is.na(banerjee_expenditures)] <- 0
banerjee_expenditures_non_imputed <- raw_data_banerjee_2015$bizexpense_1
banerjee_assets <- raw_data_banerjee_2015$bizassets_1
banerjee_consumption <- raw_data_banerjee_2015$total_exp_mo_1  
banerjee_consumerdurables <- raw_data_banerjee_2015$durables_exp_mo_1
banerjee_foodspend <- raw_data_banerjee_2015$food_exp_mo_1
banerjee_eduspend <- raw_data_banerjee_2015$educ_exp_mo_1
banerjee_healthspend <- raw_data_banerjee_2015$health_exp_mo_1  
banerjee_existingbusiness <- as.numeric(raw_data_banerjee_2015$any_old_biz)-1
banerjee_womenempower <- raw_data_banerjee_2015$women_emp_index_1
banerjee_income <- raw_data_banerjee_2015$income_index_1
banerjee_cluster_id <- raw_data_banerjee_2015$areaid
banerjee_temptation <- raw_data_banerjee_2015$temptation_exp_mo_1
banerjee_district <-  rep(1, length(banerjee_treatment))
banerjee_survey_weights <- raw_data_banerjee_2015$w1
banerjee_urban <- rep(1, length(banerjee_treatment))
banerjee_grouploan <- rep(1, length(banerjee_treatment))
banerjee_indicator <- rep(4, length(banerjee_treatment))
banerjee_hhid <- raw_data_banerjee_2015$hhid

banerjee_takeup <- rep(0, length(raw_data_banerjee_2015$spandana_1))
raw_data_banerjee_2015$spandana_1 <- as.numeric(raw_data_banerjee_2015$spandana_1)
banerjee_takeup <- raw_data_banerjee_2015$spandana_1 - 1
banerjee_loan_amount <- raw_data_banerjee_2015$spandana_amt_1

banerjee_other_MFI_loan_indicator <- rep(0, length(raw_data_banerjee_2015$othermfi_1))
raw_data_banerjee_2015$othermfi_1 <- as.numeric(raw_data_banerjee_2015$othermfi_1)
banerjee_other_MFI_loan_indicator[(raw_data_banerjee_2015$othermfi_1) > 0] <- 1
banerjee_other_MFI_loan_amount <- raw_data_banerjee_2015$othermfi_amt_1

banerjee_bank_loan_indicator <- rep(0, length(raw_data_banerjee_2015$anybank_1))
raw_data_banerjee_2015$anybank_1 <- as.numeric(raw_data_banerjee_2015$anybank_1)
banerjee_bank_loan_indicator<- (raw_data_banerjee_2015$anybank_1)- 1
banerjee_bank_loan_amount <- raw_data_banerjee_2015$bank_amt_1

banerjee_informal_loan_indicator <- rep(0, length(raw_data_banerjee_2015$anyinformal_1))
raw_data_banerjee_2015$anyinformal_1 <- as.numeric(raw_data_banerjee_2015$anyinformal_1)
banerjee_informal_loan_indicator[(raw_data_banerjee_2015$anyinformal_1) > 0] <- 1

banerjee_informal_loan_amount <- raw_data_banerjee_2015$informal_amt_1

banerjee_targetwomen <- 1
banerjee_promotion <- 0
banerjee_collateralised <- 0 
banerjee_averageloansize <- 200 #in USD (but not PPP)
banerjee_APR <- 24
banerjee_timegap <- 14 #months 
banerjee_currentmarket <- 3  #scale is integers 0-1-2-3 marking no/little/some/full existing credit market activity 
banerjee_village_rand <- 1
banerjee_loansize_percentincome <- 22 #percent of income in country per year
banerjee_loanterm <- 12 #months




# Load Crepon et al 2015 #

raw_data_crepon_2015 <- read.dta("//bbkinghome/rmeager/winprofile/mydocs/Research work/bayesian meta analysis/microcredit analysis/microcredit rct data/crepon et al 2015/Data&Code_AEJApp_MicrocreditMorocco/Input/Microcredit_EL_mini_anonym_processed_stata12file.dta")
# This processing consolidates all the little variables which is good
# I HAD doctored the original file in Stata to remove the imputed zeroes for missing values but that seems to mess with the outcome of the simple regression

crepon_dataframe_with_priorbiz <- read.dta("//bbkinghome/rmeager/winprofile/mydocs/Research work/bayesian meta analysis/microcredit analysis/microcredit rct data/crepon et al 2015/Data&Code_AEJApp_MicrocreditMorocco/Input/Microcredit_MS_anonym_stata12file.dta")
raw_data_crepon_2015 <- merge(raw_data_crepon_2015,crepon_dataframe_with_priorbiz,by="ident")


crepon_treatment <- raw_data_crepon_2015$treatment
crepon_expenditures <- raw_data_crepon_2015$expense_total
crepon_expenditures[is.na(crepon_expenditures)] <- 0
crepon_expenditures_non_imputed <- raw_data_crepon_2015$expense_total


crepon_profit <- raw_data_crepon_2015$profit_total
crepon_profit[is.na(crepon_profit)] <- 0
crepon_profit_non_imputed <- raw_data_crepon_2015$profit_total

# now this is a controversial choice for revenues
# there is no clearly defined revenue variable in the data set
# in my view the most sensible option is to construct it by subtracting expense_total from profit_total, I think....
# but this produces negatives so i cannot do it
# crepon_revenues <- crepon_profit - crepon_expenditures
# crepon_revenues_non_imputed <- crepon_profit_non_imputed - crepon_expenditures_non_imputed

# here is the OLD way I used to do it, which is more like "total productive output" than revenues
# I must return to this now 
crepon_revenues_source_dataframe <- data.frame(raw_data_crepon_2015$sale_agri, raw_data_crepon_2015$sale_livestock, raw_data_crepon_2015$sale_business, raw_data_crepon_2015$sale_livestockanim, raw_data_crepon_2015$sale_livestockprod)
crepon_revenues <- rowSums(crepon_revenues_source_dataframe, na.rm = TRUE)
crepon_revenues_is_missing <- which(rowSums(is.na(crepon_revenues_source_dataframe))==ncol(crepon_revenues_source_dataframe))
crepon_revenues_non_imputed <- rowSums(crepon_revenues_source_dataframe, na.rm = FALSE)
crepon_revenues_non_imputed[crepon_revenues_is_missing] <- NA

crepon_assets <- raw_data_crepon_2015$assets_total
crepon_healthspend <- raw_data_crepon_2015$cons_health
crepon_eduspend <- raw_data_crepon_2015$cons_schooling
crepon_foodspend <- raw_data_crepon_2015$cons_food
crepon_consumption <- raw_data_crepon_2015$consumption #monthly
crepon_consumerdurables <- raw_data_crepon_2015$cons_durables #monthly
crepon_womenempower <- raw_data_crepon_2015$women_index
crepon_income <- raw_data_crepon_2015$income
crepon_urban <- rep(0, length(crepon_treatment))
crepon_grouploan <- rep(1, length(crepon_treatment))
crepon_savings_source_dataframe <- data.frame(raw_data_crepon_2015$savings_agri, raw_data_crepon_2015$savings_livestock)
designates_missing_crepon_savings <- max(crepon_savings_source_dataframe, na.rm=TRUE)+100
crepon_savings_source_dataframe[!rowSums(!is.na(crepon_savings_source_dataframe)),] <- rep(designates_missing_crepon_savings,dim(crepon_savings_source_dataframe)[2])
crepon_savings <- rowSums(crepon_savings_source_dataframe, na.rm = TRUE)  
crepon_savings[crepon_savings== (dim(crepon_savings_source_dataframe)[2])*designates_missing_crepon_savings ] <- NA
crepon_temptation <- raw_data_crepon_2015$cons_temp #monthly


crepon_existingbusiness <- raw_data_crepon_2015$me_m7c
crepon_indicator <- rep(5, length(crepon_treatment))
crepon_village_id <- raw_data_crepon_2015$id3
crepon_cluster_id <- raw_data_crepon_2015$id3
crepon_district <- raw_data_crepon_2015$id1

crepon_takeup <- raw_data_crepon_2015$client #Authors confirmed they used this as their measure of takeup
crepon_loan_amount <- raw_data_crepon_2015$loansamt_alamana

crepon_other_formal_loan_indicator <- raw_data_crepon_2015$borrowed_oformal
crepon_other_formal_loan_amount <- raw_data_crepon_2015$loansamt_oformal

crepon_informal_loan_indicator <- raw_data_crepon_2015$borrowed_informal
crepon_informal_loan_amount <- raw_data_crepon_2015$loansamt_informal

crepon_other_MFI_loan_indicator <- raw_data_crepon_2015$borrowed_oamc
crepon_other_MFI_loan_amount <- raw_data_crepon_2015$loansamt_oamc


crepon_targetwomen <- 0
crepon_promotion <- 1
crepon_collateralised <- 0
crepon_averageloansize <- 1188 #in USD (but not PPP)
crepon_APR <- 13.5
crepon_timegap <-  24 #months
crepon_currentmarket <-   0 #scale is integers 0-1-2-3 marking no/little/some/full existing credit market activity 
crepon_village_rand <- 1 
crepon_loansize_percentincome <- 21 #percent of income in country per year
crepon_loanterm <- 16 #months



# Load Karlan and Zinman 2010 #

# This data set gave me a lot of trouble. I still cannot run their do file. Many variables which are supposed to exist - called in do file - are not in the data set.

raw_data_karlan_2010 <- read.dta("//bbkinghome/rmeager/winprofile/mydocs/Research work/bayesian meta analysis/microcredit analysis/microcredit rct data/karlan and zinman 2010/1200138sdataset_clean.dta")
# Actually I may need to run their file first to process the data into useable bits
# if I take only the subsample where sample == 1 I do their analysis. I am doing it here. I may change my mind later.
raw_data_karlan_2010 <- raw_data_karlan_2010[raw_data_karlan_2010$sample ==1,]
karlan_treatment <- as.numeric(raw_data_karlan_2010$css_loandecision_raw)
karlan_takeup <- as.numeric(raw_data_karlan_2010$css_loandecision)-1
karlan_profit_source_dataframe <- data.frame(raw_data_karlan_2010$fu_profit_1, raw_data_karlan_2010$fu_profit_2, raw_data_karlan_2010$fu_profit_3, raw_data_karlan_2010$fu_profit_4, raw_data_karlan_2010$fu_profit_5,raw_data_karlan_2010$fu_profit_6, raw_data_karlan_2010$fu_profit_7, raw_data_karlan_2010$fu_profit_8)
karlan_profit <- rowSums(karlan_profit_source_dataframe, na.rm = TRUE)
karlan_profit_is_missing <- which(rowSums(is.na(karlan_profit_source_dataframe))==ncol(karlan_profit_source_dataframe))
karlan_profit_non_imputed <- rowSums(karlan_profit_source_dataframe, na.rm = TRUE)
karlan_profit_non_imputed[karlan_profit_is_missing] <- NA



karlan_expenditures_source_dataframe <-  data.frame(raw_data_karlan_2010$fu_expenses_inventory_1, raw_data_karlan_2010$fu_expenses_inventory_2, raw_data_karlan_2010$fu_expenses_inventory_3, raw_data_karlan_2010$fu_expenses_inventory_4, raw_data_karlan_2010$fu_expenses_inventory_5,raw_data_karlan_2010$fu_expenses_inventory_6, raw_data_karlan_2010$fu_expenses_inventory_7, raw_data_karlan_2010$fu_expenses_inventory_8,
                                                    raw_data_karlan_2010$fu_expenses_bills_1, raw_data_karlan_2010$fu_expenses_bills_2, raw_data_karlan_2010$fu_expenses_bills_3, raw_data_karlan_2010$fu_expenses_bills_4, raw_data_karlan_2010$fu_expenses_bills_5,raw_data_karlan_2010$fu_expenses_bills_6, raw_data_karlan_2010$fu_expenses_bills_7, raw_data_karlan_2010$fu_expenses_bills_8,
                                                    raw_data_karlan_2010$fu_expenses_wages_1, raw_data_karlan_2010$fu_expenses_wages_2, raw_data_karlan_2010$fu_expenses_wages_3, raw_data_karlan_2010$fu_expenses_wages_4, raw_data_karlan_2010$fu_expenses_wages_5,raw_data_karlan_2010$fu_expenses_wages_6, raw_data_karlan_2010$fu_expenses_wages_7, raw_data_karlan_2010$fu_expenses_wages_8,
                                                    raw_data_karlan_2010$fu_expenses_rent_equip_1, raw_data_karlan_2010$fu_expenses_rent_equip_2, raw_data_karlan_2010$fu_expenses_rent_equip_3, raw_data_karlan_2010$fu_expenses_rent_equip_4, raw_data_karlan_2010$fu_expenses_rent_equip_5,raw_data_karlan_2010$fu_expenses_rent_equip_6, raw_data_karlan_2010$fu_expenses_rent_equip_7, raw_data_karlan_2010$fu_expenses_rent_equip_8,
                                                    raw_data_karlan_2010$fu_expenses_rent_building_1, raw_data_karlan_2010$fu_expenses_rent_building_2, raw_data_karlan_2010$fu_expenses_rent_building_3, raw_data_karlan_2010$fu_expenses_rent_building_4, raw_data_karlan_2010$fu_expenses_rent_building_5,raw_data_karlan_2010$fu_expenses_rent_building_6, raw_data_karlan_2010$fu_expenses_rent_building_7, raw_data_karlan_2010$fu_expenses_rent_building_8,
                                                    raw_data_karlan_2010$fu_expenses_taxes_1, raw_data_karlan_2010$fu_expenses_taxes_2, raw_data_karlan_2010$fu_expenses_taxes_3, raw_data_karlan_2010$fu_expenses_taxes_4, raw_data_karlan_2010$fu_expenses_taxes_5,raw_data_karlan_2010$fu_expenses_taxes_6, raw_data_karlan_2010$fu_expenses_taxes_7, raw_data_karlan_2010$fu_expenses_taxes_8,
                                                    raw_data_karlan_2010$fu_expenses_maint_1, raw_data_karlan_2010$fu_expenses_maint_2, raw_data_karlan_2010$fu_expenses_maint_3, raw_data_karlan_2010$fu_expenses_maint_4, raw_data_karlan_2010$fu_expenses_maint_5,raw_data_karlan_2010$fu_expenses_maint_6, raw_data_karlan_2010$fu_expenses_maint_7, raw_data_karlan_2010$fu_expenses_maint_8,
                                                    raw_data_karlan_2010$fu_expenses_transport_1, raw_data_karlan_2010$fu_expenses_transport_2, raw_data_karlan_2010$fu_expenses_transport_3, raw_data_karlan_2010$fu_expenses_transport_4, raw_data_karlan_2010$fu_expenses_transport_5,raw_data_karlan_2010$fu_expenses_transport_6, raw_data_karlan_2010$fu_expenses_transport_7, raw_data_karlan_2010$fu_expenses_transport_8,
                                                    raw_data_karlan_2010$fu_expenses__other_1, raw_data_karlan_2010$fu_expenses__other_2, raw_data_karlan_2010$fu_expenses__other_3, raw_data_karlan_2010$fu_expenses__other_4, raw_data_karlan_2010$fu_expenses__other_5,raw_data_karlan_2010$fu_expenses__other_6, raw_data_karlan_2010$fu_expenses__other_7, raw_data_karlan_2010$fu_expenses__other_8
)

karlan_expenditures <- rowSums(karlan_expenditures_source_dataframe, na.rm = TRUE) 
karlan_expenditures_is_missing <- which(rowSums(is.na(karlan_expenditures_source_dataframe))==ncol(karlan_expenditures_source_dataframe))
karlan_expenditures_non_imputed <- rowSums(karlan_expenditures_source_dataframe, na.rm = TRUE)
karlan_expenditures_non_imputed[karlan_expenditures_is_missing] <- "NA"
karlan_expenditures_non_imputed <- destring(karlan_expenditures_non_imputed)

karlan_gross_sales_source_dataframe <- data.frame(raw_data_karlan_2010$fu_gross_sales_1, raw_data_karlan_2010$fu_gross_sales_2, raw_data_karlan_2010$fu_gross_sales_3, raw_data_karlan_2010$fu_gross_sales_4, raw_data_karlan_2010$fu_gross_sales_5,raw_data_karlan_2010$fu_gross_sales_6, raw_data_karlan_2010$fu_gross_sales_7, raw_data_karlan_2010$fu_gross_sales_8)
karlan_gross_sales <- rowSums(karlan_gross_sales_source_dataframe, na.rm = TRUE)
karlan_gross_sales_is_missing <- which(rowSums(is.na(karlan_gross_sales_source_dataframe))==ncol(karlan_gross_sales_source_dataframe))
karlan_gross_sales_non_imputed <- rowSums(karlan_gross_sales_source_dataframe, na.rm = TRUE)
karlan_gross_sales_non_imputed[karlan_gross_sales_is_missing] <- NA

karlan_revenues <- karlan_gross_sales # this used to be profit - expenditures in old version 
karlan_revenues_non_imputed <- karlan_gross_sales_non_imputed


karlan_income <- raw_data_karlan_2010$fu_hh_total_inc
karlan_savings <- as.numeric(raw_data_karlan_2010$fu_savings) #beware the factor variable translated... this was a weird bins indicator
karlan_assets_source_dataframe <-  data.frame(raw_data_karlan_2010$fu_assets_mkt_val_land,
                                              raw_data_karlan_2010$fu_assets_mkt_val_buildings,
                                              raw_data_karlan_2010$fu_assets_mkt_val_cell_phones,
                                              raw_data_karlan_2010$fu_assets_mkt_val_stoves,
                                              raw_data_karlan_2010$fu_assets_mkt_val_fans,
                                              raw_data_karlan_2010$fu_assets_mkt_val_ac,
                                              raw_data_karlan_2010$fu_assets_mkt_val_tv_sets,
                                              raw_data_karlan_2010$fu_assets_mkt_val_vhs_dvd,
                                              raw_data_karlan_2010$fu_assets_mkt_val_cars,
                                              raw_data_karlan_2010$fu_assets_mkt_val_motor_tri_bi,
                                              raw_data_karlan_2010$fu_assets_mkt_val_jeepneys,
                                              raw_data_karlan_2010$fu_assets_mkt_val_h2o_pur_dev,
                                              raw_data_karlan_2010$fu_assets_mkt_val_fridge,
                                              raw_data_karlan_2010$fu_assets_mkt_val_computers,
                                              raw_data_karlan_2010$fu_assets_mkt_val_bus_equip,
                                              raw_data_karlan_2010$fu_assets_mkt_val_livestock,
                                              raw_data_karlan_2010$fu_assets_mkt_val_furniture,
                                              raw_data_karlan_2010$fu_assets_mkt_val_other1,
                                              raw_data_karlan_2010$fu_assets_mkt_val_other2
)  

designates_missing_karlan_assets <- max(karlan_assets_source_dataframe, na.rm=TRUE)+100
karlan_assets_source_dataframe[!rowSums(!is.na(karlan_assets_source_dataframe)),] <- rep(designates_missing_karlan_assets,dim(karlan_assets_source_dataframe)[2])
karlan_assets <- rowSums(karlan_assets_source_dataframe, na.rm = TRUE)  
karlan_assets[karlan_assets== (dim(karlan_assets_source_dataframe)[2])*designates_missing_karlan_assets ] <- NA



karlan_foodspend <- raw_data_karlan_2010$fu_hh_food_cost
karlan_existingbusiness <- rep(1, length(raw_data_karlan_2010$css_loandecision))
karlan_grouploan <- rep(0, length(raw_data_karlan_2010$css_loandecision))
karlan_individualloan <- rep(1, length(raw_data_karlan_2010$css_loandecision))
karlan_urban <- rep(1, length(raw_data_karlan_2010$css_loandecision))
karlan_indicator <- rep(6, length(karlan_treatment))
karlan_district <- raw_data_karlan_2010$branch # again i am using branch but that is a great way to subdistrict people
# karlan_takeup <- this I do not see !!!!!  
# consumption??
# temptation goods???

karlan_targetwomen <- 0
karlan_promotion <- 0
karlan_collateralised <- 0
karlan_averageloansize <- 400 #in USD (but not PPP)
karlan_APR <- 63
karlan_timegap <-  12 #months
karlan_currentmarket <-  1  #scale is integers 0-1-2-3 marking no/little/some/full existing credit market activity 
karlan_village_rand <- 0 
karlan_loansize_percentincome <- (400/1660)*100 #percent of income in country per year
# karlan_loanterm <- NOPE it's not in there pal



# Load Tarozzi et al 2015 #

#currency units are 2006 Birr unless otherwise noted

raw_data_tarozzi_2015 <- read.dta("//bbkinghome/rmeager/winprofile/mydocs/Research work/bayesian meta analysis/microcredit analysis/microcredit rct data/tarozzi et al 2015/TarozziEtAlReplicationFiles/data_stata12.dta")
# This contains another treatment and their interaction which I shall drop here
raw_data_tarozzi_2015 <- raw_data_tarozzi_2015[raw_data_tarozzi_2015$time == "Endline",]
raw_data_tarozzi_2015 <- raw_data_tarozzi_2015[raw_data_tarozzi_2015$patypen != 'Both',]
raw_data_tarozzi_2015 <- raw_data_tarozzi_2015[raw_data_tarozzi_2015$patypen != 'FP',]
tarozzi_treatment <- raw_data_tarozzi_2015$t_D
tarozzi_revenues <- raw_data_tarozzi_2015$revenues
tarozzi_revenues[is.na(tarozzi_revenues)] <- 0
tarozzi_revenues_non_imputed <- raw_data_tarozzi_2015$revenues
tarozzi_expenditures <- raw_data_tarozzi_2015$costs
tarozzi_expenditures[is.na(tarozzi_expenditures)] <- 0
tarozzi_expenditures_non_imputed <- raw_data_tarozzi_2015$costs
tarozzi_profit <- raw_data_tarozzi_2015$net #this also functions as business income
tarozzi_profit[is.na(tarozzi_profit)] <- 0
tarozzi_profit_non_imputed <- raw_data_tarozzi_2015$net #this also functions as business incom

tarozzi_assets <- raw_data_tarozzi_2015$assets
tarozzi_income_source_dataframe <- data.frame(raw_data_tarozzi_2015$inc_wage, raw_data_tarozzi_2015$inc_other, raw_data_tarozzi_2015$inc_otherall)
tarozzi_income <- rowSums(tarozzi_income_source_dataframe, na.rm = TRUE)
tarozzi_income_is_missing <- which(rowSums(is.na(tarozzi_income_source_dataframe))==ncol(tarozzi_income_source_dataframe))
tarozzi_income_non_imputed <- rowSums(tarozzi_income_source_dataframe, na.rm = TRUE)
tarozzi_income_non_imputed[tarozzi_income_is_missing] <- NA
tarozzi_healthspend <- raw_data_tarozzi_2015$health_exp
tarozzi_existingbusiness <- raw_data_tarozzi_2015$oldbusiness # double checked with authors that I interpret this right
tarozzi_indicator <- rep(7, length(tarozzi_treatment))
tarozzi_cluster_id <- raw_data_tarozzi_2015$pa


tarozzi_takeup <- raw_data_tarozzi_2015$loans_mf
tarozzi_loan_amount <- raw_data_tarozzi_2015$loanstot_mf

tarozzi_informal_loan_indicator <- raw_data_tarozzi_2015$loans_informal
tarozzi_informal_loan_amount <- raw_data_tarozzi_2015$loanstot_informal

tarozzi_moneylender_loan_indicator <- raw_data_tarozzi_2015$loans_ml

tarozzi_formal_loan_indicator <- raw_data_tarozzi_2015$loans_formal
tarozzi_formal_loan_amount <- raw_data_tarozzi_2015$loanstot_formal

tarozzi_rosca_loan_indicator <- raw_data_tarozzi_2015$loans_rca
tarozzi_rosca_loan_amount <- raw_data_tarozzi_2015$loanstot_rca

tarozzi_ngo_loan_indicator <- raw_data_tarozzi_2015$loans_ngo
tarozzi_ngo_loan_amount <- raw_data_tarozzi_2015$loanstot_ngo

tarozzi_urban <- 0*tarozzi_treatment
tarozzi_targetwomen <- 0 # BUT THIS WAS NOT ACTUALLY DONE BY OFFICERS...
tarozzi_promotion <- 0
tarozzi_grouploan <- rep(1, length(tarozzi_treatment))
tarozzi_collateralised <- 0 # BUT IN REALITY PEOPLE DID GET ASKED FOR COLLATERAL GAAAAAH
tarozzi_averageloansize <- 150 #in USD (but not PPP)
tarozzi_APR <- 12
tarozzi_timegap <- 36 #months
tarozzi_currentmarket <- 1 #scale is integers 0-1-2-3 marking no/little/some/full existing credit market activity 
tarozzi_village_rand <- 1
tarozzi_loansize_percentincome <- 118 #percent of income in country per year
tarozzi_loanterm <- 12 #months


### NOW DROP THE ORIGINALS TO FREE UP RAM ### 

rm(raw_data_angelucci_2015, raw_data_attanasio_2015, raw_data_augsberg_2015, raw_data_banerjee_2015, raw_data_crepon_2015, raw_data_karlan_2010, raw_data_tarozzi_2015)


### IMPORT THE REGRESSION RESULTS AND UNITS INFORMATION FROM THE PAPERS ###

# general variables #

study_names <- c("Angelucci", "Attanasio", "Augsberg", "Banerjee", "Crepon", "Karlan", "Tarozzi") 
study_country <- c("Mexico", "Mongolia", "Bosnia", "India", "Morocco", "Philippines", "Ethiopia")

reported_coefficients_years <- c(2012, 2009, 2010, 2007, 2009, 2008, 2006)
reported_coefficients_currency <- c("Mexican Pesos", "Mongolian Togrog", "Bosnia Herzegovina Convertible Mark", "Indian Rupee", "Moroccan Currency", "Philippines Currency", "Ethiopian Birr")

USD_PPP_conversion <- c(1/9.18, 1/513.24, 1/0.88, 1/11.09, 1/4.31, 1/17.42, 1/2.29  ) #this is to reflect the NEW years, the endline years 
USD_convert_to_2009_dollars <- 100*c(1/106.121, 1/100, 1/101.653, 1/97.101, 1/100, 1/100.065, 1/94.729)

# profit import and standardiser code # 

time_period_measurement_profit <- c("2 weeks", "1 year", " 1 year", "1 month", "1 year", " 1 month",  "1 year")
fortnight_converter_profit <- c(1, 1/26, 1/26, 14/30, 1/26, 14/30, 1/26)

reported_coefficients_profit <- c(0, -4789, 671.9, 354, 2005, 2482.57, 513)
reported_coefficients_profit_sds <- c(39,5302, 541.3, 314, 1210, 2114.02, 431)

the_standardiser_USD_PPP <- USD_PPP_conversion*USD_convert_to_2009_dollars 
the_profit_standardiser_USD_PPP_per_fortnight <- USD_PPP_conversion*USD_convert_to_2009_dollars*fortnight_converter_profit

standardised_reported_coefficients_profit <- the_profit_standardiser_USD_PPP_per_fortnight*reported_coefficients_profit
standardised_reported_coefficients_profit_sds <- the_profit_standardiser_USD_PPP_per_fortnight*reported_coefficients_profit_sds


# now bind it into a data frame
study_meta_data_summary_table_profit <- data.frame(study_names, study_country, reported_coefficients_years, reported_coefficients_currency,time_period_measurement_profit, 
                                            fortnight_converter_profit, USD_PPP_conversion,USD_convert_to_2009_dollars)

# print a nice table
pdf("//bbkinghome/rmeager/winprofile/mydocs/Research work/bayesian meta analysis/microcredit analysis/study meta data.pdf", height=4, width=20)
grid.table(study_meta_data_summary_table_profit)
dev.off()

# revenues #

reported_coefficients_revenues <- c( 121, NA , NA, NA, NA, NA, NA )
reported_coefficients_revenues_sds <- c( 52, NA , NA, NA, NA, NA, NA    )


time_period_measurement_revenues <- c("2 weeks", "1 year", " 1 year", "1 month",   "1 year", " 1 month",  "1 year")
                                            
fortnight_converter_revenues <- c(1, 1/26, 1/26, 14/30, 1/26, 14/30, 1/26)
                                         
the_revenues_standardiser_USD_PPP_per_fortnight <- USD_PPP_conversion*USD_convert_to_2009_dollars*fortnight_converter_revenues


standardised_reported_coefficients_revenues <- the_revenues_standardiser_USD_PPP_per_fortnight*reported_coefficients_revenues
standardised_reported_coefficients_revenues_sds <- the_revenues_standardiser_USD_PPP_per_fortnight*reported_coefficients_revenues_sds


# expenditures # 

reported_coefficients_expenditures <- c( 119, NA , NA, NA, NA, NA, NA)
reported_coefficients_expenditures_sds <- c(47 , NA , NA, NA, NA, NA, NA  )


time_period_measurement_expenditures <- c("2 weeks", "1 year", " 1 year", "1 month", "1 year", " 1 month",  "1 year")
                                              
fortnight_converter_expenditures <- c(1,  1/26,  1/26, 14/30, 1/26, 14/30, 1/26)

the_expenditures_standardiser_USD_PPP_per_fortnight <- USD_PPP_conversion*USD_convert_to_2009_dollars*fortnight_converter_expenditures


standardised_reported_coefficients_expenditures <- the_expenditures_standardiser_USD_PPP_per_fortnight*reported_coefficients_expenditures
standardised_reported_coefficients_expenditures_sds <- the_expenditures_standardiser_USD_PPP_per_fortnight*reported_coefficients_expenditures_sds

                                                
# consumption # 
 
reported_coefficients_consumption <- c(NA, NA, NA, NA, NA)
reported_coefficients_consumption_sds <- c(NA, NA , NA, NA, NA)

time_period_measurement_consumption <- c("1 week", "1 month",  " 1 year", "1 month", " 1 month")                                         

fortnight_converter_consumption <- c(2, 14/30,   1/26, 14/30,  14/30)
the_consumption_standardiser_USD_PPP_per_fortnight <- USD_PPP_conversion*USD_convert_to_2009_dollars*fortnight_converter_consumption

the_consumption_standardiser_USD_PPP_per_fortnight <- the_consumption_standardiser_USD_PPP_per_fortnight[1:5]

standardised_reported_coefficients_consumption <- the_consumption_standardiser_USD_PPP_per_fortnight*reported_coefficients_consumption
standardised_reported_coefficients_consumption_sds <- the_consumption_standardiser_USD_PPP_per_fortnight*reported_coefficients_consumption_sds

# temptation # 
 
reported_coefficients_temptation <- c(NA, NA, NA, NA, NA)
reported_coefficients_temptation_sds <- c(NA, NA , NA, NA, NA)

time_period_measurement_temptation <- c("1 week", "1 month",  " 1 year", "1 month", " 1 month")                                         

fortnight_converter_temptation <- c(2, 14/30,   1/26, 14/30,  14/30)
the_temptation_standardiser_USD_PPP_per_fortnight <- USD_PPP_conversion*USD_convert_to_2009_dollars*fortnight_converter_consumption

the_temptation_standardiser_USD_PPP_per_fortnight <- the_temptation_standardiser_USD_PPP_per_fortnight[1:5]

standardised_reported_coefficients_temptation <- the_temptation_standardiser_USD_PPP_per_fortnight*reported_coefficients_temptation
standardised_reported_coefficients_temptation_sds <- the_temptation_standardiser_USD_PPP_per_fortnight*reported_coefficients_temptation_sds


# consumer durables # 

# consumption # 

reported_coefficients_consumerdurables <- c(NA, NA, NA, NA)
reported_coefficients_consumerdurables_sds <- c(NA, NA , NA, NA)

time_period_measurement_consumerdurables <- c("1 month",  " 1 year", "1 month", " 1 month")                                         

fortnight_converter_consumerdurables <- c(14/30,   1/26, 14/30,  14/30)
the_consumerdurables_standardiser_USD_PPP_per_fortnight <- USD_PPP_conversion*USD_convert_to_2009_dollars*fortnight_converter_consumerdurables

the_consumerdurables_standardiser_USD_PPP_per_fortnight <- the_consumerdurables_standardiser_USD_PPP_per_fortnight[2:5]

standardised_reported_coefficients_consumerdurables <- the_consumerdurables_standardiser_USD_PPP_per_fortnight*reported_coefficients_consumerdurables
standardised_reported_coefficients_consumerdurables_sds <- the_consumerdurables_standardiser_USD_PPP_per_fortnight*reported_coefficients_consumerdurables_sds


## table of contextual variables ## 

# here's the set of contextual variables
targetwomen <- c(angelucci_targetwomen, attanasio_targetwomen, augsberg_targetwomen,banerjee_targetwomen, crepon_targetwomen, karlan_targetwomen, tarozzi_targetwomen)
village_rand <- c(angelucci_village_rand, attanasio_village_rand, augsberg_village_rand,banerjee_village_rand, crepon_village_rand, karlan_village_rand, tarozzi_village_rand)
individual_rand <- 1-village_rand
APR <- c(angelucci_APR, attanasio_APR, augsberg_APR,banerjee_APR, crepon_APR, karlan_APR, tarozzi_APR)
market_saturation <- c(angelucci_currentmarket, attanasio_currentmarket, augsberg_currentmarket,banerjee_currentmarket, crepon_currentmarket, karlan_currentmarket, tarozzi_currentmarket)
promotion <- c(angelucci_promotion, attanasio_promotion, augsberg_promotion,banerjee_promotion, crepon_promotion, karlan_promotion, tarozzi_promotion)
collateralised <- c(angelucci_collateralised, attanasio_collateralised, augsberg_collateralised,banerjee_collateralised, crepon_collateralised, karlan_collateralised, tarozzi_collateralised)
loansize_percentincome <- c(angelucci_loansize_percentincome, attanasio_loansize_percentincome, augsberg_loansize_percentincome,banerjee_loansize_percentincome, crepon_loansize_percentincome, karlan_loansize_percentincome, tarozzi_loansize_percentincome)

X_data <- cbind(individual_rand, targetwomen, APR, market_saturation, promotion, collateralised, loansize_percentincome)
study_names <- c("Mexico (Angelucci)", "Mongolia (Attanasio)", "Bosnia (Augsberg)", "India (Banerjee)", "Morocco (Crepon)", "Philippines (Karlan)", "Ethiopia (Tarozzi)")
dataframe_X <- data.frame(X_data)
rownames(dataframe_X) <- study_names
colnames(dataframe_X) <- c("Rand unit", 
                           "Target Women",
                           "APR",
                           "Market Saturation (1-4)",
                           "Promotion",
                           "Collateralisation",
                           "Loan size")

contextuals_table <- xtable(dataframe_X)
sink("contextual_vars.txt")
print.xtable(contextuals_table)
sink()




#now save every object we created 

save.image(file = "microcredit_project_data.RData")

# and now save it to the dist TE project 

save.image(file = "//bbkinghome/rmeager/winprofile/mydocs/Research work/aggregating distributional effects/data/microcredit_project_data.RData")

