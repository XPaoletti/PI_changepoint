# PI_changepoint
Codes for the stat Med paper on a new joint model using changepoint for dose finding 

Maria-Athina Altzerinakou and Xavier Paoletti



This file contains a description of the “2+2_crm_jm.R” R code that was used for the simulation study presented in the article “Change-point joint model for identification of plateau of activity in early phase trials”. This code contains the 16 functions required to run the proposed design.  All 16 functions should be executed before running the code. A description of each function has been added within the R code. 

The following packages should be installed: 
•	nlme 
•	mnormt
•	plyr
•	tensor
•	numDeriv
•	matrixcalc
•	lme4
•	pbkrtest
•	car

List of functions
•	Sigma.is
•	L.is
•	makedata
•	databegin
•	updown
•	data2plcrmjm
•	sequentialentry0
•	escalate
•	loglikcrm
•	nextdose
•	sequentialentry
•	crm
•	loglik
•	checkconv
•	nextdose1
•	jmproc



After executing all functions then the function jmproc can be used to run a simulation study. To execute jmproc the arguments below should be assigned to certain values. 

jmproc(meanpars, varpars, gampars, nt, n1, d, dose, nvisit, maxiter, target, censor, pl, meandiff, stopping_rule, inter_miss, n_min_jm, se_max_jm, fix_param, seed)  



Arguments in jmproc

•	meanpars → sets the parameters for the longitudinal and the survival models (in the order longitudinal: intercept, time, dose survival: intercept, cycle, dose)
•	varpars → sets the variance parameters (in the order: residual variance and random effect variance)
•	gampars → sets the gamma parameter
•	nt → sets the number of treatment cycles (should always be set to 6)
•	n1 → sets the number of subjects in the trial 
•	d → sets the total number of doses to be tested
•	dose → sets the vector of the dose levels
•	nvisit → sets the maximum number of repeated measurements (should always be set to 18)
•	maxiter → sets the maximum number of iterations in nlminb
•	target → sets the target toxicity level 
•	censor → sets the uniform censoring (TRUE/FALSE)
•	pl → indicates the dose where the plateau begins. In case of no plateau set pl=1
•	meandiff → sets the maximum accepted difference so that one or more doses can be considered equivalent to the MTD in terms of activity
•	stopping_rule → sets the sequential boundaries to test for excessive toxicity
•	inter_miss → sets the probability of intermittent missing visits (per visit)
•	n_min_jm → sets the minimum acceptable number of subjects before applying the joint modeling
•	se_max_jm → sets the maximum acceptable standard error for parameter estimation in the joint modeling 
•	fix_param → sets the dose slope of the survival model to a fixed value
•	seed → sets the seed to run the simulation study

To illustrate how the algorithm works we have provided an example. This example can be found in the file “Example”. Within this file we provide a list of the results that we obtain from this study. A description of these results can be found at the end of the jmproc function.

In order to reproduce this simulation study one can use the arguments below for the jmproc function: 
jmproc(meanpars=c(550, 0.5, -75, 5.6, 0.3, -3.4), varpars=c(3, 1), gampars=c(0.2), nt=6, n1=40, d=10, dose=seq(1.0, 1.54, 0.06), nvisit=18, maxiter=10000, target=0.4, censor=T, pl=5, meandiff=20,  stopping_rule=c(0, 0, 3, 4, 4, 5, 5, 5, 6, 6, 7, 7, 8, 8, 8, rep(1000,25)), inter_miss = 0.07, n_min_jm =16, se_max_jm=20, fix_param =-3.297350, seed=120)

These arguments are the ones we used for scenario 1.2 in Table 1 in the associated article.  

Important note 1:  The code was created uniquely for the design we propose in the article. In order to be adapted to different types of designs, number of cycles, and number of parameters, changes should be made within the functions.
Important note 2: Based on the number of doses, the joint modeling algorithm will try to find the model that best fits the data. Among the candidate models certain will produce an error as they cannot be fitted. These errors, like the ones shown below, will be ignored by the algorithm, and the next candidate model will be evaluated. 
Examples of error: 
Error in intervals.lme(lmedat) : 
  cannot get confidence intervals on var-cov components: Non-positive definite approximate variance-covariance
 Consider 'which = "fixed"'

Error in solve.default(estimates[dimE[1] - (p:1), dimE[2] - (p:1), drop = FALSE]) : 
  system is computationally singular: reciprocal condition number = 2.71053e-17
In addition: Warning message:
In log(Psi[ids, i]) : NaNs produced



Altzerinakou Maria-Athina
Novartis Pharma AG, Rueil-Malmaison, France

Email: m.altzerinakou@gmail.com

Xavier Paoletti
Institut Curie  & University of Versailles St Quentin
