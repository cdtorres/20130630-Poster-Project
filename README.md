20130630-Poster-Project
=======================

Simulations of clinical trials using certain Bayesian adaptive methods

simulatetrials.R defines two functions:

clinicaltrial() simulates a single clinical trial with adaptive randomization, as well as stopping early for futility/efficacy

update.evaluate() updates posterior probabilities and randomization probabilities

simulatetrials() calls clinicaltrial() multiple times (which in turn calls update.evaluate() multiple times), returning a data frame

the function simulatetrials() requires package 'multicore'

'multicore' only runs in non-Windows environments.

work.R initializes some variables and calls the simulatetrials() function, while determining how long it takes to run

additionally, it provides some summary statistics on the results with the function simsum()


