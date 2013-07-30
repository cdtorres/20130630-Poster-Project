20130630-Poster-Project
=======================

Simulations of clinical trials using certain Bayesian adaptive methods

simulatetrials.R defines five functions:

clinicaltrial() simulates a single clinical trial with adaptive randomization, as well as stopping early for futility/efficacy

update.evaluate() updates posterior probabilities and randomization probabilities

simulatetrials() calls clinicaltrial() multiple times (which in turn calls update.evaluate() multiple times), returning a data frame

simsum() gives a summary of the simulatetrials() object

get.thresholds() gets thresholds required for a certain set of parameters, to get the desired operating characteristics



the function simulatetrials() requires package 'multicore'

'multicore' only runs in non-Windows environments.

initialize.R initializes some variables, while execute.R calls the simulatetrials() or get.thresholds() function, as well as determining how long it takes to run



