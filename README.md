20130630-Poster-Project
=======================

Simulations of clinical trials using certain Bayesian adaptive methods

simulatetrials.R defines two functions:

clinicaltrial() simulates a single clinical trial

simulatetrials() calls clinicaltrial() multiple times, returning a data frame

these functions require packages 'multicore' and 'pracma'

'multicore' only runs in non-Windows environments.

work.R sets the seed and calls the simulatetrials() command, while determining how long it takes to run

additionally, it provides some summary statistics on the results




Some notes:

simulatetrials(theta_a = .2, theta_b = .2, prior = .2, B = 500, how_often = 1) takes 25 minutes

2/3 of the trials here didn't stop early
