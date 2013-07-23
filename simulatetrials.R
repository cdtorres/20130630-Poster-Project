# 2013/07 Cesar Torres
# There are 4 functions defined in this file: clinicaltrial(), simulatetrials(), update.evaluate(),
# and simsum(). update.evaluate() requires the package 'dblquad', and simulatetrials() requires
# the package 'multicore'. Please note that 'multicore' does not run on Windows systems.



# clinicaltrial() simulates a single clinical trial
# it takes in a 'seed' value for the purposes of replication
# theta_a and theta_b are the 'true' effects of placebo and treatment, respectively
# prior is the prior mean for the Beta distribution for both placebo and treatment group
# in Beta(a, b), b is the second_parameter
# N is the maximum number of patients
# integral_tolerance is for the two-dimensional integral when evalualing P(theta_a + delta < theta_b | data)
# (making integral_tolerance small greatly increases the time that a simulation takes to run)
# how_often is the number of patients we treat before evaluating P(theta_a + delta < theta_b | data),
# after the initial 10 patients and initial evaluation
clinicaltrial <-function(seed, theta_a, theta_b, prior, second_parameter, N, efficacy_threshold,
                         futility_threshold, integral_tolerance, how_often, delta)
{
  
  efficacious = F
  futile = F
  stopped_early = F#indicator for if the trial stopped before all N patients were treated
  
  #This is a vector of values that we're interested in. The seven values in order are:
  #1 number of patients treated so far
  #2 probability of being assigned to treatment group a (the placebo group)
  #3 current value of a_1, given the prior and the data obtained thus far
  #4 current value of a_2, given the prior and the data obtained thus far
  #5 current value of b_1, given the prior and the data obtained thus far
  #6 current value of b_2, given the prior and the data obtained thus far
  #7 P(theta_a + delta < theta_b | data)
  v = rep(NA, 7)
  
  #number of patients treated so far
  v[1] = 0
  
  #group a is control, group b is treatment
  #For the control,   theta_a is distributed Beta(a_1, a_2)
  #For the treatment, theta_b is distributed Beta(b_1, b_2)
  
  #The higher the second parameter the more informative the prior is. (This needs to be > 0.)
  
  #theta_a and theta_b are the "true" probabilities of recovering with placebo and treatment, respectively
  #under the null hypothesis, they're equal
  
  v[4] = second_parameter#initialization of a_2, given priors
  v[6] = second_parameter#initialization of b_2, given priors
  
  v[3] = v[4]*prior/(1 - prior)#initialization of a_1, given priors
  v[5] = v[6]*prior/(1 - prior)#initialization of b_1, given priors
  
  #initial randomization probability, set to 1/2
  #probability of being assigned to treatment group a (the control group)
  v[2] = 1/2
  
  set.seed(seed)#for the purposes of replication

  #treat 10 patients before updating the randomization probability
  #or evaluating P(theta_a + delta < theta_b | data)
  v = update.evaluate(theta_a, theta_b, delta, integral_tolerance, how_many = 10, ov = v)
  
  #remember, v[7] is P(theta_a + delta < theta_b | data)
  efficacious = (v[7] > efficacy_threshold)
  futile      = (v[7] < futility_threshold)
  
  #loop while we haven't treated all patients and while we haven't established efficacy/futility
  while(v[1] < N & !efficacious & !futile)
  {    
    v = update.evaluate(theta_a, theta_b, delta, integral_tolerance, how_many = how_often, ov = v)
    
    efficacious = (v[7] > efficacy_threshold)
    futile      = (v[7] < futility_threshold)
  }
  #Did we stop before treating all of the patients?
  stopped_early = (v[1] < N)
  
  #best estimate for theta_a, given the prior and the data we obtained before stopping
  theta_a_hat = v[3]/(v[3] + v[4])
  #best estimate for theta_b, given the prior and the data we obtained before stopping
  theta_b_hat = v[5]/(v[5] + v[6])
  
  #remember, v[1] is the number of patients treated
  #and v[7] is P(theta_a + delta < theta_b | data) given all the data up to the point we stopped
  return(c(theta_a_hat, theta_b_hat, efficacious, futile, stopped_early, v[1], v[7]))
}


# simulatetrials() simulates B clinical trials
# For now, I use only 3 cores because my laptop isn't _that_ powerful.
simulatetrials <- function(theta_a, theta_b, prior, B, delta, second_parameter = 5, N = 100,
                           efficacy_threshold = 0.95, futility_threshold = 0.05,
                           integral_tolerance = 1e-5, how_often = 5, cores = 3)
{
  require(multicore)#required for mclapply()
  seeds = sample(1:(B*10), B)#for the purposes of replication
  cat("The seed is set.\n")
  #apply the clinicaltrial() function to 'seeds', passing in all the other necessary variables
  dat = mclapply(seeds, clinicaltrial, theta_a, theta_b, prior, second_parameter, N, efficacy_threshold,
                 futility_threshold, integral_tolerance, how_often, delta, mc.cores = cores)
  cat("mclapply successfully excecuted.\n")
  #format the data into something I am more comfortable with
  x = as.data.frame(matrix(unlist(dat), ncol=7, byrow=TRUE))
  cat("Data frame created.\n")
  colnames(x) = c("placebo", "treatment", "efficacy", "futility", "early", "n", "probability")
  cat("Columns renamed.\n")
  return(x)
}

# This is where the magic happens (where distributions are updated and P(theta_a + delta < theta_b | data)
# is evaluated).
# pass in the true theta_a and theta_b, delta, integral tolerance,
# the number of patients before next evaluation, and old values
update.evaluate <- function(theta_a, theta_b, delta, integral_tolerance, how_many, ov)
{
  require(pracma)#required for dblquad()
  nv = rep(NA, length(ov))#the updated values that will be returned by this function (new values)
  
  #simulate the random assignments
  assigned_to_a = rbinom(1, how_many, ov[2])
  assigned_to_b = how_many - assigned_to_a
  nv[1] = ov[1] + how_many
  
  #of those assigned to treatment a, simulate the number of "successes" (recovery)
  sum_a = rbinom(1, assigned_to_a, theta_a)
  sum_b = rbinom(1, assigned_to_b, theta_b)
  
  #update parameters a_1, a_2, b_1, b_2
  nv[3] = ov[3] + sum_a
  nv[4] = ov[4] + assigned_to_a - sum_a
  nv[5] = ov[5] + sum_b
  nv[6] = ov[6] + assigned_to_b - sum_b
  
  #estimate the treatment effects
  theta_a_hat = nv[3]/(nv[3] + nv[4])
  theta_b_hat = nv[5]/(nv[5] + nv[6])
  
  #update the probability of being assigned to treatment group a (the control group)
  nv[2] = theta_a_hat/(theta_a_hat + theta_b_hat)
  
  #calculate the probability that theta_a + delta < theta_b, given our current distributions for them
  #integrand is joint density for theta_a and theta_b, defined on theta_a + delta < theta_b
  #theta_a is denoted as x, while theta_b is denoted as y
  
  integrand = function(x, y)
    (1/beta(nv[3], nv[4]))*(1/beta(nv[5], nv[6]))*x^(nv[3] - 1)*(1 - x)^(nv[4] - 1)*y^(nv[5] - 1)*(1 - y)^(nv[6] - 1)*(x + delta < y)
  
  nv[7] = dblquad(integrand, 0, 1, 0, 1, tol = integral_tolerance)
  
  return(nv)
}

# This simply gives a certain kind of summary of my simulations.
simsum <- function(x)
{
  est_placebo =   mean(x[,1])#average estimated placebo effect
  est_treatment = mean(x[,2])#average estimated treatment effect
  n_efficacy =    sum( x[,3])#number of trials that ended with perceived efficacy
  n_futility =    sum( x[,4])#number of trials that ended with perceived futility
  n_early =       sum( x[,5])#number of trials that stopped early
  avg_patients =  mean(x[,6])#average number of patients per trial
  mean_p =        mean(x[,7])#mean P(theta_a < theta_b | data)
  
  thesummary = data.frame(rbind(c(est_placebo, est_treatment, n_efficacy, n_futility, n_early,
                                  avg_patients, mean_p)))
  colnames(thesummary) = c("est_placebo", "est_treatment", "n_efficacy", "n_futility", "n_early",
                           "avg_patients", "mean_p")
  return(thesummary)
}