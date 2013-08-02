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
clinicaltrial <-function(seed, theta_a, theta_b, var_a, var_b, prior, second_parameter, N,
                         efficacy_threshold, futility_threshold, integral_tolerance, how_often, delta,
                         type, initial)
{
  if(type == 'binary')
  {
    var_a = NULL
    var_b = NULL
  }
  efficacious = F
  futile = F
  stopped_early = F#indicator for if the trial stopped before all N patients were treated
  
  #This is a vector of values that we're interested in. The seven values in order are:
  #1 number of patients treated so far
  #2 probability of being assigned to treatment group a (the placebo group)
  #3 current value of a_1, given the prior and the data obtained thus far (m_a for continuous)
  #4 current value of a_2, given the prior and the data obtained thus far (v_a for continuous)
  #5 current value of b_1, given the prior and the data obtained thus far (m_b for continuous)
  #6 current value of b_2, given the prior and the data obtained thus far (v_b for continuous)
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
  
  v[4] = second_parameter#initialization of a_2 (v_a), given priors
  v[6] = second_parameter#initialization of b_2 (v_b), given priors
  
  if(type == 'binary')
  {
    v[3] = v[4]*prior/(1 - prior)#initialization of a_1, given priors
    v[5] = v[6]*prior/(1 - prior)#initialization of b_1, given priors
  }
  else if(type == 'continuous')
  {
    v[3] = prior
    v[5] = prior
  }
  
  #initial randomization probability, set to 1/2
  #probability of being assigned to treatment group a (the control group)
  v[2] = 1/2
  
  set.seed(seed)#for the purposes of replication

  #treat 10 patients before updating the randomization probability
  #or evaluating P(theta_a + delta < theta_b | data)
  if(initial <= N)
    v = update.evaluate(theta_a = theta_a, theta_b = theta_b, var_a = var_a, var_b = var_b, delta = delta,
                        integral_tolerance = integral_tolerance, how_many = initial, ov = v, type = type,
                        N = N)
  else#you can't treat more than the maximum amount of patients!
    v = update.evaluate(theta_a = theta_a, theta_b = theta_b, var_a = var_a, var_b = var_b, delta = delta,
                        integral_tolerance = integral_tolerance, how_many = N, ov = v, type = type, N = N)
  
  #remember, v[7] is P(theta_a + delta < theta_b | data)
  efficacious = (v[7] > efficacy_threshold)
  futile      = (v[7] < futility_threshold)
  
  #loop while we haven't treated all patients and while we haven't established efficacy/futility
  while(v[1] < N & !efficacious & !futile)
  {    
    if(v[1] + how_often <= N)
      v = update.evaluate(theta_a = theta_a, theta_b = theta_b, var_a = var_a, var_b = var_b,
                          delta = delta, integral_tolerance = integral_tolerance, how_many = how_often,
                          ov = v, type = type, N = N)
    else#you can't treat more than the maximum amount of patients!
      v = update.evaluate(theta_a = theta_a, theta_b = theta_b, var_a = var_a, var_b = var_b,
                          delta = delta, integral_tolerance = integral_tolerance, how_many = (N - v[1]),
                          ov = v, type = type, N = N)
    
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
simulatetrials <- function(theta_a, theta_b, var_a = NULL, var_b = NULL, prior, B, delta, type,
                           second_parameter = 5, N = 100, efficacy_threshold = 0.95,
                           futility_threshold = 0.05, integral_tolerance = 1e-5, how_often = 5, cores = 3,
                           gseed = 619, initial = 10)
{
  set.seed(gseed)
  if(!(type %in% c('binary', 'continuous')))
    stop("Please put in 'binary' or 'continuous' for the variable 'type'.")
  require(multicore)#required for mclapply()
  seeds = sample(1:(B*10), B)#for the purposes of replication
  #cat("The seed is set.\n")
  #apply the clinicaltrial() function to 'seeds', passing in all the other necessary variables
  if(type == 'binary')
  {
    dat = mclapply(seeds, clinicaltrial, theta_a, theta_b, NULL, NULL, prior, second_parameter, N,
                   efficacy_threshold, futility_threshold, integral_tolerance, how_often, delta, type,
                   initial, mc.cores = cores)
  }
  else if (type == 'continuous')
  {
    dat = mclapply(seeds, clinicaltrial, theta_a, theta_b, var_a, var_b, prior, second_parameter, N,
                   efficacy_threshold, futility_threshold, integral_tolerance, how_often, delta, type,
                   initial, mc.cores = cores)
  }
  #cat("mclapply successfully excecuted.\n")
  #format the data into something I am more comfortable with
  x = as.data.frame(matrix(unlist(dat), ncol=7, byrow=TRUE))
  #cat("Data frame created.\n")
  colnames(x) = c("placebo", "treatment", "efficacy", "futility", "early", "n", "probability")
  #cat("Columns renamed.\n")
  return(x)
}

# This is where the magic happens (where distributions are updated and P(theta_a + delta < theta_b | data)
# is evaluated).
# pass in the true theta_a and theta_b, delta, integral tolerance,
# the number of patients before next evaluation, and old values
update.evaluate <- function(theta_a, theta_b, var_a, var_b, delta, integral_tolerance,
                            how_many, ov, type, N)
{
  #require(pracma)#required for dblquad()
  nv = rep(NA, length(ov))#the updated values that will be returned by this function (new values)
  
  #simulate the random assignments (this is the same regardless of outcome type)
  assigned_to_a = rbinom(1, how_many, ov[2])
  assigned_to_b = how_many - assigned_to_a
  nv[1] = ov[1] + how_many
  
  if(type == 'binary')
  {
    #of those assigned to treatment a, simulate the number of "successes" (recovery)
    sum_a = rbinom(1, assigned_to_a, theta_a)
    sum_b = rbinom(1, assigned_to_b, theta_b)
  }
  else if(type == 'continuous')
  {
    sum_a = sum(rnorm(assigned_to_a, theta_a, var_a))
    sum_b = sum(rnorm(assigned_to_b, theta_b, var_b))
  }
  
  if(type == 'binary')
  {
    #update parameters a_1, a_2, b_1, b_2
    nv[3] = ov[3] + sum_a
    nv[4] = ov[4] + assigned_to_a - sum_a
    nv[5] = ov[5] + sum_b
    nv[6] = ov[6] + assigned_to_b - sum_b
  }
  else if(type == 'continuous')
  {
    #update parameters m_a, v_a, m_b, v_b
    nv[4] = (1/ov[4] + assigned_to_a/var_a)^-1
    nv[3] = (ov[3]/ov[4] + sum_a/var_a)*nv[4]
    nv[6] = (1/ov[6] + assigned_to_b/var_b)^-1
    nv[5] = (ov[5]/ov[6] + sum_b/var_b)*nv[6]
  }
  
  if(type == 'binary')
  {
    #estimate the treatment effects
    theta_a_hat = nv[3]/(nv[3] + nv[4])
    theta_b_hat = nv[5]/(nv[5] + nv[6])
  }
  else if(type == 'continuous')
  {
    #estimate the treatment effects
    theta_a_hat = nv[3]
    theta_b_hat = nv[5]
  }
  
  #calculate the probability that theta_a + delta < theta_b, given our current distributions for them
  #integrand is joint density for theta_a and theta_b, defined on theta_a + delta < theta_b
  #theta_a is denoted as x, while theta_b is denoted as y
  
  #integrand = function(x, y)
  #  (1/beta(nv[3], nv[4]))*(1/beta(nv[5], nv[6]))*x^(nv[3] - 1)*(1 - x)^(nv[4] - 1)*y^(nv[5] - 1)*(1 - y)^(nv[6] - 1)*(x + delta < y)
  
  #nv[7] = dblquad(integrand, 0, 1, 0, 1, tol = integral_tolerance)
  
  if(type == 'binary')
  {
    integrand = function(y)
    {
      pbeta((y - delta), nv[3], nv[4])*dbeta(y, nv[5], nv[6])
    }
  }
  else if(type == 'continuous')
  {
    integrand = function(y)
    {
      pnorm((y - delta), nv[3], nv[4])*dnorm(y, nv[5], nv[6])
    }
  }
  
  if(type == 'binary')
  {
    nv[7] = integrate(integrand, 0, 1)$value
  }
  else if(type == 'continuous')
  {
    #We're calculating P(X + delta < Y | data) here. By defining W = X - Y and standardizing this new
    #variable, and calling it Z, we can find the probability of P(Z < z | data) instead, using the z
    #defined below
    z = (nv[5] - nv[3] - delta)/sqrt(nv[4] + nv[6])
    nv[7] = pnorm(z)
  }
  
  #update the probability of being assigned to treatment group a (the control group)
  if(type == 'binary')
  {
    #nv[2] = theta_a_hat/(theta_a_hat + theta_b_hat)
    r_integrand = function(y)
    {
      pbeta(y, nv[3], nv[4])*dbeta(y, nv[5], nv[6])
    }
    r_integral = integrate(r_integrand, 0, 1)$value
    c = nv[1]/(2*N)
    r_1 = (1 - r_integral)^c
    r_2 = (r_integral)^c
    nv[2] = r_1/(r_1 + r_2)
  }
  else if(type == 'continuous')
  {
    
    #going off of page 156
    #We're calculating P(X < Y | data) here. By defining W = X - Y and standardizing this new variable,
    #and call it Z, we can find the probability of P(Z < r_z | data) instead, using the r_z defined below
    r_z = (nv[5] - nv[3])/sqrt(nv[4] + nv[6])
    r_integral = pnorm(r_z)
    c = nv[1]/(2*N)
    r_1 = (1 - r_integral)^c
    r_2 = (r_integral)^c
    nv[2] = r_1/(r_1 + r_2)
  }
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







# Pass in the true theta_a and theta_b (also var_a and var_b if continuous), prior, number of simulations
# per evaluation, delta, type of data, second parameter for the prior information, maximum number of
# patients per trial, starting efficacy and futility thresholds, integral tolerance (for binary data),
# how often to update the posterior and randomization probabilities (in terms of number of patients),
# number of cores to use when running this function, number of patients in each trial before starting to
# adapt, desired errors, max outer while loops, and an indicator as to whether report the results after
# each loop.
get.thresholds <- function(theta_a, theta_b, var_a = NULL, var_b = NULL, prior, B, delta, type,
                           second_parameter = 5, N = 100, efficacy_threshold = 1,
                           futility_threshold = 0, integral_tolerance = 1e-5, how_often = 5, cores = 3,
                           initial = 10, desired_type_1_error = 0.05, desired_type_2_error = 0.05,
                           maxloops = 5, report = F)
{
  type_1_error = 1#initial type 1 error
  type_2_error = 1#initial type 2 error
  
  #This is an indicator for whether something changed in the previous outer loop. If it is set to true,
  #then it goes through the outer loop one more time to make sure that everything has "settled" and/or
  #"converged".
  something_changed = T
  loops = 0#number of loops
  
  #As long as the type_1 or type_2 error isn't where we want it to be, or if something changed in the
  #last loop, we want to loop through again. However, on the off-chance that this condition is always
  #true for some reason, we don't want it to loop forever, so we set a maximum number of loops. If when
  #the function finishes evaluating, it took the maximum number of loops to finish, we may want to
  #investigate more as to why this happened.
  while((type_1_error > desired_type_1_error | type_2_error > desired_type_2_error |
           something_changed == T) & loops < maxloops)
  {
    #The increment is 0.1 for the first loop, 0.01 for the second loop, 0.001 for the third loop, ...
    increment = 10^(-(loops + 1))
    something_changed = F
    #Get the type 1 error, by simulating under the null hypothesis.
    x = simulatetrials(theta_a = theta_a, theta_b = theta_a, var_a = var_a, var_b = var_b,
                       prior = prior, second_parameter = second_parameter, B = B,
                       how_often = how_often, delta = delta, type = type,
                       efficacy_threshold = efficacy_threshold,
                       futility_threshold = futility_threshold, initial = initial)
    type_1_error = as.numeric(simsum(x)[3]/B)#initial type 1 error
    if(type_1_error > desired_type_1_error)
    {
      #If the type 1 error is too high, we increase the efficacy threshold until the type 1 error is
      #low enough. But we can't let the efficacy threshold be more than 1, because that makes no sense!
      while(type_1_error > desired_type_1_error & efficacy_threshold <= 1)
      {
        efficacy_threshold = efficacy_threshold + increment
        x = simulatetrials(theta_a = theta_a, theta_b = theta_a, var_a = var_a, var_b = var_b,
                           prior = prior, second_parameter = second_parameter, B = B,
                           how_often = how_often, delta = delta, type = type,
                           efficacy_threshold = efficacy_threshold,
                           futility_threshold = futility_threshold, initial = initial)
        type_1_error = as.numeric(simsum(x)[3]/B)
        something_changed = T#indicate that there was a change
      }
    }
    else
    {
      #If the type 1 error is already at an acceptable level, we check to see if we can decrease the
      #efficacy threshold while still maintaining an acceptable type 1 error.
      
      #This flag is to indicate that we successfully decreased the efficacy threshold while maintaing
      #an acceptable type 1 error. If we were successful, we want to loop through and check one more
      #time. If we weren't successful, the flag will remain at F, and we will stop looping through
      #because we found the lowest efficacy threshold such that the type 1 error is still at an
      #acceptable level.
      something_changed_2 = T
      while(type_1_error < desired_type_1_error & efficacy_threshold >= 0 & something_changed_2 == T)
      {
        something_changed_2 = F
        
        #We store our results in some temporary variables.
        potential_efficacy_threshold = efficacy_threshold - increment
        potential_x = simulatetrials(theta_a = theta_a, theta_b = theta_a, var_a = var_a, var_b = var_b,
                                     prior = prior, second_parameter = second_parameter, B = B,
                                     how_often = how_often, delta = delta, type = type,
                                     efficacy_threshold = potential_efficacy_threshold,
                                     futility_threshold = futility_threshold, initial = initial)
        potential_type_1_error = as.numeric(simsum(potential_x)[3]/B)
        #This is where we check if the type 1 error is still at an acceptable level.
        #If we so, we overwrite our regular variables with the temporary ones, and indicate that we need
        #to loop through at least once more, for both the inner and outer while loops.
        if(potential_type_1_error <= desired_type_1_error & potential_efficacy_threshold <= 1)
        {
          x = potential_x
          efficacy_threshold = potential_efficacy_threshold
          type_1_error = potential_type_1_error
          something_changed   = T
          something_changed_2 = T
        }
      }
    }
    
    #Now we are going to work with the type 2 error
    #Get type 2 error, by simulating under the alternative hypothesis
    
    x = simulatetrials(theta_a = theta_a, theta_b = theta_b, var_a = var_a, var_b = var_b,
                       prior = prior, second_parameter = second_parameter, B = B,
                       how_often = how_often, delta = delta, type = type,
                       efficacy_threshold = efficacy_threshold,
                       futility_threshold = futility_threshold, initial = initial)
    type_2_error = as.numeric(simsum(x)[4]/B)
    
    
    if(type_2_error > desired_type_2_error)
    {
      #If the type 2 error is too high, we decrease the futility threshold until the type 2 error is
      #at a more acceptable level, or until the futility threshold reaches 0.
      while(type_2_error > desired_type_2_error & futility_threshold >= 0)
      {
        futility_threshold = futility_threshold - increment
        x = simulatetrials(theta_a = theta_a, theta_b = theta_b, var_a = var_a, var_b = var_b,
                           prior = prior, second_parameter = second_parameter, B = B,
                           how_often = how_often, delta = delta, type = type,
                           efficacy_threshold = efficacy_threshold,
                           futility_threshold = futility_threshold, initial = initial)
        type_2_error = as.numeric(simsum(x)[4]/B)
        something_changed = T#indicate that something has changed
      }
    }
    else
    {
      #If the type 2 error is already at an acceptable level, we check to see if we can increase the
      #futility threshold while still maintaining an acceptable type 2 error.
      
      #This flag is to indicate that we successfully increased the futility threshold while maintaing
      #an acceptable type 2 error. If we were successful, we want to loop through and check one more
      #time. If we weren't successful, the flag will remain at F, and we will stop looping through
      #because we found the highest futility threshold such that the type 2 error is still at an
      #acceptable level.
      something_changed_2 = T
      while(type_2_error < desired_type_2_error & futility_threshold <= 1 & something_changed_2 == T)
      {
        something_changed_2 = F
        #We store our results in some temporary variables.
        potential_futility_threshold = futility_threshold + increment
        potential_x = simulatetrials(theta_a = theta_a, theta_b = theta_b, var_a = var_a, var_b = var_b,
                                     prior = prior, second_parameter = second_parameter, B = B,
                                     how_often = how_often, delta = delta, type = type,
                                     efficacy_threshold = efficacy_threshold,
                                     futility_threshold = potential_futility_threshold, initial = initial)
        potential_type_2_error = as.numeric(simsum(potential_x)[4]/B)
        #This is where we check if the type 2 error is still at an acceptable level.
        #If we so, we overwrite our regular variables with the temporary ones, and indicate that we need
        #to loop through at least once more, for both the inner and outer while loops.
        if(potential_type_2_error <= desired_type_2_error)
        {
          x = potential_x
          futility_threshold = potential_futility_threshold
          type_2_error = potential_type_2_error
          something_changed   = T
          something_changed_2 = T
        }
      }
    }
    
    if(report)
    {
      #Just printing out some stuff to state what the results where after each outer loop.
      cat("Loop ")
      cat(loops + 1, "efficacy: ")
      cat(efficacy_threshold, ", ")
      cat("type 1 error: ")
      cat(type_1_error, ", futility: ")
      cat(futility_threshold, ", ")
      cat("type 2 error: ")
      cat(type_2_error, ".\n")
    }
    
    #I want the while loop to get at least to increments of 0.01, and I force it in this manner.
    if(loops == 0)
      something_changed = T
    
    loops = loops + 1
  }
  
  
  #Get final type 1 error
  x_null = simulatetrials(theta_a = theta_a, theta_b = theta_a, var_a = var_a, var_b = var_b,
                          prior = prior, second_parameter = second_parameter, B = B,
                          how_often = how_often, delta = delta, type = type,
                          efficacy_threshold = efficacy_threshold,
                          futility_threshold = futility_threshold, initial = initial)
  final_type_1_error = as.numeric(simsum(x_null)[3]/B)
  
  #Get final type 2 error
  x_alt = simulatetrials(theta_a = theta_a, theta_b = theta_b, var_a = var_a, var_b = var_b,
                         prior = prior, second_parameter = second_parameter, B = B,
                         how_often = how_often, delta = delta, type = type,
                         efficacy_threshold = efficacy_threshold,
                         futility_threshold = futility_threshold, initial = initial)
  final_type_2_error = as.numeric(simsum(x_alt)[4]/B)
  
  #Store the sumsum(x_alt) results in vector form.
  simsum_x_alt = c(simsum(x_alt)[1,1], simsum(x_alt)[1,2], simsum(x_alt)[1,3], simsum(x_alt)[1,4],
                   simsum(x_alt)[1,5], simsum(x_alt)[1,6], simsum(x_alt)[1,7])
  
  #Create data frame that will be returned.
  values_of_interest = c(final_type_1_error, final_type_2_error, efficacy_threshold, futility_threshold,
                         theta_a, simsum_x_alt[1], theta_b, simsum_x_alt[2:7], loops)
  values_of_interest = as.data.frame(rbind(values_of_interest))
  colnames(values_of_interest) = c("type_1_error", "type_2_error", "efficacy_threshold",
                                   "futility_threshold", "placebo", "est_placebo",
                                   "treatment", "est_treatment", "n_efficacy", "n_futility", "n_early",
                                   "avg_patients", "mean_p", "loops")
  return(values_of_interest)#Return data frame
}


