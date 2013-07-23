#theta_a is the true rate for the placebo group
#theta_b is the true rate for the treatment group
#prior is the prior mean
#for a Beta(a, b) distribution, b is the second_parameter
#N is the maximum number of patients (for the meantime, let's keep this in multiples of 10)
#integral_tolerance is for the two-dimensional integral when evalualing P(theta_a < theta_b | data)
#making integral_tolerance small greatly increases the time that a simulation takes to run
clinicaltrial <-function(seed, theta_a, theta_b, prior, second_parameter, N, efficacy_threshold,
                         futility_threshold, integral_tolerance, how_often, delta)
{
  #require(cubature)
  require(pracma)
  
  efficacious = F
  futile = F
  stopped_early = F
  
  v = rep(NA, 7)
  #number of patients treated so far
  #n = 0
  v[1] = 0
  
  #group a is control, group b is treatment
  #For the control,   theta_a is distributed Beta(a_1, a_2)
  #For the treatment, theta_b is distributed Beta(b_1, b_2)
  
  #The higher the second parameter the more informative the prior is. (This needs to be > 0.)
  
  #theta_a and theta_b are the "true" probabilities of recovering with placebo and treatment, respectively
  #under the null hypothesis, they're equal
  
  #a_2 = second_parameter
  #b_2 = second_parameter
  v[4] = second_parameter
  v[6] = second_parameter
  
  #a_1 = a_2*prior/(1 - prior)
  #b_1 = b_2*prior/(1 - prior)
  v[3] = v[4]*prior/(1 - prior)
  v[5] = v[6]*prior/(1 - prior)
  
  #initial randomization probabilities
  #probability of being assigned to treatment group a (the control group)
  #p_a = 1/2
  #p_b = 1 - p_a
  v[2] = 1/2
  
  set.seed(seed)
#   #simulate the random assignments
#   assigned_to_a = rbinom(1, 10, p_a)
#   assigned_to_b = 10 - assigned_to_a
#   n = n + 10
#   
#   #of those assigned to treatment a, simulate the number of "successes" (recovery)
#   sum_a = rbinom(1, assigned_to_a, theta_a)
#   sum_b = rbinom(1, assigned_to_b, theta_b)
#   
#   #update the parameters (find posterior distributions of theta_a and theta_b)
#   a_1 = a_1 + sum_a
#   a_2 = a_2 + assigned_to_a - sum_a
#   b_1 = b_1 + sum_b
#   b_2 = b_2 + assigned_to_b - sum_b
#   
#   #estimate the treatment effects
#   theta_a_hat = a_1/(a_1 + a_2)
#   theta_b_hat = b_1/(b_1 + b_2)
#   
#   #update the probability of being assigned to treatment group a (the control group)
#   p_a = theta_a_hat/(theta_a_hat + theta_b_hat)
#   
#   #calculate the probability that theta_a < theta_b, given our current distributions for them
#   #joint density for theta_a and theta_b, defined on theta_a < theta_b
#   #theta_a is denoted as x[1], while theta_b is denoted as x[2]
#   
#   #integrand = function(x)
#   #  (1/beta(a_1, a_2))*(1/beta(b_1, b_2))*x[1]^(a_1 - 1)*(1 - x[1])^(a_2 - 1)*x[2]^(b_1 - 1)*(1 - x[2])^(b_2 - 1)*(x[1] < x[2])
#   integrand = function(x, y)
#     (1/beta(a_1, a_2))*(1/beta(b_1, b_2))*x^(a_1 - 1)*(1 - x)^(a_2 - 1)*y^(b_1 - 1)*(1 - y)^(b_2 - 1)*(x + delta < y)
#   
#   #a_less_than_b = adaptIntegrate(integrand, c(0, 0), c(1, 1), tol = integral_tolerance)
#   a_less_than_b = dblquad(integrand, 0, 1, 0, 1, tol = integral_tolerance)

  v = update.evaluate(theta_a, theta_b, delta, integral_tolerance, how_many = 10, ov = v)
  
  #efficacious = (a_less_than_b$integral > efficacy_threshold)
  #futile      = (a_less_than_b$integral < futility_threshold)
  #efficacious = (a_less_than_b > efficacy_threshold)
  #futile      = (a_less_than_b < futility_threshold)
  
  efficacious = (v[7] > efficacy_threshold)
  futile      = (v[7] < futility_threshold)
  
  while(v[1] < N & !efficacious & !futile)
  {
#     #simulate the random assignments
#     assigned_to_a = rbinom(1, how_often, p_a)
#     assigned_to_b = how_often - assigned_to_a
#     n = n + how_often
#     
#     #of those assigned to treatment a, simulate the number of "successes" (recovery)
#     sum_a = rbinom(1, assigned_to_a, theta_a)
#     sum_b = rbinom(1, assigned_to_b, theta_b)
#     
#     #update the parameters (find posterior distributions of theta_a and theta_b)
#     a_1 = a_1 + sum_a
#     a_2 = a_2 + assigned_to_a - sum_a
#     b_1 = b_1 + sum_b
#     b_2 = b_2 + assigned_to_b - sum_b
#     
#     #estimate the treatment effects
#     theta_a_hat = a_1/(a_1 + a_2)
#     theta_b_hat = b_1/(b_1 + b_2)
#     
#     #update the probability of being assigned to treatment group a (the control group)
#     p_a = theta_a_hat/(theta_a_hat + theta_b_hat)
#     
#     #calculate the probability that theta_a < theta_b, given our current distributions for them
#     #joint density for theta_a and theta_b, defined on theta_a < theta_b
#     #theta_a is denoted as x[1], while theta_b is denoted as x[2]
#     
#     #integrand = function(x)
#     #  (1/beta(a_1, a_2))*(1/beta(b_1, b_2))*x[1]^(a_1 - 1)*(1 - x[1])^(a_2 - 1)*x[2]^(b_1 - 1)*(1 - x[2])^(b_2 - 1)*(x[1] < x[2])
#     integrand = function(x, y)
#       (1/beta(a_1, a_2))*(1/beta(b_1, b_2))*x^(a_1 - 1)*(1 - x)^(a_2 - 1)*y^(b_1 - 1)*(1 - y)^(b_2 - 1)*(x + delta < y)
#     
#     #a_less_than_b = adaptIntegrate(integrand, c(0, 0), c(1, 1), tol = integral_tolerance)
#     a_less_than_b = dblquad(integrand, 0, 1, 0, 1, tol = integral_tolerance)
    
    v = update.evaluate(theta_a, theta_b, delta, integral_tolerance, how_many = how_often, ov = v)
    
    #efficacious = (a_less_than_b$integral > efficacy_threshold)
    #futile      = (a_less_than_b$integral < futility_threshold)
    #efficacious = (a_less_than_b > efficacy_threshold)
    #futile      = (a_less_than_b < futility_threshold)
    efficacious = (v[7] > efficacy_threshold)
    futile      = (v[7] < futility_threshold)
  }
  stopped_early = (v[1] < N)
  
  theta_a_hat = v[3]/(v[3] + v[4])
  theta_b_hat = v[5]/(v[5] + v[6])
  
  #return(c(theta_a_hat, theta_b_hat, efficacious, futile, stopped_early, n, a_less_than_b$integral))
  #return(c(theta_a_hat, theta_b_hat, efficacious, futile, stopped_early, n, a_less_than_b))
  return(c(theta_a_hat, theta_b_hat, efficacious, futile, stopped_early, v[1], v[7]))
}



simulatetrials <- function(theta_a, theta_b, prior, B, delta, second_parameter = 5, N = 100,
                           efficacy_threshold = 0.95, futility_threshold = 0.05,
                           integral_tolerance = 1e-5, how_often = 5)
{
  require(multicore)
  seeds = sample(1:(B*10), B)
  cat("The seed is set.\n")
  dat = mclapply(seeds, clinicaltrial, theta_a, theta_b, prior, second_parameter, N, efficacy_threshold,
                 futility_threshold, integral_tolerance, how_often, delta, mc.cores = 3)
  cat("mclapply successfully excecuted.\n")
  x = as.data.frame(matrix(unlist(dat), ncol=7, byrow=TRUE))
  cat("Data frame created.\n")
  colnames(x) = c("placebo", "treatment", "efficacy", "futility", "early", "n", "probability")
  cat("Columns renamed.\n")
  return(x)
}

#pass in the true theta_a and theta_b, delta, integral tolerance,
#the number of patients before next evaluation, and old variables
update.evaluate <- function(theta_a, theta_b, delta, integral_tolerance, how_many, ov)
{
  nv = rep(NA, length(ov))#the updated values that will be returned by this function (new variables)
  
  #simulate the random assignments
  assigned_to_a = rbinom(1, how_many, ov[2])
  assigned_to_b = how_many - assigned_to_a
  nv[1] = ov[1] + how_many
  
  #of those assigned to treatment a, simulate the number of "successes" (recovery)
  sum_a = rbinom(1, assigned_to_a, theta_a)
  sum_b = rbinom(1, assigned_to_b, theta_b)
  
  #update the parameters (find posterior distributions of theta_a and theta_b)
  nv[3] = ov[3] + sum_a
  nv[4] = ov[4] + assigned_to_a - sum_a
  nv[5] = ov[5] + sum_b
  nv[6] = ov[6] + assigned_to_b - sum_b
  
  #estimate the treatment effects
  theta_a_hat = nv[3]/(nv[3] + nv[4])
  theta_b_hat = nv[5]/(nv[5] + nv[6])
  
  #update the probability of being assigned to treatment group a (the control group)
  nv[2] = theta_a_hat/(theta_a_hat + theta_b_hat)
  
  #calculate the probability that theta_a < theta_b, given our current distributions for them
  #joint density for theta_a and theta_b, defined on theta_a < theta_b
  #theta_a is denoted as x[1], while theta_b is denoted as x[2]
  
  #integrand = function(x)
  #  (1/beta(nv[3], nv[4]))*(1/beta(nv[5], nv[6]))*x[1]^(nv[3] - 1)*(1 - x[1])^(nv[4] - 1)*x[2]^(nv[5] - 1)*(1 - x[2])^(nv[6] - 1)*(x[1] < x[2])
  integrand = function(x, y)
    (1/beta(nv[3], nv[4]))*(1/beta(nv[5], nv[6]))*x^(nv[3] - 1)*(1 - x)^(nv[4] - 1)*y^(nv[5] - 1)*(1 - y)^(nv[6] - 1)*(x + delta < y)
  
  #a_less_than_b = adaptIntegrate(integrand, c(0, 0), c(1, 1), tol = integral_tolerance)
  nv[7] = dblquad(integrand, 0, 1, 0, 1, tol = integral_tolerance)
  
  return(nv)
}

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