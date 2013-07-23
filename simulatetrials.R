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
  
  #number of patients treated so far
  n = 0
  
  #group a is control, group b is treatment
  #For the control,   theta_a is distributed Beta(a_1, a_2)
  #For the treatment, theta_b is distributed Beta(b_1, b_2)
  
  #The higher the second parameter the more informative the prior is. (This needs to be > 0.)
  
  #theta_a and theta_b are the "true" probabilities of recovering with placebo and treatment, respectively
  #under the null hypothesis, they're equal
  
  a_2 = second_parameter
  b_2 = second_parameter
  
  a_1 = a_2*prior/(1 - prior)
  b_1 = b_2*prior/(1 - prior)
  
  #initial randomization probabilities
  #probability of being assigned to treatment group a (the control group)
  p_a = 1/2
  #p_b = 1 - p_a
  
  set.seed(seed)
  #simulate the random assignments
  assigned_to_a = rbinom(1, 10, p_a)
  assigned_to_b = 10 - assigned_to_a
  n = n + 10
  
  #of those assigned to treatment a, simulate the number of "successes" (recovery)
  sum_a = rbinom(1, assigned_to_a, theta_a)
  sum_b = rbinom(1, assigned_to_b, theta_b)
  
  #update the parameters (find posterior distributions of theta_a and theta_b)
  a_1 = a_1 + sum_a
  a_2 = a_2 + assigned_to_a - sum_a
  b_1 = b_1 + sum_b
  b_2 = b_2 + assigned_to_b - sum_b
  
  #estimate the treatment effects
  theta_a_hat = a_1/(a_1 + a_2)
  theta_b_hat = b_1/(b_1 + b_2)
  
  #update the probability of being assigned to treatment group a (the control group)
  p_a = theta_a_hat/(theta_a_hat + theta_b_hat)
  
  #calculate the probability that theta_a < theta_b, given our current distributions for them
  #joint density for theta_a and theta_b, defined on theta_a < theta_b
  #theta_a is denoted as x[1], while theta_b is denoted as x[2]
  
  #integrand = function(x)
  #  (1/beta(a_1, a_2))*(1/beta(b_1, b_2))*x[1]^(a_1 - 1)*(1 - x[1])^(a_2 - 1)*x[2]^(b_1 - 1)*(1 - x[2])^(b_2 - 1)*(x[1] < x[2])
  integrand = function(x, y)
    (1/beta(a_1, a_2))*(1/beta(b_1, b_2))*x^(a_1 - 1)*(1 - x)^(a_2 - 1)*y^(b_1 - 1)*(1 - y)^(b_2 - 1)*(x + delta < y)
  
  #a_less_than_b = adaptIntegrate(integrand, c(0, 0), c(1, 1), tol = integral_tolerance)
  a_less_than_b = dblquad(integrand, 0, 1, 0, 1, tol = integral_tolerance)
  
  #efficacious = (a_less_than_b$integral > efficacy_threshold)
  #futile      = (a_less_than_b$integral < futility_threshold)
  efficacious = (a_less_than_b > efficacy_threshold)
  futile      = (a_less_than_b < futility_threshold)
  
  while(n < N & !efficacious & !futile)
  {
    #simulate the random assignments
    assigned_to_a = rbinom(1, how_often, p_a)
    assigned_to_b = how_often - assigned_to_a
    n = n + how_often
    
    #of those assigned to treatment a, simulate the number of "successes" (recovery)
    sum_a = rbinom(1, assigned_to_a, theta_a)
    sum_b = rbinom(1, assigned_to_b, theta_b)
    
    #update the parameters (find posterior distributions of theta_a and theta_b)
    a_1 = a_1 + sum_a
    a_2 = a_2 + assigned_to_a - sum_a
    b_1 = b_1 + sum_b
    b_2 = b_2 + assigned_to_b - sum_b
    
    #estimate the treatment effects
    theta_a_hat = a_1/(a_1 + a_2)
    theta_b_hat = b_1/(b_1 + b_2)
    
    #update the probability of being assigned to treatment group a (the control group)
    p_a = theta_a_hat/(theta_a_hat + theta_b_hat)
    
    #calculate the probability that theta_a < theta_b, given our current distributions for them
    #joint density for theta_a and theta_b, defined on theta_a < theta_b
    #theta_a is denoted as x[1], while theta_b is denoted as x[2]
    
    #integrand = function(x)
    #  (1/beta(a_1, a_2))*(1/beta(b_1, b_2))*x[1]^(a_1 - 1)*(1 - x[1])^(a_2 - 1)*x[2]^(b_1 - 1)*(1 - x[2])^(b_2 - 1)*(x[1] < x[2])
    integrand = function(x, y)
      (1/beta(a_1, a_2))*(1/beta(b_1, b_2))*x^(a_1 - 1)*(1 - x)^(a_2 - 1)*y^(b_1 - 1)*(1 - y)^(b_2 - 1)*(x + delta < y)
    
    #a_less_than_b = adaptIntegrate(integrand, c(0, 0), c(1, 1), tol = integral_tolerance)
    a_less_than_b = dblquad(integrand, 0, 1, 0, 1, tol = integral_tolerance)
    
    #efficacious = (a_less_than_b$integral > efficacy_threshold)
    #futile      = (a_less_than_b$integral < futility_threshold)
    efficacious = (a_less_than_b > efficacy_threshold)
    futile      = (a_less_than_b < futility_threshold)
  }
  stopped_early = (n < N)
  
  #return(c(theta_a_hat, theta_b_hat, efficacious, futile, stopped_early, n, a_less_than_b$integral))
  return(c(theta_a_hat, theta_b_hat, efficacious, futile, stopped_early, n, a_less_than_b))
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
