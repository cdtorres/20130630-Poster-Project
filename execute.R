#indicators for if I want to simulate trials, or if I want to get thresholds
onesim = T
thresholds = T

if(type == 'binary')
{
  efficacy_threshold = 0.44
  futility_threshold = 0.03666
}
if(type == 'continuous')
{
  efficacy_threshold = 0.368
  futility_threshold = 0.101
}

a = Sys.time()#to time my code

if(thresholds)
{
  #Get thresholds required for good operating characteristics for these particular settings
  
  y = get.thresholds(theta_a = theta_a, theta_b = theta_b, var_a = var_a, var_b = var_b, prior = prior,
                     second_parameter = second_parameter, B = B, how_often = how_often, delta = delta,
                     type = type, desired_type_1_error = 0.05, desired_type_2_error = 0.20, report = T,
                     adapt_r = adapt_r)
  y
}

if(onesim & thresholds)
{
  efficacy_threshold = as.numeric(y[3])
  futility_threshold = as.numeric(y[4])
}

if(onesim)
{
  #set.seed(619)#for the purposes of replication
  x = simulatetrials(theta_a = theta_a, theta_b = theta_b, var_a = var_a, var_b = var_b, prior = prior,
                     second_parameter = second_parameter, B = B, how_often = how_often, delta = delta,
                     type = type, efficacy_threshold = efficacy_threshold,
                     futility_threshold = futility_threshold, adapt_r = adapt_r)
  x_null = simulatetrials(theta_a = theta_a, theta_b = theta_a, var_a = var_a, var_b = var_b, prior = prior,
                          second_parameter = second_parameter, B = B, how_often = how_often, delta = delta,
                          type = type, efficacy_threshold = efficacy_threshold,
                          futility_threshold = futility_threshold, adapt_r = adapt_r)
}

b = Sys.time()#to time my code

if(F)
{
  #dim(unique(as.matrix(x)))[1]#number of unique simulations
  #table(x[,6])#table of the number of patients used in the clinical trials
  b - a#amount of time my code took to run
  simsum(x)
  simsum(x_null)
  with(x, c(mean(placebo_0.025), mean(placebo_0.975)))
  with(x, c(mean(treatment_0.025), mean(treatment_0.975)))
}

