# true arm effects
theta_a = 0.2
theta_b = 1.1
var_a   = .5
var_b   = .5
type = 'continuous'
##################################################################
# note: the Beta distribution with mean 0.5 and second parameter 1
# is in fact the Uniform(0,1) distribution
if(type == 'continuous')
{
  prior = 0
  second_parameter = 1
}
if(type != 'continuous')
{
  prior = 0.5
  second_parameter = 1
}
##################################################################
# number of simulations
B = 10000
# after the initial 10 patients, how often do we want to update the posterior distributions?
how_often = 1
delta = 1


a = Sys.time()#to time my code
set.seed(619)#for the purposes of replication
x = simulatetrials(theta_a = theta_a, theta_b = theta_b, var_a = var_a, var_b = var_b, prior = prior,
                   second_parameter = second_parameter, B = B, how_often = how_often, delta = delta,
                   type = type)
b = Sys.time()#to time my code
b - a#amount of time my code took to run


dim(unique(as.matrix(x)))[1]#number of unique simulations
table(x[,6])#table of the number of patients used in the clinical trials
simsum(x)
