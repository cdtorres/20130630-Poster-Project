type = 'binary'
##################################################################
# note: the Beta distribution with mean 0.5 and second parameter 1
# is in fact the Uniform(0,1) distribution
if(type == 'continuous')
{
  # true arm effects
  theta_a = 0
  theta_b = 1.5
  var_a   = 1
  var_b   = 1
  
  prior = 0
  second_parameter = 1
  delta = 1.5
}
if(type != 'continuous')
{
  # true arm effects
  theta_a = 0.2
  theta_b = 0.6
  var_a   = NA
  var_b   = NA
  
  prior = 0.5
  second_parameter = 1
  delta = .4
}
##################################################################
# number of simulations
B = 1000
# after the initial 10 patients, how often do we want to update the posterior distributions?
how_often = 1
adapt_r = T
