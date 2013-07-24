#x = simulatetrials(theta_a = .2, theta_b = .2, prior = .2, B = 1000, how_often = 1, delta = .2)
#took 12 minutes 22 seconds. 946/1000 trials stopped early.

theta_a = 0.2
theta_b = 0.2
#note: the Beta distribution with mean 0.5 and second parameter 1 is in fact the Uniform(0,1) distribution
prior = 0.5
second_parameter = 1
B = 10000
how_often = 1
delta = 0.2
a = Sys.time()#to time my code
set.seed(619)#for the purposes of replication
x = simulatetrials(theta_a = .2, theta_b = .2, prior = .2, B = 1000, how_often = 1, delta = .2)
#x = simulatetrials(theta_a = theta_a, theta_b = theta_b, prior = prior,
#                   second_parameter = second_parameter, B = B, how_often = how_often, delta = delta)
b = Sys.time()#to time my code
b - a#amount of time my code took to run


dim(unique(as.matrix(x)))[1]#number of unique simulations
table(x[,6])#table of the number of patients used in the clinical trials
simsum(x)
