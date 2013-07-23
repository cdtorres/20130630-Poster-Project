#simulatetrials(theta_a = .2, theta_b = .2, prior = .2, B = 1000, how_often = 1, delta = .2)
#took 12 minutes 22 seconds. 946/1000 trials stopped early.

B = 1000
how_often = 1
set.seed(619)#for the purposes of replication
a = Sys.time()#to time my code
x = simulatetrials(theta_a = .2, theta_b = .2, prior = .2, B = B, how_often = how_often, delta = .2)
b = Sys.time()#to time my code
b - a#amount of time my code took to run


dim(unique(as.matrix(x)))[1]#number of unique simulations
table(x[,6])#table of the number of patients used in the clinical trials
simsum(x)
