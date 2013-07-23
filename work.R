#simulatetrials(theta_a = .2, theta_b = .2, prior = .2, B = 1000, how_often = 1, delta = .2)
#took almost 12 minutes. 946/1000 trials stopped early.

set.seed(619)
a = Sys.time()
x = simulatetrials(theta_a = .2, theta_b = .2, prior = .2, B = 1000, how_often = 1, delta = .2)
b = Sys.time()
b - a
#dat = as.data.frame(matrix(unlist(x), ncol=7, byrow=TRUE))
#colnames(dat) = c("placebo", "treatment", "efficacy", "futility", "early", "n", "probability")
mean(x[,1])#average estimated placebo effect
mean(x[,2])#average estimated treatment effect
sum( x[,3])#number of trials that ended with perceived efficacy
sum( x[,4])#number of trials that ended with perceived futility
sum( x[,5])#number of trials that stopped early
mean(x[,6])#average number of patients per trial
mean(x[,7])#mean P(theta_a < theta_b | data)

dim(unique(as.matrix(x)))[1]
table(x[,6])