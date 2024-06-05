# install and load required libraries 
install.packages('R2OpenBUGS')
install.packages('coda')
install.packages('mcmcplots')
install.packages("rjags")
install.packages("runjags")
want = c("parallel","compute.es")
have = want %in% rownames(installed.packages())
if ( any(!have) ) { install.packages( want[!have] ) }

library(R2OpenBUGS) 
library(coda)
library(mcmcplots) 
library(runjags)
library(rjags) 
testjags()
library(parallel)
library(compute.es)

# create a dataset as a matrix 
C_PCR <- c(rep(1, 31), rep(0, 23)) 
length(C_PCR)
N_PCR <- c(rep(1, 21), rep(0, 10), rep(1, 11), rep(0,12))
length(N_PCR)
S_PCR <- c(rep(1,15),rep(0,6),rep(1,5),rep(0,5),rep(1,5),rep(0,6),rep(1,1),rep(0,11))
length(S_PCR)
HIS <- c(rep(1,6), rep(0,9), rep(1,2), rep(0,4), rep(1,0), rep(0,5), rep(1,0),rep(0,5),  rep(1,0),rep(0,5), rep(1,0),rep(0,6), rep(1,0), rep(0,1), rep(1,0), rep(0,11))

dat <- data.frame(C_PCR,N_PCR,S_PCR,HIS)
m.ca <- as.matrix(dat)

dump("m.ca") 

N <- 54
bov.beis <-
  structure(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
              1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 
              1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 
              0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 
              0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
              0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 
              1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), .Dim = c(54L, 
                                                                                 3L), .Dimnames = list(NULL, c("C_PCR", "N_PCR", "S_PCR")))
ones <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)

################################################################################################
## Bayesian latent-class model code for three diagnostic tests 
################################################################################################

#######################################################
##Definition of the variables in the model
#######################################################

hw_4t_1p <- c("var p[N], q[N,16], pr[N], L[N], checks[N,32]; 

# N    <- observations (N = 54 cattle)
# p    <- individual samples
# q    <- different combinations of test results
# prc  <- prevalence
# s    <- test sensitivities
# c    <- test specificities
# covs <- conditional dependency between tests sensitivities
# covc <- conditional dependency between tests specificities
# bov.beis <- dataset name

model {
  
  for(i in 1:N){
    
    q[i,1]<-prc*(s1*s2*s3*s4+covs12+covs13+covs23)+(1-prc)*((1-c1)*(1-c2)*(1-c3)*(1-c4)+covc12+covc13+covc23);
    q[i,2]<-prc*(s1*s2*s3*(1-s4)+covs12-covs13-covs23)+(1-prc)*((1-c1)*(1-c2)*(1-c3)*c4+covc12-covc13-covc23);
    q[i,3]<-prc*(s1*s2*(1-s3)*s4-covs12+covs13-covs23)+(1-prc)*((1-c1)*(1-c2)*c3*(1-c4)-covc12+covc13-covc23);
    q[i,4]<-prc*(s1*s2*(1-s3)*(1-s4)-covs12-covs13+covs23)+(1-prc)*((1-c1)*(1-c2)*c3*c4-covc12-covc13+covc23);
    q[i,5]<-prc*(s1*(1-s2)*s3*s4-covs12-covs13+covs23)+(1-prc)*((1-c1)*c2*(1-c3)*c4-covc12-covc13+covc23);
    q[i,6]<-prc*(s1*(1-s2)*s3*(1-s4)-covs12+covs13-covs23)+(1-prc)*((1-c1)*c2*(1-c3)*c4-covc12+covc13-covc23);
    q[i,7]<-prc*(s1*(1-s2)*(1-s3)*s4+covs12-covs13-covs23)+(1-prc)*((1-c1)*c2*c3*(1-c4)+covc12-covc13-covc23);
    q[i,8]<-prc*(s1*(1-s2)*(1-s3)*(1-s4)+covs12+covs13+covs23)+(1-prc)*((1-c1)*c2*c3*c4+covc12+covc13+covc23);
    q[i,9]<-prc*((1-s1)*s2*s3*s4+covs12+covs13+covs23)+(1-prc)*(c1*(1-c2)*(1-c3)*(1-c4)+covc12+covc13+covc23);
    q[i,10]<-prc*((1-s1)*s2*s3*(1-s4)+covs12-covs13-covs23)+(1-prc)*(c1*(1-c2)*(1-c3)*c4+covc12-covc13-covc23);
    q[i,11]<-prc*((1-s1)*s2*(1-s3)*s4-covs12+covs13-covs23)+(1-prc)*(c1*(1-c2)*c3*(1-c4)-covc12+covc13-covc23);
    q[i,12]<-prc*((1-s1)*s2*(1-s3)*(1-s4)-covs12-covs13+covs23)+(1-prc)*(c1*(1-c2)*c3*c4-covc12-covc13+covc23);
    q[i,13]<-prc*((1-s1)*(1-s2)*s3*s4-covs12-covs13+covs23)+(1-prc)*(c1*c2*(1-c3)*(1-c4)-covc12-covc13+covc23);
    q[i,14]<-prc*((1-s1)*(1-s2)*s3*(1-s4)-covs12+covs13-covs23)+(1-prc)*(c1*c2*(1-c3)*c4-covc12+covc13-covc23);
    q[i,15]<-prc*((1-s1)*(1-s2)*(1-s3)*s4+covs12-covs13-covs23)+(1-prc)*(c1*c2*c3*(1-c4)+covc12-covc13-covc23);
    q[i,16]<-prc*((1-s1)*(1-s2)*(1-s3)*(1-s4)+covs12+covs13+covs23)+(1-prc)*(c1*c2*c3*c4+covc12+covc13+covc23);

    #######################################################
    ## Check and correct potential errors of probabilities exceeding (0,1) bounds 
    #######################################################
    
    checks[i,1]<-   s1*s2*s3*s4+covs12+covs13+covs23;
    checks[i,2]<-   (1-c1)*(1-c2)*(1-c3)*(1-c4)+covc12+covc13+covc23;
    checks[i,3]<-  s1*s2*s3*(1-s4)+covs12-covs13-covs23;
    checks[i,4]<-  (1-c1)*(1-c2)*(1-c3)*c4+covc12-covc13-covc23;
    checks[i,5]<-  s1*s2*(1-s3)*s4-covs12+covs13-covs23;
    checks[i,6]<-  (1-c1)*(1-c2)*c3*(1-c4)-covc12+covc13-covc23;
    checks[i,7]<-  s1*s2*(1-s3)*(1-s4)-covs12-covs13+covs23;
    checks[i,8]<-  (1-c1)*(1-c2)*c3*c4-covc12-covc13+covc23;
    checks[i,9]<-  s1*(1-s2)*s3*s4-covs12-covs13+covs23;
    checks[i,10]<- (1-c1)*c2*(1-c3)*c4-covc12-covc13+covc23;
    checks[i,11]<- s1*(1-s2)*s3*(1-s4)-covs12+covs13-covs23;
    checks[i,12]<- (1-c1)*c2*(1-c3)*c4-covc12+covc13-covc23;
    checks[i,13]<- s1*(1-s2)*(1-s3)*s4+covs12-covs13-covs23;
    checks[i,14]<- (1-c1)*c2*c3*(1-c4)+covc12-covc13-covc23;
    checks[i,15]<- s1*(1-s2)*(1-s3)*(1-s4)+covs12+covs13+covs23;
    checks[i,16]<- (1-c1)*c2*c3*c4+covc12+covc13+covc23;
    checks[i,17]<- (1-s1)*s2*s3*s4+covs12+covs13+covs23;
    checks[i,18]<- c1*(1-c2)*(1-c3)*(1-c4)+covc12+covc13+covc23;
    checks[i,19]<- (1-s1)*s2*s3*(1-s4)+covs12-covs13-covs23;
    checks[i,20]<- c1*(1-c2)*(1-c3)*c4+covc12-covc13-covc23;
    checks[i,21]<- (1-s1)*s2*(1-s3)*s4-covs12+covs13-covs23);
    checks[i,22]<- c1*(1-c2)*c3*(1-c4)-covc12+covc13-covc23;
    checks[i,23]<- (1-s1)*s2*(1-s3)*(1-s4)-covs12-covs13+covs23;
    checks[i,24]<- c1*(1-c2)*c3*c4-covc12-covc13+covc23;
    checks[i,25]<- (1-s1)*(1-s2)*s3*s4-covs12-covs13+covs23;
    checks[i,26]<- c1*c2*(1-c3)*(1-c4)-covc12-covc13+covc23;
    checks[i,27]<- (1-s1)*(1-s2)*s3*(1-s4)-covs12+covs13-covs23;
    checks[i,28]<- c1*c2*(1-c3)*c4-covc12+covc13-covc23;
    checks[i,29]<- (1-s1)*(1-s2)*(1-s3)*s4+covs12-covs13-covs23;
    checks[i,30]<- c1*c2*c3*(1-c4)+covc12-covc13-covc23;
    checks[i,31]<- (1-s1)*(1-s2)*(1-s3)*(1-s4)+covs12+covs13+covs23;
    checks[i,32]<- c1*c2*c3*c4+covc12+covc13+covc23;
    
  
    valid[i]<- step(1-q[i,1])*step(q[i,1])*
      step(1-q[i,2])*step(q[i,2])*
      step(1-q[i,3])*step(q[i,3])* 
      step(1-q[i,4])*step(q[i,4])*
      step(1-q[i,5])*step(q[i,5])*
      step(1-q[i,6])*step(q[i,6])*
      step(1-q[i,7])*step(q[i,7])*
      step(1-q[i,8])*step(q[i,8])*
      step(1-q[i,9])*step(q[i,9])*
      step(1-q[i,10])*step(q[i,10])*
      step(1-q[i,11])*step(q[i,11])*
      step(1-q[i,12])*step(q[i,12])*
      step(1-q[i,13])*step(q[i,13])*
      step(1-q[i,14])*step(q[i,14])*
      step(1-q[i,15])*step(q[i,15])*
      step(1-q[i,16])*step(q[i,16])*
      step(1-checks[i,1])*step(checks[i,1])*
      step(1-checks[i,2])*step(checks[i,2])*
      step(1-checks[i,3])*step(checks[i,3])*
      step(1-checks[i,4])*step(checks[i,4])*
      step(1-checks[i,5])*step(checks[i,5])*
      step(1-checks[i,6])*step(checks[i,6])*
      step(1-checks[i,7])*step(checks[i,7])*
      step(1-checks[i,8])*step(checks[i,8])*
      step(1-checks[i,9])*step(checks[i,9])*
      step(1-checks[i,10])*step(checks[i,10])*
      step(1-checks[i,11])*step(checks[i,11])*
      step(1-checks[i,12])*step(checks[i,12])*
      step(1-checks[i,13])*step(checks[i,13])*
      step(1-checks[i,14])*step(checks[i,14])*
      step(1-checks[i,15])*step(checks[i,15])*
      step(1-checks[i,16])*step(checks[i,16])*
      step(1-checks[i,17])*step(checks[i,17])*
      step(1-checks[i,18])*step(checks[i,18])*
      step(1-checks[i,19])*step(checks[i,19])*
      step(1-checks[i,20])*step(checks[i,20])*
      step(1-checks[i,21])*step(checks[i,21])*
      step(1-checks[i,22])*step(checks[i,22])*
      step(1-checks[i,23])*step(checks[i,23])*
      step(1-checks[i,24])*step(checks[i,24])*
      step(1-checks[i,25])*step(checks[i,25])*
      step(1-checks[i,26])*step(checks[i,26])*
      step(1-checks[i,27])*step(checks[i,27])*
      step(1-checks[i,28])*step(checks[i,28])*
      step(1-checks[i,29])*step(checks[i,29])*
      step(1-checks[i,30])*step(checks[i,30])*
      step(1-checks[i,31])*step(checks[i,31])*
      step(1-checks[i,32])*step(checks[i,32]);
    
    #######################################################
    ## Contribution to the likelihood for each observation
    #######################################################
    
    L[i]<- equals(valid[i],1)*(
      equals(bov.beis [i,1],1)*equals(bov.beis[i,2],1)*equals(bov.beis [i,3],1)*equals(bov.beis [i,4],1)*q[i,1]
      + equals(bov.beis [i,1],1)*equals(bov.beis[i,2],1)*equals(bov.beis [i,3],0)*equals(bov.beis [i,4],0)*q[i,2]
      + equals(bov.beis [i,1],1)*equals(bov.beis[i,2],1)*equals(bov.beis [i,3],0)*equals(bov.beis [i,4],1)*q[i,3]
      + equals(bov.beis [i,1],1)*equals(bov.beis[i,2],1)*equals(bov.beis [i,3],0)*equals(bov.beis [i,4],0)*q[i,4]
      + equals(bov.beis [i,1],1)*equals(bov.beis[i,2],0)*equals(bov.beis [i,3],1)*equals(bov.beis [i,4],1)*q[i,5]
      + equals(bov.beis [i,1],1)*equals(bov.beis[i,2],0)*equals(bov.beis [i,3],1)*equals(bov.beis [i,4],0)*q[i,6]
      + equals(bov.beis [i,1],1)*equals(bov.beis[i,2],0)*equals(bov.beis [i,3],0)*equals(bov.beis [i,4],1)*q[i,7]
      + equals(bov.beis [i,1],1)*equals(bov.beis[i,2],0)*equals(bov.beis [i,3],0)*equals(bov.beis [i,4],0)*q[i,8]
      + equals(bov.beis [i,1],0)*equals(bov.beis[i,2],1)*equals(bov.beis [i,3],1)*equals(bov.beis [i,4],1)*q[i,9]
      + equals(bov.beis [i,1],0)*equals(bov.beis[i,2],1)*equals(bov.beis [i,3],1)*equals(bov.beis [i,4],0)*q[i,10]
      + equals(bov.beis [i,1],0)*equals(bov.beis[i,2],1)*equals(bov.beis [i,3],0)*equals(bov.beis [i,4],1)*q[i,11]
      + equals(bov.beis [i,1],0)*equals(bov.beis[i,2],1)*equals(bov.beis [i,3],0)*equals(bov.beis [i,4],0)*q[i,12]
      + equals(bov.beis [i,1],0)*equals(bov.beis[i,2],0)*equals(bov.beis [i,3],1)*equals(bov.beis [i,4],1)*q[i,13]
      + equals(bov.beis [i,1],0)*equals(bov.beis[i,2],0)*equals(bov.beis [i,3],1)*equals(bov.beis [i,4],0)*q[i,14]
      + equals(bov.beis [i,1],0)*equals(bov.beis[i,2],0)*equals(bov.beis [i,3],0)*equals(bov.beis [i,4],1)*q[i,15]
      + equals(bov.beis [i,1],0)*equals(bov.beis[i,2],0)*equals(bov.beis [i,3],0)*equals(bov.beis [i,4],0)*q[i,16]
    ) +(1-equals(valid[i],1)) *(1e-14);
    
    
    #######################################################
    ## Trick to ensure the probabilities are always less than 1
    #######################################################
    
    p[i] <- L[i] / 1;## divided by a constant just to ensure all p's <1
    ones[i] ~ dbern(p[i]);     
    
  }
  
  #######################################################
  ## Definition of model priors which should be modified for sensitivity analysis.
  ## Assuming a uniform distribution for the covariance terms
  ## Assuming a beta distribution for prevalence and test accuracy estimates.
  ## Parameters of the beta prior distribution obtained by Betabuster
  ## https://betabuster.software.informer.com/1.0/
  #######################################################
  
  covs12 ~ dbeta(1,1);
  covs13 ~ dbeta(1,1);
  covs23  ~ dbeta(1,1);
  covc12 ~ dbeta(1,1);
  covc13 ~ dbeta(1,1);
  covc23 ~ dbeta(1,1);
  
  prc ~ dbeta(1.532, 1.532);  ##Mode 0.50, Min 0.10 (95% confidence), Quantiles: 5%, 95%	(0.1, 0.9) 
c1  ~ dbeta(1,1);                                 
c2  ~ dbeta(1,1);   
c3  ~ dbeta(1,1);  
c4  ~ dbeta(42.111, 1.839);   ## Sp Mode 0.98, Min 0.90 (95% confidence),   Quantiles: 5%, 95%	(0.9, 0.993) 
s1  ~ dbeta(1,1);               
s2  ~ dbeta(1,1);        
s3  ~ dbeta(1,1);  
s4  ~ dbeta(1,1);  
  logL<-sum(log(p[1:N])); 
  
}") 

# Set different starting values for the three different chains to assess convergence

inits1 = list(".RNG.name" ="base::Mersenne-Twister",
              ".RNG.seed" = 100022)
inits2 = list(".RNG.name" ="base::Mersenne-Twister",
              ".RNG.seed" = 300022)
inits3 = list(".RNG.name" ="base::Mersenne-Twister",
              ".RNG.seed" = 500022)


# Run the model with the package runjags 
results <- run.jags(hw_4t_1p,
                    data=list(N=N, m.ca=m.ca, ones=ones),
                    inits=list(inits1, inits2, inits3), 
                    n.chains = 3,
                    burnin = 10000,
                    sample = 100000,
                    adapt = 1000,
                    monitor = c('prc','c1','c2','c3','c4','s1','s2','s3','s4', 'covs12','covs13','covs23','covc12','covc13','covc23', 'dic', 'deviance'))

# Model checking and diagnostics:

results 

plot(results,
     vars = list("prc", "c1","c2", "c3"),
     layout = c(5,4),
     plot.type = c("trace", "histogram", "autocorr", "ecdf"))


plot(results,
     vars = list("s1", "s2", "s3"),
     layout = c(5,4),
     plot.type = c("trace", "histogram", "autocorr", "ecdf"))

plot(results, plot.type='density', vars=c('c1', 'c2','c3','s1','s2','s3'))

round(results$summary$quantiles,3)*100

#Extract DIC
extract.runjags(results, what = 'dic') #To obtain mean deviance and penalized deviance