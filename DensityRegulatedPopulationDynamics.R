
## R and OpenBUGS code for "Density regulation amplifies environmentally induced population fluctuations"


##  R code for simulating data from the stochastic Gompertz model

dataGompertz<-function(n, r, k, sdr, init=k){ 
  # n is the sample size, sdr the variance of the noise
  # init is the initial population size, set to 100
  y<-rep(0,n); epsillon<-rep(0,n)
  gr<-rep(0,n)
  y[1]<-init
  for(t in 2: n){
    epsillon[t]<-rnorm(1, 0, sdr)
    gr[t]<-exp(r*(1-log(y[t-1])/k)+ epsillon[t])
    y[t]<-y[t-1]*gr[t]
  }
  y
}

############################################################################################################################################
# We will fit the stochastic Gompertz and stochastic Ricker models to the simulated data replicates with a Bayesian approach using MCMC 
# throug OpenBUGS. We will use the package BRugs, provides an interface netween R and OpenBUGS. We start by defining two functions for the 
# Gomprtz model and the Ricker model.
############################################################################################################################################

gompertzModel=function() {
  for(t in 2: n){
    y[t]~dnorm(m[t],tau.y)
    m[t]<- r + beta*y[t-1]
    ypred[t]~dnorm(m[t],tau.y)
    #ypred is drawn from the PPD at time t 
    squerr[t]<-pow((ypred[t]-y[t]),2)
  }
  k<-r/(1-beta)
  r~dnorm(0,1)%_%I(0,)
  beta~dnorm(0,1)
  tau.y~dgamma(1,1)
  sigma2.y<-1/tau.y
  rmse<-mean(squerr[2:n])
  # rmse is the root mean squared error
}


rickerModel=function() {
  for(t in 2: n){
    m[t]<- y[t-1] + r*(1-exp(y[t-1])/K)
    y[t]~dnorm(m[t],tau.y)
    ypred[t]~dnorm(m[t],tau.y)
    #ypred is drawn from the PPD at time t 
    squerr[t]<-pow(ypred[t]-y[t],2)
  }
  K~dgamma(1,1)
  r~dnorm(0,1)%_%I(0,)
  tau.y~dgamma(1,1)
  sigma2.y<-1/tau.y
  rmse<-mean(squerr[2:n])
  # rmse is the root mean squared error
} 


############################################### DATA SIMULATION AND ANALYSIS ################################################################

nsim=50 # number of replications for each combination of parameters

# We'll simulate 300 observations starting from k and drop the first 200 samples to ensure that the last n=100 observations come from the
# stationary distribution


## Root Mean Squred Errors under Gompertz (RMSE1) and Ricker (RMSE2)

RMSE1=rep(0, nsim) 

RMSE2<-rep(0, nsim) 


## Deviance Information Criteria under Gompertz (DIC1) and Ricker (DIC2)

DIC1=rep(0, nsim) 

DIC2=rep(0, nsim) 


# Stationary variance (vs), estimated environmental variance from Gompertz(sigma2y1) and from Ricker(sigma2y2)

vs = NULL

sigma2y1 = NULL

sigma2y2 = NULL


#Fitting the Gompertz model (model1) and the Ricker model (model2) to the simulated data

# Simulation setup

m=300; r=1; k=4; sd=1

beta=(1-r/k) # corresponding value of the AR(1) parameter

## k values to be considered: k=1.334; k=2; k=4



# Using the R library BRugs for fitting OpenBUGS from within R

library(BRugs)

for(i in 1:nsim){ 
  
  #simulate a dataset
  
  simData<- dataGompertz(m, 1, k, 1)
  
  ## vs<-var(log(simData)) 
  
  vs[i]<-var(log(simData)) 
  
  # Formatting the data for BUGS
  
  DataToBUGS=list(y=log(simData[201:300]), n=100)
  
  writeModel(gompertzModel, "model1.bug")
  
  writeModel(rickerModel, "model2.bug")
  
  bugsData(DataToBUGS, "Data1.bug")
  
  
  thing1=BRugsFit("model1.bug", "Data1.bug", numChains = 1, parametersToSave=c("r", "beta", "tau.y", "ypred"), nBurnin = 4000, nIter = 6000, nThin = 10,  DIC = TRUE, working.directory = NULL, digits = 3)
  
  thing2=BRugsFit("model2.bug", "Data1.bug", numChains = 1, parametersToSave=c("r",  "K", "tau.y", "ypred"), nBurnin = 4000, nIter = 6000, nThin = 4,    DIC = TRUE, working.directory = NULL, digits = 3)
  
  ypred1=thing1$Stats[4:102,1]
  
  ypred2=thing2$Stats[4:102,1]
  
  sigma2y1[i]=1/thing1$Stats[3,1]
  
  sigma2y2[i]=1/thing2$Stats[3,1]
  
  
  ## Root Mean Squared Errors under Gompertz model (RMSE1) and under the Ricker model (RMSE2)
  
  RMSE1[i]<-sqrt(mean( log(simData[2:n])-ypred1)^2)
  
  RMSE2[i]<-sqrt(mean( log(simData[2:n])-ypred2)^2)
  
  ## Deviance Information Criteria under Gompertz model (DIC1) and the Ricker model (DIC2) 
  
  DIC1[i]=thing1$DIC[3,3]
  
  DIC2[i]=thing2$DIC[3,3]
  
  ## Saving the results to the file res_sim2.txt
  
  sink("C:\\Users\\cmutshinda\\Desktop\\res_sim2.txt", append = TRUE)
  
  cat(vs[i], sigma2y1[i], sigma2y2[i], RMSE1[i], RMSE2[i], DIC1[i], DIC2[i], "\n", sep=" ")
  
  sink()
  
}