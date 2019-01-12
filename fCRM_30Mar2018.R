## -------------------------------------------------------------------------- ##
## -                    Implement the fractional CRM design                 - ##
## -                                                                        - ##
## - Created by YANG Zhao and Guosheng Yin, 22 Oct, 2017                    - ##
## - Modified by YANG Zhao and Guosheng Yin, 27 Dec, 2017                   - ##
## - Maintained by Guosheng Yin, gyin@hku.hk                                - ##
## -------------------------------------------------------------------------- ##

nSimulation <- 10   ## No. of simulated studies
target <- 0.3        ## The target toxicity probability
prob <- 0.1          ## The DLTs probability of DLTs at the first half of tau
tau <- 3             ## Assessment window
a <- 1               ## Accrual window
ndose <- 6           ## No. of dose levels
set.seed(6)

## Fixed scenarios
p.true <- matrix(NA, ncol = 6, nrow = 13)
p.true[1,] <- c(0.30, 0.38, 0.45, 0.60, 0.68, 0.75)
p.true[2,] <- c(0.17, 0.30, 0.43, 0.55, 0.65, 0.80)
p.true[3,] <- c(0.05, 0.13, 0.30, 0.38, 0.65, 0.85)
p.true[4,] <- c(0.10, 0.16, 0.30, 0.42, 0.60, 0.70)
p.true[5,] <- c(0.01, 0.05, 0.17, 0.30, 0.40, 0.45)
p.true[6,] <- c(0.08, 0.12, 0.15, 0.30, 0.45, 0.65)
p.true[7,] <- c(0.06, 0.08, 0.10, 0.20, 0.30, 0.50)
p.true[8,] <- c(0.01, 0.04, 0.08, 0.10, 0.15, 0.30)
p.true[9,] <- c(0.08, 0.10, 0.12, 0.15, 0.31, 0.55)
p.true[10,] <- c(0.05, 0.10, 0.15, 0.28, 0.40, 0.58)
p.true[11,]<- c(0.15, 0.29, 0.36, 0.43, 0.52, 0.59)
p.true[12,]<- c(0.06, 0.12, 0.18, 0.27, 0.37, 0.45)
p.true[13,]<- c(0.45, 0.50, 0.58, 0.65, 0.73, 0.85)


## Load the required packages 
library(tidyverse)
library(foreach)
library(iterators)
library(parallel)
library(doParallel)
library(dfcrm)
library(survival)
library(purrr)
library(magrittr)
library(listviewer)
library(Iso)
library(rstan)
library(R2jags)

## Prepare the simulated stidies 
cl <- makeCluster(detectCores())
registerDoParallel(cl)

## Summary
SummarySim <- function(x, ndose = 6, Sec = nrow(p.true), nSim = nSimulation ) {
  Result <- NULL
  for (sec in 1:Sec) {
    Res <- matrix(0, ncol = 10 + ndose, nrow = 4)
    colnames(Res) <- c(paste0("Dose", 1:ndose, sep = ""),
                       "# overdose",
                       "Risk of high toxicity",
                       "# DLTs",
                       "Trial duration",
                       "PCT for MTD Selection",
                       "% patients treated at MTD",
                       "% of trials select overdose",
                       "% of patients allocated to overdose",
                       "% of patients experiencing DLTs",
                       "# treated patients")
    row.names(Res) <- c(paste("Scenario ", sec, sep = ""),
                        "Selection %", "# patients", "# DLTs")
    doseselect <- ntrts <- ndlts <- noverdose <- 0
    nstop <- ntotaldlt <- trialperiod <- peroverdose <- 0
    pernoverdose <- perdlts <- 0
    for (nsim in 1:nSim) {
      OC <- x[[sec]][[nsim]]
      dMTD <- which.min(abs(0.3 - OC$p.true))
      ## Scenario
      Res[1, 1:ndose] <- OC$p.true
      doseselect <- doseselect + OC$DoseSelect/nSim * 100
      ntrts <- ntrts + OC$nTRTs/nSim
      ndlts <- ndlts + OC$nDLTs/nSim
      noverdose <- noverdose + OC$nOverDose/nSim
      nstop <- nstop + OC$nStop/nSim * 100
      ntotaldlt <- ntotaldlt + OC$nTotalDLT/nSim
      trialperiod <- trialperiod + OC$TrialPeriod/nSim
      peroverdose <- peroverdose + (OC$MTD > dMTD)/nSim * 100
      pernoverdose<- pernoverdose + (OC$nOverDose/36)/nSim
      perdlts <- perdlts + (sum(OC$nDLTs)/36)/nSim
    }
    Res[2, 1:ndose] <- round(doseselect,2)    ## 1-6 columns
    Res[3, 1:ndose] <- round(ntrts,2)         ## 1-6 columns
    Res[4, 1:ndose] <- round(ndlts,2)         ## 1-6 columns
    Res[2, 1 + ndose] <- round(noverdose,2)   ## 7th column
    Res[2, 2 + ndose] <- round(nstop,2)       ## 8th column
    Res[2, 3 + ndose] <- round(ntotaldlt,2)   ## 9th column
    Res[2, 4 + ndose] <- round(trialperiod,2) ## 10th column
    Res[2, 5 + ndose] <- round(doseselect[dMTD],2) ## 11th column
    Res[2, 6 + ndose] <- round(ntrts[dMTD]/36 * 100, 2) ## 12 column
    Res[2, 7 + ndose] <- round(peroverdose, 2)          ## 13 column
    Res[2, 8 + ndose] <- round(pernoverdose, 2)         ## 14 column
    Res[2, 9 + ndose] <- round(perdlts, 2)              ## 15 column
    Res[2, 10+ ndose] <- round(sum(ntrts), 2)
    Result <- rbind(Result, Res)
  }
  return(Result)
}

## Fractional contribution
get.FractionalContribution <- function(time,         ## Time for DLT
                                       DLT,          ## Whether DLT or not
                                       current.time, ## Current Time
                                       assess.time,  ## Assessment period
                                       tau){         ## Follow-up period
  ## Loading the survival package for implmenting the Kaplan-Meier Method
  ## require(survival)
  km   <- survfit(Surv(time,DLT) ~ 1)  ## K-M Method for Estimating Survival Probability
  stau <- km$surv[which.min(km$surv[km$time <= tau])] ## Calc the Survival Probability = S(tau)
  U    <- current.time - assess.time + tau  ## Calc the Observed Time Period
  fDLT <- DLT                                         ##
  for(k in 1:length(DLT)){
    if(current.time < assess.time[k] && DLT[k] == 0){ ## DLT == 0 && Within the Assessment Window
      sTime   <- km$surv[which.min(km$surv[km$time <= U[k]])]
      fDLT[k] <- (sTime-stau)/sTime                   ## Calc the Fractional Contribution
    }
  }
  return(fDLT)
}


## fCRM
fCRM <- function(
  nSim, nScenario,
  p,                 ## Prior Toxicity Probability
  p.true,            ## Ture Toxicity Probability
  NCohort    = 12,   ## No. of Cohorts
  CohortSize = 3,    ## No. of Subs at each cohort
  tau,               ## Assess Window
  a,                 ## Accrual Window 
  target = 0.30,     ## Target Toxicity Probability
  p.eli  = 0.95,     ## Overtoxic Probability for Trial Termination
  p.offset= 0.00,    ## Small positive number restrict dose-finding rules
  prob   = 0.30,     ## Probability of DLT occurs in the 1st half of tau
  method  = "Weibull",# Distribution of Time-To-Event 
  Sigma2 = 2) {
  cat("The fractional CRM Scenario: ",nScenario, "Simulation", nSim, "\n")
  
  ## set.seed(Seed)
  ## ndose <- length(p)
  ## Likelihood Function
  PosteriorH <- function(alpha,DoseLevel,DLT,p,Sigma2){
    likelihood <- 1  
    for (i in 1:length(DLT)){
      phat <- p[DoseLevel[i]]^exp(alpha)
      likelihood <- likelihood*(phat^DLT[i])*(1-phat)^(1 - DLT[i])
    }
    return(likelihood*exp(-0.5*alpha*alpha/Sigma2))
  }
  
  ## posterior mean of toxicity probability
  PosteriorPi <- function(alpha, DoseLevel,DLT,p,j,Sigma2){
    p[j]^(exp(alpha))*PosteriorH(alpha,DoseLevel,DLT,p,Sigma2)
  }
  
  ## Weibull distribution for Time-To-Event
  #if (method == "Weibull"){ 
  #  alpha = log(log(1-p.true)/log(1-0.3*p.true))/log(2)
  #  beta  = tau/(-log(1-p.true))^(1/alpha)
  #}
  if (method == "Weibull") {
    alpha.Weibull <- log(log(1 - p.true)/log(1 - prob*p.true))/log(2)
    beta.Weibull  <- tau/(-log(1 - p.true))**(1/alpha.Weibull)
  } else if (method == "Exponential") { 
    lambda <- -log(1 - p.true)/tau
  } else if (method == "Gamma") {
    f <- function(aa, p = prob, q = 1 - prob) { 
      qgamma(p, aa, 1)/qgamma(p - p * q, aa, 1) - 2
    }
    alpha.Gamma <- rep(0, ndose)
    beta.Gamma  <- rep(0, ndose)
    for (k in 1:ndose) {
      alpha.Gamma[k] <- uniroot(f, c(0.01, 10), p = p.true[k], q = 1 - prob)$root
      beta.Gamma[k]  <- qgamma(p.true[k], alpha.Gamma[k], 1)/tau
    }
  } else if (method == "Log-Normal") {
    sigma <- log(2)/(qnorm(p.true) - qnorm(prob * p.true))
    mu    <- (log(1.5) * qnorm(p.true) - log(tau) * qnorm(prob * p.true))/(qnorm(p.true) - qnorm(prob * p.true))
  }
  
  DoseSelect <- rep(0, length(p))
  nDLTs      <- rep(0, length(p))
  nTRTs      <- rep(0, length(p))
  #TrialPeriod<- 0
  NStop <- DoseBest <- 0
  #pb <- txtProgressBar(min = 1, max = NTrial, style = 3)
  ## Simulations 
  #for (Trial in 1:NTrial){
  #Sys.sleep(0.01)
  #setTxtProgressBar(pb, Trial)
  ##cat("Simulations:",Trial,"\n")
  DoseCurr <- 1;
  DLT      <- NULL
  TimeY    <- NULL
  DoseLevel<- NULL
  pihat    <- rep(0, length(p))
  Stop     <- 0
  CurrentTime <- 1
  AssessTime  <- NULL
  for (i in 1:NCohort){
    ## 1st Cohort
    if (method == "Weibull") { 
      TimeY <- c(TimeY, rweibull(CohortSize, 
                                 alpha.Weibull[DoseCurr],
                                 beta.Weibull[DoseCurr]))
    } else if (method == "Exponential") { 
      TimeY <- c(TimeY, rexp(CohortSize, 
                             lambda[DoseCurr]))
    } else if (method == "Gamma")       { 
      TimeY <- c(TimeY, rgamma(CohortSize, 
                               alpha.Gamma[DoseCurr], 
                               beta.Gamma[DoseCurr]))      
    } else if (method == "Lognormal")   { 
      TimeY <- c(TimeY, rlnorm(CohortSize, 
                               sigma[DoseCurr], 
                               mu[DoseCurr]))
    }
    # 
    # TimeY       <- c(TimeY,Y)                                      ## Time-To-event
    AssessTime  <- c(AssessTime,rep(CurrentTime + tau,CohortSize))   ## Assessment Window
    DoseLevel   <- c(DoseLevel,rep(DoseCurr,CohortSize))             ## Current Dose Level
    
    ## Next Cohort
    CurrentTime <- CurrentTime + a      ## EntryDate Time for  Next Cohort
    DLT <- (TimeY < (CurrentTime - AssessTime + tau))*(CurrentTime < AssessTime) +
      (TimeY <= tau)*(CurrentTime >= AssessTime);

    Time <- TimeY*DLT + (1-DLT)*(CurrentTime - AssessTime + tau)*(CurrentTime < AssessTime) +
      (1-DLT)*TimeY*(CurrentTime >= AssessTime)    
    
    while (sum(DLT) == 0) {
      CurrentTime <- CurrentTime + a
      DLT <- (TimeY <= tau) * (AssessTime <= CurrentTime) + 
        (TimeY <= (CurrentTime - AssessTime + tau)) * 
        (AssessTime > CurrentTime)
      if (CurrentTime == max(AssessTime)) {
        break
      }
    }
    

    
    if(sum(DLT) == 0 && DoseCurr <  length(p)){ 
      DoseCurr <- DoseCurr + 1 
    }
    if(sum(DLT) == 0 && DoseCurr == length(p)){ 
      DoseCurr <- DoseCurr    
    }
    if(sum(DLT) != 0){
      ## Fractional Controbution
      km   <- survfit(Surv(Time, DLT)~1)
      stau <- km$surv[which.min(km$surv[km$time <= tau])]
      U    <- CurrentTime - AssessTime + tau
      fDLT <- DLT
      for(i in 1:length(DLT)){
        if(CurrentTime < AssessTime[i] && DLT[i] == 0){
          sTime   <- km$surv[which.min(km$surv[km$time <= U[i]])]
          fDLT[i] <- (sTime-stau)/sTime 
        }
      }
      ## Dose-Finding Algorithm
      marginal <- integrate(PosteriorH,lower = -Inf, upper = Inf, 
                            DoseLevel, fDLT, p, Sigma2)$value
      ## Early Termination
      POverTox1 <- integrate(PosteriorH, lower = -Inf, upper = log(log(target)/log(p[1])),
                             DoseLevel, fDLT, p, Sigma2)$value/marginal
      if (POverTox1 > (p.eli - p.offset)){ 
        Stop = 1; 
        break; 
      }
      ## Posterior mean of toxicity probabiity
      for (j in 1:length(p)){ 
        pihat[j] <- integrate(PosteriorPi, lower=-Inf, upper=Inf,
                              DoseLevel, fDLT, p, j, Sigma2)$value/marginal 
      }
      DoseBest <- which.min(abs(pihat-target))
      if (DoseBest > DoseCurr && DoseCurr < length(p)) { DoseCurr <- DoseCurr + 1 }
      if (DoseBest < DoseCurr && DoseCurr > 1        ) { DoseCurr <- DoseCurr - 1 }
      DoseCurr <- DoseCurr 
    }
  } ## NCohort
  
  #  if (Stop == 0){
  ## Each Subs Should be Fully Followed = AssessTime Window (tau)
  DLT    <- (TimeY <= tau) 
  ## Dose-Finding Algorithm
  marginal <- integrate(PosteriorH, lower = -Inf, upper = Inf, 
                        DoseLevel, DLT, p, Sigma2)$value
  ## Early Termination
  POverTox1 <- integrate(PosteriorH, lower=-Inf, upper=log(log(target)/log(p[1])), 
                         DoseLevel, DLT, p, Sigma2)$value/marginal
  if(POverTox1 > (p.eli - p.offset)) {
    Stop = 1;
  }
  for (j in 1:length(p)){ pihat[j] <- integrate(PosteriorPi, lower=-Inf, upper=Inf, 
                                                DoseLevel, DLT, p, j,
                                                Sigma2)$value/marginal }
  if (Stop == 0){
    DoseBest <- which.min(abs(pihat-target))
    DoseSelect[DoseBest] = 1
  }
  if ((sum(DLT)/sum(length(DLT))) > target) { NStop = NStop + 1 }
  #  }
  ## DLT should be determined after the full Assessment Window
  for (i in 1:length(p)) { nDLTs[i] = nDLTs[i] + sum(DLT[DoseLevel == i]) }    
  for (i in 1:length(p)) { nTRTs[i] = nTRTs[i] + sum(DoseLevel == i)      }
  TrialPeriod = AssessTime[which.max(AssessTime)]
  #} ## NTrial
  
  #DoseSelect <- DoseSelect/NTrial*100
  #nStop      <- NStop/NTrial*100
  nTotalDLT  <- sum(nDLTs)
  DoseSet    <- which.min(abs(p.true - target))
  nOverDose  <- sum(nTRTs[DoseSet:length(p)])-nTRTs[DoseSet]
  
  res <- list(prior.p = p,
              p.true  = p.true,
              DoseSelect = DoseSelect,   ## % of Selection
              nTRTs = nTRTs,        ## No. of TRTs at Each Dose Level
              nDLTs = nDLTs,        ## No. of DLTs at Each Dose Level
              nStop = NStop,        ## Risk of High Toxicity
              nTotalDLT = nTotalDLT,    ## Total No. of DLTs
              nOverDose = nOverDose,    ## Total No. of Subs Assigned to Over Dose Level
              TrialPeriod = TrialPeriod,
              MTD = DoseBest)  ## Average Trial Period
  return(res)
}

## fCRM simulation
time <- Sys.time()
fCRM.OC <- foreach(sec = 1:nrow(p.true)) %:%
  foreach(nsim = 1:nSimulation) %do% {
    nu <- which.min(abs(p.true[sec, ] - target))
    p <- getprior(halfwidth = 0.05, target = target, nu = nu, nlevel = ndose)
    fCRM(nSim = nsim, nScenario = sec,
         p = p, p.true = p.true[sec,],
         NCohort = 12, CohortSize = 3,
         tau = tau, a = a,
         target = target, p.eli = 0.85, prob = prob)
  }
fCRMTime <- Sys.time() - time
fCRMOC <- SummarySim(x = fCRM.OC)



