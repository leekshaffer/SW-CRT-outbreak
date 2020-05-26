##########################################
##### SWT_Analysis_PH.R ##################
##### Perform Cox proportional hazards estimation and asymptotic inference on SWT results
##### Kennedy-Shaffer and Lipsitch #######
##### Update: May 26, 2020 ###############
##### Please Cite https://doi.org/10.1101/2020.05.01.20087429
##### Based on methods described by Durovni et al. (https://doi.org/10.1016/S1473-3099(13)70187-7)
##### and Bellan et al. (https://doi.org/10.1016/S1473-3099(15)70139-8)
##### For questions, contact Lee Kennedy-Shaffer: lee_kennedyshaffer@g.harvard.edu
##### Approximate Runtime for SWT-A from paper: 3 minutes for 500 permutations for one simulation
##########################################

NumPerms <- 500  ## The number of permutations to be used for hypothesis testing

### Packages Needed: ###
# install.packages("survival")
require(survival)

# Parameter values needed from Outbreak_Simulations: ave_inc_period, trial_startday, trial_length
trial_lastday <- trial_startday + trial_length

# Function to Analyze Results
SurvAnalysis <- function(trial_nodes, crossovers, trial_inf_nodes) {
  df <- merge(trial_nodes, crossovers, by=c("Community")) #Get crossover date for each community
  Infected_Nodes <- trial_inf_nodes
  Infected_Nodes$Node <- Infected_Nodes$InfectedNode
  df <- merge(df, Infected_Nodes[,c("Community","Node","DayInfected")], by=c("Community","Node"), all.x=TRUE)
  df_neverinf <- df[is.na(df$DayInfected),]
  df_vaccinf <- df[!is.na(df$DayInfected) & df$DayInfected >= df$DayVaccinated,]
  df_continf <- df[!is.na(df$DayInfected) & df$DayInfected < df$DayVaccinated,]
  
  #Create Survival-type datasets:
  df_neverinf_format1 <- cbind(df_neverinf, event=rep(0,dim(df_neverinf)[1]),
                               start=df_neverinf$DayEnrolled + ave_inc_period - df_neverinf$DayEnrolled,
                               stop=df_neverinf$DayVaccinated-df_neverinf$DayEnrolled,
                               interv=rep(0,dim(df_neverinf)[1]))
  df_neverinf_format2 <- cbind(df_neverinf, event=rep(0,dim(df_neverinf)[1]),
                               start=df_neverinf$DayVaccinated + ave_inc_period - df_neverinf$DayEnrolled,
                               stop=rep(trial_lastday,dim(df_neverinf)[1])-df_neverinf$DayEnrolled,
                               interv=rep(1,dim(df_neverinf)[1]))
  df_vaccinf_format1 <- cbind(df_vaccinf, event=rep(0,dim(df_vaccinf)[1]),
                              start=df_vaccinf$DayEnrolled + ave_inc_period - df_vaccinf$DayEnrolled,
                              stop=df_vaccinf$DayVaccinated-df_vaccinf$DayEnrolled,
                              interv=rep(0,dim(df_vaccinf)[1]))
  df_vaccinf_format2 <- cbind(df_vaccinf, event=rep(1,dim(df_vaccinf)[1]),
                              start=df_vaccinf$DayVaccinated + ave_inc_period - df_vaccinf$DayEnrolled,
                              stop=df_vaccinf$DayInfected-df_vaccinf$DayEnrolled,
                              interv=rep(1,dim(df_vaccinf)[1]))
  df_continf_format <- cbind(df_continf, event=rep(1,dim(df_continf)[1]),
                             start=df_continf$DayEnrolled + ave_inc_period - df_continf$DayEnrolled,
                             stop=df_continf$DayInfected-df_continf$DayEnrolled,
                             interv=rep(0,dim(df_continf)[1]))
  
  SurvDF <- rbind(df_neverinf_format1,df_neverinf_format2,df_vaccinf_format1,df_vaccinf_format2,df_continf_format)
  SurvDF <- SurvDF[order(SurvDF$Node),]
  # Remove events within average incubation period of enrollment or crossover:
  SurvDF2 <- SurvDF[SurvDF$stop > SurvDF$start,]
  numev_v <- sum(SurvDF2$event*SurvDF2$interv)
  numev_c <- sum(SurvDF2$event*(1-SurvDF2$interv))
  
  if (numev_v > 0 & numev_c > 0) {
    # Run Cox model with Gamma-distributed shared frailty for community and time-varying intervention indicator:
    res2 <- try(coxph(Surv(start,stop,event)~interv+frailty(Community, distribution="gamma", sparse=FALSE),
                      data=SurvDF2), silent=TRUE)
    useres2 <- !inherits(res2, 'try-error')
    FullOut <- list(SurvDF=SurvDF,SurvDF2=SurvDF2,res2=res2)
    
    if (useres2) {
      Est <- unname(res2$coefficients['interv'])
      PVal <- unname(summary(res2)$coefficients['interv','p'])
      ErrType <- 0
    } else {
      ErrType <- 1
      Est <- NA
      PVal <- NA
    }
  } else if (numev_v > 0 & numev_c==0) {
    Est <- log(2) # Gives VE = -1 if there are events in vaccine arm but not in control arm
    PVal <- 1
    FullOut <- list(SurvDF=SurvDF,SurvDF2=SurvDF2,res2=c(numev_v=numev_v,numev_c=numev_c))
    ErrType <- 2
  } else if (numev_v==0 & numev_c > 0) {
    Est <- -100 # Gives VE = 1 if there are events in control arm but not in vaccine arm
    # To get p-value, we take the median event time in the control arm. We add one event
    # at that time in a cluster where that time is on control and one where it is on intervention.
    # We take an intervention cluster that will have the event at the shortest time after intervention
    # and the most conservative control cluster, to get a conservative p-value.
    Medtime <- median(SurvDF2$stop[SurvDF2$event==1 & SurvDF2$interv==0], na.rm=TRUE)
    DayEnrolled <- median(SurvDF2$DayEnrolled, na.rm=TRUE)
    IntClusts <- crossovers$Community[crossovers$DayVaccinated + ave_inc_period < Medtime + DayEnrolled]
    ContClusts <- crossovers$Community[crossovers$DayVaccinated > Medtime + DayEnrolled]
    ErrType <- 3
    if (length(IntClusts)==0 | length(ContClusts==0)) {
      Est <- NA
      PVal <- NA
    } else {
      pval_set <- NULL
      ClustIvaccday <- min(crossovers$DayVaccinated[crossovers$Community %in% IntClusts])
      ClustI <- min(crossovers$Community[crossovers$DayVaccinated==ClustIvaccday])
      for (ClustC in ContClusts) {
        ClustIvaccday <- crossovers$DayVaccinated[crossovers$Community==ClustI]
        ClustCvaccday <- crossovers$DayVaccinated[crossovers$Community==ClustC]
        df_adj <- rbind(SurvDF2,
                        data.frame(Community=c(ClustI,ClustI,ClustC),Node=c(100000,100000,100001),
                                   SimulationNumber=c(NA,NA,NA), TrialStatus=c(1,1,0),
                                   DayEnrolled=rep(DayEnrolled,3), DayVaccinated=c(rep(ClustIvaccday,2),ClustCvaccday),
                                   DayInfected=rep(Medtime+DayEnrolled,3),
                                   event=c(0,1,1), start=c(0,ClustIvaccday-DayEnrolled,0),
                                   stop=c(ClustIvaccday-DayEnrolled,Medtime,Medtime),interv=c(0,1,0)))
        res_adj <- try(coxph(Surv(start,stop,event)~interv+frailty(Community, distribution="gamma", sparse=FALSE),
                             data=df_adj))
        useres_adj <- !inherits(res_adj, 'try-error')
        if (useres_adj) {
          pval_set <- c(pval_set,summary(res_adj)$coefficients['interv','p'])
        }
      }
      if (is.null(pval_set)) {
        PVal <- NA
      } else  {
        PVal <- max(pval_set, na.rm=TRUE)
      }
    }
    FullOut <- list(SurvDF=SurvDF,SurvDF2=SurvDF2,res2=c(numev_v=numev_v,numev_c=numev_c))
  } else { #No events in either arm --> no information from trial
    ErrType <- 4
    Est <- NA
    PVal <- NA
    FullOut <- list(SurvDF=SurvDF,SurvDF2=SurvDF2,res2=c(numev_v=numev_v,numev_c=numev_c))
  }
  FullOut <- append(FullOut, list(Est=Est, PVal=PVal, ErrType=ErrType))
  return(FullOut)
}

PermSurvAn <- function(trial_nodes, crossovers, trial_inf_nodes) {
  PermCrossovers <- crossovers
  PermCrossovers$Community <- sample(crossovers$Community, size=dim(crossovers)[1], replace=FALSE)
  Analysis <- SurvAnalysis(trial_nodes, PermCrossovers, trial_inf_nodes)
  return(Analysis$Est)
}

st <- proc.time()
# Estimation and Asymptotic Inference:
Truth <- SurvAnalysis(fullRes_SWT$trial_nodes, fullRes_SWT$crossovers, fullRes_SWT$trial_inf_nodes)
TrueEst <- Truth$Est
AsyPVal <- Truth$PVal
ErrType <- Truth$ErrType
# Permutation Inference:
if (!is.na(TrueEst) & NumPerms > 0) {
  PermEsts <- replicate(n=NumPerms,
                        expr=PermSurvAn(fullRes_SWT$trial_nodes, fullRes_SWT$crossovers,
                                        fullRes_SWT$trial_inf_nodes),
                        simplify="vector")
  PermPVal <- mean(ifelse(abs(PermEsts) >= abs(TrueEst), 1, 0), na.rm=TRUE)
} else {
  PermPVal <- NA
}


Analysis_PH_Res <- list(Est=c(PH=TrueEst), VE=c(PH=1-exp(TrueEst)), PVal=c(PH.Asy=AsyPVal, PH.Perm=PermPVal), 
                        ErrType=ErrType)

proc.time() - st
