######################################################
### Code to Analyze SW-CRTs in Outbreak Setting    ###
### Using Methods from                             ###
### Kennedy-Shaffer et al., Stat Med, 2020         ###
### https://doi.org/10.1002/sim.8451               ###
### Update: 3/24/2020                              ###
### Contact: lee_kennedyshaffer@g.harvard.edu      ###
### See Latest Update at:                          ###
### https://github.com/leekshaffer/SW-CRT-analysis ###
######################################################


######################################################
###  check/install/load needed packages            ###
######################################################

if (!require(MASS)) {
  package.install("MASS")
  library(MASS)
}
if (!require(Synth)) {
  package.install("Synth")
  library(Synth)
}
if (!require(lme4)) {
  package.install("lme4")
  library(lme4)
}


######################################################
###  SC-SWT Functions                              ###
######################################################

### Helper Functions ###

expit <- function(x) exp(x)/(1+exp(x))

## Functions Defining Contrasts of Interests and Inverting those Contrasts ##
RDFunc <- function(a,b) a-b
RDInvApply <- function(out,effect) out-effect
logRRFunc <- function(a,b) log(a/b)
logRRInvApply <- function(out,effect) out*exp(-1*effect)
logORFunc <- function(a,b) log(a*(1-b)/((1-a)*b))
logORInvApply <- function(out,effect) expit(-1*effect+log(out/(1-out)))

## Create a vector with each cluster's first intervention period, by observation ##
StartTimeVec <- function(Periods, Clusters, Trts) {
  StartTimes <- tapply(Periods[Trts==1], Clusters[Trts==1], FUN=min)
  NumPds <- tapply(Periods, Clusters, FUN=length)
  return(rep(as.vector(StartTimes), times=NumPds))
}


### Methods for Estimation from SW-CRT ###

## NPWP Method from Thompson et al. 2018##
# Inputs:
## Periods: A vector of period numbers for each cluster
## Outcomes: A vector of cluster-level outcomes in the same order as Periods
## StartTimes: A vector of the first period on intervention for each cluster,
###            repeated for each appearance of the cluster in Outcomes
## ContrastFunc: A function that takes a "treated" value as the first argument and "control" value as the second argument
###             and returns one value that specifies the causal contrast of interest
# Outputs:
## A list with one element:
### Est.NPWP: the estimated treatment effect on the specified contrast scale

# NPWP.Effect.Est <- function(Periods, Outcomes, StartTimes, ContrastFunc) {
#   PercInt <- tapply(ifelse(StartTimes-Periods<=0,1,0), Periods, FUN=mean)
#   ContrastPds <- as.numeric(labels(PercInt[PercInt > 0 & PercInt < 1])[[1]])
#   ContrastPds <- ContrastPds[order(ContrastPds)]
#   Contrasts <- matrix(rep(NA,2*length(ContrastPds)),nrow=2, ncol=length(ContrastPds))
#   for (i in 1:length(ContrastPds)) {
#     Pd <- ContrastPds[i]
#     Ints <- Outcomes[Periods == Pd & StartTimes <= Pd]
#     Conts <- Outcomes[Periods == Pd & StartTimes > Pd]
#     NumCont <- length(Conts)
#     ContMean <- mean(Conts)
#     ContVar <- ifelse(NumCont==1,0,var(Conts))
#     NumInt <- length(Ints)
#     IntMean <- mean(Ints)
#     IntVar <- ifelse(NumInt==1,0,var(Ints))
#     Contrast <- ContrastFunc(IntMean, ContMean)
#     Weight <- ifelse(IntVar==0 & ContVar==0, (10^(-5)*(1/NumCont + 1/NumInt))^(-1), 
#                      ((((NumCont-1)*ContVar + (NumInt-1)*IntVar)/(NumCont+NumInt-2))*(1/NumCont + 1/NumInt))^(-1))
#     if (is.finite(Contrast)) {
#       Contrasts[,i] <- c(Contrast,Weight)
#     }
#   }
#   TxEst <- sum(Contrasts[1,]*Contrasts[2,], na.rm=TRUE)/sum(Contrasts[2,], na.rm=TRUE)
#   return(list(Est.NPWP=TxEst,ContrastMat=Contrasts))
# }


## Adjusted NPWP Method: Provides Cluster-Pd Specific Results AND Accounts for Zeros in RR/OR Contrasts ##
# Inputs:
## Periods: A vector of period numbers for each cluster
## Outcomes: A vector of cluster-level outcomes in the same order as Periods
## Clusters: A vector of cluster numbers/names in the same order as Periods
## StartTimes: A vector of the first period on intervention for each cluster,
###            repeated for each appearance of the cluster in Outcomes
## ContrastFunc: A function that takes a "treated" value as the first argument and "control" value as the second argument
###             and returns one value that specifies the causal contrast of interest
# Outputs:
## A list with one element:
### Est.NPWP: the estimated treatment effect on the specified contrast scale

NPWPAdj.Effect.Est <- function(Periods, Outcomes, Clusters, StartTimes, ContrastFunc, 
                               Indivs=NULL, ContrastName=NULL, 
                               EstsOnly=FALSE) {
  PercInt <- tapply(ifelse(StartTimes-Periods<=0,1,0), Periods, FUN=mean)
  ContrastPds <- as.numeric(labels(PercInt[PercInt > 0 & PercInt < 1])[[1]])
  ContrastPds <- ContrastPds[order(ContrastPds)]
  for (k in 1:length(ContrastFunc)) {
    assign(x=paste0("Contrasts_",k), value=NULL)
  }
  
  for (i in 1:length(ContrastPds)) {
    Pd <- ContrastPds[i]
    Ints <- Outcomes[Periods == Pd & StartTimes <= Pd]
    FirstPds <- ifelse(StartTimes[Periods == Pd & StartTimes <= Pd]==Pd,1,0)
    IntClusts <- Clusters[Periods == Pd & StartTimes <= Pd]
    StartPds <- StartTimes[Periods == Pd & StartTimes <= Pd]
    Conts <- Outcomes[Periods == Pd & StartTimes > Pd]
    NumCont <- length(Conts)
    ContMean <- mean(Conts)
    ContIndivs <- Indivs[Periods==Pd & StartTimes > Pd]
    ContVar <- ifelse(NumCont==1,0,var(Conts))
    NumInt <- length(Ints)
    IntVar <- ifelse(NumInt==1,0,var(Ints))
    AdjInts <- Ints
    AdjConts <- rep(ContMean,NumInt)
    OverallCases <- sum(Outcomes[Periods==Pd]*Indivs[Periods==Pd])
    OverallMean <- mean(c(Ints,Conts))
    ThompWeight <- rep(ifelse(IntVar==0 & ContVar==0, (10^(-5)*(1/NumCont + 1/NumInt))^(-1), 
                          ((((NumCont-1)*ContVar + (NumInt-1)*IntVar)/(NumCont+NumInt-2))*(1/NumCont + 1/NumInt))^(-1))/NumInt,
                       NumInt)
    
    for (k in 1:length(ContrastFunc)) {
      CFunck <- ContrastFunc[[k]]
      if (all.equal(CFunck,logRRFunc)==TRUE | all.equal(CFunck,logORFunc)==TRUE) {
        if (sum(AdjInts==0) > 0 | ContMean==0) {
          for (j in 1:NumInt) {
            if (xor(AdjInts[j] == 0,ContMean == 0)) {
              AdjInts[j] <- (Ints[j]*Indivs[Clusters==IntClusts[j] & Periods==Pd] + 1/2)/(Indivs[Clusters==IntClusts[j] & Periods==Pd] + 1)
              AdjConts[j] <- mean((ContIndivs*Conts+1/2)/(ContIndivs+1))
            }
          }
        }
      }
    
    ContrastVec <- CFunck(AdjInts, AdjConts)
    ContrastVec_Zero <- ifelse(is.finite(ContrastVec),ContrastVec,0)
    NonZeroWt <- ifelse(is.finite(ContrastVec),1,0)
    assign(x=paste0("Contrasts_",k), value=rbind(get(paste0("Contrasts_",k)),
                                                  data.frame(Cluster=IntClusts, Period=rep(Pd,NumInt), 
                                                             Contrast=ContrastVec_Zero, UseCP=NonZeroWt,
                                                             StartPd=StartPds, 
                                                             ThompWeight=ThompWeight, FirstPd=FirstPds,
                                                             OverallWeight=rep(OverallCases,NumInt))))
    }
  }
  
  outlist2 <- NULL
  for (k in 1:length(ContrastFunc)) {
    Contrasts <- get(paste0("Contrasts_",k))
    TxEst.TW <- sum(Contrasts$Contrast*Contrasts$ThompWeight*Contrasts$UseCP)/sum(Contrasts$ThompWeight*Contrasts$UseCP)
    TxEst.Mean <- sum(Contrasts$Contrast*Contrasts$UseCP)/sum(Contrasts$UseCP)
    TxEst.MeanEx1 <- sum(Contrasts$Contrast[Contrasts$FirstPd==0]*Contrasts$UseCP[Contrasts$FirstPd==0])/sum(Contrasts$UseCP[Contrasts$FirstPd==0])
    TxEst.OverallWt <- sum(Contrasts$Contrast*Contrasts$UseCP*Contrasts$OverallWeight)/sum(Contrasts$UseCP*Contrasts$OverallWeight)
    TxEst.OverallWtEx1 <- sum(Contrasts$Contrast[Contrasts$FirstPd==0]*Contrasts$UseCP[Contrasts$FirstPd==0]*Contrasts$OverallWeight[Contrasts$FirstPd==0])/sum(Contrasts$UseCP[Contrasts$FirstPd==0]*Contrasts$OverallWeight[Contrasts$FirstPd==0])
    if (EstsOnly) {
      outlist <- list(list(Est.NPWP.TW=TxEst.TW,
                           Est.NPWP.Mean=TxEst.Mean,
                           Est.NPWP.MeanEx1=TxEst.MeanEx1,
                           Est.NPWP.OW=TxEst.OverallWt,
                           Est.NPWP.OWEx1=TxEst.OverallWtEx1))
    } else {
      outlist <- list(list(Est.NPWP.TW=TxEst.TW,
                           Est.NPWP.Mean=TxEst.Mean,
                           Est.NPWP.MeanEx1=TxEst.MeanEx1,
                           Est.NPWP.OW=TxEst.OverallWt,
                           Est.NPWP.OWEx1=TxEst.OverallWtEx1,
                           ContrastMat=Contrasts))
    }
    names(outlist) <- ContrastName[k]
    outlist2 <- append(outlist2, outlist)
  }
  if (length(ContrastFunc)==1) {
    outlist2 <- unlist(outlist2)
  }
  return(outlist2)
}


## Crossover Method from Kennedy-Shaffer et al. 2019 ##
# Inputs:
## Periods: A vector of period numbers for each cluster
## Outcomes: A vector of cluster-level outcomes in the same order as Periods
## Clusters: A vector of cluster numbers/names in the same order as Periods
## StartTimes: A vector of the first period on intervention for each cluster,
###            repeated for each appearance of the cluster in Outcomes
## ContrastFunc: A function that takes a "treated" value as the first argument and "control" value as the second argument
###             and returns one value that specifies the causal contrast of interest
## CtrlType: one  of "Ctrl","Both","CtWt","BoWt", indicating whether only clusters
###          on control in both periods should be used as comparisons ("Ctrl"/"CtWt")
###          or cluster on intervention in both periods should be used as well ("Both"/"BoWt")
###          and whether the weights should be equal across periods ("Ctrl"/"Both")
###          or proportional to the harmonic mean of the number of intervention and comparison clusters  ("CtWt"/"BoWt").
###          Multiple values may be specified as a vector to get multiple effect estimates.
###          In the article, "Ctrl" is CO-1, "CtWt" is CO-2, and "Both" is CO-3.
# Outputs:
## A list with one element for each value in CtrlType:
### Est.CO.[CtrlType]: the estimated treatment effect for the specified CtrlType

CO.Effect.Est <- function(Periods, Outcomes, Clusters, StartTimes, ContrastFunc, 
                          Indivs=NULL, ContrastName=NULL, 
                          CtrlType=c("Ctrl","Both","CtWt","BoWt"), FwdOffset=0, BwdOffset=0, 
                          EstsOnly=FALSE) {
  Pds <- sort(unique(Periods))
  CPType <- ifelse(Periods - StartTimes >= FwdOffset,"Int",ifelse(Periods - StartTimes < -1*BwdOffset,"Cont","None"))
  XODiff <- FwdOffset + BwdOffset + 1
  PercInt <- tapply(ifelse(CPType=="Int",1,0), Periods, FUN=mean)
  ContrastPds <- as.numeric(labels(PercInt[PercInt > 0])[[1]])
  ContrastPds <- sort(ContrastPds)
  ContrastPds <- ContrastPds[ContrastPds > Pds[XODiff]]
  StartPdsVec <- StartTimes[Periods==Pds[1]]
  names(StartPdsVec) <- as.character(Clusters[Periods==Pds[1]])
  
  for (k in 1:length(ContrastFunc)) {
    assign(x=paste0("Contrasts_",k), NULL)
    assign(x=paste0("CPContrasts_",k), NULL)
  }
  
  for (i in 1:length(ContrastPds)) {
    Pdi <- ContrastPds[i]
    Clusts1 <- intersect(Clusters[Periods == Pdi & CPType=="Int"],Clusters[Periods == Pdi-XODiff & CPType=="Cont"])
    Clusts2 <- intersect(Clusters[Periods == Pdi & CPType=="Int"],Clusters[Periods == Pdi-XODiff & CPType=="Int"])
    Clusts0 <- intersect(Clusters[Periods == Pdi & CPType=="Cont"],Clusters[Periods == Pdi-XODiff & CPType=="Cont"])
    n1 <- length(Clusts1)
    n0 <- length(Clusts0)
    n2 <- length(Clusts2)
    
    AdjOutcomes <- Outcomes[Periods==Pdi]
    names(AdjOutcomes) <- Clusters[Periods==Pdi]
    AdjOutcomes <- AdjOutcomes[as.character(c(Clusts1,Clusts2,Clusts0))]
    AdjPrevOuts <- Outcomes[Periods==Pdi-XODiff]
    names(AdjPrevOuts) <- Clusters[Periods==Pdi-XODiff]
    AdjPrevOuts <- AdjPrevOuts[as.character(c(Clusts1,Clusts2,Clusts0))]
    
    for (k in 1:length(ContrastFunc)) {
      CFunck <- ContrastFunc[[k]]
      if (all.equal(CFunck,logRRFunc)==TRUE | all.equal(CFunck,logORFunc)==TRUE) {
        if (sum(AdjOutcomes==0) > 0 | sum(AdjPrevOuts==0) > 0) {
          if (is.null(Indivs)) {
            warning("Need Indivs to Handle Zero Processing. Results May Be Non-Finite")
          } else {
            for (i in as.character(c(Clusts1,Clusts2,Clusts0))) {
              if (AdjOutcomes[i]==0 | AdjPrevOuts[i]==0) {
                AdjOutcomes[i] <- (AdjOutcomes[i]*Indivs[Periods == Pdi & Clusters==as.numeric(i)] + 1/2)/(Indivs[Periods == Pdi & Clusters==as.numeric(i)] + 1/2)
                AdjPrevOuts[i] <- (AdjPrevOuts[i]*Indivs[Periods == Pdi-XODiff & Clusters==as.numeric(i)] + 1/2)/(Indivs[Periods == Pdi-XODiff & Clusters==as.numeric(i)] + 1/2)
              }
            }
          }
        }
      }

    ContrastVec <- CFunck(AdjOutcomes, AdjPrevOuts)
    TypeVec <- ifelse(names(ContrastVec) %in% as.character(Clusts1),1,
                      ifelse(names(ContrastVec) %in% as.character(Clusts0),0,
                             ifelse(names(ContrastVec) %in% as.character(Clusts2),2,NA)))
    StartPds <- unname(StartPdsVec[as.character(names(ContrastVec))])
    # ContrastDF <- data.frame(Clusters=names(ContrastVec),Types=TypeVec,Contrasts=unlist(ContrastVec))
    # CPContrasts <- get(paste0("CPContrasts_",k))
    # CPContrasts[[Pdi]] <- ContrastDF
    assign(x=paste0("CPContrasts_",k), value=rbind(get(paste0("CPContrasts_",k)),
                                                   data.frame(Cluster=names(ContrastVec),
                                                              Period=rep(Pdi,length(ContrastVec)),
                                                              Contrast=unname(ContrastVec),
                                                              StartPd=StartPds,
                                                              Type=TypeVec)))
    
    Contrast.Ctrl <- ifelse(n0 > 0, mean(ContrastVec[TypeVec==1]) - mean(ContrastVec[TypeVec==0]), NA)
    Contrast.Both <- ifelse(n0 + n2 > 0, mean(ContrastVec[TypeVec==1]) - mean(ContrastVec[TypeVec==0 | TypeVec==2]), NA)
    Weight.CtWt <- ifelse(n0 > 0, (1/n1 + 1/n0)^(-1), NA)
    Weight.BoWt <- ifelse(n0 + n2 > 0, (1/n1 + 1/(n0+n2))^(-1), NA)
    
    assign(x=paste0("Contrasts_",k),
           value=rbind(get(paste0("Contrasts_",k)),
                       data.frame(Cont.Ctrl=Contrast.Ctrl, Cont.Both=Contrast.Both,
                                  Wt.CtWt=Weight.CtWt, Wt.BoWt=Weight.BoWt)))
    }
  }
  
  outlist2 <- NULL
  for (k in 1:length(ContrastFunc)) {
    Contrasts <- get(paste0("Contrasts_",k))
    if (EstsOnly) {
      outlist <- NULL
    } else {
      outlist <- list(ContrastMat=get(paste0("CPContrasts_",k)),
                    PdConts=Contrasts)
    }
    if ("Ctrl" %in% CtrlType) {
      TxEff <- mean(Contrasts$Cont.Ctrl, na.rm=TRUE)
      outlist <- append(outlist, list(Est.CO.Ctrl=TxEff))
    }
    if ("Both" %in% CtrlType) {
      TxEff <- mean(Contrasts$Cont.Both, na.rm=TRUE)
      outlist <- append(outlist, list(Est.CO.Both=TxEff))
    }
    if ("CtWt" %in% CtrlType) {
      TxEff <- sum(Contrasts$Cont.Ctrl*Contrasts$Wt.CtWt, na.rm=TRUE)/sum(Contrasts$Wt.CtWt, na.rm=TRUE)
      outlist <- append(outlist, list(Est.CO.CtWt=TxEff))
    }
    if ("BoWt" %in% CtrlType) {
      TxEff <- sum(Contrasts$Cont.Both*Contrasts$Wt.BoWt, na.rm=TRUE)/sum(Contrasts$Wt.BoWt, na.rm=TRUE)
      outlist <- append(outlist, list(Est.CO.BoWt=TxEff))
    }
    outlist <- list(outlist)
    names(outlist) <- ContrastName[k]
    outlist2 <- append(outlist2, outlist)
  }
  
  return(outlist2)
}


## Mixed Effects Model Method from Hussey & Hughes 2007 ##
# Inputs:
## Periods: A vector of period numbers for each cluster
## Outcomes: A vector of cluster-level outcomes in the same order as Periods
## Clusters: A vector of cluster numbers/names in the same order as Periods
## StartTimes: A vector of the first period on intervention for each cluster,
###            repeated for each appearance of the cluster in Outcomes
## Indivs: For binary outcomes, a vector of the number of individuals measured in each
###        cluster-period, in the same order as Periods; or, a single number if the number
###        of individuals measured is the same in each cluster-period.
###        For non-binary outcomes, leave as NULL and specify each individual's
###        outcome, cluster, period, and start time in the respective input vector.
## family: a GLM family, see family or glm. E.g., binomial or gaussian
## link: a GLM link, see family or glm. E.g., "logit" or "identity"
# Outputs:
## A list with two elements:
### Est.MEM: the estimated treatment effect
### PVal.MEM: the asymptotic p-value for the null hypothesis of no treatment effect

MEM.Effect.Est <- function(Periods, Outcomes, Clusters, StartTimes, Indivs=NULL,
                           family, link, EstsOnly=FALSE) {
  Trts <- ifelse(Periods >= StartTimes, 1, 0)
  
  if (!is.null(Indivs)) {
    if (length(Indivs)==1) {
      Indivs = rep(Indivs, length(Periods))
    }
    NewTrts <- rep(Trts, times=Indivs)
    NewClusts <- rep(Clusters, times=Indivs)
    NewPds <- rep(Periods, times=Indivs)
    NewOuts <- NULL
    for (i in 1:length(Outcomes)) {
      NewOuts <- append(NewOuts, rep(1,round(Outcomes[i]*Indivs[i], digits=0)))
      NewOuts <- append(NewOuts, rep(0,Indivs[i]-round(Outcomes[i]*Indivs[i], digits=0)))
    }
  } else {
    NewTrts <- Trts
    NewClusts <- Clusters
    NewPds <- Periods
    NewOuts <- Outcomes
  }
  PdFactor <- relevel(factor(NewPds), ref=1)
  Indexes <- NewClusts*max(NewPds)+NewPds
  fam <- family(link=link)
  capture.output(res <- tryCatch({glmer(NewOuts~NewTrts+PdFactor+(1|NewClusts), 
                                        family=fam)},
                                      error=function(err) {1}), file="/dev/null")
  if (is.numeric(res)) {
    Est <- NA
    PVal <- NA
  } else {
    Est <- as.numeric(summary(res)$coefficients['NewTrts','Estimate'])
    SE <- as.numeric(summary(res)$coefficients['NewTrts','Std. Error'])
    PVal <- as.numeric(summary(res)$coefficients['NewTrts','Pr(>|z|)'])
  }
  if (EstsOnly) {
    return(list(Est.MEM=Est))
  } else {
    return(list(Est.MEM=Est, PVal.MEM=PVal))
  }
}

## Mixed Effects Model with Cluster-Period Interaction Method from Hooper et al. 2016 ##
# Inputs:
## Periods: A vector of period numbers for each cluster
## Outcomes: A vector of cluster-level outcomes in the same order as Periods
## Clusters: A vector of cluster numbers/names in the same order as Periods
## StartTimes: A vector of the first period on intervention for each cluster,
###            repeated for each appearance of the cluster in Outcomes
## Indivs: For binary outcomes, a vector of the number of individuals measured in each
###        cluster-period, in the same order as Periods; or, a single number if the number
###        of individuals measured is the same in each cluster-period
###        For non-binary outcomes, leave as NULL and specify each individual's
###        outcome, cluster, period, and start time in the respective input vector.
## family: a GLM family, see family or glm. E.g., binomial or gaussian
## link: a GLM link, see family or glm. E.g., "logit" or "identity"
# Outputs:
## A list with two elements:
### Est.CPI: the estimated treatment effect
### PVal.CPI: the asymptotic p-value for the null hypothesis of no treatment effect

CPI.Effect.Est <- function(Periods, Outcomes, Clusters, StartTimes, Indivs=NULL,
                           family, link, EstsOnly=FALSE) {
  Trts <- ifelse(Periods >= StartTimes, 1, 0)
  
  if (!is.null(Indivs)) {
    if (length(Indivs)==1) {
      Indivs = rep(Indivs, length(Periods))
      }
    NewTrts <- rep(Trts, times=Indivs)
    NewClusts <- rep(Clusters, times=Indivs)
    NewPds <- rep(Periods, times=Indivs)
    NewOuts <- NULL
    for (i in 1:length(Outcomes)) {
      NewOuts <- append(NewOuts, rep(1,round(Outcomes[i]*Indivs[i], digits=0)))
      NewOuts <- append(NewOuts, rep(0,Indivs[i]-round(Outcomes[i]*Indivs[i], digits=0)))
    }
  } else {
    NewTrts <- Trts
    NewClusts <- Clusters
    NewPds <- Periods
    NewOuts <- Outcomes
  }
  PdFactor <- relevel(factor(NewPds), ref=1)
  Indexes <- NewClusts*max(NewPds)+NewPds
  fam <- family(link=link)
  capture.output(res <- tryCatch({glmer(NewOuts~NewTrts+PdFactor+(1|NewClusts)+(1|Indexes), 
                                        family=fam)},
                                 error=function(err) {1}), file="/dev/null")
  if (is.numeric(res)) {
    Est <- NA
    PVal <- NA
  } else {
    Est <- as.numeric(summary(res)$coefficients['NewTrts','Estimate'])
    SE <- as.numeric(summary(res)$coefficients['NewTrts','Std. Error'])
    PVal <- as.numeric(summary(res)$coefficients['NewTrts','Pr(>|z|)'])
  }
  if (EstsOnly) {
    return(list(Est.CPI=Est))
  } else {
    return(list(Est.CPI=Est, PVal.CPI=PVal))
  }
}

## Synthetic Control Method from Kennedy-Shaffer et al. 2019 ##
# Inputs:
## Periods: A vector of period numbers for each cluster
## Outcomes: A vector of cluster-level outcomes in the same order as Periods
## Clusters: A vector of cluster numbers/names in the same order as Periods
## StartTimes: A vector of the first period on intervention for each cluster,
###            repeated for each appearance of the cluster in Outcomes
## ContrastFunc: A function that takes a "treated" value as the first argument and "control" value as the second argument
###             and returns one value that specifies the causal contrast of interest
# Outputs:
## A list with three elements:
### Est.SCSWT1: the estimated treatment effect for equal weighting (SC-1 in the article)
### Est.SCSWT2: the estimated treatment effect for weighting by inverse-MSPE within equivalence groups (SC-2 in the article)
### ContrastMat: the matrix of cluster-period-specific treatment effect estimates, as well as the weights used in SC-1 and SC-2.
####             Can be used for alternate weighting schemes

SCSWT.Effect.Est <- function(Periods, Outcomes, Clusters, StartTimes, ContrastFunc, 
                             Indivs=NULL, ContrastName=NULL, EstsOnly=FALSE) {
  
  Index <- Clusters*max(Periods)+Periods
  DataMat <- matrix(data=Outcomes[order(Index)], nrow=length(unique(Periods)), 
                    ncol=length(unique(Clusters)), byrow=FALSE,
                    dimnames=list(sort(unique(Periods)),sort(unique(Clusters))))
  ## Note: DataMat here only works if all clusters have same periods ##
  
  PercInt <- tapply(ifelse(StartTimes-Periods<=0,1,0), Periods, FUN=mean)
  ContrastPds <- as.numeric(labels(PercInt[PercInt > 0 & PercInt < 1])[[1]])
  ContrastPds <- ContrastPds[order(ContrastPds)]
  
  for (k in 1:length(ContrastFunc)) {
    assign(x=paste0("Contrasts_",k), value=NULL)
  }
  
  for (i in 1:length(ContrastPds)) {
    Pd <- ContrastPds[i]
    IntClusts <- Clusters[Periods == Pd & StartTimes <= Pd]
    IntIndivs <- Indivs[Periods==Pd & StartTimes <= Pd]
    FirstPds <- ifelse(StartTimes[Periods == Pd & StartTimes <= Pd]==Pd,1,0)
    ContClusts <- Clusters[Periods == Pd & StartTimes > Pd]
    NumInts <- length(IntClusts)
    NumConts <- length(ContClusts)
    ClustMinVec <- rep(min(NumInts,NumConts),NumInts)
    OverallCases <- sum(Outcomes[Periods==Pd]*Indivs[Periods==Pd])
    OverallMean <- mean(Outcomes[Periods == Pd])
    
    if (Pd == min(Periods)) {
      SCOutcomes <- rep(mean(Outcomes[Periods == Pd & Clusters %in% ContClusts]),NumInts)
      MSPEs <- rep(10^(-8), NumInts)
      StartTimeForDF <- rep(1, NumInts)
    } else {
      PrevPds <- sort(unique(Periods[Periods < Pd]))
      DataMati <- DataMat[as.character(PrevPds),,drop=FALSE]
    
      MSPEs <- rep(10^(-8), NumInts) ## Truncates MSPE at minimum of 10^(-8) to prevent inverse of zero issue
      StartTimeForDF <- rep(NA, NumInts)
      if (length(ContClusts) == 1) {
        SCOutcomes <- rep(Outcomes[Periods == Pd & Clusters == ContClusts], NumInts)
        for (j in 1:NumInts) {
          StartTimej <- min(StartTimes[Clusters==IntClusts[j]])
          StartTimeForDF[j] <- StartTimej
          PrevPdsj <- PrevPds[PrevPds < StartTimej]
          MSPEs[j] <- max(MSPEs[j], mean((DataMati[as.character(PrevPdsj),as.character(ContClusts)]-DataMati[as.character(PrevPdsj),as.character(IntClusts[j])])^2))
        }
      } else {
        SCOutcomes <- rep(0, NumInts)
        for (j in 1:NumInts) {
          StartTimej <- min(StartTimes[Clusters==IntClusts[j]])
          StartTimeForDF[j] <- StartTimej
          PrevPdsj <- PrevPds[PrevPds < StartTimej]
          X0ij <- DataMati[as.character(PrevPdsj),as.character(ContClusts),drop=FALSE]
          X1ij <- DataMati[as.character(PrevPdsj),as.character(IntClusts[j]),drop=FALSE]
          capture.output(SynthOut <- tryCatch({synth(X1=X1ij, X0=X0ij, Z0=X0ij, Z1=X1ij)},
                                              error=function(err) {1}), file="/dev/null")
          if (is.numeric(SynthOut)) {
            SCOutcomes[j] <- mean(Outcomes[Periods == Pd & Clusters %in% ContClusts])
            MSPEs[j] <- max(MSPEs[j],mean((X1ij - 
                                apply(X=X0ij,MARGIN=1,FUN=mean))^2))
          } else {
            MSPEs[j] <- max(MSPEs[j],as.numeric(SynthOut$loss.v))
            SCOutcomes[j] <- as.numeric(SynthOut$solution.w) %*% Outcomes[Periods == Pd & Clusters %in% ContClusts]
          }
        }
      }
    }
    TargOutcomes <- Outcomes[Periods == Pd & StartTimes <= Pd]
    for (k in 1:length(ContrastFunc)) {
      CFunck <- ContrastFunc[[k]]
      Contrastsk <- get(paste0("Contrasts_",k))
      if (all.equal(CFunck,logRRFunc)==TRUE | all.equal(CFunck,logORFunc)==TRUE) {
        if (xor(sum(TargOutcomes==0) > 0, sum(SCOutcomes==0) > 0)) {
          if (is.null(Indivs)) {
            warning("Need Indivs to Handle Zero Processing. Results May Be Non-Finite")
          } else {
            for (j in 1:length(TargOutcomes)) {
              if (xor(TargOutcomes[j]==0,SCOutcomes[j]==0)) {
                TargOutcomes[j] <- (TargOutcomes[j]*IntIndivs[j] + 1/2)/(IntIndivs[j] + 1)
                SCOutcomes[j] <- (SCOutcomes[j]*IntIndivs[j]+ 1/2)/(IntIndivs[j] + 1)
              }
            }
          }
        }
      }
      ContrastVec <- CFunck(TargOutcomes,SCOutcomes)
      ContrastVec_Zero <- ifelse(is.finite(ContrastVec),ContrastVec,0)
      ZeroWt <- ifelse(is.finite(ContrastVec),1,0)
      Contrastsk <- rbind(Contrastsk, data.frame(Cluster=IntClusts,
                                                 Period=rep(Pd, NumInts),
                                                 Contrast=ContrastVec_Zero,
                                                 UseCP=ZeroWt, OverallWeight=rep(OverallCases,NumInts),
                                                 StartPd=StartTimeForDF,
                                                 MSPE=MSPEs,
                                                 ClustMin=ClustMinVec,
                                                 FirstPd=FirstPds))
      assign(x=paste0("Contrasts_",k), value=Contrastsk)
    }
  }
  
  outlist2 <- NULL
  for (k in 1:length(ContrastFunc)) {
    Contrasts <- get(paste0("Contrasts_",k))
    Contrasts$PdMin <- rep(0,length(Contrasts$Contrast))
    Contrasts$SumMSPE <- rep(0,length(Contrasts$Contrast))
    for (i in 1:length(Contrasts$Contrast)) {
      Clusteri <- Contrasts$Cluster[i]
      ClustStart <- Contrasts$StartPd[i]
      Contrasts$SumMSPE[i] <- sum(Contrasts$UseCP[Contrasts$StartPd==ClustStart]/Contrasts$MSPE[Contrasts$StartPd==ClustStart])
      Contrasts$PdMin[i] <- min(sum(ContrastPds >= ClustStart),
                                length(unique(Periods[Periods < ClustStart])))
    }
    Contrasts$Wt2 <- ifelse(Contrasts$SumMSPE==0,0,(Contrasts$MSPE*Contrasts$SumMSPE)^(-1)*Contrasts$UseCP)
    
    TxEst1 <- sum(Contrasts$Contrast*Contrasts$UseCP)/sum(Contrasts$UseCP)
    TxEst2 <- mean(Contrasts$Contrast*Contrasts$Wt2)
    TxEst3 <- sum(Contrasts$Contrast*Contrasts$UseCP*Contrasts$OverallWeight)/sum(Contrasts$UseCP*Contrasts$OverallWeight)
    TxEst.1.Ex1 <- sum(Contrasts$Contrast[Contrasts$FirstPd==0]*Contrasts$UseCP[Contrasts$FirstPd==0])/sum(Contrasts$UseCP[Contrasts$FirstPd==0])
    TxEst.2.Ex1 <- sum(Contrasts$Contrast[Contrasts$FirstPd==0]*Contrasts$Wt2[Contrasts$FirstPd==0])/sum(Contrasts$Wt2[Contrasts$FirstPd==0])
    TxEst.3.Ex1 <- sum(Contrasts$Contrast[Contrasts$FirstPd==0]*Contrasts$UseCP[Contrasts$FirstPd==0]*Contrasts$OverallWeight[Contrasts$FirstPd==0])/sum(Contrasts$UseCP[Contrasts$FirstPd==0]*Contrasts$OverallWeight[Contrasts$FirstPd==0])
    if(EstsOnly) {
      outlist <- list(list(Est.SCSWT1=TxEst1,Est.SCSWT2=TxEst2,Est.SCSWT3=TxEst3,
                           Est.SCSWT1.Ex1=TxEst.1.Ex1,Est.SCSWT2.Ex1=TxEst.2.Ex1,Est.SCSWT3.Ex1=TxEst.3.Ex1))
    } else {
      outlist <- list(list(Est.SCSWT1=TxEst1,Est.SCSWT2=TxEst2,Est.SCSWT3=TxEst3,
                           Est.SCSWT1.Ex1=TxEst.1.Ex1,Est.SCSWT2.Ex1=TxEst.2.Ex1,Est.SCSWT3.Ex1=TxEst.3.Ex1,
                         ContrastMat=Contrasts))
    }
    names(outlist) <- ContrastName[k]
    outlist2 <- append(outlist2,outlist)
  }
  if (length(ContrastFunc)==1) {
    outlist2 <- unlist(outlist2)
  }
  return(outlist2)
}



## Crossover Synthetic Control Method from Kennedy-Shaffer et al. 2019 ##
# Inputs:
## Periods: A vector of period numbers for each cluster
## Outcomes: A vector of cluster-level outcomes in the same order as Periods
## Clusters: A vector of cluster numbers/names in the same order as Periods
## StartTimes: A vector of the first period on intervention for each cluster,
###            repeated for each appearance of the cluster in Outcomes
## ContrastFunc: A function that takes a "treated" value as the first argument and "control" value as the second argument
###             and returns one value that specifies the causal contrast of interest
# Outputs:
## A list with three elements:
### Est.COSC1: the estimated treatment effect for equal weighting (COSC-1 in the article)
### Est.COSC2: the estimated treatment effect for weighting by inverse-MSPE within equivalence groups (COSC-2 in the article)
### ContrastMat: the matrix of cluster-period-specific treatment effect estimates, as well as the weights used in SC-1 and SC-2.
####             Can be used for alternate weighting schemes

# COSC.Effect.Est <- function(Periods, Outcomes, Clusters, StartTimes, ContrastFunc) {
#   Pds <- sort(unique(Periods))
#   Clusts <- sort(unique(Clusters))
#   StartTimesOrd <- StartTimes[order(Clusters)]
#   StartTimesByClust <- StartTimesOrd[seq(length(Pds),length(StartTimes),by=length(Pds))]
#   names(StartTimesByClust) <- as.character(Clusts)
#   COPds <- Pds[-1]
#   Index <- Clusters*max(Periods)+Periods
#   DataMat <- matrix(data=Outcomes[order(Index)], nrow=length(Pds), 
#                     ncol=length(Clusts), byrow=FALSE,
#                     dimnames=list(Pds,Clusts))
#   ## Note: DataMat here only works if all clusters have same periods ##
#   COMat <- ContrastFunc(DataMat[2:length(Pds),], DataMat[1:length(COPds),])
#   TypeMat <- COMat
#   for (i in 1:(length(COPds))) {
#     Pdi <- COPds[i]
#     TypeMat[i,] <- rep(0,dim(TypeMat)[2]) + (Pdi >= StartTimesByClust) + (Pdi > StartTimesByClust)
#   }
#   ContrastPds <- COPds[apply(TypeMat,1,function(x) sum(x==0)) > 0 & apply(TypeMat,1,function(x) sum(x==1)) > 0]
#   Contrasts <- data.frame(Period = NULL, Cluster = NULL,Contrast = NULL, MSPE = NULL,
#                           ClustMin = NULL)
#   
#   for (i in 1:length(ContrastPds)) {
#     Pd <- ContrastPds[i]
#     COMatPd <- COMat[i,]
#     TypeMatPd <- TypeMat[i,]
#     n1 <- sum(TypeMatPd==1)
#     n0 <- sum(TypeMatPd==0)
#     IntClusts <- as.character(names(TypeMatPd[TypeMatPd==1]))
#     ContClusts <- as.character(names(TypeMatPd[TypeMatPd==0]))
#     ClustMinVec <- rep(min(n1,n0),n1)
#     
#     if (Pd == min(COPds)) {
#       SCOutcomes <- rep(mean(COMatPd[TypeMatPd==0]),n1)
#       MSPEs <- rep(10^(-8), n1)
#     } else {
#       PrevPds <- COPds[COPds < Pd]
#       COMati <- COMat[as.character(PrevPds),,drop=FALSE]
#       
#       MSPEs <- rep(10^(-8), n1) ## Truncates MSPE at minimum of 10^(-8) to prevent inverse of zero issue
#       if (n0 == 1) {
#         SCOutcomes <- rep(COMatPd[TypeMatPd==0], n1)
#         for (j in 1:n1) {
#           MSPEs[j] <- max(MSPEs[j], 
#                           mean((COMati[as.character(PrevPds),TypeMatPd==0]-COMati[as.character(PrevPds),IntClusts[j]])^2))
#         }
#       } else {
#         SCOutcomes <- rep(0, n1)
#         X0i <- COMati[,TypeMatPd==0,drop=FALSE]
#         for (j in 1:n1) {
#           X1ij <- COMati[,IntClusts[j],drop=FALSE]
#           capture.output(SynthOut <- tryCatch({synth(X1=X1ij, X0=X0i, Z0=X0i, Z1=X1ij)},
#                                               error=function(err) {1}), file="/dev/null")
#           if (is.numeric(SynthOut)) {
#             SCOutcomes[j] <- mean(COMatPd[TypeMatPd==0])
#             MSPEs[j] <- max(MSPEs[j],mean((X1ij - 
#                                              apply(X=X0i,MARGIN=1,FUN=mean))^2))
#           } else {
#             MSPEs[j] <- max(MSPEs[j],as.numeric(SynthOut$loss.v))
#             SCOutcomes[j] <- as.numeric(SynthOut$solution.w) %*% COMatPd[TypeMatPd==0]
#           }
#         }
#       }
#     }
#     ContrastVec <- as.numeric(COMatPd[TypeMatPd==1] - SCOutcomes)
#     Contrasts <- data.frame(Period = append(Contrasts$Period, rep(Pd, n1)),
#                             Cluster = append(Contrasts$Cluster, as.numeric(IntClusts)),
#                             Contrast = append(Contrasts$Contrast, ContrastVec),
#                             MSPE = append(Contrasts$MSPE, MSPEs),
#                             ClustMin = append(Contrasts$ClustMin, ClustMinVec))
#   }
#   Contrasts$SumMSPE <- rep(0,length(Contrasts$Contrast))
#   for (i in 1:length(Contrasts$Contrast)) {
#     Periodi <- Contrasts$Period[i]
#     Contrasts$SumMSPE[i] <- sum(1/Contrasts$MSPE[Contrasts$Period==Periodi])
#   }
#   Contrasts$Wt2 <- (Contrasts$MSPE*Contrasts$SumMSPE)^(-1)
# 
#   TxEst1 <- mean(Contrasts$Contrast)
#   TxEst2 <- sum(Contrasts$Contrast*Contrasts$Wt2)/sum(Contrasts$Wt2)
#   
#   return(list(Est.COSC1=TxEst1,Est.COSC2=TxEst2,ContrastMat=Contrasts))
# }


## Generic Ensemble Method from Kennedy-Shaffer et al. 2019 ##
# Inputs:
## Estimates: A named vector or list of treatment effect estimates that includes those to be used in the ensemble
## Prefix: The prefix that denotes the names of Estimates. For the above methods, "Est."
## Types: A vector of the names (excluding Prefix) of the Estimates to be used in the Method
## Weights: A vector of the weights to be given to the Types. Can sum to 1 or not; if not, weights will be scaled to sum to 1
# Outputs:
## The value of the ensemble estimate with given Weights to each of the given Types.

Ens.Effect.Est <- function(Estimates, Prefix, Types, Weights) {
  if (sum(Weights) != 1) {
    print("Warning: Weights do not sum to 1. Adjusting to maintain proportions given.")
    Weights = Weights/sum(Weights)
  }
  Vals <- Estimates[paste0(Prefix,Types)]
  out <- sum(Vals*Weights)
  return(out)
}

## Specific Ensemble Method from Kennedy-Shaffer et al. 2019 ##
# Inputs:
## Estimates: A named vector or list of treatment effect estimates that includes those to be used in the ensemble
###           Must include values for "SCSWT2" and "CO.CtWt"
## Prefix: The prefix that denotes the names of Estimates. For the above methods, "Est."
# Output:
## A list with one element:
### Est.ENS: the estimator for the ENS method described in the article,
####         an unweighted mean of SC-2 and CO-2
Ens1.Effect.Est <- function(Estimates, Prefix) {
  outEst <- Ens.Effect.Est(Estimates, Prefix, Types=c("SCSWT2","CO.CtWt"), Weights=c(1/2,1/2))
  return(list(Est.ENS=outEst))
}


### Methods for Inference for Above Estimators ###

## Permutated Effect Estimates ##
# Inputs:
## Periods: A vector of period numbers for each cluster
## Outcomes: A vector of cluster-level outcomes in the same order as Periods
## Clusters: A vector of cluster numbers/names in the same order as Periods
## StartTimes: A vector of the first period on intervention for each cluster,
###            repeated for each appearance of the cluster in Outcomes
## Indivs: If "MEM" or "CPI" is included in Type and the outcomes are binary, the number of
###        individuals represented by the cluster-level outcomes in Outcomes.
###        Same format as in MEM.Effect.Est and CPI.Effect.Est
## family: a GLM family, see family or glm. E.g., binomial or gaussian
## link: a GLM link, see family or glm. E.g., "logit" or "identity".
###      For methods other than MEM/CPI, the ContrastFunc will be imputed from family and link.
## NumPerms: the number of permutations to perform. Must be an integer >= 1
## Type: a vector of character strings specifying the methods to use.
###      Each method is specified by the name following "Est." in the outputs of above methods. E.g., "NPWP" or "SCSWT1"
# Note: there is random permutation generation in this method. Set a seed before running for replicability.
# Outputs:
## A list with one element for each Type:
### PermEsts.[Type]: A vector of [NumPerms] effect estimates by the [Type] Method,
###                  one for each permutation conducted, with each permutation
###                  under the null hypothesis of no treatment effect.

Perm.Effects <- function(Periods, Outcomes, Clusters, StartTimes, Indivs=NULL,
                         family, link, NumPerms, 
                         Type=c("MEM","CPI","NPWP",
                                "SCSWT1","SCSWT2",
                                "CO.Ctrl","CO.Both","CO.CtWt","CO.BoWt",
                                "ENS"),
                         FwdOffset=0,BwdOffset=0) {
  Starts <- as.numeric(tapply(StartTimes, Clusters, FUN=min))
  NumPds <- as.numeric(tapply(Periods, Clusters, FUN=length))
  
  if (NumPerms > 0) {
    RandStartsMat <- replicate(n=NumPerms, 
                               expr=rep(sample(x=Starts, size=length(Starts), replace=FALSE), 
                                        times=NumPds))
  }
  
  ContrastFunc <- NULL
  ContrastName <- NULL
  Numks <- length(link)
  for (k in 1:Numks) {
    if (link[k]=="logit") {
      ContrastFunc <- c(ContrastFunc,logORFunc)
      ContrastName <- c(ContrastName,"logOR")
    } else if (link[k]=="log") {
      ContrastFunc <- c(ContrastFunc,logRRFunc)
      ContrastName <- c(ContrastName,"logRR")
    } else if (link[k]=="identity") {
      ContrastFunc <- c(ContrastFunc,RDFunc)
      ContrastName <- c(ContrastName,"RD")
    } else {
      print("Warning: Contrast function unknown for this family. Using difference for non-parametric methods.")
      ContrastFunc <- c(ContrastFunc,RDFunc)
      ContrastName <- c(ContrastName,"RD")
    }
  }
  if (length(ContrastFunc)==1) {
    ContrastFunc <- unlist(ContrastFunc)
  }
  
  for (k in 1:Numks) {
    if (NumPerms==0) {
      assign(x=paste0("outdf_",k), value=data.frame(perm=0))
      outlist0 <- NULL
    } else {
      assign(x=paste0("outdf_",k), value=data.frame(perm=(1:NumPerms)))
    }
  }
  
  if ("MEM" %in% Type) {
    MEM.AsyPVals <- NULL
    for (k in 1:Numks) {
      outdf <- get(paste0("outdf_",k))
      if (NumPerms > 0) {
        Res.MEM <- apply(X=RandStartsMat,
                              MARGIN=2,
                              FUN=function(x) unlist(MEM.Effect.Est(Periods, Outcomes, Clusters, x, 
                                                             Indivs,family,link[k],EstsOnly=TRUE)))
        outdf$PermEsts.MEM <- Res.MEM
        
      } else if (NumPerms == 0) {
        Res.MEM <- MEM.Effect.Est(Periods, Outcomes, Clusters, StartTimes, 
                                  Indivs, family, link[k], EstsOnly=FALSE)
        MEM.AsyPVals[[ContrastName[k]]] <- Res.MEM$PVal.MEM
        outdf$PermEsts.MEM <- Res.MEM$Est.MEM
      }
      assign(x=paste0("outdf_",k), value=outdf)
    }
    if (NumPerms == 0) {
      outlist0 <- append(outlist0, list(MEM.AsyPVals=MEM.AsyPVals))
    }
  }
  
  if ("CPI" %in% Type) {
    CPI.AsyPVals <- NULL
    for (k in 1:Numks) {
      outdf <- get(paste0("outdf_",k))
      if (NumPerms > 0) {
        Res.CPI <- apply(X=RandStartsMat,
                         MARGIN=2,
                         FUN=function(x) unlist(CPI.Effect.Est(Periods, Outcomes, Clusters, x, 
                                                               Indivs,family,link[k],EstsOnly=TRUE)))
        outdf$PermEsts.CPI <- Res.CPI
        
      } else if (NumPerms == 0) {
        Res.CPI <- CPI.Effect.Est(Periods, Outcomes, Clusters, StartTimes, 
                                  Indivs, family, link[k], EstsOnly=FALSE)
        CPI.AsyPVals[[ContrastName[k]]] <- Res.CPI$PVal.CPI
        outdf$PermEsts.CPI <- Res.CPI$Est.CPI
      }
      assign(x=paste0("outdf_",k), value=outdf)
    }
    if (NumPerms == 0) {
      outlist0 <- append(outlist0, list(CPI.AsyPVals=CPI.AsyPVals))
    }
  }

  if ("NPWP" %in% Type) {
    NPWP.Vars <- c("NPWP.TW","NPWP.Mean","NPWP.MeanEx1","NPWP.OW","NPWP.OWEx1")
    invar <- paste0("Est.",NPWP.Vars)
    outvar <- paste0("PermEsts.",NPWP.Vars)
    
    if (NumPerms > 0) {
      FullRes.NPWP <- apply(X=RandStartsMat,
                            MARGIN=2,
                            FUN=function(x) unlist(NPWPAdj.Effect.Est(Periods, Outcomes, Clusters, StartTimes=x, 
                                                            ContrastFunc=ContrastFunc, Indivs, ContrastName,
                                                            EstsOnly=TRUE)))
      for (k in 1:Numks) {
        outdf <- get(paste0("outdf_",k))
        outdf[outvar] <- t(FullRes.NPWP[paste0(ContrastName[k],".",invar),])
        assign(x=paste0("outdf_",k), value=outdf)
      }
    } else if (NumPerms == 0) {
      NPWP.ContrastMats <- NULL
      FullRes.NPWP <- NPWPAdj.Effect.Est(Periods, Outcomes, Clusters, StartTimes,
                                         ContrastFunc, Indivs, ContrastName, EstsOnly=FALSE)
      for (k in 1:Numks) {
        outdf <- get(paste0("outdf_",k))
        Resk <- FullRes.NPWP[[ContrastName[k]]]
        outdf[outvar] <- unlist(Resk)[invar]
        assign(x=paste0("outdf_",k), value=outdf)
        NPWP.ContrastMats[[ContrastName[k]]] <- Resk$ContrastMat
      }
      outlist0 <- append(outlist0, list(NPWP.ContrastMats=NPWP.ContrastMats))
    }
  }
  
  if ("SC" %in% Type) {
    AllSCSWTs <- c("SCSWT1","SCSWT2","SCSWT3","SCSWT1.Ex1","SCSWT2.Ex1","SCSWT3.Ex1")
    invar <- paste0("Est.",AllSCSWTs)
    outvar <- paste0("PermEsts.",AllSCSWTs)
    
    if (NumPerms > 0) {
      FullRes.SC <- apply(X=RandStartsMat,
                            MARGIN=2,
                            FUN=function(x) unlist(SCSWT.Effect.Est(Periods, Outcomes, Clusters, StartTimes=x, 
                                                                      ContrastFunc=ContrastFunc, Indivs, ContrastName,
                                                                      EstsOnly=TRUE)))
      for (k in 1:Numks) {
        outdf <- get(paste0("outdf_",k))
        outdf[outvar] <- t(FullRes.SC[paste0(ContrastName[k],".",invar),])
        assign(x=paste0("outdf_",k), value=outdf)
      }
    } else if (NumPerms == 0) {
      SC.ContrastMats <- NULL
      FullRes.SC <- SCSWT.Effect.Est(Periods, Outcomes, Clusters, StartTimes,
                                         ContrastFunc, Indivs, ContrastName, EstsOnly=FALSE)
      for (k in 1:Numks) {
        outdf <- get(paste0("outdf_",k))
        Resk <- FullRes.SC[[ContrastName[k]]]
        outdf[outvar] <- unlist(Resk)[invar]
        assign(x=paste0("outdf_",k), value=outdf)
        SC.ContrastMats[[ContrastName[k]]] <- Resk$ContrastMat
      }
      outlist0 <- append(outlist0, list(SC.ContrastMats=SC.ContrastMats))
    }
  }
  
  if (sum(c("CO.Ctrl","CO.Both","CO.CtWt","CO.BoWt") %in% Type) > 0) {
    types <- c("Ctrl","Both","CtWt","BoWt")
    Usetypes <- types[paste0("CO.",types) %in% Type]
    invar <- paste0("Est.CO.",Usetypes)
    outvar <- paste0("PermEsts.CO.",Usetypes)

    if (NumPerms > 0) {
      FullRes.CO <- apply(X=RandStartsMat,
                            MARGIN=2,
                            FUN=function(x) unlist(CO.Effect.Est(Periods, Outcomes, Clusters, StartTimes=x, 
                                                                 ContrastFunc=ContrastFunc, Indivs, ContrastName,
                                                                 CtrlType=Usetypes, FwdOffset, BwdOffset, 
                                                                 EstsOnly=TRUE)))
      for (k in 1:Numks) {
        outdf <- get(paste0("outdf_",k))
        outdf[outvar] <- t(FullRes.CO[paste0(ContrastName[k],".",invar),])
        assign(x=paste0("outdf_",k), value=outdf)
      }
    } else if (NumPerms == 0) {
      CO.ContrastMats <- NULL
      FullRes.CO <- CO.Effect.Est(Periods, Outcomes, Clusters, StartTimes, 
                                  ContrastFunc, Indivs, ContrastName, 
                                  CtrlType=Usetypes, FwdOffset, BwdOffset, EstsOnly=FALSE)
      for (k in 1:Numks) {
        outdf <- get(paste0("outdf_",k))
        Resk <- FullRes.CO[[ContrastName[k]]]
        outdf[outvar] <- unlist(Resk)[invar]
        assign(x=paste0("outdf_",k), value=outdf)
        CO.ContrastMats[[ContrastName[k]]] <- Resk$ContrastMat
      }
      outlist0 <- append(outlist0, list(CO.ContrastMats=CO.ContrastMats))
    }
  }

  
  if ("ENS" %in% Type) {
    for (k in 1:Numks) {
      outdf <- get(paste0("outdf_",k))
      outlistint <- apply(X=t(as.matrix(outdf)), MARGIN=2,
                          FUN=function(x) Ens1.Effect.Est(x, Prefix="PermEsts.")$Est.ENS)
      outdf$PermEsts.ENS <- outlistint
      assign(x=paste0("outdf_",k), value=outdf)
    }
  }
  
  if (Numks==1) {
    return(outdf_1)
  } else {
    outlist <- list()
    for (k in 1:Numks) {
      outlist <- append(outlist, list(get(paste0("outdf_",k))))
    }
    names(outlist) <- ContrastName
    if (NumPerms==0) {
      outlist <- append(outlist, outlist0)
    }
    return(outlist)
  }
}

## Hypothesis Testing of H_0: beta = 0 ##
# Inputs:
## Periods: A vector of period numbers for each cluster
## Outcomes: A vector of cluster-level outcomes in the same order as Periods
## Clusters: A vector of cluster numbers/names in the same order as Periods
## StartTimes: A vector of the first period on intervention for each cluster,
###            repeated for each appearance of the cluster in Outcomes
## Indivs: If "MEM" or "CPI" is included in Type and the outcomes are binary, the number of
###        individuals represented by the cluster-level outcomes in Outcomes.
###        Same format as in MEM.Effect.Est and CPI.Effect.Est
## family: a GLM family, see family or glm. E.g., binomial or gaussian
## link: a GLM link, see family or glm. E.g., "logit" or "identity".
###      For methods other than MEM/CPI, the ContrastFunc will be imputed from family and link.
## NumPerms: the number of permutations to perform. Must be an integer >= 1
## Type: a vector of character strings specifying the methods to use.
###      Each method is specified by the name following "Est." in the outputs of above methods. E.g., "NPWP" or "SCSWT1"
## Estimate: a named vector or list of the observed (non-permuted) effect estimates,
###          with names in the form Est.[Type].
## Alternative: a character string specifying the type of alternative hypothesis:
###             "Both" to test H_a: beta != 0; "Greater" to test H_a: beta > 0; "Less" to test H_a: beta < 0
# Note: there is random permutation generation in this method (through calling Perm.Effects).
##      Set a seed before running for replicability.
# Outputs:
## A named vector with one element for each Type:
### PVal.[Type]: the p-value for the specified permutation-based hypothesis test
####             for the method [Type] using [NumPerms] permutations of the data
####             and the observed effect estimate [Estimate]

Null.Test <- function(Periods, Outcomes, Clusters, StartTimes, 
                      Indivs=NA, family, link, NumPerms, 
                      Type=c("MEM","CPI","NPWP",
                             "SC",
                             "CO.Ctrl","CO.Both","CO.CtWt","CO.BoWt","ENS"),
                      Estimate, Alternative="Both", Alpha=.05, 
                      FwdOffset=0, BwdOffset=0) {
  AllTypes <- Type
  if ("NPWP" %in% AllTypes) {
    AllTypes <- AllTypes[AllTypes != "NPWP"]
    AllTypes <- c(AllTypes,"NPWP.TW","NPWP.Mean","NPWP.MeanEx1","NPWP.OW","NPWP.OWEx1")
  }
  if ("SC" %in% AllTypes) {
    AllTypes <- AllTypes[AllTypes != "SC"]
    AllTypes <- unique(c(AllTypes,"SCSWT1","SCSWT2","SCSWT3","SCSWT1.Ex1","SCSWT2.Ex1","SCSWT3.Ex1"))
  }
  
  Perm.Ests <- Perm.Effects(Periods, Outcomes, Clusters, StartTimes, Indivs, 
                            family, link, NumPerms, Type, 
                            FwdOffset, BwdOffset)
  
  if (Alternative=="Both") {
    PValFunc <- function(Truth,Perms) mean(abs(Truth) <= abs(Perms), na.rm=TRUE)
  } else if (Alternative=="Greater") {
    PValFunc <- function(Truth,Perms) mean(Truth >= Perms, na.rm=TRUE)
  } else if (Alternative=="Less") {
    PValFunc <- function(Truth,Perms) mean(Perms <= Truth, na.rm=TRUE)
  } else {
    stop("Please enter 'Both', 'Greater', or 'Less' for the variable Alternative")
  }
  
  PVal <- NULL
  Res <- NULL
  if (length(link) == 1) {
    for (i in 1:length(AllTypes)) {
      outstr1 <- paste0("PVal.",AllTypes[i])
      Perms <- Perm.Ests[,paste0("PermEsts.",AllTypes[i])]
      Truth <- unlist(Estimate[paste0("Est.",AllTypes[i])])
      PVali <- PValFunc(Truth,Perms)
      namesprior1 <- names(PVal)
      PVal <- append(PVal, PVali)
      names(PVal) <- c(namesprior1, outstr1)
      Resi <- ifelse(PVali <= Alpha, "Reject", "Accept")
      Res <- append(Res, Resi)
    }
  } else {
    for (k in 1:length(link)) {
      Estimatek <- Estimate[[k]]
      Permsk <- Perm.Ests[[k]]
      
      outstr <- paste0("PVal.",AllTypes)
      instr <- paste0("PermEsts.",AllTypes)
      Perms <- Permsk[,instr]
      Truth <- Estimatek[1,instr]
      PValk <- mapply(function(x,y) {PValFunc(x,y)}, Truth, Perms)
      Resk <- ifelse(PValk <= Alpha, "Reject", "Accept")
      
      PVal <- append(PVal, list(PValk))
      Res <- append(Res, list(Resk))
    }
    names(PVal) <- names(Estimate)[1:length(link)]
    names(Res) <- names(Estimate)[1:length(link)]
  }
  return(list(PVal=PVal, Res=Res))
}


### Full Analysis with Estimation ###
### and Hypothesis Tests for Specified Null Values ###
### for Selected Methods ###

# Inputs:
## Periods: A vector of period numbers for each cluster
## Outcomes: A vector of cluster-level outcomes in the same order as Periods
## Clusters: A vector of cluster numbers/names in the same order as Periods
## StartTimes: A vector of the first period on intervention for each cluster,
###            repeated for each appearance of the cluster in Outcomes
## Treatments: A vector of treatment indicators (0 for control, 1 for intervention)
###            in the same order as Periods.
#### Note: can specify either StartTimes or Treatments
## Indivs: If "MEM" or "CPI" is included in Type and the outcomes are binary, the number of
###        individuals represented by the cluster-level outcomes in Outcomes.
###        Same format as in MEM.Effect.Est and CPI.Effect.Est
## family: a GLM family, see family or glm. E.g., binomial or gaussian
## link: a GLM link, see family or glm. E.g., "logit" or "identity".
###      For methods other than MEM/CPI, the ContrastFunc will be imputed from family and link.
## NumPerms: the number of permutations to perform for inference. Must be an integer >= 1
## Type: a vector of character strings specifying the methods to use.
###      Each method is specified by the name following "Est." in the outputs of above methods. E.g., "NPWP" or "SCSWT1"
## NullVals: a vector of null hypothesis values to test. By default, only tests 0 to get the p-value.
###          Can specify additional NullVals to test whether the 1-alpha CIs include those values.
## Alternative: a character string specifying the type of alternative hypothesis:
###             "Both" to test H_a: beta != 0; "Greater" to test H_a: beta > 0; "Less" to test H_a: beta < 0
## Alpha: the significance level at which to perform hypothesis tests and CI inclusion checks
# Note: there is random permutation generation in this method (through calling Perm.Effects).
##      Set a seed before running for replicability.
# Outputs:
## A data frame with one row for each component of NullVals. Each row has the elements:
### NullVal: the null hypothesis value being tested
### For each Type:
#### Est.[Type]: the treatment effect estimate for the [Type] method
#### PVal.[Type]: the p-value for H_0: beta = NullVal for the [Type] estimate
#### Res.[Type]: "Accept" if PVal.[Type] > alpha; otherwise, "Reject".
#####            "Accept" indicates failing to reject the null hypothesis in the test or 
#####             inclusion of NullVal in the (1-alpha) CI
#####            "Reject" indicates rejecting the null hypothesis in the test or 
#####             exclusion of NullVal from the (1-alpha) CI
#### AsyPVal.[Type]: only for MEM and CPI types; the same as PVal 
#####                but using asymptotic inference instead of permutation-based exact inference
#### AsyRes.[Type]: only for MEM and CPI types; the same as Res 
#####                but using asymptotic inference instead of permutation-based exact inference

SWT.Permutation.Analysis <- function(Periods, Outcomes, Clusters, 
                                     StartTimes=NA, Treatments=NA, 
                                     Indivs=NA,
                                     family, link, 
                                     NumPerms, 
                                     Type=c("MEM","CPI","NPWP",
                                            "SC",
                                            "CO.Ctrl","CO.Both","CO.CtWt","CO.BoWt",
                                            "ENS"), 
                                     NullVals=0,
                                     Alternative="Both", Alpha=.05, 
                                     FwdOffset=0, BwdOffset=0) {
  if (sum(is.na(StartTimes)) > 0) {
    if (sum(is.na(Treatments)) > 0) {
      stop("Must Specify Either StartTimes or Vector of Treatment Indicators")
    } else {
      StartTimes <- StartTimeVec(Periods=Periods, Clusters=Clusters, Trts=Treatments)
    }
  }
  
  # if (link=="logit") {
  #   ContrastFunc <- logORFunc
  #   ContrastInvApply <- logORInvApply
  # } else if (link=="log") {
  #   ContrastFunc <- logRRFunc
  #   ContrastInvApply <- logRRInvApply
  # } else if (link=="identity") {
  #   ContrastFunc <- RDFunc
  #   ContrastInvApply <- RDInvApply
  # } else {
  #   print("Warning: Contrast function unknown for this family. Using difference for non-parametric methods.")
  #   ContrastFunc <- RDFunc
  #   ContrastInvApply <- RDInvApply
  # }
  
  OutDF <- NULL
  for (m in 1:length(NullVals)) {
    # if (NullVals[m] == 0) {
    #   AdjOutcomes <- Outcomes
    # } else {
    #   AdjOutcomes <- ifelse(Periods >= StartTimes, pmax(0,ContrastInvApply(Outcomes,NullVals[m])), Outcomes)
    # }
    AdjOutcomes <- Outcomes
    
    Ests <- Perm.Effects(Periods, AdjOutcomes, Clusters, StartTimes, Indivs, family, link,
                         NumPerms=0, Type, FwdOffset, BwdOffset)
    
    if (NumPerms == 0) {
      return(Ests)
    } else {
      Test.Res <- Null.Test(Periods, AdjOutcomes, Clusters, StartTimes,
                            Indivs, family, link, NumPerms, 
                            Type,
                            Estimate=Ests, Alternative, Alpha,
                            FwdOffset, BwdOffset)
      return(append(Ests, list(PVal=Test.Res$PVal, Res=Test.Res$Res)))
    }
  }
}




