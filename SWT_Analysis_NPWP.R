##########################################
##### SWT_Analysis_NPWP.R ################
##### Perform NPWP estimation and permutation inference on SWT results
##### Kennedy-Shaffer and Lipsitch #######
##### Update: May 26, 2020 ###############
##### Please Cite https://doi.org/10.1101/2020.05.01.20087429
##### Based on methods described by Thompson et al. (http://doi.org/10.1002/sim.7668)
##### and used by Kennedy-Shaffer et al. (http://doi.org/10.1002/sim.8451)
##### For questions, contact Lee Kennedy-Shaffer: lee_kennedyshaffer@g.harvard.edu
##### Approximate Runtime for SWT-A from paper: 3 seconds for 500 permutations of one simulation
##########################################

NumPerms <- 500 ## The number of permutations to be used for hypothesis testing

SimpleContrasts <- function(Cont,Int,ContIndivs,IntIndivs) {
  RD <- Int - Cont
  RD <- ifelse(is.finite(RD),RD,NA)
  if (!is.finite(Cont) | !is.finite(Int)) {
    RR <- NA
  } else {
    if (Cont > 0 & Int > 0) {
      RR <- Int/Cont
    } else if (Cont==0| Int==0) {
      IntAdj <- (Int*IntIndivs+0.5)/(IntIndivs+1)
      ContAdj <- (Cont*ContIndivs+0.5)/(ContIndivs+1)
      RR <- IntAdj/ContAdj
    } else {
      RR <- NA
    }
  }
  lRR <- log(RR)
  return(c("RD"=RD,"RR"=RR,"lRR"=lRR))
}

NPPermConts <- function(Data,OutMat,Permute) {
  if (Permute==TRUE) {
    PermStarts <- sample(Data[,"StartPd"], size=length(Data[,"StartPd"]), replace=FALSE)
  } else {
    PermStarts <- Data[,"StartPd"]
  }
  Indivs <- Data[,"indivs"]
  cols <- colnames(OutMat)
  NumPds <- dim(OutMat)[2]
  NumClusts <- dim(OutMat)[1]
  PdMat <- matrix(rep(1:NumPds,NumClusts),nrow=NumClusts,ncol=NumPds,byrow=TRUE)
  colnames(PdMat) <- cols
  TypeMat <- apply(PdMat, 2, function(x) ifelse(x>PermStarts,2,
                                                ifelse(x==PermStarts,1,0)))
  
  ## NP Calcs:
  Avgs <- apply(OutMat*ifelse(TypeMat==0,1,NA), 2, mean, na.rm=TRUE)
  Avgs <- rbind(Avgs, apply(OutMat*ifelse(TypeMat>0,1,NA), 2,
                            mean, na.rm=TRUE))
  Avgs <- rbind(Avgs, apply(OutMat*ifelse(TypeMat==2,1,NA), 2,
                            mean, na.rm=TRUE))
  Avgs <- rbind(Avgs, apply(Data[,cols]*ifelse(TypeMat==0,1,NA), 2, 
                            sum, na.rm=TRUE))
  Avgs <- rbind(Avgs, apply(Data[,cols], 2, sum, na.rm=TRUE))
  Avgs <- rbind(Avgs, apply(Data[,cols]*ifelse(TypeMat==1,NA,1), 2, 
                            sum, na.rm=TRUE))
  Avgs <- rbind(Avgs, apply(TypeMat, 2, function(x) sum(x==0)))
  Avgs <- rbind(Avgs, apply(TypeMat, 2, function(x) sum(x>0)))
  Avgs <- rbind(Avgs, apply(TypeMat, 2, function(x) sum(x==2)))
  Avgs <- rbind(Avgs, apply(Indivs*ifelse(TypeMat==0,1,NA), 2, sum, na.rm=TRUE))
  Avgs <- rbind(Avgs, apply(Indivs*ifelse(TypeMat>0,1,NA), 2, sum, na.rm=TRUE))
  Avgs <- rbind(Avgs, apply(Indivs*ifelse(TypeMat==2,1,NA), 2, sum, na.rm=TRUE))
  rownames(Avgs) <- c("Cont","Int","IntEx1",
                      "ContCases","TotCases","TotCasesEx1",
                      "ContClusts","IntClusts","IntClustsEx1",
                      "ContIndivs","IntIndivs","IntIndivsEx1")
  
  Res <- rbind(Avgs, apply(Avgs, 2,
                           function(x) SimpleContrasts(x["Cont"],x["Int"],
                                                       x["ContIndivs"],x["IntIndivs"])))
  Res <- rbind(Res, apply(Avgs, 2,
                          function(x) SimpleContrasts(x["Cont"],x["IntEx1"],
                                                      x["ContIndivs"],x["IntIndivsEx1"])))
  rownames(Res) <- c(rownames(Avgs),c("NP.RD","NP.RR","NP.lRR",
                                      "NP.RD.Ex1","NP.RR.Ex1","NP.lRR.Ex1"))
  return(Res)
}

NPEsts <- function(ResOut) {
  NPEsts <- NULL
  for (conttype in c("RD","RR","lRR")) {
    for (ex1 in c("","Ex1")) {
      row0ex1 <- ifelse(ex1=="Ex1",".Ex1","")
      row0 <- paste0("NP.",conttype,row0ex1)
      innames <- names(NPEsts)
      NPEsts <- c(NPEsts,
               mean(ResOut[row0,], na.rm=TRUE),
               sum(ResOut[row0,]*ResOut["IntClusts",], 
                   na.rm=TRUE)/sum(ifelse(is.finite(ResOut[row0,]),1,0)*ResOut[paste0("IntClusts",ex1),], 
                                   na.rm=TRUE),
               sum(ResOut[row0,]*ResOut["ContCases",], 
                   na.rm=TRUE)/sum(ifelse(is.finite(ResOut[row0,]),1,0)*ResOut["ContCases",], 
                                   na.rm=TRUE),
               sum(ResOut[row0,]*ResOut["TotCases",], 
                   na.rm=TRUE)/sum(ifelse(is.finite(ResOut[row0,]),1,0)*ResOut[paste0("TotCases",ex1),], 
                                   na.rm=TRUE))
      names(NPEsts) <- c(innames, paste0(conttype,".NP.",c("MeanPd","MeanCP","CoCaWt","ToCaWt"),row0ex1))
    }
  }
  return(NPEsts)
}


RDNames <- paste0("RD.NP.",c("MeanPd","MeanCP","CoCaWt","ToCaWt"))
RRNames <- paste0("RR.NP.",c("MeanPd","MeanCP","CoCaWt","ToCaWt"))
lRRNames <- paste0("lRR.NP.",c("MeanPd","MeanCP","CoCaWt","ToCaWt"))
RDNames <- c(RDNames,paste0(RDNames,".Ex1"))
RRNames <- c(RRNames,paste0(RRNames,".Ex1"))
lRRNames <- c(lRRNames,paste0(lRRNames,".Ex1"))

EstByPeriod <- NPPermConts(fullRes_SWT$Data,fullRes_SWT$OutMat,FALSE)
TrueEst <- NPEsts(EstByPeriod)
VE <- TrueEst
VE[c(RRNames,lRRNames)] <- 1-exp(VE[c(RRNames,lRRNames)])
if (NumPerms > 0) {
  PermEsts <- replicate(n=NumPerms, expr=NPEsts(NPPermConts(fullRes_SWT$Data,fullRes_SWT$OutMat,TRUE)))
  PermCompRD <- apply(X=PermEsts[RDNames,], MARGIN=2,
                    FUN=function(x) ifelse(abs(x)>=abs(TrueEst[RDNames]),1,0))
  PermCompRR <- apply(X=PermEsts[RRNames,], MARGIN=2,
                      FUN=function(x) ifelse(abs(log(x))>=abs(log(TrueEst[RRNames])),1,0))
  PermComplRR <- apply(X=PermEsts[lRRNames,], MARGIN=2,
                      FUN=function(x) ifelse(abs(x)>=abs(TrueEst[lRRNames]),1,0))
  PermComp <- rbind(PermCompRD,PermCompRR,PermComplRR)
  PVal <- apply(X=PermComp, MARGIN=1, mean, na.rm=TRUE)
} else {
  PVal <- NA
}


Analysis_NPWP_Res <- list(Est=TrueEst, PVal=PVal, VE=VE)
Analysis_NPWP_ByPeriod <- EstByPeriod[c("NP.RD","NP.RR","NP.lRR",
                                        "NP.RD.Ex1","NP.RR.Ex1","NP.lRR.Ex1"),]
  