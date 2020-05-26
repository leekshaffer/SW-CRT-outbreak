##########################################
##### SWT_Analysis_SC.R ##################
##### Perform Synthetic Control estimation and permutation inference on SWT results
##### Kennedy-Shaffer and Lipsitch #######
##### Update: May 26, 2020 ###############
##### Please Cite https://doi.org/10.1101/2020.05.01.20087429
##### Based on methods described by Kennedy-Shaffer et al. (http://doi.org/10.1002/sim.8451)
##### For questions, contact Lee Kennedy-Shaffer: lee_kennedyshaffer@g.harvard.edu
##### Approximate Runtime for SWT-A from paper: 1.5 minutes per permutation for one simulation
##########################################

NumPerms <- 500  ## The number of permutations to be used for hypothesis testing

### Packages Needed: ###
# install.packages("Synth")
require(Synth)

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

SCPermConts <- function(Data,OutMat,Permute) {
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
  
  ## SC Calcs:
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
  
  ContrastPds <- (1:NumPds)[apply(TypeMat, 2, function(x) mean(x>0)>0 & mean(x==0)>0)]
  SCOutMat <- matrix(NA, nrow=dim(TypeMat)[1], ncol=dim(TypeMat)[2])
  colnames(SCOutMat) <- cols
  MSPEMat <- SCOutMat
  for (i in 1:length(ContrastPds)) {
    Pd <- ContrastPds[i]
    Types <- TypeMat[,paste0("Pd_",Pd)]
    Ints <- which(Types>0)
    NumInts <- length(Ints)
    Conts <- which(Types==0)
    NumConts <- length(Conts)
    Outs <- OutMat[,paste0("Pd_",Pd)]
    if (Pd==1) {
      SCOutMat[Ints,"Pd_1"] <- rep(mean(Outs[Conts], na.rm=TRUE), NumInts)
    } else {
      Mati <- OutMat[,paste0("Pd_",(1:(Pd-1))),drop=FALSE]
      Typei <- TypeMat[,paste0("Pd_",(1:(Pd-1))),drop=FALSE]
      for (k in 1:NumInts) {
        j <- Ints[k]
        LastContj <- max(which(Typei[j,]==0))
        Outsj <- Mati[j,1:LastContj,drop=FALSE]
        Repeat <- 0
        if (k > 1) {
          for (ell in 1:(k-1)) {
            LastContl <- max(which(Typei[Ints[ell],]==0))
            if (LastContl==LastContj) {
              Outsl <- Mati[Ints[ell],1:LastContj,drop=FALSE]
              if (!is.na(sum(Outsj==Outsl)) & sum(Outsj==Outsl)==length(Outsj)) {
                Repeat <- 1
                SCOutMat[j,paste0("Pd_",Pd)] <- SCOutMat[Ints[ell],paste0("Pd_",Pd)]
                MSPEMat[j,paste0("Pd_",Pd)] <- MSPEMat[Ints[ell],paste0("Pd_",Pd)]
              }
            }
          }
        }
        if (Repeat==0) {
          Contsj <- Mati[Conts,1:LastContj,drop=FALSE]
          UsePreds <- apply(Contsj,2,sd)>0
          ContsjUse <- Contsj[,UsePreds,drop=FALSE]
          OutsjUse <- Outsj[,UsePreds,drop=FALSE]
          if (NumConts==1) {
            SCOutMat[j,paste0("Pd_",Pd)] <- Outs[Conts]
            MSPEMat[j,paste0("Pd_",Pd)] <- sum((Outsj-Contsj)^2)
          } else {
            capture.output(SynthOut <- tryCatch({synth(X1=t(OutsjUse), X0=t(ContsjUse), 
                                                       Z1=t(OutsjUse), Z0=t(ContsjUse))},
                                                error=function(err) {1}), file="/dev/null")
            if (is.numeric(SynthOut)) { ## if synth does not converge, use simple average
              SCOutMat[j,paste0("Pd_",Pd)] <- mean(Outs[Types==0], na.rm=TRUE)
              MSPEMat[j,paste0("Pd_",Pd)] <- mean((Outsj - apply(Contsj,2,mean,na.rm=TRUE))^2)
            } else { ## if synth does converge, use SC fit
              SCOutMat[j,paste0("Pd_",Pd)] <- SynthOut$solution.w[,1] %*% Outs[Types==0]
              MSPEMat[j,paste0("Pd_",Pd)] <- as.numeric(SynthOut$loss.v)
            }
          }
        }
      }
    }
  }
  
  MSPEMat <- ifelse(MSPEMat<10^(-10),min(MSPEMat[MSPEMat>=10^(-10)], na.rm=TRUE),MSPEMat)
  MSPEwtMat <- MSPEMat
  for (c in 1:dim(MSPEwtMat)[2]) {
    sums <- aggregate(MSPEMat[,c], by=list(PermStarts),
                      FUN=function(x) sum(1/x, na.rm=TRUE))
    for (r in 1:dim(MSPEwtMat)[1]) {
      MSPEwtMat[r,c] <- (1/MSPEMat[r,c])/sums$x[sums$Group.1==PermStarts[r]]
    }
  }
  WtSums <- apply(MSPEwtMat, 2, sum, na.rm=TRUE)
  WtSumsEx1 <- apply(MSPEwtMat*ifelse(TypeMat==2,1,NA), 2, sum, na.rm=TRUE)
  MSPEwtMat2 <- t(apply(MSPEwtMat, 1, function(x) x/WtSums))
  MSPEwtMatEx1 <- t(apply(MSPEwtMat*ifelse(TypeMat==2,1,NA), 1, function(x) x/WtSumsEx1))
  MSPEwtMatEx1 <- ifelse(is.finite(MSPEwtMatEx1),MSPEwtMatEx1,NA)
  
  SCAvgs <- apply(SCOutMat, 2, mean, na.rm=TRUE)
  SCAvgs <- rbind(SCAvgs, apply(SCOutMat*ifelse(TypeMat==2,1,NA), 
                                2, mean, na.rm=TRUE))
  SCAvgs <- rbind(SCAvgs, apply(SCOutMat*MSPEwtMat2, 2, sum, na.rm=TRUE))
  SCAvgs <- rbind(SCAvgs, apply(Outs*MSPEwtMat2, 2, sum, na.rm=TRUE))
  SCAvgs <- rbind(SCAvgs, apply(SCOutMat*MSPEwtMatEx1, 2, sum, na.rm=TRUE))
  SCAvgs <- rbind(SCAvgs, apply(Outs*MSPEwtMatEx1, 2, sum, na.rm=TRUE))
  rownames(SCAvgs) <- c("SCCont","SCContEx1","SCContWt","SCIntWt","SCContWtEx1","SCIntWtEx1")
  SCAvgs[,apply(SCOutMat, 2, function(x) sum(is.finite(x))==0)] <- NA
  SCAvgs[c("SCContEx1","SCContWtEx1","SCIntWtEx1"),
         apply(SCOutMat*ifelse(TypeMat==2,1,NA), 2, function(x) sum(is.finite(x))==0)] <- NA
  CombAvgs <- rbind(Avgs, SCAvgs)
  
  SCRes <- apply(CombAvgs, 2, 
                 function(x) SimpleContrasts(x["SCCont"],x["Int"],
                                             x["ContIndivs"],x["IntIndivs"]))
  SCRes <- rbind(SCRes, apply(CombAvgs, 2,
                              function(x) SimpleContrasts(x["SCContEx1"],x["IntEx1"],
                                                          x["ContIndivs"],x["IntIndivsEx1"])))
  SCRes <- rbind(SCRes, apply(CombAvgs, 2,
                              function(x) SimpleContrasts(x["SCContWt"],x["SCIntWt"],
                                                          x["ContIndivs"],x["IntIndivs"])))
  SCRes <- rbind(SCRes, apply(CombAvgs, 2,
                              function(x) SimpleContrasts(x["SCContWtEx1"],x["SCIntWtEx1"],
                                                          x["ContIndivs"],x["IntIndivsEx1"])))
  rownames(SCRes) <- c("SC.RD","SC.RR","SC.lRR","SC.RD.Ex1","SC.RR.Ex1","SC.lRR.Ex1",
                       "SCWt.RD","SCWt.RR","SCWt.lRR","SCWt.RD.Ex1","SCWt.RR.Ex1","SCWt.lRR.Ex1")
  return(rbind(CombAvgs,SCRes))
}

SCEsts <- function(ResOut) {
  SCEsts <- NULL
  for (antype in c("SC","SCWt")) {
    for (conttype in c("RD","RR","lRR")) {
      for (ex1 in c("","Ex1")) {
        row0ex1 <- ifelse(ex1=="Ex1",".Ex1","")
        row0 <- paste0(antype,".",conttype,row0ex1)
        innames <- names(SCEsts)
        SCEsts <- c(SCEsts,
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
        names(SCEsts) <- c(innames, paste0(conttype,".",antype,".",c("MeanPd","MeanCP","CoCaWt","ToCaWt"),row0ex1))
      }
    }
  }
  return(SCEsts)
}


RDNames <- paste0("RD.",rep(c("SC","SCWt"),each=4),".",rep(c("MeanPd","MeanCP","CoCaWt","ToCaWt"),2))
RRNames <- paste0("RR.",rep(c("SC","SCWt"),each=4),".",rep(c("MeanPd","MeanCP","CoCaWt","ToCaWt"),2))
lRRNames <- paste0("lRR.",rep(c("SC","SCWt"),each=4),".",rep(c("MeanPd","MeanCP","CoCaWt","ToCaWt"),2))
RDNames <- c(RDNames,paste0(RDNames,".Ex1"))
RRNames <- c(RRNames,paste0(RRNames,".Ex1"))
lRRNames <- c(lRRNames,paste0(lRRNames,".Ex1"))

EstByPeriod <- SCPermConts(fullRes_SWT$Data,fullRes_SWT$OutMat,FALSE)
TrueEst <- SCEsts(EstByPeriod)
VE <- TrueEst
VE[c(RRNames,lRRNames)] <- 1-exp(VE[c(RRNames,lRRNames)])
if (NumPerms > 0) {
  PermEsts <- NULL
  for (n in 1:NumPerms) {
    print(paste0("Permutation Number: ",n))
    PermEst <- SCEsts(SCPermConts(fullRes_SWT$Data,fullRes_SWT$OutMat,TRUE))
    PermEsts <- rbind(PermEsts,c(PermNo=n,PermEst))
  }
  PermCompRD <- t(apply(X=PermEsts[,RDNames], MARGIN=1,
                    FUN=function(x) ifelse(abs(x)>=abs(TrueEst[RDNames]),1,0)))
  PermCompRR <- t(apply(X=PermEsts[,RRNames], MARGIN=1,
                      FUN=function(x) ifelse(abs(log(x))>=abs(log(TrueEst[RRNames])),1,0)))
  PermComplRR <- t(apply(X=PermEsts[,lRRNames], MARGIN=1,
                      FUN=function(x) ifelse(abs(x)>=abs(TrueEst[lRRNames]),1,0)))
  PermComp <- cbind(PermCompRD,PermCompRR,PermComplRR)
  PVal <- apply(X=PermComp, MARGIN=2, mean, na.rm=TRUE)
} else {
  PVal <- NA
}
  
  
Analysis_SC_Res <- list(Est=TrueEst, PVal=PVal, VE=VE)
Analysis_SC_ByPeriod <- EstByPeriod[c("SC.RD","SC.RR","SC.lRR","SC.RD.Ex1","SC.RR.Ex1","SC.lRR.Ex1",
                                      "SCWt.RD","SCWt.RR","SCWt.lRR","SCWt.RD.Ex1","SCWt.RR.Ex1","SCWt.lRR.Ex1"),]
