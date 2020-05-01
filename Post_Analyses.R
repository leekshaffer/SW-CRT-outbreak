
# basefolder <- "/Users/leekennedy-shaffer/Documents/Harvard/Research/ID Time-Varying Estimand/Code/Outbreak_Sims/"
basefolder <- "/home/lk143/TimeVarID/OB_Sims"
infolder <- paste0(basefolder,"/FullComp")
res2folder <- paste0(basefolder,"/FullRes2")
NPSCfolder <- paste0(basefolder,"/NPSCRes")
Compfolder <- paste0(basefolder,"/Compiled")
load(paste0(basefolder,"/Varying.Rda"))

scen <- as.numeric(substring(Sys.getenv('SLURM_JOB_NAME'), first=10, last=10))
param <- substring(Sys.getenv('SLURM_JOB_NAME'), first=13, last=14)
simno <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
PerArray <- 250
nums <- (PerArray*(simno-1)+1):(PerArray*simno)
VaryingRow <- 12*(scen-1)+as.numeric(param)

######################################################
###  function to easily use R packages on cluster  ###
######################################################
package.install = function(pack) {
  local({r <- getOption("repos");r["CRAN"] <- "http://cran.r-project.org"; options(repos=r)})
  
  # name of package to install / load
  pack = pack
  
  if (pack %in% rownames(installed.packages())) {
    library(pack, character.only=T)
  } else {
    if (pack %in% rownames(installed.packages(lib.loc='/home/lk143/apps/R_3.5.1/library'))) {
      library(pack, lib.loc='/home/lk143/apps/R_3.5.1/library', character.only=T)
    } else {
      install.packages(pack, lib='/home/lk143/apps/R_3.5.1/library')
      library(pack, lib.loc='/home/lk143/apps/R_3.5.1/library', character.only=T)
    }
  }
}


### Packages Needed: ###
package.install("Synth")
package.install("abind")

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

NPSCContrasts <- function(mat,NumPds,NumClusts) {
  cols <- paste0("Pd_",(1:NumPds))
  PdMat <- matrix(rep(1:NumPds,NumClusts),nrow=NumClusts,ncol=NumPds,byrow=TRUE)
  colnames(PdMat) <- cols
  TypeMat <- apply(PdMat, 2, function(x) ifelse(x>mat[,"StartPd"],2,
                                                ifelse(x==mat[,"StartPd"],1,0)))
  Indivs <- mat[,"indivs"]
  OutMat <- apply(mat[,cols], 2, function(x) x/Indivs)
  
  ## NP Calcs:
  Avgs <- apply(OutMat*ifelse(TypeMat==0,1,NA), 2, mean, na.rm=TRUE)
  Avgs <- rbind(Avgs, apply(OutMat*ifelse(TypeMat>0,1,NA), 2,
                            mean, na.rm=TRUE))
  Avgs <- rbind(Avgs, apply(OutMat*ifelse(TypeMat==2,1,NA), 2,
                            mean, na.rm=TRUE))
  Avgs <- rbind(Avgs, apply(mat[,cols]*ifelse(TypeMat==0,1,NA), 2, 
                            sum, na.rm=TRUE))
  Avgs <- rbind(Avgs, apply(mat[,cols], 2, sum, na.rm=TRUE))
  Avgs <- rbind(Avgs, apply(mat[,cols]*ifelse(TypeMat==1,NA,1), 2, 
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
  
  ## SC Calcs:
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
              if (sum(Outsj==Outsl)==length(Outsj)) {
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
    sums <- aggregate(MSPEMat[,c], by=list(mat[,"StartPd"]),
                      FUN=function(x) sum(1/x, na.rm=TRUE))
    for (r in 1:dim(MSPEwtMat)[1]) {
      MSPEwtMat[r,c] <- (1/MSPEMat[r,c])/sums$x[sums$Group.1==mat[r,"StartPd"]]
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
  CombAvgs <- rbind(Res, SCAvgs)

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
  
  CombRes <- rbind(CombAvgs, SCRes)
  
  fullRes <- list(Data=mat,TypeMat=TypeMat,OutMat=OutMat,
                  SCOutMat=SCOutMat,MSPEMat=MSPEMat,
                  CombRes=CombRes)
  
  return(list(CombRes=CombRes,fullRes=fullRes))
}

st <- proc.time()
# for (scen in c(2)) {
# for (scen in 1:2) {
  PdLen <- Varying$PdLength[VaryingRow]
  NumPds <- 28/PdLen+40/Varying$XOperPd[VaryingRow]
  NumClusts <- 40
  # for (param in c("01")) {
  # for (param in c("01","02","03","04","05","06","07","08","09","10","11","12")) {
    name <- paste0("S",scen,"_P",param)
    name2 <- paste0(name,"_",simno)
    load(paste0(infolder,"/FullComp_",name,".Rda"))
    ResArray <- NULL
    TE.RD <- NULL
    TE.RR <- NULL
    TE.lRR <- NULL
    for (num in nums) {
      print(paste0("Starting Sim Number ",num))
      if (file.exists(paste0(res2folder,"/FullRes2_",name,"_",num,".Rda"))) {
        load(paste0(res2folder,"/FullRes2_",name,"_",num,".Rda"))
      } else {
        ResOut0 <- NPSCContrasts(out2$FullArray[,,num], NumPds, NumClusts)
        fullRes <- ResOut0$fullRes
        save(fullRes,
             file=paste0(res2folder,"/FullRes2_",name,"_",num,".Rda"))
      }
      ResOut <- fullRes$CombRes
      ResArray <- abind(ResArray, ResOut, along=3)
      for (conttype in c("RD","RR","lRR")) {
        matname <- paste0("TE.",conttype,".num")
        assign(matname, c("Sim"=num))
        for (antype in c("NP","SC","SCWt")) {
          for (ex1 in c("","Ex1")) {
            row0ex1 <- ifelse(ex1=="Ex1",".Ex1","")
            row0 <- paste0(antype,".",conttype,row0ex1)
            innames <- names(get(matname))
            dat <- c(get(matname),
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
            names(dat) <- c(innames, paste0(antype,".",c("MeanPd","MeanCP","CoCaWt","ToCaWt"),row0ex1))
            assign(matname, dat)
          }
        }
        assign(paste0("TE.",conttype), rbind(get(paste0("TE.",conttype)), get(matname)))
      }
    }
    save(ResArray, file=paste0(NPSCfolder,"/ResArray_",name2,".Rda"))
    save(TE.RD, file=paste0(NPSCfolder,"/TERD_",name2,".Rda"))
    save(TE.RR, file=paste0(NPSCfolder,"/TERR_",name2,".Rda"))
    save(TE.lRR, file=paste0(NPSCfolder,"/TElRR_",name2,".Rda"))
    # fullResComp <- ResArray
    # save(fullResComp, file=paste0(Compfolder,"/FullRes2Comp_",name,".Rda"))
    # save(TE.RD, file=paste0(Compfolder,"/TERD_",name,".Rda"))
    # save(TE.RR, file=paste0(Compfolder,"/TERR_",name,".Rda"))
    # save(TE.lRR, file=paste0(Compfolder,"/TElRR_",name,".Rda"))
  # }
# }
proc.time() - st

