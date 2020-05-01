
basefolder <- "/home/lk143/TimeVarID/OB_Sims"
res2folder <- paste0(basefolder,"/FullRes2")
NPSCfolder <- paste0(basefolder,"/NPSCRes")
Outfolder <- paste0(basefolder,"/Compiled")

nums <- 1:500
NumPerms <- 500
scens <- as.numeric(substring(Sys.getenv('SLURM_JOB_NAME'), first=16, last=16))
params <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))


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
# package.install("Synth")
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

NPPermConts <- function(Data,OutMat) {
  PermStarts <- sample(Data[,"StartPd"], size=length(Data[,"StartPd"]), replace=FALSE)
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

for (scen in scens) {
  for (p in params) {
    param <- formatC(p, width=2, format="d", flag="0")
    name <- paste0("S",scen,"_P",param)
    errs <- NULL
    NPPVals <- NULL
    for (num in 1:500) {
      if (file.exists(paste0(res2folder,"/FullRes2_",name,"_",num,".Rda"))) {
        load(file=paste0(res2folder,"/FullRes2_",name,"_",num,".Rda"))
        TrueEsts <- NPEsts(fullRes$CombRes)
        PermEsts <- replicate(n=NumPerms, expr=NPEsts(NPPermConts(fullRes$Data,fullRes$OutMat)))
        PermCompRD <- apply(X=PermEsts[RDNames,], MARGIN=2,
                          FUN=function(x) ifelse(abs(x)>=abs(TrueEsts[RDNames]),1,0))
        PermCompRR <- apply(X=PermEsts[RRNames,], MARGIN=2,
                            FUN=function(x) ifelse(abs(log(x))>=abs(log(TrueEsts[RRNames])),1,0))
        PermComplRR <- apply(X=PermEsts[lRRNames,], MARGIN=2,
                            FUN=function(x) ifelse(abs(x)>=abs(TrueEsts[lRRNames]),1,0))
        PermComp <- rbind(PermCompRD,PermCompRR,PermComplRR)
        PVals <- c(Sim=num,apply(X=PermComp, MARGIN=1, mean, na.rm=TRUE))
        NPPVals <- rbind(NPPVals, PVals)
      } else {
        errs <- c(errs, num)
      }
    }
    save(NPPVals, file=paste0(Outfolder,"/PVals_NP_",name,".Rda"))
    print(paste0("Errors for Simulation: ",name))
    print(paste(errs, collapse=","))
  }
}