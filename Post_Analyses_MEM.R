
# basefolder <- "/Users/leekennedy-shaffer/Documents/Harvard/Research/ID Time-Varying Estimand/Code/Outbreak_Sims/"
basefolder <- "/home/lk143/TimeVarID/OB_Sims"
infolder <- paste0(basefolder,"/FullComp")
res2folder <- paste0(basefolder,"/FullRes2")
MEMfolder <- paste0(basefolder,"/MEMRes")
load(paste0(basefolder,"/Varying.Rda"))

scen <- as.numeric(substring(Sys.getenv('SLURM_JOB_NAME'), first=10, last=10))
param <- substring(Sys.getenv('SLURM_JOB_NAME'), first=13, last=14)
simno <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
PerArray <- 50
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
package.install("Matrix")
package.install("lme4")

PdLen <- Varying$PdLength[VaryingRow]
NumPds <- 28/PdLen+40/Varying$XOperPd[VaryingRow]
NumClusts <- 40
name <- paste0("S",scen,"_P",param)
name2 <- paste0(name,"_",simno)
load(paste0(infolder,"/FullComp_",name,".Rda"))
MEMests <- NULL
MEMpvals <- NULL
st <- proc.time()
for (num in nums) {
  OutMat <- out2$FullArray[,,num]
  Outs <- NULL
  Clust <- NULL
  Pd <- NULL
  Trt <- NULL
  OutsEx1 <- NULL
  ClustEx1 <- NULL
  PdEx1 <- NULL
  TrtEx1 <- NULL
  CPIndex <- NULL
  CPIndexEx1 <- NULL
  index0 <- 0
  index1 <- 0
  for (i in 1:(dim(OutMat)[1])) {
    comm <- unname(OutMat[i,"Community"])
    SP <- unname(OutMat[i,"StartPd"])
    enrolled <- unname(OutMat[i,"indivs"])
    for (j in 1:NumPds) {
      if (enrolled > 0) {
        newcases <- unname(OutMat[i,paste0("Pd_",j)])
        index0 <- index0 + 1
        Clust <- c(Clust, rep(comm,enrolled))
        Pd <- c(Pd, rep(j,enrolled))
        Trt <- c(Trt, rep(ifelse(j>=SP,1,0),enrolled))
        Outs <- c(Outs, c(rep(1,newcases),rep(0,enrolled-newcases)))
        CPIndex <- c(CPIndex, rep(index0,enrolled))
        if (j != SP) {
          index1 <- index1 + 1
          ClustEx1 <- c(ClustEx1, rep(comm,enrolled))
          PdEx1 <- c(PdEx1, rep(j,enrolled))
          TrtEx1 <- c(TrtEx1, rep(ifelse(j>SP,1,0),enrolled))
          OutsEx1 <- c(OutsEx1, c(rep(1,newcases),rep(0,enrolled-newcases)))
          CPIndexEx1 <- c(CPIndexEx1, rep(index1,enrolled))
        }
        enrolled <- enrolled - newcases
      }
    }
  }
  PdFactor <- relevel(factor(Pd), ref=1)
  PdFactorEx1 <- relevel(factor(PdEx1), ref=1)
  MEMests.num <- c(num,rep(NA,4))
  names(MEMests.num) <- c("Sim","MEM","CPI","MEM.Ex1","CPI.Ex1")
  MEMpvals.num <- MEMests.num
  st <- proc.time()
  capture.output(MEM <- tryCatch({glmer(Outs~Trt+PdFactor+(1|Clust),
                                        family=binomial(link="logit"))},
                                 error=function(err) {1}), file="/dev/null")
  if (!is.numeric(MEM)) {
    MEMests.num["MEM"] <- summary(MEM)$coefficients['Trt','Estimate']
    MEMpvals.num["MEM"] <- summary(MEM)$coefficients['Trt','Pr(>|z|)']
  }
  capture.output(CPI <- tryCatch({glmer(Outs~Trt+PdFactor+(1|Clust)+(1|CPIndex),
                                        family=binomial(link="logit"))},
                                 error=function(err) {1}), file="/dev/null")
  if (!is.numeric(CPI)) {
    MEMests.num["CPI"] <- summary(CPI)$coefficients['Trt','Estimate']
    MEMpvals.num["CPI"] <- summary(CPI)$coefficients['Trt','Pr(>|z|)']
  }
  capture.output(MEMEx1 <- tryCatch({glmer(OutsEx1~TrtEx1+PdFactorEx1+(1|ClustEx1),
                                           family=binomial(link="logit"))},
                                    error=function(err) {1}), file="/dev/null")
  if (!is.numeric(MEMEx1)) {
    MEMests.num["MEM.Ex1"] <- summary(MEMEx1)$coefficients['TrtEx1','Estimate']
    MEMpvals.num["MEM.Ex1"] <- summary(MEMEx1)$coefficients['TrtEx1','Pr(>|z|)']
  }
  capture.output(CPIEx1 <- tryCatch({glmer(OutsEx1~TrtEx1+PdFactorEx1+(1|ClustEx1)+(1|CPIndexEx1),
                                           family=binomial(link="logit"))},
                                    error=function(err) {1}), file="/dev/null")
  if (!is.numeric(CPIEx1)) {
    MEMests.num["CPI.Ex1"] <- summary(CPIEx1)$coefficients['TrtEx1','Estimate']
    MEMpvals.num["CPI.Ex1"] <- summary(CPIEx1)$coefficients['TrtEx1','Pr(>|z|)']
  }
  MEMests <- rbind(MEMests,MEMests.num)
  MEMpvals <- rbind(MEMpvals,MEMpvals.num)
  FullRes <- list(MEM=MEM,CPI=CPI,MEM.Ex1=MEMEx1,CPI.Ex1=CPIEx1)
  save(FullRes, file=paste0(MEMfolder,"/FullOuts/FullMEM_",name,"_",num,".Rda"))
  rm(FullRes,MEM,CPI,MEMEx1,CPIEx1,MEMests.num,MEMpvals.num)
}
proc.time() - st
MEMRes <- list(MEM.Ests=MEMests,MEM.Pvals=MEMpvals)
save(MEMRes, file=paste0(MEMfolder,"/MEMRes_",name2,".Rda"))


