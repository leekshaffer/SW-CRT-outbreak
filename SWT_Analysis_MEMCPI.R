##########################################
##### SWT_Analysis_MEMCPI.R ##############
##### Perform MEM and CPI estimation and asymptotic inference on SWT results
##### Kennedy-Shaffer and Lipsitch #######
##### Update: May 26, 2020 ###############
##### Please Cite https://doi.org/10.1101/2020.05.01.20087429
##### Based on methods described in Kennedy-Shaffer et al. (http://doi.org/10.1002/sim.8451)
##### For questions, contact Lee Kennedy-Shaffer: lee_kennedyshaffer@g.harvard.edu
##### Approximate Runtime for SWT-A from paper: 7.5 minutes for one simulation for all 4 methods
##########################################

### Packages Needed: ###
# install.packages("Matrix","lme4")
require("Matrix")
require("lme4")

NumPds <- dim(fullRes_SWT$OutMat)[2]
NumClusts <- dim(fullRes_SWT$OutMat)[1]

MEMests <- rep(NA,4)
names(MEMests) <- c("MEM","CPI","MEM.Ex1","CPI.Ex1")
MEMpvals <- MEMests

Data <- fullRes_SWT$Data
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
for (i in 1:NumClusts) {
  comm <- unname(Data[i,"Community"])
  SP <- unname(Data[i,"StartPd"])
  enrolled <- unname(Data[i,"indivs"])
  for (j in 1:NumPds) {
    if (enrolled > 0) {
      newcases <- unname(Data[i,paste0("Pd_",j)])
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

capture.output(MEM <- tryCatch({glmer(Outs~Trt+PdFactor+(1|Clust),
                                      family=binomial(link="logit"))},
                               error=function(err) {1}), file="/dev/null")
if (!is.numeric(MEM)) {
  MEMests["MEM"] <- summary(MEM)$coefficients['Trt','Estimate']
  MEMpvals["MEM"] <- summary(MEM)$coefficients['Trt','Pr(>|z|)']
}
capture.output(CPI <- tryCatch({glmer(Outs~Trt+PdFactor+(1|Clust)+(1|CPIndex),
                                      family=binomial(link="logit"))},
                               error=function(err) {1}), file="/dev/null")
if (!is.numeric(CPI)) {
  MEMests["CPI"] <- summary(CPI)$coefficients['Trt','Estimate']
  MEMpvals["CPI"] <- summary(CPI)$coefficients['Trt','Pr(>|z|)']
}
capture.output(MEMEx1 <- tryCatch({glmer(OutsEx1~TrtEx1+PdFactorEx1+(1|ClustEx1),
                                         family=binomial(link="logit"))},
                                  error=function(err) {1}), file="/dev/null")
if (!is.numeric(MEMEx1)) {
  MEMests["MEM.Ex1"] <- summary(MEMEx1)$coefficients['TrtEx1','Estimate']
  MEMpvals["MEM.Ex1"] <- summary(MEMEx1)$coefficients['TrtEx1','Pr(>|z|)']
}
capture.output(CPIEx1 <- tryCatch({glmer(OutsEx1~TrtEx1+PdFactorEx1+(1|ClustEx1)+(1|CPIndexEx1),
                                         family=binomial(link="logit"))},
                                  error=function(err) {1}), file="/dev/null")
if (!is.numeric(CPIEx1)) {
  MEMests["CPI.Ex1"] <- summary(CPIEx1)$coefficients['TrtEx1','Estimate']
  MEMpvals["CPI.Ex1"] <- summary(CPIEx1)$coefficients['TrtEx1','Pr(>|z|)']
}


Analysis_MEMCPI_Res <- list(Est=MEMests, PVal=MEMpvals, VE=1-exp(MEMests))


