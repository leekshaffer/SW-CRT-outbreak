# SW-CRT-outbreak
Methods and Simulations for Stepped Wedge Cluster-Randomized Trials in Infectious Disease Outbreaks
Accompanies:
Kennedy-Shaffer L and Lipsitch M. Statistical properties of stepped wedge cluster-randomized trials in infectious disease outbreaks. <i>American Journal of Epidemiology<i> 2020; 189: 1324-1332. DOI:10.1093/aje/kwaa141 https://doi.org/10.1093/aje/kwaa141.

Last Commit: May 26, 2020

Outbreak_Simulations.R runs the simulated outbreaks with three trial types: IRT, CRT, and SWT. As written, it generates one network graph and simulates an outbreak and trial for each of the three types. Analyses for IRT and CRT are performed. Multiple simulations can be run within this file, but then the user must add code to save fullRes_SWT from each simulation for later or analysis or add analysis code within the file, e.g., by sourcing the SWT_Analysis files.

From fullRes_SWT, the following files can be run to analyze the SWT:
1. SWT_Analysis_PH.R (requires several parameters that are set in Outbreak_Simulations.R)
2. SWT_Analysis_MEMCPI.R ***Note: this can be slow to run with large numbers of clusters, periods, or enrollees per cluster
3. SWT_Analysis_NPWP.R
4. SWT_Analysis_SC.R ***Note: this is very slow to run with a high number of permutations
Permuation Inference can be run in PH, NPWP, and SC. Default is 500 permutations, but this can be changed with NumPerms. If NumPerms is set to 0, estimation will proceed by p-values will be NA.

Figure_Generation.R demonstrates how to construct figures like those shown in the manuscript. Each figure requires a specific combination of results of the simulations + analyses. A description of the required data frame/matrix/array is given within this file.

Parameters.Rda gives the parameters for the simulation settings used in the manuscript. Each is run for 1,000 simulations.

These programs are a work in progress, as we work to improve usability, error-catching, and speed of analysis. If you find errors, please contact Lee Kennedy-Shaffer at lee_kennedyshaffer (at) g (dot) harvard (dot) edu. 


