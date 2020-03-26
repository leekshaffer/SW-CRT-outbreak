# SW-CRT compared to IRT & CRT in Outbreak Setting

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
package.install("Rlab")
package.install("igraph")
package.install("deSolve")
package.install("reshape2")
package.install("ggplot2")
package.install("caTools")
package.install("devtools")
package.install("survival")
package.install("coxme")
package.install("frailtypack")
package.install("lme4")
package.install("NetSurv")


#devtools::install_github("jdwilson4/NetSurv")
#library(NetSurv, lib.loc='/home/lk143/apps/R_3.5.1/library')

# Set Simulation Parameters
basefolder <- "/home/lk143/TimeVarID/OB_Sims"
source(paste0(basefolder,"/SW-CRT Analysis Methods.R"))
load(paste0(basefolder,"/Varying.Rda"))
scen <- as.numeric(substring(Sys.getenv('SLURM_JOB_NAME'), first=10, last=10))
param <- as.numeric(substring(Sys.getenv('SLURM_JOB_NAME'), first=13, last=14))
simno <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
VaryingRow <- 12*(scen-1)+param

outname <- paste0(substring(Sys.getenv('SLURM_JOB_NAME'), first=9, last=14),"_",simno)

# Set Folders
scenfolder <- basefolder
resfolder <- paste0(scenfolder,"/FullRes")

# Number of simulated trials
nsim<-20

# Number of Permutations for Permutation Tests. If NumP==0, gives point esimates & asy. inference only for SW-CRTs
NumP <- 0

# Population structure parameters:
# Average size of one community
ave_community_size<-100
# Range of community sizes (sizes are uniformly distributed on this range)
community_size_range<-40
# Number of communities
num_communities<-40
# Probability of an edge between two nodes in the same community
rate_within<-0.15
# Probability of an edge between two nodes in different communities
rate_between<-0

# Disease characteristics:
# Per-time-step hazard of infection for a susceptible nodes from an infectious
# neighbour
beta<-Varying$Beta[VaryingRow]
# Expected number of importations to the population over two years
num_introductions<-80
# Leaky multiplicative efficacy of vaccine
direct_VE<-Varying$VE[VaryingRow]
# Gamma-distribution parameters of incubation and infectious period
incperiod_shape<-5.807 ## Parameters from Lauer et al. 2020
incperiod_rate<-1/0.948
infperiod_shape<-1.13
infperiod_rate<-0.226
ave_inc_period <- ceiling(incperiod_shape/incperiod_rate)
FixInc <- Varying$FixInc[VaryingRow]
# First day of trial enrollment, relative to start of epidemic
trial_startday<-56
# Days of follow-up
trial_length<-308

# General CRT Parameters:
# Target community enrollment proportion
cluster_coverage<-0.5

# Parallel-Arm CRT Parameters:
# Number of clusters targeted for enrollment
# Must be less than or equal to the number of communities
num_clusters_enrolled_per_day<-40
if (num_clusters_enrolled_per_day > num_communities) {stop("Enrolling too many communities!")}
# Number of days over which subjects are enrolled
enrollment_period<-1

# SW-CRT Parameters:
## Number of days between crossover events
step_interval <- Varying$PdLength[VaryingRow]
## First crossover event: this should be equal to trial_startday or a multiple of step_interval after
first_crossover <- 84
## Number of clusters enrolled per crossover event
num_clusters_per_step <- Varying$XOperPd[VaryingRow]
## Number of crossover events
num_steps <- num_communities/num_clusters_per_step
if (num_steps*num_clusters_per_step > num_clusters_enrolled_per_day*enrollment_period) {stop("Crossing over too many communities!")}
if (first_crossover + (num_steps-1)*step_interval >  trial_startday + trial_length) {stop("Final enrollment is after trial ends!")}


# Calculate R0
R0 <- (1 - (infperiod_rate/(infperiod_rate+beta))^infperiod_shape) *
  (((ave_community_size-1)*(1-rate_within)*rate_within + 
      (num_communities-1)*ave_community_size*(1-rate_between)*rate_between + 
      ((ave_community_size-1)*rate_within + (num_communities-1)*ave_community_size*rate_between)^2)/
     ((ave_community_size-1)*rate_within + (num_communities-1)*ave_community_size*rate_between)-1)
# cat("R0: ",R0,"\n")
# cat("Average number of individuals enrolled: ",num_clusters_enrolled_per_day*enrollment_period*cluster_coverage*ave_community_size)

list <- structure(NA,class="result")
"[<-.result" <- function(x,...,value) {
  args <- as.list(match.call())
  args <- args[-c(1:2,length(args))]
  length(value) <- length(args)
  for(i in seq(along=args)) {
    a <- args[[i]]
    if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
  }
  x
}

make_network <- function(ave_community_size, community_size_range, 
                         num_communities, rate_within, rate_between) {
  # Function to make a network of closely-connected communities that are more sparsely
  # connected with each other. I use a stochastic block model.
  # Inputs:
  # Ave_community_size is the average number of people in a community
  # Community_size_range is the range of community sizes. Currently size of a community
  # is being chosen from a uniform distribution on (ave-range/2, ave+range/2)
  # Num_communities is the number of communities in the study population
  # rate_within is the probability of an edge between any two nodes within the same community
  # rate_between is the probability of an edge between any two nodes in different communities
  
  # Create the network, and assign all members a community number
  # The community number will be output by the epidemic function and used in the 
  # gamma frailty model and calculation of the ICC
  community_sizes <- ave_community_size + round(runif(num_communities,-community_size_range/2,community_size_range/2))
  studypop_size <- sum(community_sizes)
  # All communities have the same connectedness, and all communities are equally
  # connected to each other
  within_rates <- diag(nrow=num_communities,ncol=num_communities,x=rate_within)
  between_rates <- matrix(rate_between,nrow=num_communities,ncol=num_communities) -
    diag(nrow=num_communities,ncol=num_communities,x=rate_between)
  rates<-within_rates+between_rates
  g <- sample_sbm(studypop_size,rates,community_sizes)
  # Give the nodes a name so that igraph remembers them
  V(g)$name<-1:studypop_size
  V(g)$community<-rep(1:num_communities,community_sizes)
  # Trial status will track whether a node is not in the trial (NA), in the control arm (0) or
  # in the vaccine arm (1)
  V(g)$trialstatus<-NA
  V(g)$enrollmentday<-NA
  
  return(g)
  
}

network_epidemic<-function(g,beta,num_introductions,VE,
                           incperiod_shape,incperiod_rate,infperiod_shape,infperiod_rate,
                           bTrial,bCluster,
                           trial_startday,trial_length,
                           num_enrolled_per_day,enrollment_period,cluster_coverage,simnum) {
  # Inputs:
  # g - the graph to run the epidemic on
  # beta - Every infectious individual contacts all their neighbours in a time step
  # and infects each susceptible with hazard beta. So beta represents the hazard of infection from one
  # contact between an infectious and a susceptible.
  # num_introductions - how many separate introductions we expect on average from the main epidemic. This is used
  # to calibrate the external force of infection
  # VE - direct leaky efficacy of the vaccine
  # bTrial - whether we are running a trial or not
  # bCluster - indicator of whether we are running a SW-CRT (2), a cRCT (1) or an iRCT (0)
  # trial_startday - first day of trial enrollment
  # trial_length - end of follow-up of trial partcipants, counting from the first day of enrollment
  # num_enrolled_per_day - number of individuals/clusters enrolled per day in parallel-arm CRT
  # enrollment_period - length of enrollment period
  # cluster_coverage - The proportion of each cluster we expect to enroll
  
  list <- structure(NA,class="result")
  "[<-.result" <- function(x,...,value) {
    args <- as.list(match.call())
    args <- args[-c(1:2,length(args))]
    length(value) <- length(args)
    for(i in seq(along=args)) {
      a <- args[[i]]
      if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
    }
    x
  }
  
  # Recover and spread functions
  recover<-function(e_nodes,i_nodes,r_nodes,infperiod_shape,infperiod_rate) {
    # Input is a list of the exposed nodes, 
    # with number of days since infection and total incubation/latent
    # period, and equivalently for the infectious nodes.
    # For each of these nodes, we will add it to newinfectious if the number of days exposed has
    # reached the total length of the incubation period, and equivalently for the infectious nodes.
    
    # Advance infectious nodes first, otherwise you will doubly advance any nodes switching from
    # exposed to infectious at this time step
    indices_to_remove <- i_nodes[2,]>=i_nodes[3,]
    newremoved<-as.vector(i_nodes[1,indices_to_remove])
    
    # Add one day to the length of each infected individual's time infected
    i_nodes[2,] <- i_nodes[2,]+1
    
    # Remove any recovered from i_nodes and add to r_nodes
    i_nodes <- i_nodes[,!(i_nodes[1,] %in% newremoved),drop=FALSE]
    r_nodes <- c(r_nodes,newremoved)
    
    # Now advance exposed nodes
    indices_to_remove <- e_nodes[2,]>=e_nodes[3,]
    newinfectious<-as.vector(e_nodes[1,indices_to_remove])
    
    # Add one day to the length of each infected individual's time infected
    e_nodes[2,] <- e_nodes[2,]+1
    
    # Remove any progressing from e_nodes and add to i_nodes
    e_nodes <- e_nodes[,!(e_nodes[1,] %in% newinfectious),drop=FALSE]
    inf_periods <- rgamma(length(newinfectious),infperiod_shape,infperiod_rate)
    i_nodes <- cbind(i_nodes,rbind(newinfectious,rep(0,length(newinfectious)),inf_periods))
    
    list(e_nodes, i_nodes, r_nodes, sort(newinfectious))
  }
  
  spread<-function(g, s_nodes, v_nodes, e_nodes, i_nodes, beta, VE,
                   incperiod_shape, incperiod_rate,
                   connected_nodes,external_inf_F,source_num_inf){
    # Spread will create new infected nodes from two sources: infectious nodes within the the study
    # population, and external pressure from the source population
    # Inputs:
    # g is the graph, used to find neighbours of infected nodes
    # s_nodes, e_nodes and i_nodes are susceptible, exposed and infected nodes
    # beta is the hazard of infection for one contact
    # incperiod_shape and rate are used to assign each newly exposed node a latent/incubation period
    # length, currently drawn from a gamma distribution
    # connected_nodes is a list of nodes that are connected the the source population
    # external_inf_F is a constant of proportionality that defines infectious pressure from source population 
    # to an individual
    # source_num_inf is the number of infectious individuals in the source population
    
    # Process: go through list of i_nodes, and choose a random number of its susceptible neighbours to
    # be infected, according to beta and choose a random number of its susceptible vaccinated neighbours to
    # be infected, according to beta and VE
    # Then go through list of nodes that are connected to source population and infect each susceptible
    # one with probability 1-exp(-FI), where I is the number/proportion of infectious, and F is a constant
    
    if (ncol(i_nodes)>0) {
      
      # Make a beta vector
      betavec <- rep(beta,ncol(i_nodes))
      beta_v <- betavec*(1-VE)
      
      # Get a list of all neighbours of all infected nodes
      potential_contacts<-lapply(i_nodes[1,],function(x) neighbors(g,x))
      susc_contacts<-lapply(potential_contacts,function(x,susceptibles) intersect(x,susceptibles),susceptibles=s_nodes)
      num_neighbours_susc<-rapply(susc_contacts,length)
      # Sample from each group of neighbours in turn
      # First choose how many neighbours each node infects
      num_contacts_susc<-rbinom(length(num_neighbours_susc),num_neighbours_susc,1-exp(-betavec))
      # Then sample from the neighbours
      # If one node gets picked twice by different nodes, just discard the duplicate.
      # In the very rare case that each i_nodes makes a number of new infectees equal to the number
      # of infectious nodes, mapply will make a matrix and unlist won't work. Therefore, the c() around
      # it ensures we turn the matrix into a vector. Unique then removes duplicates.
      infectees_susc<-c(unlist(mapply(function(x,y) x[sample.int(length(x),y)],x=susc_contacts,y=num_contacts_susc)))
      infectees_susc<-unique(infectees_susc)
      
      if (length(v_nodes)>0) {
        # Same as above but with vaccinated susceptible nodes
        vacc_contacts<-lapply(potential_contacts,function(x,vacc) intersect(x,vacc),vacc=v_nodes)
        num_neighbours_vacc<-rapply(vacc_contacts,length)
        num_contacts_vacc<-rbinom(length(num_neighbours_vacc),num_neighbours_vacc,1-exp(-beta_v))
        infectees_vacc<-c(unlist(mapply(function(x,y) x[sample.int(length(x),y)],x=vacc_contacts,y=num_contacts_vacc)))
        infectees_vacc<-unique(infectees_vacc)
      } else {
        infectees_vacc<-c()
      }
      
      infectees <- c(infectees_susc, infectees_vacc)
      
    } else {
      infectees_susc <- c()
      infectees_vacc <- c()
      infectees <- c()
    }
    # Pick out the nodes connected to the source that are still susceptible 
    # and haven't just been infected
    target_cnodes_susc <- setdiff(intersect(connected_nodes,s_nodes),infectees_susc)
    target_cnodes_vacc <- setdiff(intersect(connected_nodes,v_nodes),infectees_vacc)
    
    # Make a vector to represent external infection hazard for each individual
    communities_s <- V(g)[target_cnodes_susc]$community
    communities_v <- V(g)[target_cnodes_vacc]$community
    comm_sizes_s <- sapply(1:num_communities,function(x) sum(communities_s==x))
    comm_sizes_v <- sapply(1:num_communities,function(x) sum(communities_v==x))
    
    # Hazard of infection
    extFs_s<-rep(extF,comm_sizes_s)
    extFs_v<-rep(extF,comm_sizes_v)
    
    # Probability of infection
    prob_inf_fromsource <- 1 - exp(-mean(extFs_s)*source_num_inf)
    prob_inf_fromsource_v <- 1 - exp(-(1-VE)*mean(extFs_v)*source_num_inf)
    
    # Choose a number of individuals to be infected, then sample those individuals
    if (length(target_cnodes_susc)>0) {
      num_conn_inf_susc <- rbinom(1,length(target_cnodes_susc),prob_inf_fromsource)
      conn_inf_susc <- target_cnodes_susc[sample.int(length(target_cnodes_susc),num_conn_inf_susc,prob=extFs_s)]
    } else {
      conn_inf_susc <- c()
    }
    
    # Same as above, but for vaccinated individuals
    if (length(target_cnodes_vacc)>0) {
      num_conn_inf_vacc <- rbinom(1,length(target_cnodes_vacc),prob_inf_fromsource_v)
      conn_inf_vacc <- target_cnodes_vacc[sample.int(length(target_cnodes_vacc),num_conn_inf_vacc,prob=extFs_v)]
    } else {
      conn_inf_vacc <- c()
    }
    
    newinfected_susc <- c(infectees_susc,conn_inf_susc)
    newinfected_vacc <- c(infectees_vacc,conn_inf_vacc)
    newinfected <- c(newinfected_susc, newinfected_vacc)
    newinfected <- unique(newinfected)
    
    if (length(newinfected)>0) {
      
      # Give each newly exposed node an incubation/latent period
      if (FixInc==1) {
        inc_periods <- rep(ceiling(incperiod_shape/incperiod_rate),length(newinfected))
      } else {
        inc_periods <- rgamma(length(newinfected),incperiod_shape,incperiod_rate)
      }
      # Add them to e_nodes and remove from s_nodes and v_nodes
      e_nodes <- cbind(e_nodes,rbind(newinfected,rep(0,length(newinfected)),inc_periods))
      s_nodes<-setdiff(s_nodes,newinfected_susc)
      v_nodes <- setdiff(v_nodes,newinfected_vacc)
      
    }
    
    list(s_nodes, v_nodes, e_nodes)
  }
  
  #### RUN THE EPIDEMIC IN THE SOURCE POPULATION ####
  # This is to define external infectious pressure to the network
  
  model <- function(t, y, parms) {
    with(as.list(c(y,parms)), {
      
      beta <- betahat * (1 - a2/(1 + exp(-a1 * (t - atau))))
      
      dS <- -beta * S * (I1+I2+I3) / (S+E1+E2+E3+I1+I2+I3+R)
      dE1 <- beta * S * (I1+I2+I3) / (S+E1+E2+E3+I1+I2+I3+R) - sigma * 3 * E1
      dE2 <- sigma * 3 * E1 - sigma * 3 * E2
      dE3 <- sigma * 3 * E2 - sigma * 3 * E3
      dI1 <- sigma * 3 * E3 - gamma * 3 * I1
      dI2 <- gamma * 3 * I1 - gamma * 3 * I2
      dI3 <- gamma * 3 * I2 - gamma * 3 * I3
      dR <- gamma * 3 * I3
      list(c(dS,dE1,dE2,dE3,dI1,dI2,dI3,dR))
    })
  }
  N <- 50000
  y<- c(S=N-1,E1=0,E2=0,E3=0,I1=1,I2=0,I3=0,R=0)
  times<-seq(0,730,1)
  parms<-c(betahat=0.94,a1=0.19,a2=0.6,atau=27.79,sigma=0.14,gamma=0.33)
  out<-as.data.frame(lsoda(y,times,model,parms))
  
  #### RUN THE EPIDEMIC IN THE STUDY POPULATION ####
  
  # Define how the study population is linked to the source population
  # Connect all individuals to source population at same hazard
  # Constant of proportionality varies by community
  studypop_size<-length(V(g))
  connected_to_source <- V(g)$name
  
  # Calibrate extF to the number of introductions, given the progression of the epidemic in the source population
  num_communities <- max(V(g)$community)
  comm_sizes <- sapply(1:num_communities,function(x) length(V(g)[community==x]))
  sumsqrt <- sum(sqrt(comm_sizes))
  extF <- -log(1-num_introductions/(sqrt(comm_sizes)*sumsqrt))/trapz(times,out$I1+out$I2+out$I3)
  
  # Number of timesteps to run the epidemic - only need to go until the end of the trial
  num_timesteps <- trial_startday + trial_length + enrollment_period - 1
  
  if (bTrial) {
    
    # Parameters to do with trial recruitment
    # Enrollment per day is number of clusters enrolled per day
    enrollment_schedule <- rep(num_enrolled_per_day,enrollment_period)
    enroll_endday <- trial_startday+enrollment_period-1
    
    non_trial_clusters <- 1:max(V(g)$community)
    trial_clusters <- NULL
    
    if (bCluster > 0) {
      vax_clusters <- NULL
    }
    
    if (bCluster == 2) {
      community_crossovers <- data.frame(Community=NULL,DayVaccinated=NULL)
      crossover_days <- seq(from=first_crossover, by=step_interval, length.out=num_steps)
    }
  }
  
  # Initialize the S, E, I, and R nodes. I seed the epidemic from an SIR curve in a source population,
  # so initially all nodes in the study population are susceptible
  # i_nodes and e_nodes are matrices. The first row is the  identity of the node. The second row
  # is the number of days since infection/infectiousness. The third row is the total incubation/infectious period, 
  # drawn from a distribution when it becomes infected/infectious.
  # We are only going to consider a vaccine with effects on susceptibility, so only need one
  # vaccinated class
  e_nodes<-matrix(nrow=3,ncol=0)
  i_nodes<-matrix(nrow=3,ncol=0)
  v_nodes<-c()
  s_nodes <- as.vector(V(g))
  r_nodes <- c()
  
  # Initialize results.
  # Results will be, for each newly-infected node, the identity of the node, the day it was infected,
  # the community of which it is a member and its trial status at the time of infection. 
  # This is enough information to run a Cox PH with gamma frailty.
  # Make a data frame the size of the study pop and fill it in, then trim at the end
  results<-data.frame("SimulationNumber"=rep(NA,studypop_size),
                      "InfectedNode"=rep(NA,studypop_size),
                      "DayInfected"=rep(NA,studypop_size),
                      "Community"=rep(NA,studypop_size),
                      "TrialStatus"=rep(NA,studypop_size),
                      "DayEnrolled"=rep(NA,studypop_size))
  numinfectious<-0
  
  for (t in 1:num_timesteps) {

    # I'm recovering first, so I need to ensure that everyone has at least one chance to infect.
    # I do this by initializing an infectious node with 0 days since infection, seeing whether they
    # recover, then advancing them one day along their infectious period.
    
    if (bTrial) {
      
      # Recruit and randomize if during the enrollment period
      if ((t>=trial_startday) && (t<=enroll_endday)) {
        
        num_to_enroll <- enrollment_schedule[t-trial_startday+1]
        
        if (bCluster == 0) {
          # Individually-randomized trial, stratifying on community
          
          # Need to choose from the clusters not already enrolled
          if (length(non_trial_clusters)==num_to_enroll) {
            new_clusters <- non_trial_clusters
          } else {
            new_clusters <- sample(x=non_trial_clusters,size=num_to_enroll)
          }
          
          # From the chosen clusters, choose a fraction of the non-infectious individual. That fraction is defined in the inputs
          # I will then vaccinate half of each chosen sample
          # For each new cluster, I sample from that cluster a proportion of the whole cluster,
          # but only from the susceptible or exposed individuals. If I'm trying to sample more than are available (because
          # there are lots of infectious/recovered individuals), just sample all of them.
          # These are the people who are recruited - I then assign half and half to vaccine or control
          new_recruits <- lapply(new_clusters,
                                 function(x) sample(intersect(V(g)[community==x]$name,c(e_nodes[1,],s_nodes)),
                                                    min(round(cluster_coverage*length(V(g)[community==x]$name)),
                                                        length(intersect(V(g)[community==x]$name,c(e_nodes[1,],s_nodes))))))
          new_vacc <- unlist(lapply(new_recruits,
                                    function(x) sample(x,round(length(x)/2))))
          new_controls <- setdiff(unlist(new_recruits),new_vacc)
          
          
          non_trial_clusters<-setdiff(non_trial_clusters,new_clusters)
          
        } else if (bCluster == 1) {
          
          # We try and enroll as many from the cluster as you can. I have set an
          # enrollment rate rather than cluster size, e.g. 70% enrolled in each cluster.
          # It means that every simulated trial would have slightly different numbers enrolled 
          # (=coverage*ave_community_size)
          
          # Need to choose from the clusters not already enrolled
          if (length(non_trial_clusters)==num_to_enroll) {
            new_clusters <- non_trial_clusters
          } else {
            new_clusters <- sample(x=non_trial_clusters,size=num_to_enroll)
          }
          new_clusters_v <- sample(new_clusters,num_to_enroll/2)
          new_clusters_c <- setdiff(new_clusters,new_clusters_v)
          # From the chosen clusters, a fraction of the non-infectious individual. That fraction is defined in the inputs
          # This looks complicated: for each new cluster, I sample from that cluster a proportion of the whole cluster,
          # but only from the susceptible or exposed individuals. If I'm trying to sample more than are available (because
          # there are lots of infectious/recovered individuals), just sample all of them.
          new_vacc <- unlist(lapply(new_clusters_v,
                                    function(x) sample(intersect(V(g)[community==x]$name,c(e_nodes[1,],s_nodes)),
                                                       min(round(cluster_coverage*length(V(g)[community==x]$name)),
                                                           length(intersect(V(g)[community==x]$name,c(e_nodes[1,],s_nodes)))))))
          new_controls <- unlist(lapply(new_clusters_c,
                                        function(x) sample(intersect(V(g)[community==x]$name,c(e_nodes[1,],s_nodes)),
                                                           min(round(cluster_coverage*length(V(g)[community==x]$name)),
                                                               length(intersect(V(g)[community==x]$name,c(e_nodes[1,],s_nodes)))))))
          
          non_trial_clusters<-setdiff(non_trial_clusters,new_clusters)
          vax_clusters <- c(vax_clusters,new_clusters_v)
        } else {
          # Need to choose from the clusters not already enrolled
          if (length(non_trial_clusters)==num_to_enroll) {
            new_clusters <- non_trial_clusters
          } else {
            new_clusters <- sample(x=non_trial_clusters,size=num_to_enroll)
          }
          new_controls <- unlist(lapply(new_clusters,
                                        function(x) sample(intersect(V(g)[community==x]$name,c(e_nodes[1,],s_nodes)),
                                                           min(round(cluster_coverage*length(V(g)[community==x]$name)),
                                                               length(intersect(V(g)[community==x]$name,c(e_nodes[1,],s_nodes)))))))
          non_trial_clusters<-setdiff(non_trial_clusters,new_clusters)
        }
        
        trial_clusters <- c(trial_clusters,new_clusters)
        
        V(g)[name %in% new_controls]$trialstatus<-0
        V(g)[name %in% new_controls]$enrollmentday<-t
        if (bCluster < 2) {
          V(g)[name %in% new_vacc]$trialstatus<-1
          V(g)[name %in% new_vacc]$enrollmentday<-t
          V(g)[name %in% new_vacc]$treatmentday<-t
          
          # Move the vaccinated susceptibles to from s_nodes to v_nodes
          vacc_susc <- intersect(s_nodes,new_vacc)
          s_nodes <- setdiff(s_nodes,vacc_susc)
          v_nodes <- c(v_nodes,vacc_susc)
        }
        
      }
      if (bCluster==2) {
        if (t %in% crossover_days) {
          cont_clusters <- setdiff(trial_clusters,vax_clusters)
          if (length(cont_clusters)==num_clusters_per_step) {
            new_clusters_v = cont_clusters
          } else {
            new_clusters_v <- sample(x=cont_clusters,size=num_clusters_per_step)
          }
          vax_clusters <- c(vax_clusters,new_clusters_v)
          new_vacc <- unlist(lapply(new_clusters_v,
                                    function(x) intersect(V(g)[community==x]$name[!is.na(V(g)[community==x]$trialstatus)],c(e_nodes[1,],s_nodes))))
          V(g)[name %in% new_vacc]$trialstatus<-1
          community_crossovers <- rbind(community_crossovers,data.frame(Community=new_clusters_v,DayVaccinated=rep(t,length(new_clusters_v))))
          # Move the vaccinated susceptibles from s_nodes to v_nodes
          vacc_susc <- intersect(s_nodes,new_vacc)
          s_nodes <- setdiff(s_nodes,vacc_susc)
          v_nodes <- c(v_nodes,vacc_susc)
        }
      }
    }
    
    # Only need to recover if there are any infected or exposed
    if ((ncol(i_nodes)>0)||(ncol(e_nodes)>0)) {
      list[e_nodes,i_nodes,r_nodes,newinfectious]<-
        recover(e_nodes,i_nodes,r_nodes,infperiod_shape,infperiod_rate)
      
    } else {
      newinfectious <- c()
    }
    
    list[s_nodes,v_nodes,e_nodes]<-
      spread(g,s_nodes,v_nodes,e_nodes,i_nodes,
             beta,VE,incperiod_shape,incperiod_rate,
             connected_to_source,extF,out$I1[t]+out$I2[t]+out$I3[t])
    
    numnewinfectious<-length(newinfectious)
    if (numnewinfectious>0) {
      
      newcommunities <- V(g)[name %in% newinfectious]$community
      
      # Update results
      results$SimulationNumber[(numinfectious+1):(numinfectious+numnewinfectious)]<-rep(simnum,numnewinfectious)
      results$InfectedNode[(numinfectious+1):(numinfectious+numnewinfectious)]<-newinfectious
      results$DayInfected[(numinfectious+1):(numinfectious+numnewinfectious)]<-rep(t,numnewinfectious)
      results$Community[(numinfectious+1):(numinfectious+numnewinfectious)]<-newcommunities
      results$TrialStatus[(numinfectious+1):(numinfectious+numnewinfectious)]<-V(g)[name %in% newinfectious]$trialstatus
      results$DayEnrolled[(numinfectious+1):(numinfectious+numnewinfectious)]<-V(g)[name %in% newinfectious]$enrollmentday
      
      numinfectious <- numinfectious+numnewinfectious
      
    }
    
  }
  
  trial_nodes <- V(g)[!is.na(V(g)$trialstatus)]$name
  trial_nodes_info<-data.frame("SimulationNumber"=rep(simnum,length(trial_nodes)),
                               "Node"=trial_nodes,
                               "Community"=V(g)[trial_nodes]$community,
                               "TrialStatus"=V(g)[trial_nodes]$trialstatus,
                               "DayEnrolled"=V(g)[trial_nodes]$enrollmentday)
  
  # Tidy up results
  if (numinfectious>0) {
    results<-results[1:numinfectious,]
  } else {
    results<-results[1,]
    results$SimulationNumber[1]<-simnum
  }
  
  if (bTrial) {
    if (bCluster==2) {
      list(results,trial_nodes_info,community_crossovers)
    } else {
      list(results,trial_nodes_info)
    }
  } else {
    list(results,trial_nodes_info)
  }
  
  
}

analyse_data <- function(results,trial_nodes,trial_startday,trial_length,ave_inc_period,
                         bCluster,
                         crossovers=NULL,trial_inf_nodes=NULL,trial_uninf_nodes=NULL) {
  
  list <- structure(NA,class="result")
  "[<-.result" <- function(x,...,value) {
    args <- as.list(match.call())
    args <- args[-c(1:2,length(args))]
    length(value) <- length(args)
    for(i in seq(along=args)) {
      a <- args[[i]]
      if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
    }
    x
  }
  
  coxmodel <- function(data,VEpointest) {
    
    survmodel<-try(coxph(Surv(DayInfected,eventstatus)~TrialStatus+strata(Community),data),silent=T)
    usesurvmod <- !inherits(survmodel, 'try-error')
    
    if (usesurvmod && vcov(survmodel)>=0){
      # If no error was thrown and the variance is positive, use the results of the model
      
      vaccEffEst <- 1-exp(survmodel$coefficient + c(0, 1.96, -1.96)*sqrt(survmodel$var))
      zval <- survmodel$coefficient/sqrt(survmodel$var)
      pval <- pnorm(zval, lower.tail = vaccEffEst[1]>0)*2
      
    } else {
      
      vaccEffEst<-c(VEpointest,NA,NA)
      pval <- NA
      
    }
    
    list(vaccEffEst,pval)
    
  }
  
  clustermodels <- function(data,VEpointest) {
    
    # Gaussian shared frailty, using coxme
    mod <- try(coxme(Surv(DayInfected, eventstatus) ~ TrialStatus + (1|Community), data = data),silent=T)
    usemod <- !inherits(mod, 'try-error')
    
    if (usemod && is.na(vcov(mod))) {
      vaccEffEst_gaussian_coxme <- c(VEpointest,NA,NA)
      pval_gaussian_coxme<-NA
    } else {
      if (usemod && vcov(mod)>=0){
        # If the model converges and variance is positive, use results of model
        
        vaccEffEst_gaussian_coxme <- 1-exp(mod$coefficients + c(0, 1.96, -1.96)*sqrt(vcov(mod)))
        zval <- mod$coefficients/sqrt(vcov(mod))
        pval_gaussian_coxme <- pnorm(zval, lower.tail = vaccEffEst_gaussian_coxme[1]>0)*2
      } else {
        vaccEffEst_gaussian_coxme <- c(VEpointest,NA,NA)
        pval_gaussian_coxme<-NA
      }
    }
    
    # Try with a gaussian-distributed random effect.
    mod2 <- try(coxph(Surv(DayInfected, eventstatus) ~ TrialStatus + frailty(Community,distribution="gaussian",sparse=FALSE), data = data), silent=T)
    usemod2 <- !inherits(mod2, 'try-error')
    
    if (usemod2 && is.na(vcov(mod2))) {
      vaccEffEst_gaussian_coxph <- c(VEpointest,NA,NA)
      pval_gaussian_coxph<-NA
    } else {
      if (usemod2 && vcov(mod2)>=0){
        
        vaccEffEst_gaussian_coxph <- 1-exp(mod2$coefficients[1] + c(0, 1.96, -1.96)*sqrt(vcov(mod2)[1]))
        zval2 <- mod2$coefficients[1]/sqrt(vcov(mod2)[1])
        pval_gaussian_coxph <- pnorm(zval2, lower.tail = vaccEffEst_gaussian_coxph[1]>0)*2
      } else {
        vaccEffEst_gaussian_coxph <- c(VEpointest,NA,NA)
        pval_gaussian_coxph<-NA
      }
    }
    
    # Try with a gamma-distributed random effect, using coxph and frailty
    mod3 <- try(coxph(Surv(DayInfected, eventstatus) ~ TrialStatus + frailty(Community,distribution="gamma",sparse=FALSE), data = data), silent=T)
    usemod3 <- !inherits(mod3, 'try-error')
    
    if (usemod3 && is.na(vcov(mod3))) {
      vaccEffEst_gamma_coxph <- c(VEpointest,NA,NA)
      pval_gamma_coxph<-NA
    } else {
      if (usemod3 && vcov(mod3)>=0){
        
        vaccEffEst_gamma_coxph <- 1-exp(mod3$coefficients[1] + c(0, 1.96, -1.96)*sqrt(vcov(mod3)[1]))
        zval3 <- mod3$coefficients[1]/sqrt(vcov(mod3)[1])
        pval_gamma_coxph <- pnorm(zval3, lower.tail = vaccEffEst_gamma_coxph[1]>0)*2
      } else {
        vaccEffEst_gamma_coxph <- c(VEpointest,NA,NA)
        pval_gamma_coxph<-NA
      }
    }
    
    # Generalized estimating equations/robust variance
    mod5 <- try(coxph(Surv(DayInfected, eventstatus) ~ TrialStatus + cluster(Community), data = data), silent=T)
    usemod5 <- !inherits(mod5, 'try-error')
    
    if (usemod5 && is.na(vcov(mod5))) {
      vaccEffEst_gee <- c(VEpointest,NA,NA)
      pval_gee<-NA
    } else {
      if (usemod5 && vcov(mod5)>=0){
        
        vaccEffEst_gee <- 1-exp(mod5$coefficients[1] + c(0, 1.96, -1.96)*sqrt(vcov(mod5)[1]))
        zval5 <- mod5$coefficients[1]/sqrt(vcov(mod5)[1])
        pval_gee <- pnorm(zval5, lower.tail = vaccEffEst_gee[1]>0)*2
      } else {
        vaccEffEst_gee <- c(VEpointest,NA,NA)
        pval_gee<-NA
      }
    }
    
    list(vaccEffEst_gaussian_coxme,pval_gaussian_coxme,
         vaccEffEst_gaussian_coxph,pval_gaussian_coxph,
         vaccEffEst_gamma_coxph,pval_gamma_coxph,
         vaccEffEst_gee,pval_gee)
    
  }
  
  results$DayInfected <- results$DayInfected - results$DayEnrolled
  
  # Get a list of nodes that were enrolled in the trial but never infected
  noninf<-setdiff(trial_nodes$Node,results$InfectedNode)
  # Get list of nodes that became infectious while they were in the trial
  # This is the step that excludes those who were R at time of enrollment
  results_analysis<-results[!is.na(results$TrialStatus),]
  # Get a list of nodes who were infected after their follow-up time was over
  # (i.e. those enrolled at the beginning but infected right at the end)
  censored <- results_analysis[results_analysis$DayInfected>trial_length,]
  results_analysis<-results_analysis[results_analysis$DayInfected<=trial_length,]
  # Assign them eventstatus=1 for the Cox analysis
  results_analysis$eventstatus<-rep(1,nrow(results_analysis))
  # Make data frame for those who were never infected (i.e. censored by end of study)
  noninfdf<-data.frame(InfectedNode=noninf,DayInfected=rep(trial_length,length(noninf)),
                       Community=trial_nodes$Community[trial_nodes$Node %in% noninf],
                       TrialStatus=trial_nodes$TrialStatus[trial_nodes$Node %in% noninf],
                       eventstatus=rep(0,length(noninf)),
                       DayEnrolled=trial_nodes$DayEnrolled[trial_nodes$Node %in% noninf])
  if (nrow(censored)>0) {
    censored$DayInfected<-trial_length
    censored$eventstatus<-0
  }
  # Remove column with simulation number so the columns match up
  results_analysis$SimulationNumber<-NULL
  results_analysis<-rbind(results_analysis,noninfdf,censored)
  
  # Finally, exclude any cases who were infected during the first n days of follow-up
  # This tries to rid of those who were already latently infected when enrolled
  results_analysis<-results_analysis[results_analysis$DayInfected>ave_inc_period,]
  
  numevents_vacc <- nrow(results_analysis[(results_analysis$eventstatus==1) & (results_analysis$TrialStatus==1),])
  numevents_cont <- nrow(results_analysis[(results_analysis$eventstatus==1) & (results_analysis$TrialStatus==0),])
  
  total_vacc_pt <- sum(results_analysis$DayInfected[results_analysis$TrialStatus==1])
  total_cont_pt <- sum(results_analysis$DayInfected[results_analysis$TrialStatus==0])
  VE_pointest <- 1 - (numevents_vacc/total_vacc_pt)/(numevents_cont/total_cont_pt)
  
  sample_size <- nrow(results_analysis)
  
  # Calculate ICC
  # Number of events in each cluster
  events_by_cluster<-aggregate(results_analysis$eventstatus,
                               by=list(Community=results_analysis$Community),FUN=sum)
  # Size of each cluster
  cluster_sizes<-aggregate(results_analysis$InfectedNode,by=list(Community=results_analysis$Community),FUN=length)
  # Overall trial size
  N <- sum(cluster_sizes$x)
  # Overall number of clusters
  K <- nrow(cluster_sizes)
  n0 <- 1/(K-1) * (N - sum(cluster_sizes$x^2)/N)
  n01 <- 1/(K-1) * ((K-1)*n0 - sum(cluster_sizes$x^2)/N)
  MSB <- 1/(K-1) * sum((events_by_cluster$x-mean(events_by_cluster$x))^2/cluster_sizes$x)
  MSW <- 1/(N-K-1) * sum(events_by_cluster$x-events_by_cluster$x^2/cluster_sizes$x)
  ICC <- (MSB - MSW) / (MSB + (n01-1) * MSW)
  deff <- 1+(mean(cluster_sizes$x)-1)*ICC
  
  # Proportion of clusters that have zero cases
  prop_zeros <- sum(events_by_cluster$x==0)/K
  
  if (bCluster == 0) {
    # Analysis for iRCT
    if ((numevents_vacc>0)&&(numevents_cont>0)) {
      # If we have events in both arms, can try a Cox PH. It can still be singular, so if it
      # throws an error, the trial has failed (not enough events)
      list[vaccEffEst,pval] <- coxmodel(results_analysis,VE_pointest)
      
    } else if ((numevents_vacc>0)&&(numevents_cont==0)) {
      # If there are no events in the control arm but
      # events in the vaccine arm, VE estimate is -1 and p-value is 1
      vaccEffEst<-c(-1,-1,-1)
      pval <- 1
      
    } else if ((numevents_vacc==0)&&(numevents_cont>0)) {
      # If there are no events in the vaccine arm and events in the control arm, VE is 1, but
      # for the p-value we add one event to both arms and do a Cox regression on that data
      # I give both events the median time among control events.
      newevent_v_rownum <- min(which((results_analysis$eventstatus==0)&(results_analysis$TrialStatus==1)))
      newevent_c_rownum <- min(which((results_analysis$eventstatus==0)&(results_analysis$TrialStatus==0)))
      
      eventtime <- median(results_analysis$DayInfected[results_analysis$eventstatus==1])
      
      results_analysis$DayInfected[newevent_v_rownum] <- eventtime
      results_analysis$eventstatus[newevent_v_rownum] <- 1
      
      results_analysis$DayInfected[newevent_c_rownum] <- eventtime
      results_analysis$eventstatus[newevent_c_rownum] <- 1
      
      list[,pval] <- coxmodel(results_analysis,VE_pointest)
      vaccEffEst<-1
      
    } else {
      # If no events are observed in either arm, the trial has failed and no result can be obtained
      vaccEffEst<-c(NA,NA,NA)
      pval <- NA
      
    }
    
    list(vaccEffEst,pval,numevents_vacc,numevents_cont,sample_size)
    
  } else if (bCluster==1) {
    
    if ((numevents_vacc>0)&&(numevents_cont>0)) {
      # Run the models for the clustered data
      list[vaccEffEst_gaussian_coxme,pval_gaussian_coxme,
           vaccEffEst_gaussian_coxph,pval_gaussian_coxph,
           vaccEffEst_gamma_coxph,pval_gamma_coxph,
           vaccEffEst_gee,pval_gee] <- 
        clustermodels(results_analysis,VE_pointest)
      
    } else if ((numevents_vacc>0)&&(numevents_cont==0)) {
      # If there are no events in the control arm but
      # events in the vaccine arm, VE estimate is -1 and p-value is 1
      vaccEffEst_gaussian_coxme<-c(-1,-1,-1)
      pval_gaussian_coxme<-1
      vaccEffEst_gaussian_coxph<-c(-1,-1,-1)
      pval_gaussian_coxph<-1
      vaccEffEst_gamma_coxph<-c(-1,-1,-1)
      pval_gamma_coxph<-1
      vaccEffEst_gee<-c(-1,-1,-1)
      pval_gee<-1
      
    } else if ((numevents_vacc==0)&&(numevents_cont>0)) {
      # If there are no events in the vaccine arm and events in the control arm, VE is 1
      # For p-value, add one event to both arms. Give it the median event time among control events.
      # Cluster for event in vaccine arm doesn't matter
      # Choose cluster for event in control arm to be most conservative
      
      communities <- results_analysis$Community[((results_analysis$eventstatus==1)&(results_analysis$TrialStatus==0))]
      freq_table <- sort(table(communities),decreasing=TRUE)
      community <- as.numeric(names(freq_table[1]))
      
      newevent_v_rownum <- min(which((results_analysis$eventstatus==0)&(results_analysis$TrialStatus==1)))
      newevent_c_rownum <- min(which((results_analysis$eventstatus==0)&(results_analysis$TrialStatus==0)&(results_analysis$Community==community)))
      
      eventtime <- median(results_analysis$DayInfected[results_analysis$eventstatus==1])
      
      results_analysis$DayInfected[newevent_v_rownum] <- eventtime
      results_analysis$eventstatus[newevent_v_rownum] <- 1
      
      results_analysis$DayInfected[newevent_c_rownum] <- eventtime
      results_analysis$eventstatus[newevent_c_rownum] <- 1
      
      list[,pval_gaussian_coxme,
           ,pval_gaussian_coxph,
           ,pval_gamma_coxph,
           ,pval_gee] <- 
        clustermodels(results_analysis,VE_pointest)
      vaccEffEst_gaussian_coxme<-1
      vaccEffEst_gaussian_coxph<-1
      vaccEffEst_gamma_coxph<-1
      vaccEffEst_gee<-1
      
    } else {
      # If no events are observed in either arm, the trial has failed and no result can be obtained
      vaccEffEst_gaussian_coxme<-c(NA,NA,NA)
      pval_gaussian_coxme<-NA
      vaccEffEst_gaussian_coxph<-c(NA,NA,NA)
      pval_gaussian_coxph<-NA
      vaccEffEst_gamma_coxph<-c(NA,NA,NA)
      pval_gamma_coxph<-NA
      vaccEffEst_gee<-c(NA,NA,NA)
      pval_gee<-NA
    }
    
    list(vaccEffEst_gaussian_coxme,pval_gaussian_coxme,
         vaccEffEst_gaussian_coxph,pval_gaussian_coxph,
         vaccEffEst_gamma_coxph,pval_gamma_coxph,
         vaccEffEst_gee,pval_gee,
         numevents_vacc,numevents_cont,sample_size,ICC,deff,prop_zeros)
  } else {
    ## Analysis for SW-CRT
    numpds <- ceiling(trial_length/step_interval)
    numclusts <- num_clusters_per_step*num_steps
    df <- data.frame(Period=rep(1:numpds,numclusts))
    df$StartDay <- trial_startday+(df$Period-1)*step_interval
    df$EndDay <- df$StartDay+step_interval
    df$EndDay <- ifelse(df$Period==numpds,df$EndDay+1,df$EndDay) #Captures last day of trial if that's the end of a period
    df$Cluster <- rep(crossovers$Community,each=numpds)
    StartVaccinated <- rep(crossovers$DayVaccinated,each=numpds)
    df$StartPd <- (StartVaccinated - trial_startday)/step_interval+1
    df$Infections <- apply(df, 1, 
                           function(x) sum(trial_inf_nodes$DayInfected[trial_inf_nodes$Community==x["Cluster"]] >= x["StartDay"] & trial_inf_nodes$DayInfected[trial_inf_nodes$Community==x["Cluster"]] < x["EndDay"]))
    df$InfAdj <- apply(df, 1, 
                           function(x) sum(trial_inf_nodes$DayInfected[trial_inf_nodes$Community==x["Cluster"]] >= x["StartDay"]+ave_inc_period & trial_inf_nodes$DayInfected[trial_inf_nodes$Community==x["Cluster"]] < x["EndDay"]+ave_inc_period))
    df$Enrolled <- apply(df, 1, 
                         function(x) sum(trial_inf_nodes$Community==x["Cluster"])+sum(trial_uninf_nodes$Community==x["Cluster"]))
    df$Indivs <- df$Enrolled
    BadClusts <- unique(df$Cluster[df$Indivs==0])
    df <- df[!(df$Cluster %in% BadClusts),] ## Removes All data from clusters that have 0 enrolled
    df$Outcomes <- df$Infections/df$Indivs
    df$OutcomesAdj <- df$InfAdj/df$Indivs
    
    events_vacc <- sum(df$Infections[df$Period >= df$StartPd])
    events_cont <- sum(df$Infections[df$Period < df$StartPd])
    analysed_trialsize <- length(trial_inf_nodes$SimulationNumber)+length(trial_uninf_nodes$SimulationNumber)
    
    SWTRes <- SWT.Permutation.Analysis(df$Period, df$Outcomes, df$Cluster,
                                       StartTimes=df$StartPd, Indivs=df$Indivs,
                                       family=binomial, link=c("identity","log"),
                                       NumPerms=NumP, Type=c("NPWP","SC","CO.Ctrl","CO.Both","CO.CtWt","CO.BoWt",
                                                             "ENS"),
                                       NullVals=0, Alternative="Both", Alpha=0.05, 
                                       FwdOffset=1, BwdOffset=0)
    # SWTResMEM <- SWT.Permutation.Analysis(df$Period, df$Outcomes, df$Cluster,
    #                                     StartTimes=df$StartPd, Indivs=df$Indivs,
    #                                     family=binomial, link="logit",
    #                                     NumPerms=NumP, Type=c("MEM","CPI"),
    #                                     NullVals=0, Alternative="Both", Alpha=0.05)
    # SWTResMEMAdj <- SWT.Permutation.Analysis(df$Period, df$OutcomesAdj, df$Cluster,
    #                                       StartTimes=df$StartPd, Indivs=df$Indivs,
    #                                       family=binomial, link="logit",
    #                                       NumPerms=NumP, Type=c("MEM","CPI"),
    #                                       NullVals=0, Alternative="Both", Alpha=0.05)
    
    Ests.Ident <- unlist(SWTRes$RD)[-1]
    Ests.Log <- unlist(SWTRes$logRR)[-1]
    # Ests.MEM <- unlist(SWTResMEM)[-1]
    # Ests.MEMAdj <- unlist(SWTResMEMAdj)[-1]
    
    VEs.Ident <- ifelse(!is.finite(Ests.Ident),NA,Ests.Ident)
    VEs.Log <- ifelse(!is.finite(Ests.Log),NA,1-exp(Ests.Log))
    # VEs.MEM <- ifelse(!is.finite(Ests.MEM),NA,1-exp(Ests.MEM))
    # VEs.MEMAdj <- ifelse(!is.finite(Ests.MEMAdj),NA,1-exp(Ests.MEMAdj))
    
    if (NumP > 0) {
      PVals.Ident <- SWTRes$PVal$RD
      PVals.Log <- SWTRes$PVal$logRR
      # PVals.MEM <- SWTResMem$PVal
      # PVals.MEMAdj <- SWTResMemAdj$PVal
      list(PVals.Log,VEs.Log,PVals.Ident,VEs.Ident,
           # PVals.MEM,VEs.MEM,PVals.MEMAdj,VEs.MEMAdj,
           events_vacc,events_cont,analysed_trialsize,
           SWTRes$NPWP.ContrastMats,SWTRes$CO.ContrastMats,SWTRes$SC.ContrastMats
           # SWTResMEM$MEM.AsyPVals, SWTResMEM$CPI.AsyPVals,
           # SWTResMEMAdj$MEM.AsyPVals, SWTResMEMAdj$CPI.AsyPVals
           )
    } else {
      list(VEs.Log,VEs.Ident,
           # VEs.MEM,VEs.MEMAdj,
           events_vacc,events_cont,analysed_trialsize,
           SWTRes$NPWP.ContrastMats,SWTRes$CO.ContrastMats,SWTRes$SC.ContrastMats)
    }
  }
}

# Initialize vectors/data frames to store results
pvals<-data.frame('iRCT'=rep(NA,nsim),
                  'cRCT_gaussian_coxme'=rep(NA,nsim),
                  'cRCT_gaussian_coxph'=rep(NA,nsim),
                  'cRCT_gamma_coxph'=rep(NA,nsim),
                  'cRCT_gee'=rep(NA,nsim),
                  'SWT_log_NPWP.TW'=rep(NA,nsim),
                  'SWT_log_NPWP.Mean'=rep(NA,nsim),
                  'SWT_log_NPWP.MeanEx1'=rep(NA,nsim),
                  'SWT_log_NPWP.OW'=rep(NA,nsim),
                  'SWT_log_NPWP.OWEx1'=rep(NA,nsim),
                  'SWT_log_SCSWT1'=rep(NA,nsim),
                  'SWT_log_SCSWT2'=rep(NA,nsim),
                  'SWT_log_SCSWT3'=rep(NA,nsim),
                  'SWT_log_SCSWT1.Ex1'=rep(NA,nsim),
                  'SWT_log_SCSWT2.Ex1'=rep(NA,nsim),
                  'SWT_log_SCSWT3.Ex1'=rep(NA,nsim),
                  'SWT_log_CO.Ctrl'=rep(NA,nsim),
                  'SWT_log_CO.Both'=rep(NA,nsim),
                  'SWT_log_CO.CtWt'=rep(NA,nsim),
                  'SWT_log_CO.BoWt'=rep(NA,nsim),
                  'SWT_log_ENS'=rep(NA,nsim),
                  # 'SWT_logit_MEM'=rep(NA,nsim),
                  # 'SWT_logit_CPI'=rep(NA,nsim),
                  # 'SWT_logit_MEMAdj'=rep(NA,nsim),
                  # 'SWT_logit_CPIAdj'=rep(NA,nsim),
                  'SWT_ident_NPWP.TW'=rep(NA,nsim),
                  'SWT_ident_NPWP.Mean'=rep(NA,nsim),
                  'SWT_ident_NPWP.MeanEx1'=rep(NA,nsim),
                  'SWT_ident_NPWP.OW'=rep(NA,nsim),
                  'SWT_ident_NPWP.OWEx1'=rep(NA,nsim),
                  'SWT_ident_SCSWT1'=rep(NA,nsim),
                  'SWT_ident_SCSWT2'=rep(NA,nsim),
                  'SWT_ident_SCSWT3'=rep(NA,nsim),
                  'SWT_ident_SCSWT1.Ex1'=rep(NA,nsim),
                  'SWT_ident_SCSWT2.Ex1'=rep(NA,nsim),
                  'SWT_ident_SCSWT3.Ex1'=rep(NA,nsim),
                  'SWT_ident_CO.Ctrl'=rep(NA,nsim),
                  'SWT_ident_CO.Both'=rep(NA,nsim),
                  'SWT_ident_CO.CtWt'=rep(NA,nsim),
                  'SWT_ident_CO.BoWt'=rep(NA,nsim),
                  'SWT_ident_ENS'=rep(NA,nsim)
                  # 'SWT_asyMEM_logit'=rep(NA,nsim),
                  # 'SWT_asyCPI_logit'=rep(NA,nsim),
                  # 'SWT_asyMEMAdj_logit'=rep(NA,nsim),
                  # 'SWT_asyCPIAdj_logit'=rep(NA,nsim)
                  )
VEs<-pvals[,!(names(pvals) %in% c('SWT_asyMEM_logit','SWT_asyCPI_logit',
                                  'SWT_asyMEMAdj_logit','SWT_asyCPIAdj_logit'))]
ICCs<-rep(NA,nsim)
deffs<-rep(NA,nsim)
props_zeros <- rep(NA,nsim)
numevents_cont<-matrix(NA,nrow=3,ncol=nsim)
numevents_vacc<-matrix(NA,nrow=3,ncol=nsim)
ss<-matrix(NA,nrow=3,ncol=nsim)

st <- proc.time()

for (sim in 1:nsim) {
  
  set.seed(Varying$Seed[VaryingRow]+simno*20+sim-1)
  
  g<-make_network(ave_community_size, community_size_range, num_communities,rate_within, rate_between)
  
  list[results,trial_nodes]<-
    network_epidemic(g,beta,num_introductions,direct_VE,
                     incperiod_shape,incperiod_rate,infperiod_shape,infperiod_rate,1,0,
                     trial_startday,trial_length,num_clusters_enrolled_per_day,
                     enrollment_period,cluster_coverage,sim)
  
  list[VE,pval,events_vacc,events_cont,analysed_trialsize]<-
      analyse_data(results,trial_nodes,trial_startday,trial_length,ave_inc_period,0)
  pvals$iRCT[sim]<-pval
  VEs$iRCT[sim]<-VE[1]
  numevents_vacc[1,sim]<-events_vacc
  numevents_cont[1,sim]<-events_cont
  ss[1,sim]<-analysed_trialsize
  
  fullRes <- list(results=results, trial_nodes=trial_nodes)
  
  save(fullRes, file=paste0(resfolder,"/FullRes_",outname,"_IRT_",sim,".Rda"))
  
  list[results,trial_nodes]<-
    network_epidemic(g,beta,num_introductions,direct_VE,
                     incperiod_shape,incperiod_rate,infperiod_shape,infperiod_rate,1,1,
                     trial_startday,trial_length,num_clusters_enrolled_per_day,
                     enrollment_period,cluster_coverage,sim)
  list[VE_gaussian_coxme,pval_gaussian_coxme,
         VE_gaussian_coxph,pval_gaussian_coxph,
         VE_gamma_coxph,pval_gamma_coxph,
         VE_gee,pval_gee,
         events_vacc,events_cont,analysed_trialsize,ICC,deff,prop_zeros]<-
      analyse_data(results,trial_nodes,trial_startday,trial_length,ave_inc_period,1)
  
  pvals$cRCT_gaussian_coxme[sim]<-pval_gaussian_coxme
  VEs$cRCT_gaussian_coxme[sim]<-VE_gaussian_coxme[1]
  
  pvals$cRCT_gaussian_coxph[sim]<-pval_gaussian_coxph
  VEs$cRCT_gaussian_coxph[sim]<-VE_gaussian_coxph[1]
  
  pvals$cRCT_gamma_coxph[sim]<-pval_gamma_coxph
  VEs$cRCT_gamma_coxph[sim]<-VE_gamma_coxph[1]
  
  pvals$cRCT_gee[sim]<-pval_gee
  VEs$cRCT_gee[sim]<-VE_gee[1]
  
  numevents_vacc[2,sim]<-events_vacc
  numevents_cont[2,sim]<-events_cont
  ss[2,sim]<-analysed_trialsize
  
  ICCs[sim]<-ICC
  deffs[sim]<-deff
  props_zeros[sim]<-prop_zeros
  
  fullRes <- list(results=results, trial_nodes=trial_nodes)
  
  save(fullRes, file=paste0(resfolder,"/FullRes_",outname,"_CRT_",sim,".Rda"))
  
  list[results,trial_nodes,community_crossovers]<-
    network_epidemic(g,beta,num_introductions,direct_VE,
                     incperiod_shape,incperiod_rate,infperiod_shape,infperiod_rate,1,2,
                     trial_startday,trial_length,num_clusters_enrolled_per_day,
                     enrollment_period,cluster_coverage,sim)
  
  trial_infected <- results[!is.na(results$TrialStatus),]
  trial_uninfected <- trial_nodes[!(trial_nodes$Node %in% results$InfectedNode),]
  
  if (NumP > 0) {
    list[PVals.Log,VEs.Log,PVals.Ident,VEs.Ident,
         # PVals.MEM,VEs.MEM,PVals.MEMAdj,VEs.MEMAdj,
         events_vacc,events_cont,analysed_trialsize,
         NPWP.ContrastMats,CO.ContrastMats,
         SC.ContrastMats
         # MEM.AsyPVals,CPI.AsyPVals,
         # MEMAdj.AsyPVals,CPIAdj.AsyPVals
         ]<-analyse_data(results,trial_nodes,trial_startday,
                                        trial_length,ave_inc_period,2,
                                        community_crossovers,trial_infected,trial_uninfected)
  } else {
    list[VEs.Log,VEs.Ident,
         # VEs.MEM,VEs.MEMAdj,
         events_vacc,events_cont,analysed_trialsize,
         NPWP.ContrastMats,CO.ContrastMats,
         SC.ContrastMats
         # MEM.AsyPVals,CPI.AsyPVals,
         # MEMAdj.AsyPVals,CPIAdj.AsyPVals
         ]<-analyse_data(results,trial_nodes,trial_startday,
                                               trial_length,ave_inc_period,2,
                                               community_crossovers,trial_infected,trial_uninfected)
  }
  
  numevents_vacc[3,sim]<-events_vacc
  numevents_cont[3,sim]<-events_cont
  ss[3,sim]<-analysed_trialsize
  
  if (NumP > 0) {
    pvals[sim,grep("^SWT_log_",colnames(pvals))] <- PVals.Log
    pvals[sim,grep("^SWT_ident_",colnames(pvals))] <- PVals.Ident
    # pvals[sim,"SWT_asyMEM_logit"] <- MEM.AsyPVals
    # pvals[sim,"SWT_asyCPI_logit"] <- CPI.AsyPVals
    # pvals[sim,"SWT_asyMEMAdj_logit"] <- MEMAdj.AsyPVals
    # pvals[sim,"SWT_asyCPIAdj_logit"] <- CPIAdj.AsyPVals
  }
  VEs[sim,grep("^SWT_log_",colnames(VEs))] <- VEs.Log
  VEs[sim,grep("^SWT_ident_",colnames(VEs))] <- VEs.Ident
  
  fullRes <- list(results=results, trial_nodes=trial_nodes, crossovers=community_crossovers,
                  trial_inf_nodes=trial_infected, trial_uninf_nodes=trial_uninfected, 
                  NPWP.Contrasts=NPWP.ContrastMats,CO.Contrasts=CO.ContrastMats,SC.Contrasts=SC.ContrastMats)
  
  save(fullRes, file=paste0(resfolder,"/FullRes_",outname,"_SWT_",sim,".Rda"))
  
  print(paste0("Simulation ",sim))
  
}

proc.time() - st

OverallRes <- list(beta=beta,pvals=pvals,VEs=VEs,ICCs=ICCs,deffs=deffs,props_zeros=props_zeros,
                   numevents_cont=numevents_cont,numevents_vacc=numevents_vacc,ss=ss)
save(OverallRes, file=paste0(scenfolder,"/OverallRes/OverallRes_",outname,".Rda"))


