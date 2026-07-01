## Master template script used for simulation-based power analysis for To et al. (2016) birth-death trees 
##

####################################################################

# Install necessary packages

#install.packages('ape')
#install.packages('phytools')
#install.packages('treedater') # install from source for alternative optimization schemes
#install.packages("phylopomp",repos="https://kingaa.github.io") # install from repo for most up-to-date version
#install.packages('matrixStats')
#install.packages('tryCatchLog')
#install.packages('future')
#install.packages('doFuture')
#install.packages('castor')
#install.packages('diptest')
#install.packages('stats')

# Load packages

library(ape)
library(treedater) # install from source for alternative optimization schemes
library(phytools)
library(phylopomp) # epidemiological sims
library(matrixStats)
library(tryCatchLog)
library(future)
library(doFuture)
library(castor)
library(diptest)
library(stats)

stopifnot(packageVersion("ape")>="5.8")
stopifnot(packageVersion("treedater")>="1.0.2")
stopifnot(packageVersion("phytools")>="2.3.0")
stopifnot(packageVersion("phylopomp")>="0.14.8.0")
stopifnot(packageVersion("matrixStats")>="1.5.0")
stopifnot(packageVersion("tryCatchLog")>="1.3.1")
stopifnot(packageVersion("future")>="1.34.0")
stopifnot(packageVersion("doFuture")>="1.0.1")
stopifnot(packageVersion("castor")>="1.8.3")
stopifnot(packageVersion("diptest")>="0.77.1")
stopifnot(packageVersion("stats")>="4.3.3")

startTime <- Sys.time() # Runtime

####################################################################
######################## Set-up ####################################
####################################################################

# Functions for simulation analyses
source('sim_funcs.R')

####################################################################
################## Run-modifying variables #########################
####################################################################
###
### Modify these variables as desired for each simulation run

### Birth-death runs ###

tree_type <- 'RMC_750_3_25'
scenario <- 'kscale_highshape'

## Clock simulation - parameter values (a - shape, b - scale)
a <- 6.25 # Gamma shape and scale that yield same mean and sd as lognormal mean = 1 and sd = 0.4, simulated in To et al. (2016)
b <- 0.16*0.006 # scaled to desired values using constant used in To et al. (2016)
# a <- a/2 # If low-shape

## Effect sizes for power analysis
klist <- c(1,1.1,1.2,1.3,1.4,1.5,2,5,10,25)

#### ENSURE TO FIRST SET UP FILE PATHS AS BELOW ###

# Set paths
timepath <- paste0('./to_etal_2016_sims/',tree_type,'_time/')
ratepath <- paste0('./to_etal_2016_sims/',tree_type,'/',scenario,'/',tree_type,'_rate/')
subpath <- paste0('./to_etal_2016_sims/',tree_type,'/',scenario,'/',tree_type,'_subs/')
timeinfpath <-  paste0('./to_etal_2016_sims/',tree_type,'/',scenario,'/',tree_type,'_time_inferred/')
rateinfpath <-  paste0('./to_etal_2016_sims/',tree_type,'/',scenario,'/',tree_type,'_rate_inferred/')
powerpath <-  paste0('./to_etal_2016_sims/',tree_type,'/',scenario,'/PowerAnalysis/')
errorpath <-  paste0('./to_etal_2016_sims/',tree_type,'/',scenario,'/Errors/')

# Create paths (time should already exist)
dir.create(ratepath, recursive = TRUE)
dir.create(subpath, recursive = TRUE)
dir.create(timeinfpath, recursive = TRUE)
dir.create(rateinfpath, recursive = TRUE)
dir.create(powerpath, recursive = TRUE)
dir.create(errorpath, recursive = TRUE)

####################################################################
################## Other Initializations ###########################
####################################################################

## Initialize number of cores dependent on system
ncores = floor(0.75*(parallel::detectCores()))

# Number of independent trees to conduct power analysis
power_n <- 100

file_list1 <- list.files(path = timepath, pattern = '.nwk')
file_list <- file_list1[1:power_n]

# Determine array length
v_t <- read.tree(paste0(timepath,file_list[1]))
tlen <- length(v_t$edge.length)
for (i in 2:length(file_list)) {
  v_t <- read.tree(paste0(timepath,file_list[i]))
  tlen2 <- length(v_t$edge.length)
  if (tlen2 > tlen) {
    tlen <- tlen2
  }
}

# Sequence length
nsites <- 1000

# Power array initialization
pwr_calc_array <- array(NaN,dim=c(length(klist),7,2))

pwr_arr_mean <- c()
pwr_arr_ratio <- c()
pwr_arr_tstat <- c()
pwr_arr_cv <- c()
pwr_arr_KST <- c()
pwr_arr_MWWT <- c()
pwr_arr_DT <- c()

# Effect size counter
kind = 1 


####################################################################
############ Parallelized Power Analysis ###########################
####################################################################

#### Running in parallel
registerDoFuture()
plan(multisession, workers = ncores)

for (k in klist) {

y <- foreach (
  j = c(1:power_n), # Number of trees
  .combine = 'cbind',
  .options.future = list(seed = TRUE),
  .errorhandling = 'pass'
) %dofuture% {
  
  ## Libraries re-downloaded for each parallel computation
  library(ape)
  library(treedater)
  library(phytools)
  #library(phylopomp) # epidemiological sims
  library(matrixStats)
  library(tryCatchLog)
  library(diptest)
  library(stats)

  error_file <- paste0(errorpath,"error_dump_",k,"_",j,".rds")

  v_t <- read.tree(paste0(timepath,tree_type,'_time_',j,'.nwk'))

  tree_n = length(v_t$tip.label) + v_t$Nnode - 1
  stopifnot(tree_n == length(v_t$edge.length))
  
  # Get deme edges
  d <- edg_sep(v_t)
  deme1_edges <- d[2][[1]]
  deme0_edges <- d[1][[1]]
  


  ############################################ 
  ##### Simulate rates and substitutions #####
  ############################################
  
  # Reproducibility of each parallel run
  set.seed(12*j)

  ### Uncorrelated Gamma-Poisson model for rate variation ###

  # Rate tree
  v_r <- v_t
  v_r$edge.length <- rgamma(tree_n, shape = a, scale = b) # Uncorrelated non-additive

  # Simulate structured rate variation
  if (k != 1) {
     v_r$edge.length[deme1_edges] <- rgamma(length(deme1_edges), shape = a, scale = b*k)
  }

  # Substitution tree
  v_p1 <- v_t
  v_p1$edge.length <- rpois(tree_n, v_r$edge.length * v_t$edge.length * nsites) / nsites # rate*time*nsites
  
  
  # Write trees
  write.tree(v_r,paste0(ratepath,tree_type,"_k",k,"_",j,"_rate.nwk"))
  write.tree(v_p1,paste0(subpath,tree_type,"_k",k,"_",j,"_subs.nwk"))


  ################################################### 
  ##### Date tree and calculate test statistics #####
  ###################################################

  ## try-catch-log for error-logging
  tryCatchLog({
   
   ### Time-calibrate tree using clock (treedater) ###

   # st <- Sys.time() - assess clock runtime
    v_est <- sim_date(v_t,
                      v_p1,
                      clock_type = 'uncorr',
                      sl = nsites,
                      ncpus = 1) 
    # et <- Sys.time()

    # Inferred rate tree
    v_est_r <- v_est
    v_est_r$edge.length <- v_est$omegas
    
    # Save trees 
    save(v_est,file=paste0(timeinfpath,tree_type,'_k',k,'_',j,'_time_inf.RData'))
    write.tree(v_est,paste0(timeinfpath,tree_type,'_k',k,'_',j,'_time_inf.nwk'))
    write.tree(v_est_r,paste0(rateinfpath,tree_type,'_k',k,'_',j,'_rate_inf.nwk'))



    ### Ancestral state reconstruction - Mk model ###
 
    t_bh <- v_p1
    t_bh <- ladderize(t_bh,right=FALSE)
    
    names_bh <- t_bh$tip.label
    t_1_labs <- grep('_1_',names_bh)
    state_arr <- rep(1,length(names_bh))
    state_arr[t_1_labs] <- 2
    asr <- asr_mk_model(t_bh, 
                    tip_states = state_arr, 
                    rate_model = "ER", 
                    Nstates = 2, 
                    include_ancestral_likelihoods = TRUE) 
    acl2 <- asr$ancestral_likelihoods
    node_labs <- rownames(acl2)
    acl <- asr$ancestral_states
    acl <- data.frame(acl)
    rownames(acl) <- node_labs
    acl$acl[acl$acl == 1] <- 0
    acl$acl[acl$acl == 2] <- 1

    names_0_repl_ind <- which(unname(acl==0))
    names_1_repl_ind <- which(unname(acl==1))

    t_bh_anc <- t_bh
    t_bh_anc$node.label[names_0_repl_ind] <- gsub('_1_', '_0_', t_bh_anc$node.label[names_0_repl_ind])
    t_bh_anc$node.label[names_1_repl_ind] <- gsub('_0_', '_1_', t_bh_anc$node.label[names_1_repl_ind])

    # Obtain empirical deme branches
    d_asr <- edg_sep(t_bh_anc)
    deme1_edges_asr <- d_asr[2][[1]]
    deme0_edges_asr <- d_asr[1][[1]]

    v_p_deme1_edges <- deme1_edges_asr
    v_p_deme0_edges <- deme0_edges_asr



    ### Parametric bootstrap ###

    pbn <- 10000 # Number of replicates for bootstrap
    nullrates <- array(dim = c(tlen,pbn))
    for (i in 1:pbn) {
        
        # Per-branch null distribution of rates
	nullrates[1:length(v_t$edge.length),i] <- rgamma(length(v_t$edge.length), shape = v_est$r, scale = v_est$theta)/nsites

    }
      


    ### Test Statistics ###

    ## If no rate variation inferred ##
    strict_inf <- sum((v_est$coef_of_variation)*(mean(v_est$omegas)) < (0.001/2/3))    
    
    ### Location-based statistics ###

    # Difference of mean rates
    rate1_sim_mean <- mean(v_est_r$edge.length[v_p_deme1_edges])
    rate0_sim_mean <- mean(v_est_r$edge.length[v_p_deme0_edges])
    simdiff_mean <- rate1_sim_mean - rate0_sim_mean
      
    nulldiff_mean <- colMeans(nullrates[v_p_deme1_edges,]) - colMeans(nullrates[v_p_deme0_edges,])

    alpha_val_mean_025_rt <- min(sort(nulldiff_mean,decreasing=TRUE)[1:(0.025*10000)])
    alpha_val_mean_025_lft <- max(sort(nulldiff_mean,decreasing=FALSE)[1:(0.025*10000)])

    acc_rej_mean <- sum((simdiff_mean >= alpha_val_mean_025_rt) || (simdiff_mean <= alpha_val_mean_025_lft))
      
    # Ratio of mean rates
    rate1_sim_mean <- mean(v_est_r$edge.length[v_p_deme1_edges])
    rate0_sim_mean <- mean(v_est_r$edge.length[v_p_deme0_edges])
    simdiff_ratio <- rate1_sim_mean/rate0_sim_mean

    nulldiff_ratio <- colMeans(nullrates[v_p_deme1_edges,])/colMeans(nullrates[v_p_deme0_edges,])

    alpha_val_ratio_025_rt <- min(sort(nulldiff_ratio,decreasing=TRUE)[1:(0.025*10000)])
    alpha_val_ratio_025_lft <- max(sort(nulldiff_ratio,decreasing=FALSE)[1:(0.025*10000)])

    acc_rej_ratio <- sum((simdiff_ratio >= alpha_val_ratio_025_rt) || (simdiff_ratio <= alpha_val_ratio_025_lft))

    # Difference in CV of mean rates
    rate1_sim_mean <- mean(v_est_r$edge.length[v_p_deme1_edges])
    rate1_sim_sd <- sd(v_est_r$edge.length[v_p_deme1_edges])
    rate0_sim_mean <- mean(v_est_r$edge.length[v_p_deme0_edges])
    rate0_sim_sd <- sd(v_est_r$edge.length[v_p_deme0_edges])
    simdiff_cv <- (rate1_sim_sd/rate1_sim_mean)/(rate0_sim_sd/rate0_sim_mean)

    nulldiff_cv <- (colSds(nullrates[v_p_deme1_edges,])/colMeans(nullrates[v_p_deme1_edges,]))/(colSds(nullrates[v_p_deme0_edges,])/colMeans(nullrates[v_p_deme0_edges,]))

    alpha_val_cv_025_rt <- min(sort(nulldiff_cv,decreasing=TRUE)[1:(0.025*10000)])
    alpha_val_cv_025_lft <- max(sort(nulldiff_cv,decreasing=FALSE)[1:(0.025*10000)])

    acc_rej_cv <- sum((simdiff_cv >= alpha_val_cv_025_rt) || (simdiff_cv <= alpha_val_cv_025_lft))

    # Welch's t-statistic 
    rate1_sim_mean <- mean(v_est_r$edge.length[v_p_deme1_edges])
    rate0_sim_mean <- mean(v_est_r$edge.length[v_p_deme0_edges])
    simdiff_mean <- rate1_sim_mean - rate0_sim_mean

    rate1_sim_ts_sd <- var(v_est_r$edge.length[v_p_deme1_edges])/length(v_p_deme1_edges)
    rate0_sim_ts_sd <- var(v_est_r$edge.length[v_p_deme0_edges])/length(v_p_deme0_edges)
    simdiff_ts_sd <- sqrt(rate1_sim_ts_sd + rate0_sim_ts_sd)

    simdiff_ts <- simdiff_mean/simdiff_ts_sd

    nulldiff_mean <- colMeans(nullrates[v_p_deme1_edges,]) - colMeans(nullrates[v_p_deme0_edges,])
 
    nulldiff_ts_sd_1 <- colVars(nullrates[v_p_deme1_edges,])/length(v_p_deme1_edges)
    nulldiff_ts_sd_0 <- colVars(nullrates[v_p_deme0_edges,])/length(v_p_deme0_edges)
    nulldiff_ts_sd <- sqrt(nulldiff_ts_sd_1 + nulldiff_ts_sd_0)

    nulldiff_ts <- nulldiff_mean/nulldiff_ts_sd

    alpha_val_ts_025_rt <- min(sort(nulldiff_ts,decreasing=TRUE)[1:(0.025*10000)])
    alpha_val_ts_025_lft <- max(sort(nulldiff_ts,decreasing=FALSE)[1:(0.025*10000)])

    acc_rej_tstat <- sum((simdiff_ts >= alpha_val_ts_025_rt) || (simdiff_ts <= alpha_val_ts_025_lft))


    ### Distribution-based statistics ###

    # Hartigan's dip test
    DT <- dip.test(v_est_r$edge.length)
    acc_rej_DT <- sum(DT$p.value < 0.05)

    # Two-sample Kolmogorov-Smirnov test
    KST <- ks.test(v_est_r$edge.length[v_p_deme1_edges],v_est_r$edge.length[v_p_deme0_edges],simulate.p.value = TRUE)
    acc_rej_KST <- sum(KST$p.value < 0.05)
    
    # Two-sample Mann-Whitney-Wilcoxon test
    MWWT <- wilcox.test(v_est_r$edge.length[v_p_deme1_edges],v_est_r$edge.length[v_p_deme0_edges])
    acc_rej_MWWT <- sum(MWWT$p.value < 0.05)


    # Return all values
    return(c(j, acc_rej_mean, acc_rej_ratio, acc_rej_tstat, acc_rej_cv,
             acc_rej_KST, acc_rej_MWWT, acc_rej_DT, strict_inf))
      

    # Error handling with custom function
    }, execution.context.msg = "",
       write.error.dump.file = FALSE,
       silent.warnings = TRUE,
       silent.messages = TRUE,
       include.full.call.stack = FALSE,
    error = function(e) {
       error_handler(e, error_file)
    })

  }
  
  # All return values 
  acc_rej_array <- y 
  
  ################################################### 
  ############## Power Calculations #################
  ###################################################
  
  # Save accept-reject file for effect size
  save(y, acc_rej_array, nsites, a, b, file = paste0(powerpath,tree_type,"_",scenario,"_",k,"_accrejall_w2side_",format(Sys.time(), format = "%Y-%m-%d"),".RData"))
  
  # Remove NaN and error indices from consideration
  err_list <- c()
  for (i in 1:length(acc_rej_array[,1])-1) {
     err_list1 <- which(is.na(get_numbers(acc_rej_array[i,])))
     err_list <- c(err_list,err_list1)
  }
  err_list <- unique(sort(err_list))
  arr_cols <- 1:length(acc_rej_array[1,])
  
  # Indices without NaN and errors
  if (sum(err_list) > 0) {
    ok_list1 <- arr_cols[-err_list]
  } else {
    ok_list1 <- arr_cols
  }

  ok_list <- ok_list1[get_numbers(acc_rej_array[9,ok_list1]) == 0]

  ### Power calculation arrays - total accept and total runs ###

  pwr_calc_array[kind,1,1] <- sum(get_numbers(acc_rej_array[2,ok_list]))
  pwr_calc_array[kind,1,2] <- sum(get_numbers(acc_rej_array[2,ok_list1]) >= 0)
  pwr_calc_array[kind,2,1] <- sum(get_numbers(acc_rej_array[3,ok_list]))
  pwr_calc_array[kind,2,2] <- sum(get_numbers(acc_rej_array[3,ok_list1]) >= 0)
  pwr_calc_array[kind,3,1] <- sum(get_numbers(acc_rej_array[4,ok_list]))
  pwr_calc_array[kind,3,2] <- sum(get_numbers(acc_rej_array[4,ok_list1]) >= 0)
  pwr_calc_array[kind,4,1] <- sum(get_numbers(acc_rej_array[5,ok_list]))
  pwr_calc_array[kind,4,2] <- sum(get_numbers(acc_rej_array[5,ok_list1]) >= 0)
  pwr_calc_array[kind,5,1] <- sum(get_numbers(acc_rej_array[6,ok_list]))
  pwr_calc_array[kind,5,2] <- sum(get_numbers(acc_rej_array[6,ok_list1]) >= 0)
  pwr_calc_array[kind,6,1] <- sum(get_numbers(acc_rej_array[7,ok_list]))
  pwr_calc_array[kind,6,2] <- sum(get_numbers(acc_rej_array[7,ok_list1]) >= 0)
  pwr_calc_array[kind,7,1] <- sum(get_numbers(acc_rej_array[8,ok_list]))
  pwr_calc_array[kind,7,2] <- sum(get_numbers(acc_rej_array[8,ok_list1]) >= 0)

  ### Power calculation arrays - individual test statistic power values ###

  power_mean <- sum(get_numbers(acc_rej_array[2,ok_list])) / sum(get_numbers(acc_rej_array[2,ok_list1]) >= 0)
  power_ratio <- sum(get_numbers(acc_rej_array[3,ok_list])) / sum(get_numbers(acc_rej_array[3,ok_list1]) >= 0)
  power_tstat <- sum(get_numbers(acc_rej_array[4,ok_list])) / sum(get_numbers(acc_rej_array[4,ok_list1]) >= 0)
  power_cv <- sum(get_numbers(acc_rej_array[5,ok_list])) / sum(get_numbers(acc_rej_array[5,ok_list1]) >= 0)
  power_KST <- sum(get_numbers(acc_rej_array[6,ok_list])) / sum(get_numbers(acc_rej_array[6,ok_list1]) >= 0)
  power_MWWT <- sum(get_numbers(acc_rej_array[7,ok_list])) / sum(get_numbers(acc_rej_array[7,ok_list1]) >= 0)
  power_DT <- sum(get_numbers(acc_rej_array[8,ok_list])) / sum(get_numbers(acc_rej_array[8,ok_list1]) >= 0)

  # Fill in full power analysis array for each effect size
  pwr_arr_mean[kind] <- power_mean
  pwr_arr_ratio[kind] <- power_ratio
  pwr_arr_tstat[kind] <- power_tstat
  pwr_arr_cv[kind] <- power_cv
  pwr_arr_KST[kind] <- power_KST
  pwr_arr_MWWT[kind] <- power_MWWT
  pwr_arr_DT[kind] <- power_DT

  kind = kind + 1
}

# Stop parallel cluster, free up workers, print report on memory usage 
plan(sequential)
gc(verbose=TRUE,full=TRUE)


endTime <- Sys.time()
timeElapsed <- endTime - startTime



### Save power analysis file
save(
  klist,
  pwr_calc_array,
  pwr_arr_mean,
  pwr_arr_ratio,
  pwr_arr_tstat,
  pwr_arr_cv,
  pwr_arr_KST,
  pwr_arr_MWWT,
  pwr_arr_DT,
  timeElapsed,
  nsites,
  a,
  b,
  file = paste0(powerpath,tree_type,"_",scenario,
    "_",klist[1],"to",tail(klist, n = 1),"_pwrarrall_w2side_",
    format(Sys.time(), format = "%Y-%m-%d"),".RData")
)

# Print session info
sessionInfo()

