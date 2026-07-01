#!/usr/bin/env Rscript

### Detect Structured Rates ###
##
## Template script to detect structured rate heterogeneity per developed methodology
##
## Implement in CLI
##
## Required inputs: 
## Phylogram/substitution tree (simulated or empirical), and corresponding dated tree
## Tree must have tip labels, with two demes represented 
##

####################################################################

# Install necessary packages
#install.packages('ape')
#install.packages('phytools')
#install.packages('matrixStats')
#install.packages('castor')
#install.packages('diptest')
#install.packages('optparse')
#install.packages('tryCatchLog')

# Load packages
suppressMessages({
  suppressWarnings({
  library(ape)
  library(phytools)
  library(matrixStats)
  library(castor)
  library(diptest)
  library(optparse)
  library(tryCatchLog)
  })
})

####################################################################
##### Initialization ##### 

arg_opts <- list(
  make_option(c("--tp"), type = "character", help = "Directory of input substitution tree (newick)"),
  make_option(c("--td"), type = "character", help = "Directory of dated tree (newick)"),
  make_option(c("--tr"), type = "character", help = "Directory of input rate tree (newick)"),
  make_option(c("--dm1"), type = "character", help = "Deme 1 label (one word)"),
  make_option(c("--dm2"), type = "character", help = "Deme 2 label (one word)"),
  make_option(c("--nsites"), type = "integer", help = "Number of sites in alignment"),
  make_option(c("--shape"), type = "integer", help = "Gamma shape parameter inferred for tree lineage rates"),
  make_option(c("--scale"), type = "integer", help = "Gamma scale parameter inferred for tree lineage rates")
  )

arg_list <- parse_args(OptionParser(option_list = arg_opts))

### Input check
if (is.null(arg_list$tp)) stop('Input substitution tree required')
if (!file.exists(arg_list$tp)) stop('Input substitution tree not found')
if (is.null(arg_list$td)) stop('Input dated tree required')
if (!file.exists(arg_list$td)) stop('Input dated tree not found')
if (is.null(arg_list$tr)) stop('Input rate tree required')
if (!file.exists(arg_list$tr)) stop('Input rate tree not found')
if (is.null(arg_list$dm1)) stop('Deme 1 label required')
if (is.null(arg_list$dm2)) stop('Deme 2 label required')
if (is.null(arg_list$nsites)) stop('Number of sites required')
if (is.null(arg_list$shape)) stop('Gamma shape required')
if (is.null(arg_list$scale)) stop('Gamma scale required')


# Random seed for reproducibility
rand_val <- 1234
set.seed(rand_val)


####################################################################
##### Detect structured rate variation ##### 
withCallingHandlers({

  ### Function to separate tree based on deme information  ### 
  # ENSURE: DEME 0 is labeled with '_0_' in node and tip label, and DEME 1 is labeled with '_1_' in node and tip label 
  
  edg_sep <- function(phylotree_unp) {
    phylotree_unp <- ladderize(phylotree_unp,right=FALSE)
    
    # Get node labels for demes (hosts)
    h2_n <- grep('_1_',phylotree_unp$node.label)
    h2_n <- h2_n + length(phylotree_unp$tip.label)
    t2_n <- grep('_1_',phylotree_unp$tip.label)
    
    # Isolate the relevant edges
    edg_r <- list()
    r <- c()
    for (v in 1:(length(phylotree_unp$node.label)+length(phylotree_unp$tip.label))) {
      if (v %in% h2_n) { # If deme 1
        #print(v)
        r <- c()
        r1 <- which(phylotree_unp$edge[, 2] == v) # Tip
        r2 <- which(phylotree_unp$edge[, 1] == v) # Node
        for (i in 1:(length(r2))) { # Across nodes...
          # If either tip as node has a _1_ or tip as tip has a _1_
          if ((sum(h2_n %in% phylotree_unp$edge[r2[i], 2]) + sum(t2_n %in% phylotree_unp$edge[r2[i], 2])) > 0) { 
            r <- c(r,r1,r2[i]) # Counted with both node and tip
          } else {
            r <- r1
          }# If tip has a _0_ 
        } 
      }  
      if (length(r) == 0)
        r <- NA
      edg_r[[length(edg_r) + length(r)]] <- r
    }
    
    
    # Isolate the relevant edges - ensure all relevant tips are included (sometimes not leading from a node counted above...)
    edg_r2 <- list()
    r2 <- c()
    for (v2 in 1:(length(phylotree_unp$tip.label))) {
      if (v2 %in% t2_n) {
        r2 <- which(phylotree_unp$edge[,2] == v2)
        edg_r2[[length(edg_r2)+length(r2)]] <- r2
        #print(r2)
      }
    }
    
    # function to flatten edge arrays, adapted from Stack Overflow comment (https://stackoverflow.com/a/72056420)
    flattenArray <- function(x) {
      if (is.data.frame(x)) return(list(x))
      if (!is.list(x)) return(x)
      unlist(lapply(x, flattenArray), FALSE)
    }
    
    edg_rf <- flattenArray((edg_r))
    edg_rf2 <- flattenArray((edg_r2))
    # All relevant edges isolated for deme 1
    edg_rf3 <- sort(na.omit(unique(c(edg_rf,edg_rf2))))
    
    deme1_edges <- edg_rf3
    deme0_edges <- which(!((1:(length(phylotree_unp$node.label)+length(phylotree_unp$tip.label)-1)) %in% edg_rf3))
    
    d <- cbind(list(deme0_edges),list(deme1_edges))
    
    return(d)
    
  }
  

  ### Read inputs ### 
  
  # Hosts - note the naming convention
  deme1 = trimws(arg_list$dm1) # labeled as 1
  deme2 = trimws(arg_list$dm2) # labeled as 0

  ## Load substitution tree (phylogram) - simulated or empirical
  ## Also load dated tree (chronogram)
  t_bh <- read.tree(arg_list$tp)
  v_est <- read.tree(arg_list$td)
  
  # Ratogram - if treedater-dated. If not, provide rate-tree as well
  v_est_r <- read.tree(arg_list$tr)

  # Length of tree
  tlen <- length(v_est$edge.length)
  
  # Bootstrap parameters
  nsites <- arg_list$nsites
  shape <- arg_list$shape
  scale <- arg_list$scale

  
  ### Ancestral state reconstruction ### 
  
  t_bh <- ladderize(t_bh, right = FALSE)
  names_bh <- t_bh$tip.label
  t_1_labs <- grep(deme2, names_bh)
  state_arr <- rep(1, length(names_bh))
  state_arr[t_1_labs] <- 2
  asr <- asr_mk_model(
    t_bh,
    tip_states = state_arr,
    rate_model = "ER",
    Nstates = 2,
    include_ancestral_likelihoods = TRUE
  )
  acl <- asr$ancestral_states
  
  t_bh_anc <- t_bh
  nodelabs <- asr$ancestral_states
  labs <- sprintf("node_%d", 1:length(nodelabs))
  nodelabs[nodelabs == 1] <- unlist(lapply(labs[which(nodelabs == 1)] , function(x)
    paste0('t_', deme1, '_0_', x)))
  nodelabs[nodelabs == 2] <- unlist(lapply(labs[which(nodelabs == 2)] , function(x)
    paste0('t_', deme2, '_1_', x)))
  t_bh_anc$node.label <- nodelabs
  t_bh_anc$tip.label <- gsub(deme1, paste0(deme1,'_0_'), t_bh_anc$tip.label)
  t_bh_anc$tip.label <- gsub(deme2, paste0(deme2,'_1_'), t_bh_anc$tip.label)
  
  d_asr <- edg_sep(t_bh_anc)
  deme0_edges_asr <- d_asr[2][[1]]
  deme1_edges_asr <- d_asr[1][[1]]
  
  v_p_deme0_edges <- deme0_edges_asr
  v_p_deme1_edges <- deme1_edges_asr
  

  
  ### Parametric bootstrap ###
  suppressWarnings({
  pbn <- 10000 # Number of replicates for bootstrap
  nullrates <- array(dim = c(tlen, pbn))
  set.seed(rand_val)
  for (i in 1:pbn) {
    # Per-branch null distribution of rates
    nullrates[1:tlen, i] <- rgamma(tlen,shape = shape,scale = scale)/nsites
    
  }
  })
  
  
  ### Test Statistics ###
  
  ### Location-based statistics ###
  
  # Difference of mean rates
  rate1_sim_mean <- mean(v_est_r$edge.length[v_p_deme1_edges])
  rate0_sim_mean <- mean(v_est_r$edge.length[v_p_deme0_edges])
  simdiff_mean <- rate1_sim_mean - rate0_sim_mean
  
  nulldiff_mean <- colMeans(nullrates[v_p_deme1_edges, ]) - colMeans(nullrates[v_p_deme0_edges, ])
  
  alpha_val_mean_025_rt <- min(sort(nulldiff_mean, decreasing = TRUE)[1:(0.025*pbn)])
  alpha_val_mean_025_lft <- max(sort(nulldiff_mean, decreasing = FALSE)[1:(0.025*pbn)])
  
  acc_rej_mean <- sum((simdiff_mean >= alpha_val_mean_025_rt) ||
                        (simdiff_mean <= alpha_val_mean_025_lft)
  )
  
  # Ratio of mean rates
  rate1_sim_mean <- mean(v_est_r$edge.length[v_p_deme1_edges])
  rate0_sim_mean <- mean(v_est_r$edge.length[v_p_deme0_edges])
  simdiff_ratio <- rate1_sim_mean / rate0_sim_mean
  
  nulldiff_ratio <- colMeans(nullrates[v_p_deme1_edges, ]) / colMeans(nullrates[v_p_deme0_edges, ])
  
  alpha_val_ratio_025_rt <- min(sort(nulldiff_ratio, decreasing = TRUE)[1:(0.025*pbn)])
  alpha_val_ratio_025_lft <- max(sort(nulldiff_ratio, decreasing = FALSE)[1:(0.025*pbn)])
  
  acc_rej_ratio <- sum((simdiff_ratio >= alpha_val_ratio_025_rt) ||
                         (simdiff_ratio <= alpha_val_ratio_025_lft)
  )
  
  # Difference in CV of mean rates
  rate1_sim_mean <- mean(v_est_r$edge.length[v_p_deme1_edges])
  rate1_sim_sd <- sd(v_est_r$edge.length[v_p_deme1_edges])
  rate0_sim_mean <- mean(v_est_r$edge.length[v_p_deme0_edges])
  rate0_sim_sd <- sd(v_est_r$edge.length[v_p_deme0_edges])
  simdiff_cv <- (rate1_sim_sd / rate1_sim_mean) / (rate0_sim_sd / rate0_sim_mean)
  
  nulldiff_cv <- (colSds(nullrates[v_p_deme1_edges, ]) / colMeans(nullrates[v_p_deme1_edges, ])) /
    (colSds(nullrates[v_p_deme0_edges, ]) / colMeans(nullrates[v_p_deme0_edges, ]))
  
  alpha_val_cv_025_rt <- min(sort(nulldiff_cv, decreasing = TRUE)[1:(0.025*pbn)])
  alpha_val_cv_025_lft <- max(sort(nulldiff_cv, decreasing = FALSE)[1:(0.025*pbn)])
  
  acc_rej_cv <- sum((simdiff_cv >= alpha_val_cv_025_rt) ||
                      (simdiff_cv <= alpha_val_cv_025_lft))
  
  # Welch's t-statistic
  rate1_sim_mean <- mean(v_est_r$edge.length[v_p_deme1_edges])
  rate0_sim_mean <- mean(v_est_r$edge.length[v_p_deme0_edges])
  simdiff_mean <- rate1_sim_mean - rate0_sim_mean
  
  rate1_sim_ts_sd <- var(v_est_r$edge.length[v_p_deme1_edges]) / length(v_p_deme1_edges)
  rate0_sim_ts_sd <- var(v_est_r$edge.length[v_p_deme0_edges]) / length(v_p_deme0_edges)
  simdiff_ts_sd <- sqrt(rate1_sim_ts_sd + rate0_sim_ts_sd)
  
  simdiff_ts <- simdiff_mean / simdiff_ts_sd
  
  nulldiff_mean <- colMeans(nullrates[v_p_deme1_edges, ]) - colMeans(nullrates[v_p_deme0_edges, ])
  
  nulldiff_ts_sd_1 <- colVars(nullrates[v_p_deme1_edges, ]) / length(v_p_deme1_edges)
  nulldiff_ts_sd_0 <- colVars(nullrates[v_p_deme0_edges, ]) / length(v_p_deme0_edges)
  nulldiff_ts_sd <- sqrt(nulldiff_ts_sd_1 + nulldiff_ts_sd_0)
  
  nulldiff_ts <- nulldiff_mean / nulldiff_ts_sd
  
  alpha_val_ts_025_rt <- min(sort(nulldiff_ts, decreasing = TRUE)[1:(0.025*pbn)])
  alpha_val_ts_025_lft <- max(sort(nulldiff_ts, decreasing = FALSE)[1:(0.025*pbn)])
  
  acc_rej_tstat <- sum((simdiff_ts >= alpha_val_ts_025_rt) ||
                         (simdiff_ts <= alpha_val_ts_025_lft))
  
  
  
  ### Distribution-based statistics ###
  
  # Hartigan's dip test
  set.seed(rand_val)
  DT <- dip.test(v_est_r$edge.length)
  acc_rej_DT <- sum(DT$p.value < 0.05)
  
  # Two-sample Kolmogorov-Smirnov test
  set.seed(rand_val)
  KST <- ks.test(v_est_r$edge.length[v_p_deme0_edges],
                 v_est_r$edge.length[v_p_deme1_edges],
                 simulate.p.value = TRUE)
  acc_rej_KST <- sum(KST$p.value < 0.05)
  
  # Two-sample Mann-Whitney-Wilcoxon test
  set.seed(rand_val)
  MWWT <- wilcox.test(v_est_r$edge.length[v_p_deme0_edges], v_est_r$edge.length[v_p_deme1_edges])
  acc_rej_MWWT <- sum(MWWT$p.value < 0.05)
  
  
  
  
  ####################################################################
  ##### Post-Analysis #####
  
  ### P-values for all test statistics calculated ###
  p_mean <- 2 * min(sum(sort(nulldiff_mean, decreasing = TRUE) < simdiff_mean)/pbn,
                    sum(sort(nulldiff_mean, decreasing = TRUE) > simdiff_mean)/pbn)
  p_ratio <- 2 * min(sum(sort(nulldiff_ratio, decreasing = TRUE) < simdiff_ratio)/pbn,
                     sum(sort(nulldiff_ratio, decreasing = TRUE) > simdiff_ratio)/pbn)
  p_tstat <- 2 * min(sum(sort(nulldiff_ts, decreasing = TRUE) < simdiff_ts)/pbn,
                     sum(sort(nulldiff_ts, decreasing = TRUE) > simdiff_ts)/pbn)
  p_KST <- KST$p.value
  p_MW <- MWWT$p.value
  ## Poor-performing test statistics
  p_cv <- sum(sort(nulldiff_cv, decreasing = TRUE) < simdiff_cv)/pbn
  p_DT <- DT$p.value
  
  ### Print results to console ###
  
  sig <- function(pval,sigl = 0.05) {
    ifelse(pval < sigl, paste0(pval, "***"),pval)
  }
  
  cat("\nDetection of structured rate variation between",deme1, "and",deme2,"in this phylogeny\n") 
  cat("Test statistic p-values are listed, with significance at 0.05 marked\n")
  cat("\nMean-difference test statistic\n",sig(p_mean), "\n\n",sep="")
  cat("Mean-ratio test statistic\n",sig(p_ratio), "\n\n",sep="")
  cat("Welch\'s t-test statistic\n",sig(p_tstat), "\n\n",sep="")
  cat("Kolmogorov-Smirnov test statistic\n",sig(p_KST), "\n\n",sep="") 
  cat("Mann-Whitney-Wilcoxon test statistic\n",sig(p_MW), "\n\n",sep="")

### Error handling
}, warning = function(warn) {
  message("Warning thrown: ", conditionMessage(warn))
  invokeRestart("muffleWarning")
}, error = function(err) {
  message("Error thrown: ", conditionMessage(err))
})

