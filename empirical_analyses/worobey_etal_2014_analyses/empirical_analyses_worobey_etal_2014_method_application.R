## Empirical Analyses
## Worobey et al. (2014) dataset
## Analysis script - Apply method and print results
##
library(ape)
library(lubridate)
library(treedater)
library(phytools)
library(castor)
library(matrixStats)
library(diptest)

stopifnot(packageVersion("ape")>="5.8")
stopifnot(packageVersion("lubridate")>="1.9.4")
stopifnot(packageVersion("treedater")>="1.0.2")
stopifnot(packageVersion("castor")>="1.8.3")
stopifnot(packageVersion("phytools")>="2.3.0")
stopifnot(packageVersion("matrixStats")>="1.5.0")
stopifnot(packageVersion("diptest")>="0.77.1")

# Set to script directory
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

source('../../scripts/sim_funcs.R')
# Random seed for reproducibility
rand_val <- 1234
set.seed(rand_val)

sl <- 1701 # sequence length

## Load substitution tree (phylogram) - simulated or empirical
## Also load dated tree (chronogram)

load('worobey_etal_2014_trees.RData')

t_bh <- t
v_est <- t_dated

# Ratogram - if treedater-dated. If not, provide rate-tree as well
v_est_r <- t_bh
v_est_r$edge.length <- v_est$omegas

# Hosts - note the naming convention
deme1 = 'Avian' # labeled as 1
deme2 = 'Equine' # labeled as 0

# Length of tree
tlen <- length(v_est$edge.length)

# Bootstrap parameters
nsites <- sl
shape <- v_est$r
scale <- v_est$theta

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
deme0_edges_asr <- d_asr[1][[1]]
deme1_edges_asr <- d_asr[2][[1]]

v_p_deme0_edges <- deme0_edges_asr
v_p_deme1_edges <- deme1_edges_asr



### Parametric bootstrap ###
suppressWarnings({
  pbn <- 10000 # Number of replicates for bootstrap
  nullrates <- array(dim = c(tlen, pbn))
  #set.seed(rand_val)
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
#set.seed(rand_val)
DT <- dip.test(v_est_r$edge.length)
acc_rej_DT <- sum(DT$p.value < 0.05)

# Two-sample Kolmogorov-Smirnov test
#set.seed(rand_val)
KST <- ks.test(v_est_r$edge.length[v_p_deme1_edges],
               v_est_r$edge.length[v_p_deme0_edges],
               simulate.p.value = TRUE)
acc_rej_KST <- sum(KST$p.value < 0.05)

# Two-sample Mann-Whitney-Wilcoxon test
#set.seed(rand_val)
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

