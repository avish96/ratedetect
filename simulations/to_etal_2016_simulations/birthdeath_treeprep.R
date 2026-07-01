## Script used for preparing To et al. (2016) birth-death trees for power analyses
## Trees were used in To et al. (2016) and Volz and Frost (2017)
## Link to original dataset: http://www.atgc-montpellier.fr/LSD/
##

####################################################################
# Pre-requisites before running this script

# Download dataset to the right directory - ensure paths are correct

datapath = './../../data/True unrooted trees with dates/'

####################################################################

# Install necessary packages

#install.packages('ape')
#install.packages('phytools')
#install.packages('castor')

# Load packages

library(ape)
library(phytools)
library(castor)

stopifnot(packageVersion("ape")>="5.8")
stopifnot(packageVersion("phytools")>="2.3.0")
stopifnot(packageVersion("castor")>="1.8.3")

# Load required functions
source('sim_funcs.R')

# Obtain rates for rate matrix from phylopomp tree results
get_label <- function(tree, node) {
  if (node <= Ntip(tree)) {
    return(tree$tip.label[node])
  } else {
    return(tree$node.label[node - Ntip(tree)])
  }
}

deme_shifts <- c()
shifts_per_brl <- c()
pp_dir <- './../phylopomp_simulations/phylopomp_2host_trees_time' # directory with phylopomp time trees
file_list <- list.files(path = pp_dir)
for (j in 1:length(file_list)) {
  v_t <- read.tree(paste0(pp_dir,'/',file_list[j]))
  v_t <- ladderize(v_t, right = FALSE)
  edge_labels <- apply(v_t$edge, 2, function(col) {
    sapply(col, function(node) get_label(v_t, node))
  })
  edge_deme <- array(dim=c(length(edge_labels[,1]),2))
  edge_deme[grep('_0_',edge_labels[,1]),1] <- 0
  edge_deme[grep('_1_',edge_labels[,1]),1] <- 1
  edge_deme[grep('_0_',edge_labels[,2]),2] <- 0
  edge_deme[grep('_1_',edge_labels[,2]),2] <- 1
  deme_shifts[j] <- sum(abs(edge_deme[,1]-edge_deme[,2]))
  shifts_per_brl[j] <- deme_shifts[j]/(sum(v_t$edge.length/365.25))
}

pp_rate <- mean(shifts_per_brl,digits=2) 

####################################################################

# Read trees downloaded from To et al. (2016)
d1 <- read.tree(paste0(datapath,'Strict molecular clock trees/D750_11_10_out.tree'))
d2 <- read.tree(paste0(datapath,'Strict molecular clock trees/D750_3_25_out.tree'))
d3 <- read.tree(paste0(datapath,'Strict molecular clock trees/D995_11_10_out.tree'))
d4 <- read.tree(paste0(datapath,'Strict molecular clock trees/D995_3_25_out.tree'))

# Loop across scenarios
for (bdk in 1:4) { 
  ns_arr <- c()

# 100 trees per scenario
  for (j in 1:100) {
    i <- 0
    index11 <- 1
    while (i < 1) {
      if (bdk == 1) {
        v_p <- ladderize(drop.tip(d1[[j]], 'out'), right = F) # substitution tree, outgroup removed
        tree_type <- 'RMC_750_11_10'
	      timepath <- paste0('./to_etal_2016_sims/',tree_type,'_time/')
	      if (!dir.exists(timepath)) {
	        dir.create(timepath, recursive = TRUE)
	      }
      } else if (bdk == 2) {
        v_p <- ladderize(drop.tip(d2[[j]], 'out'), right = F) # substitution tree, outgroup removed
        tree_type <- 'RMC_750_3_25'
        timepath <- paste0('./to_etal_2016_sims/',tree_type,'_time/')
        if (!dir.exists(timepath)) {
          dir.create(timepath, recursive = TRUE)
        }
      } else if (bdk == 3) {
        v_p <- ladderize(drop.tip(d3[[j]], 'out'), right = F) # substitution tree, outgroup removed
        tree_type <- 'RMC_995_11_10'
        timepath <- paste0('./to_etal_2016_sims/',tree_type,'_time/')
        if (!dir.exists(timepath)) {
          dir.create(timepath, recursive = TRUE)
        }
      } else if (bdk == 4) {
        v_p <- ladderize(drop.tip(d4[[j]], 'out'), right = F) # substitution tree, outgroup removed
        tree_type <- 'RMC_995_3_25'
        timepath <- paste0('./to_etal_2016_sims/',tree_type,'_time/')
        if (!dir.exists(timepath)) {
          dir.create(timepath, recursive = TRUE)
        }
      }
      
      v_t <- v_p
      v_t$edge.length <- v_t$edge.length / 0.006 # scaling to underlying time-tree
      
      # Simulate demes and label tips and nodes
      Q = M <- matrix(c(
               -1*pp_rate, pp_rate,
                pp_rate, -1*pp_rate
                ), nrow = 2, byrow = TRUE)

      set.seed((100*(bdk-1) + j)*index11)

      allstates <- simulate_mk_model(v_t, Q)
      deme_tips <- which(v_t$edge[, 2] %in% which(allstates$tip_states == 1))
      
      all_tips = 1:length(v_t$tip.label)
      non_deme_tips <- which(!(all_tips %in% v_t$edge[deme_tips, 2]))
      v_t$tip.label[v_t$edge[deme_tips, 2]] <- unlist(lapply(v_t$tip.label[v_t$edge[deme_tips, 2]] , function(x)
        paste0('t_0_', x)))
      v_t$tip.label[non_deme_tips] <- unlist(lapply(v_t$tip.label[non_deme_tips] , function(x)
        paste0('t_1_', x)))
      
      nodelabs <- allstates$node_states
      labs <- sprintf("node_%d", 1:length(nodelabs))
      nodelabs[nodelabs == 1] <- unlist(lapply(labs[which(nodelabs == 1)] , function(x)
        paste0('t_0_', x)))
      nodelabs[nodelabs == 2] <- unlist(lapply(labs[which(nodelabs == 2)] , function(x)
        paste0('t_1_', x)))
      v_t$node.label <- nodelabs
      
      d <- edg_sep(v_t)
      deme0_edges <- d[1][[1]]
      deme1_edges <- d[2][[1]]
      
      # Save trees if at least two edges represented in each deme
      if ((length(deme0_edges) > 2) && (length(deme1_edges) > 2)) {
        write.tree(v_t,file = paste0(timepath,tree_type,'_time_', j, '.nwk'))
        i <- 1
      } else {
        index11 = index11 + 1
      }
      
      
      ### PLOTTING ###
      #cols <- rep('black', length(v_p$edge.length))
      #cols[deme1_edges] <- 'red'
      #plot(ladderize(v_p,right=F), show.tip.label = F,show.node.label = F, edge.color = cols,cex=0.3)
      
    }
  }
}

# Print session info
sessionInfo()
