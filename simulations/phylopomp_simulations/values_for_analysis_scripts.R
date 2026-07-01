## Script to calculate various essential values needed in simulation and empirical analysis scripts
##

####################################################################

#install.packages('ape')
#install.packages('phytools')
#install.packages('stats')
#install.packages('seqinr')

### Libraries to install
library(ape)
library(phytools)
library(stats)
library(seqinr)

stopifnot(packageVersion("ape")>="5.8")
stopifnot(packageVersion("phytools")>="2.3.0")
stopifnot(packageVersion("stats")>="4.3.3")
stopifnot(packageVersion("seqinr")>="4.2.36")

####################################################################

### Getting shifts/br length for phylopomp trees for use in To et al. (2016) simulations
### Code in tree-prep script for To et al. (2016) simulations
get_label <- function(tree, node) {
  if (node <= Ntip(tree)) {
    return(tree$tip.label[node])
  } else {
    return(tree$node.label[node - Ntip(tree)])
  }
}
deme_shifts <- c()
shifts_per_brl <- c()
file_list <- list.files(path = '.')
for (j in 1:length(file_list)) {
  v_t <- read.tree(file_list[j])
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

####################################################################

### Code to derive sampling weights for heterogeneous sampling through time analyses

# Set to script directory
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Number of sequences in each fasta file
l_ar <- c()
file_list <- list.files(path='../../data/H3N2_HA_GISAID/',pattern='.fasta') # directory with H3N2 fasta files
for (i in (1:length(file_list))) {
  seq <- read.fasta(paste0('../../data/H3N2_HA_GISAID/',file_list[i]),as.string=T,seqonly=T)
  l_ar[i] <- length(seq)
}

l_ar_2 <- c(0,0,1,0,4,6,11,74,1,61,104,206,230,183,867) # 2000-2014, number of sequences each year (separately manually counted for 2000 - 2008)
l_ar_3 <- c(l_ar_2,l_ar[c(-1,-2,-3,-11,-12)]) # remove 2000-2014 (since already counted), and 2022 due to anomalously high number of samples
#l_ar_3[23] <- sum(l_ar[11:12]) # replace with full 2022 length (summed from original array)

xar <- 1:24

# Sequence data for regression
sample_data <- data.frame(x=xar[11:24],y=l_ar_3[11:24]) # 2010-2024 yearly number of sequences, data frame created
poisson_model <- glm(y~x,data=sample_data,family = 'poisson') # poisson glm

# Create a basic scatterplot
# plot(sample_data$x, sample_data$y,xaxt='n',xlab='Year',ylab='Num Sequences')
# axis(1, at=1:25, labels=c('00','01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24'))

# Array you want to sample from, e.g., original 'x' values
sample_array <- seq(1:14)

# Predict values for the sample array with the fitted model to get weights
predicted_weights <- predict(poisson_model, newdata = data.frame(x = sample_array),type='response')
plot(predicted_weights)
predicted_weights[predicted_weights<0] <- 0

# Normalize weights to sum to 1
normalized_weights <- predicted_weights / max(predicted_weights)
normalized_weights <- signif(normalized_weights,digits=3) # weights used in simulation scripts


####################################################################

### Essential data on simulated trees (written for phylopomp below)

source('sim_funcs.R')

file_list <- list.files(path = '../simulation_results/simulations_SANN/phylopomp_simulations/phylopomp_2host_trees_time/',pattern='.nwk')
num_tips <- c()
hum_edges <- c()
av_edges <- c()
for (j in 1:length(file_list)) {
  v_t <- read.tree(paste0('../simulation_results/simulations_SANN/phylopomp_simulations/phylopomp_2host_trees_time/',file_list[j]))
  num_tips[j] <- Ntip(v_t)
  d <- edg_sep(v_t)
  deme1_edges <- d[2][[1]]
  deme0_edges <- d[1][[1]]
  hum_edges[j] <- length(deme1_edges)
  av_edges[j] <- length(deme0_edges)
}
hum_prop <- hum_edges/(av_edges + hum_edges)
mean(num_tips)
range(num_tips)
quantile(num_tips,probs = c(0.1, 0.9))
mean(hum_prop)
range(hum_prop)
quantile(hum_prop,probs = c(0.1, 0.9))
