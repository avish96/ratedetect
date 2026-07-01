## Empirical Analyses
## H3 dataset from GISAID
## Analysis script - tree dating
##

library(ape)
library(treedater)
library(lubridate)
library(castor)
library(phytools)
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

source('../scripts/sim_funcs.R')

# Unevenly sampled

H3_tree_uneven <- read.tree('./unevenly_sampled/RAxML_bestTree.H3AS_LBA_T_unrooted.nwk')
names <- read.csv('./unevenly_sampled/gisaid_epiflu_sequence_H3_HA_Avian_Swine_lengthdupidfilt_aligned_trimmed_edit_nohang_LBA_names.csv',header=F)
times <- read.csv('./unevenly_sampled/gisaid_epiflu_sequence_H3_HA_Avian_Swine_lengthdupidfilt_aligned_trimmed_edit_nohang_LBA_dates.csv',header=F)

times <- as.numeric(ymd(noquote(times[,1])))/365.25
names <- paste0("'",names[,1],"'") # With quotes

sts <- setNames(times,names)
sl <- 1701

start <- Sys.time()
H3_tree_uneven_dated <- dater(H3_tree_uneven,sts=sts,s=sl,clock="uncorr",ncpu = 30)
end <- Sys.time()

uneven_dating_time <- end - start

H3_avian_uneven <- drop.tip(H3_tree,H3_tree$tip.label[grep('Swine',H3_tree$tip.label)])
H3_swine_uneven <- drop.tip(H3_tree,H3_tree$tip.label[grep('Avian',H3_tree$tip.label)])
H3_avian_dated_uneven <- dater(H3_avian_uneven,sts=sts,s=sl,clock='uncorr',ncpu = 30)
H3_swine_dated_uneven <- dater(H3_swine_uneven,sts=sts,s=sl,clock='uncorr',ncpu = 30)

save(H3_tree_uneven,H3_tree_uneven_dated,H3_avian_uneven,H3_swine_uneven,H3_avian_dated_uneven,H3_swine_dated_uneven,file='./unevenly_sampled/H3_uneven_trees.RData')

# Evenly sampled (data not passage-filtered)

H3_tree_even <- read.tree('./evenly_sampled_(no_passage_filter)/RAxML_bestTree.H3AS_nopass_samp_LBA_T_unrooted.nwk')
names2 <- read.csv('./evenly_sampled_(no_passage_filter)/gisaid_epiflu_sequence_H3_HA_Avian_Swine_nopass_lengthdupidfilt_sampled_aligned_trimmed_edit_nohang_LBA_names.csv',header=F)
times2 <- read.csv('./evenly_sampled_(no_passage_filter)/gisaid_epiflu_sequence_H3_HA_Avian_Swine_nopass_lengthdupidfilt_sampled_aligned_trimmed_edit_nohang_LBA_dates.csv',header=F)

times2 <- as.numeric(ymd(noquote(times2[,1])))/365.25
names2 <- paste0("'",names2[,1],"'") # With quotes

sts2 <- setNames(times2,names2)
sl <- 1701

start2 <- Sys.time()
H3_tree_even_dated <- dater(H3_tree_even,sts=sts2,s=sl,clock="uncorr",ncpu = 30)
end2 <- Sys.time()

dating_time_even <- end2 - start2

H3_avian_even <- drop.tip(H3_tree,H3_tree$tip.label[grep('Swine',H3_tree$tip.label)])
H3_swine_even <- drop.tip(H3_tree,H3_tree$tip.label[grep('Avian',H3_tree$tip.label)])
H3_avian_dated_even <- dater(H3_avian_even,sts=sts2,s=sl,clock='uncorr',ncpu = 30)
H3_swine_dated_even <- dater(H3_swine_even,sts=sts2,s=sl,clock='uncorr',ncpu = 30)

save(H3_tree_even,H3_tree_even_dated,H3_avian_even,H3_swine_even,H3_avian_dated_even,H3_swine_dated_even,file='./evenly_sampled_(no_passage_filter)/H3_even_trees.RData')

