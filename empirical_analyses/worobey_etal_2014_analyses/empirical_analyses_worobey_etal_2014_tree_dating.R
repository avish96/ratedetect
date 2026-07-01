## Empirical Analyses
## Worobey et al. (2014) dataset
## Analysis script - tree dating
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

# Sequence lengths
sl = 1716

t <- read.tree('RAxML_result.H7fixedT_rooted.nwk')
times <- read.csv('H7_labels_dates.csv',header=F)
names <- read.csv('H7_labels.csv',header=F)

times <- as.numeric(ymd(paste0((times[,1]),'-07-01')))/365.25
names <- paste0("'",names[,1],"'") # With quotes
#names <- names[,1] # Without quotes

sts <- setNames(times,names)

t_dated <- dater(t,sts=sts,s=sl,clock="uncorr",ncpu = 1) 

###

mrca_equine <- getMRCA(t_dated, t_dated$tip.label[grep('Equine',t_dated$tip.label)])
t_equine <- extract.clade(t,mrca_equine)
mrca_avian <- getMRCA(t_dated, t_dated$tip.label[grep('Avian',t_dated$tip.label)])
t_avian <- extract.clade(t,mrca_avian)

t_dated_equine <- dater(t_equine,sts=sts,s=sl,clock="uncorr",ncpu = 1) 
t_dated_avian <- dater(t_avian,sts=sts,s=sl,clock="uncorr",ncpu = 1) 

save(t,t_dated,t_equine,t_avian,t_dated_equine,t_dated_avian,file='worobey_etal_2014_trees.RData')

