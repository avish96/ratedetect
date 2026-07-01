## Empirical Analyses
## Worobey et al. (2014) dataset
## Figure plotting script
##

library(ape)
library(lubridate)
library(treedater)
library(phytools)
library(castor)

stopifnot(packageVersion("ape")>="5.8")
stopifnot(packageVersion("lubridate")>="1.9.4")
stopifnot(packageVersion("treedater")>="1.0.2")
stopifnot(packageVersion("castor")>="1.8.3")
stopifnot(packageVersion("phytools")>="2.3.0")

source('../../scripts/sim_funcs.R')

# Set to script directory
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

##########################################################################
########################### Load analyses ################################
##########################################################################

load('../worobey_etal_2014_analyses/worobey_etal_2014_trees.RData')

t_bh <- t
t_bh <- ladderize(t_bh,right=FALSE)

names_bh <- t_bh$tip.label
t_1_labs <- grep('Equine',names_bh)
state_arr <- rep(1,length(names_bh))
state_arr[t_1_labs] <- 2
asr <- asr_mk_model(t_bh, 
                    tip_states = state_arr, 
                    rate_model = 'ER',
                    Nstates = 2, 
                    include_ancestral_likelihoods = TRUE) 
acl <- asr$ancestral_states

t_bh_anc <- t_bh
nodelabs <- asr$ancestral_states
labs <- sprintf("node_%d", 1:length(nodelabs))
nodelabs[nodelabs == 1] <- unlist(lapply(labs[which(nodelabs == 1)] , function(x)
  paste0('t_Avian_0_', x)))
nodelabs[nodelabs == 2] <- unlist(lapply(labs[which(nodelabs == 2)] , function(x)
  paste0('t_Equine_1_', x)))
t_bh_anc$node.label <- nodelabs
t_bh_anc$tip.label <- gsub('Avian', 'Avian_0_', t_bh_anc$tip.label)
t_bh_anc$tip.label <- gsub('Equine', 'Equine_1_', t_bh_anc$tip.label)

# Empirical deme branches
d_asr <- edg_sep(t_bh_anc)
deme0_edges_asr <- d_asr[1][[1]]
deme1_edges_asr <- d_asr[2][[1]]

v_p_deme0_edges <- deme0_edges_asr
v_p_deme1_edges <- deme1_edges_asr

##########################################################################
############################## FIGURES ###################################
##########################################################################

#####################################
######### FIGURE 5 A, B #############
#####################################
png("Figure5ab.png", width = 10000*.3, height = 7000*.3, res=300)
wtree <- read.tree('../worobey_etal_2014_analyses/worobey_etal_2014_original_H7_datasets/H7.fnl.2_rates.mcc.mod.nwk')

edgelabels_eq <- 1:31
edgelabels_noteq <- which(!((1:length(t_dated$edge[,1])) %in% edgelabels_eq))

par(mfrow=c(1,2))
cols <- rep('#BA65F1', length(wtree$edge.length))
cols[edgelabels_eq] <- 'orange'
max_time <- max(node.depth.edgelength(wtree))
plot(ladderize(wtree,right=F),show.tip.label = F,edge.col=cols,edge.width=1.5,
     main='H7 HA Worobey et al. (2014) + HSLC', cex.main=1,
     no.margin = FALSE, x.lim = c(0,max_time))
manual_positions <- c(0, 34, 59,109,159)
manual_labels <- rev(c(2009,1975,1950,1900,1850))
axis(1, at = manual_positions, labels = manual_labels, line = 0, cex.axis = 1)
legend('topleft', legend=c('Avian','Equine'),col=c('#BA65F1','orange'),
       pch=19, cex=1.5,bty='n')#box.lty = 1,box.lwd = 1,box.col = "black",bg = "white")

cols <- rep('#BA65F1', length(t_dated$edge.length))
cols[v_p_deme1_edges] <- 'orange'
plot(ladderize(t_dated,right=F),show.tip.label = F,edge.col=cols,edge.width=1.5,
     main='H7 HA RAxML + URC',cex.main=1,
     no.margin = FALSE, x.lim = c(0,max_time))
manual_positions <- c(0, 34, 59,109,159)
manual_labels <- rev(c(2009,1975,1950,1900,1850))
axis(1, at = manual_positions, labels = manual_labels, line = 0, cex.axis = 1)
legend('topleft', legend=c('Avian','Equine'),col=c('#BA65F1','orange'),
       pch=19, cex=1.5,bty='n')#box.lty = 1,box.lwd = 1,box.col = "black",bg = "white")
dev.off()

#####################################
######### FIGURE 5 A - Inset ########
#####################################
# Worobey reproduction
png("Figure5ab_inset.png", width = 10000*.3, height = 7000*.3, res=300)
par(mfrow=c(1,1))
h3 <- hist(rgamma(1000,shape=t_dated_avian$r*60,scale=t_dated_avian$theta/sl/60),breaks=25)
h4 <- hist(rgamma(1000,shape=t_dated_equine$r*1.25*2,scale=t_dated_equine$theta/sl/2),breaks=25)
plot(h3,xaxt='n',yaxt='n',xlim=c(0,0.005),col='#BA65F1',
     xlab='',#Substitution rate (# subs/site/yr)',
     ylab='',
     main='')#H7 substitution rate, treedater uncorrelated non-additive')
plot(h4,col='orange',add=T)
ticks_to_label <- c(0.0005,0.001, 0.002, 0.003, 0.004,0.0045)
axis(side = 1, at = ticks_to_label, labels = NA, col = "black", lwd = 1,tcl = 0)
axis(1, at = ticks_to_label[2:5], labels = format(ticks_to_label[2:5], scientific = TRUE),
     cex.axis=2)
legend(0.004,100, legend=c('Avian','Equine'),col=c('#BA65F1','orange'),pch=19, 
       cex=1,inset=0.2)
dev.off()

#####################################
######### FIGURE 5 C, D #############
#####################################
png("Figure5cd.png", width = 6000*.3, height = 7000*.3, res=300)
par(mfrow=c(2,1))
par(cex.lab=1.25)
h1 <- hist(t_dated$omegas[v_p_deme0_edges],breaks=40,plot=F)
h2 <- hist(t_dated$omegas[v_p_deme1_edges],breaks=50,plot=F)
plot(h1,xaxt='n',yaxt='n',xlim=c(0,0.009),ylim=c(0,30),col='#BA65F1',
     xlab='Substitution rate (# subs/site/yr)',
     ylab='Frequency',
     main = 'Inferred substitution rates - Full tree',
     cex.main=1.25)
#main='H7 substitution rate, treedater uncorrelated non-additive')
plot(h2,breaks=100,col='orange',add=T)
ticks_to_label <- c(0,0.001,0.002,0.003,0.004,0.005,0.006,0.0085)
axis(side = 1, at = ticks_to_label, labels = NA, col = "black", lwd = 1,tcl = 0)
axis(1, at = ticks_to_label[1:7], labels = c(0,format(ticks_to_label[2:7], scientific = TRUE)))
axis(2, at = c(0,10,20,30),labels=c(0,10,20,30))
legend('topright', legend=c('Avian','Equine'),col=c('#BA65F1','orange'),pch=19, 
       cex=1,inset=0.2)

h3 <- hist(t_dated_avian$omegas,breaks=70,plot=F)
h4 <- hist(t_dated_equine$omegas,breaks=7,plot=F)
plot(h3,xaxt='n',yaxt='n',xlim=c(0,0.009),ylim=c(0,30),col='#BA65F1',
     xlab='Substitution rate (# subs/site/yr)',
     ylab='Frequency',
     main='Inferred substitution rates - Host-specific',
     cex.main=1.25)
plot(h4,breaks=100,col='orange',add=T)
ticks_to_label <- c(0,0.001,0.002,0.003,0.004,0.005,0.008,0.0085)
axis(side = 1, at = ticks_to_label, labels = NA, col = "black", lwd = 1,tcl = 0)
axis(1, at = ticks_to_label[1:7], labels = c(0,format(ticks_to_label[2:7], scientific = TRUE)))
axis(2, at = c(0,10,20,30),labels=c(0,10,20,30))
legend('topright', legend=c('Avian','Equine'),col=c('#BA65F1','orange'),pch=19, 
       cex=1,inset=0.2)
dev.off()
