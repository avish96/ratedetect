## Empirical Analyses
## H3 dataset from GISAID
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
########################### Main analyses ################################
##########################################################################

load('../H3_analyses/evenly_sampled_(no_passage_filter)/H3_even_trees.RData')

t_bh <- H3_tree_even_dated$intree
t_bh <- ladderize(t_bh,right=FALSE)

names_bh <- t_bh$tip.label
t_1_labs <- grep('Swine',names_bh)
state_arr <- rep(1,length(names_bh))
state_arr[t_1_labs] <- 2
asr <- asr_mk_model(t_bh, 
                    tip_states = state_arr, 
                    rate_model = "ER",
                    Nstates = 2, 
                    include_ancestral_likelihoods = TRUE) 
acl <- asr$ancestral_states

t_bh_anc <- t_bh
nodelabs <- asr$ancestral_states
labs <- sprintf("node_%d", 1:length(nodelabs))
nodelabs[nodelabs == 1] <- unlist(lapply(labs[which(nodelabs == 1)] , function(x)
  paste0('t_Avian_0_', x)))
nodelabs[nodelabs == 2] <- unlist(lapply(labs[which(nodelabs == 2)] , function(x)
  paste0('t_Swine_1_', x)))
t_bh_anc$node.label <- nodelabs
t_bh_anc$tip.label <- gsub('Avian', 'Avian_0_', t_bh_anc$tip.label)
t_bh_anc$tip.label <- gsub('Swine', 'Swine_1_', t_bh_anc$tip.label)

d_asr <- edg_sep(t_bh_anc)
deme0_edges_asr <- d_asr[1][[1]]
deme1_edges_asr <- d_asr[2][[1]]

v_p_deme0_edges <- deme0_edges_asr
v_p_deme1_edges <- deme1_edges_asr

##########################################################################
############################## FIGURES ###################################
##########################################################################

#####################################
######### FIGURE 6 A ################
#####################################
par(mfrow=c(1,1))
#par(font.axis = 2, font.lab = 2)
png("Figure6a.png", width = 6000*.3, height = 7000*.3, res=300)
cols <- rep('#BA65F1', length(H3_tree_even_dated$edge.length))
cols[v_p_deme1_edges] <- 'pink'
max_time <- max(node.depth.edgelength(H3_tree_even_dated))
plot(ladderize(H3_tree_even_dated,right=F),show.tip.label = F,edge.col=cols,edge.width=1.5,
     main="H3 HA Avian-Swine phylogeny, UCGP",cex.main=1,
     no.margin = FALSE, x.lim = c(0,max_time))
manual_positions <- c(0, 7.5, 15,46)
manual_labels <- rev(c(2025,2017,2010,1979))
axis(1, at = manual_positions, labels = manual_labels, line = 0, cex.axis = 1.2) #font=2)
legend('topleft', legend=c('Avian','Swine'),col=c('#BA65F1','pink'),pch=19, cex=1.5,bty='n')
dev.off()

#####################################
######### FIGURE 6 B, C #############
#####################################
png("Figure6bc.png", width = 5000*.3, height = 7000*.3, res=300)
par(mfrow=c(2,1))
#par(font.axis = 2, font.lab = 2)
h1 <- hist(H3_tree_even_dated$omegas[v_p_deme1_edges],breaks=100,plot=F) # 750
h2 <- hist(H3_tree_even_dated$omegas[v_p_deme0_edges],breaks=50,plot=F) # 500
plot(h1,col='pink',,xaxt='n',yaxt='n',xlim=c(0,0.006),ylim=c(0,400),
     xlab='Substitution rate (# subs/site/yr)',
     ylab='Frequency',
     main='Inferred substitution rates - Full tree',
     cex.main=1)
plot(h2,col='#BA65F1',add=T)
ticks_to_label <- c(0,0.001,0.002,0.003, 0.004, 0.005, 0.006)
axis(1, at = ticks_to_label, labels = c(0,format(ticks_to_label[2:7], scientific = TRUE))) #,font=2)
axis(2, at = c(0,200,400),labels=c(0,200,400)) #,font=2)
legend('topright', legend=c('Avian','Swine'),col=c('#BA65F1','pink'),pch=19, 
       cex=1)

h3 <- hist(H3_swine_dated_even$omegas,breaks=100,plot=F)
h4 <- hist(H3_avian_dated_even$omegas,breaks=100,plot=F) #60000
plot(h3,col='pink',xaxt='n',yaxt='n',xlim=c(0,0.006),ylim=c(0,400),
     xlab='Substitution rate (# subs/site/yr)',
     ylab='Frequency',
     main='Inferred substitution rates - Host-specific',
     cex.main=1)
plot(h4,col='#BA65F1',add=T)
ticks_to_label <- c(0, 0.001,0.002,0.003, 0.004, 0.005,0.006)
axis(1, at = ticks_to_label, labels = c(0,format(ticks_to_label[2:7], scientific = TRUE))) #,font=2)
axis(2, at = c(0,200,400),labels=c(0,200,400)) #,font=2)
legend('topright', legend=c('Avian','Swine'),col=c('#BA65F1','pink'),pch=19, 
       cex=1)
dev.off()

##########################################################################
########################### Supp analyses ################################
##########################################################################

load('../H3_analyses/unevenly_sampled/H3_uneven_trees.RData')

t_bh <- H3_tree_uneven_dated$intree
t_bh <- ladderize(t_bh,right=FALSE)

names_bh <- t_bh$tip.label
t_1_labs <- grep('Swine',names_bh)
state_arr <- rep(1,length(names_bh))
state_arr[t_1_labs] <- 2
asr <- asr_mk_model(t_bh, 
                    tip_states = state_arr, 
                    rate_model = "ER",
                    Nstates = 2, 
                    include_ancestral_likelihoods = TRUE) 
acl <- asr$ancestral_states

t_bh_anc <- t_bh
nodelabs <- asr$ancestral_states
labs <- sprintf("node_%d", 1:length(nodelabs))
nodelabs[nodelabs == 1] <- unlist(lapply(labs[which(nodelabs == 1)] , function(x)
  paste0('t_Avian_0_', x)))
nodelabs[nodelabs == 2] <- unlist(lapply(labs[which(nodelabs == 2)] , function(x)
  paste0('t_Swine_1_', x)))
t_bh_anc$node.label <- nodelabs
t_bh_anc$tip.label <- gsub('Avian', 'Avian_0_', t_bh_anc$tip.label)
t_bh_anc$tip.label <- gsub('Swine', 'Swine_1_', t_bh_anc$tip.label)

d_asr <- edg_sep(t_bh_anc)
deme0_edges_asr <- d_asr[1][[1]]
deme1_edges_asr <- d_asr[2][[1]]

v_p_deme0_edges <- deme0_edges_asr
v_p_deme1_edges <- deme1_edges_asr

#####################################
########### FIGURE S15 A ############
#####################################

par(mfrow=c(1,1))
#par(font.axis = 2, font.lab = 2)
png("FigureS15a.png", width = 6000*.3, height = 7000*.3, res=300)
cols <- rep('#BA65F1', length(H3_tree_uneven_dated$edge.length))
cols[v_p_deme1_edges] <- 'pink'
max_time <- max(node.depth.edgelength(H3_tree_uneven_dated))
plot(ladderize(H3_tree_uneven_dated,right=F),show.tip.label = F,edge.col=cols,edge.width=1.5,
     main="H3 HA Avian-Swine phylogeny, UCGP",cex.main=1,
     no.margin = FALSE, x.lim = c(0,max_time))
manual_positions <- c(0, 7.5, 15,40.5)
manual_labels <- rev(c(2025,2017,2010,1985))
axis(1, at = manual_positions, labels = manual_labels, line = 0, cex.axis = 1.2) #font=2)
legend('topleft', legend=c('Avian','Swine'),col=c('#BA65F1','pink'),pch=19, cex=1.5,bty='n')
dev.off()

#####################################
######### FIGURE S15 B, C ###########
#####################################
png("FigureS15bc.png", width = 5000*.3, height = 7000*.3, res=300)
#par(font.axis = 2, font.lab = 2)
par(mfrow=c(2,1))
h1 <- hist(H3_tree_uneven_dated$omegas[v_p_deme1_edges],breaks=750,plot=F) # 750
h2 <- hist(H3_tree_uneven_dated$omegas[v_p_deme0_edges],breaks=500,plot=F) # 500
plot(h1,col='pink',,xaxt='n',yaxt='n',xlim=c(0,0.005),ylim=c(0,600),
     xlab='Substitution rate (# subs/site/yr)',
     ylab='Frequency',
     main='Inferred substitution rates - Full tree',
     cex.main=1)
plot(h2,col='#BA65F1',add=T)
ticks_to_label <- c(0, 0.001,0.002,0.003, 0.004, 0.005)
axis(1, at = ticks_to_label, labels = c(0,format(ticks_to_label[2:6], scientific = TRUE)))
axis(2, at = c(0,300,600),labels=c(0,300,600))
legend('top', legend=c('Avian','Swine'),col=c('#BA65F1','pink'),pch=19, 
       cex=1)

h3 <- hist(H3_swine_dated_uneven$omegas,breaks=200,plot=F)
h4 <- hist(H3_avian_dated_uneven$omegas,breaks=60000,plot=F) #60000
plot(h3,col='pink',xaxt='n',yaxt='n',xlim=c(0,0.005),ylim=c(0,600),
     xlab='Substitution rate (# subs/site/yr)',
     ylab='Frequency',
     main='Inferred substitution rates - Host-specific',
     cex.main=1)
plot(h4,col='#BA65F1',add=T)
ticks_to_label <- c(0, 0.001,0.002,0.003, 0.004, 0.005)
axis(1, at = ticks_to_label, labels = c(0,format(ticks_to_label[2:6], scientific = TRUE)))
axis(2, at = c(0,300,600),labels=c(0,300,600))
legend('topright', legend=c('Avian','Swine'),col=c('#BA65F1','pink'),pch=19, 
       cex=1)
dev.off()
