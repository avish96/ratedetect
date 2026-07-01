## Manuscript - Figure 1
## Example simulated phylopomp phylogeny, and simulated and inferred rates
##

library(ape)
library(castor)

stopifnot(packageVersion('ape')>='5.8')
stopifnot(packageVersion('castor')>='1.8.3')

source('sim_funcs.R')

t_t <- read.tree('../simulations_SANN/phylopomp_results/phylopomp_2host_time_1.nwk')
t_t <- ladderize(t_t,right=FALSE)

d <- edg_sep(t_t)
deme0_edges <- d[1][[1]]
deme1_edges <- d[2][[1]]

## PLOTTING
rate_tree_colors <- c()
for (i in (1:(length(t_t$node.label)+length(t_t$tip.label)))) {
  if (i %in% deme0_edges) {
    rate_tree_colors[i] <- '#BA65F1'
  } else {
    rate_tree_colors[i] <- '#17F2FD'
  }
}

plot(t_t, show.tip.label = FALSE, show.node.label = FALSE, cex=0.5, edge.color=rate_tree_colors) 
legend('topleft', legend=c('Avian','Human'),bty='n',col=c('#BA65F1','#17F2FD'),pch=19, cex=0.8)

t_p <- read.tree('../simulations_SANN/phylopomp_results/phylopomp_2host_k10_1_subs.nwk')
t_r <- read.tree('../simulations_SANN/phylopomp_results/phylopomp_2host_k10_1_rate.nwk')
t_r_inf <- read.tree('../simulations_SANN/phylopomp_results/phylopomp_2host_k10_1_rate_inf.nwk')

# Mk model ancestral state reconstruction
t_bh <- t_p
t_bh <- ladderize(t_bh,right=FALSE)

names_bh <- t_bh$tip.label
t_1_labs <- grep('_1_',names_bh)
state_arr <- rep(1,length(names_bh))
state_arr[t_1_labs] <- 2
asr <- asr_mk_model(t_bh, 
                    tip_states = state_arr, 
                    rate_model = 'ER',
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

# Empirical deme branches
d_asr <- edg_sep(t_bh_anc)
deme0_edges_asr <- d_asr[1][[1]]
deme1_edges_asr <- d_asr[2][[1]]

data_list <- list(Vector1 = t_r$edge.length[deme0_edges], 
                  Vector2 = t_r$edge.length[deme1_edges], 
                  Vector3 = t_r_inf$edge.length[deme0_edges_asr], 
                  Vector4 = t_r_inf$edge.length[deme1_edges_asr])

png('../simulations_SANN/figures_main/Figure2a.png', width = 6000*.3, height = 7000*.3, res=300)
par(mfrow=c(1,1))
phylotree_unp_yr <- t_t
phylotree_unp_yr$edge.length <- phylotree_unp_yr$edge.length/365.25
plot(phylotree_unp_yr, show.tip.label = FALSE, show.node.label = FALSE, cex=0.5, edge.color=rate_tree_colors) # edge.color(phylotree_unp,edg_rf3,c('crowngroup','stemgroup'),'red'))
legend('topleft', legend=c('Avian','Human'),bty='n',col=c('#BA65F1','#17F2FD'),pch=19, cex=1.5)
add.scale.bar(length=1,x=0.5,y=500)
dev.off()

png('../simulations_SANN/figures_main/Figure2b.png', width = 10000*.3, height = 7000*.3, res=300)
boxplot(data_list,
        #main = 'Boxplot of rate distribution in simulated vs. inferred tree',
        ylab = 'Rate (# substitutions/site/yr)',
        names = c('Low-rate','High-rate', 'Low-rate', 'High-rate'),
        #ylim = c(0,0.08),
        cex.lab=1.5,
        cex.axis=1.5)
dev.off()

