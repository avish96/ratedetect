# Avinash Subramanian
# Created April 16, 2025
# To et al. 2016, Volz and Frost 2017 birth-death trees
# Link: http://www.atgc-montpellier.fr/LSD/

# Load packages
library(ape)

stopifnot(packageVersion('ape')>='5.8')

source('sim_funcs.R')

v_t_750_3_25 <- ladderize(read.tree('../simulations_SANN/to_etal_results/RMC_750_3_25_time_1.nwk'),right=F)
v_t_750_11_10 <- ladderize(read.tree('../simulations_SANN/to_etal_results/RMC_750_11_10_time_1.nwk'),right=F)
v_t_995_3_25 <- ladderize(read.tree('../simulations_SANN/to_etal_results/RMC_995_3_25_time_1.nwk'),right=F)
v_t_995_11_10 <- ladderize(read.tree('../simulations_SANN/to_etal_results/RMC_995_11_10_time_1.nwk'),right=F)

wd <- '../simulations_SANN/figures_supp/to_etal_supp/'
png(paste0(wd,'toetal_2016_bd_trees.png'), width = 8000*0.3, height = 6000*0.3,res = 300)

layout_matrix <- matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE)

# Define the heights and widths for layouts
heights <- c(1, 1)  # Equal heights for rows
widths <- c(1, 1)   # Equal widths for columns

# Create a layout with zero inner spacing
layout(layout_matrix, widths = widths, heights = heights)

# Define outer margins, which affect space around the entire plot area
par(oma = c(1, 1, 1, 1))  # No additional outer margins
par(pty = 's')
par(mfrow = c(2, 2))
par(font.axis = 2, font.lab = 2)

par(mar = c(2, 2, 2, 2))  # Top space to reduce proximity to next plot


d <- edg_sep(v_t_750_3_25)
deme0_edges <- d[1][[1]] 
deme1_edges <- d[2][[1]]

cols <- rep('black', length(v_t_750_3_25$edge.length))
cols[deme0_edges] <- 'red'

plot(
  v_t_750_3_25,
  show.tip.label = F,
  show.node.label = F,
  edge.color = cols,
  cex = 0.3,
  main='750_3_25'
)
add.scale.bar(x = 0.660374, y = 67, len=4)

d <- edg_sep(v_t_750_11_10)
deme0_edges <- d[1][[1]]
deme1_edges <- d[2][[1]]

cols <- rep('black', length(v_t_750_11_10$edge.length))
cols[deme0_edges] <- 'red'

plot(
  v_t_750_11_10,
  show.tip.label = F,
  show.node.label = F,
  edge.color = cols,
  cex = 0.3,
  main='750_11_10'
)
add.scale.bar(x = -0.07331178, y = 97.51303, len=2)

d <- edg_sep(v_t_995_3_25)
deme0_edges <- d[1][[1]]
deme1_edges <- d[2][[1]]

cols <- rep('black', length(v_t_995_3_25$edge.length))
cols[deme0_edges] <- 'red'

plot(
  v_t_995_3_25,
  show.tip.label = F,
  show.node.label = F,
  edge.color = cols,
  cex = 0.3,
  main='995_3_25'
)
add.scale.bar(x = 0.660374, y = 67, len=2)

d <- edg_sep(v_t_995_11_10)
deme0_edges <- d[1][[1]]
deme1_edges <- d[2][[1]]

cols <- rep('black', length(v_t_995_11_10$edge.length))
cols[deme0_edges] <- 'red'

plot(
  v_t_995_11_10,
  show.tip.label = F,
  show.node.label = F,
  edge.color = cols,
  cex = 0.3,
  main='995_11_10'
)
add.scale.bar(x = -0.07331178, y = 97.51303, len=2)

dev.off()