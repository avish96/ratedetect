## Script to translate genealogies/transmission trees from phylopomp to phylogenetic trees
##

####################################################################
### Initialization

st <- Sys.time()

#install.packages("phylopomp",repos="https://kingaa.github.io") # install from repo for most up-to-date version
#install.packages('ape')
#install.packages('phytools')
#install.packages('stringr')
#install.packages('readr')

# Load libraries
library(ape)
library(phytools)
library(phylopomp)
library(stringr)
library(readr)

stopifnot(packageVersion("ape")>="5.8")
stopifnot(packageVersion("phytools")>="2.3.0")
stopifnot(packageVersion("phylopomp")>="0.14.8.0")
stopifnot(packageVersion("stringr")>="1.5.1")
stopifnot(packageVersion("readr")>="2.1.5")

# Load any required functions
source('sim_funcs.R')

wd <- getwd()
treefilepath <- './phylopomp_2host_transmissiontrees/'
treefilepath2 <- './phylopomp_2host_trees_time/'
dir.create(treefilepath2,recursive=TRUE)

####################################################################

### Make tree if not done yet
# source('phylopomp_treemaker.R') 


### Read and save phylopomp trees as newick files
file_list <- list.files(path = treefilepath, pattern = '.txt')

for (j in 1:length(file_list)) {

  v_t <- read.tree(paste0(treefilepath,file_list[j]))
  v_t <- ladderize(v_t, right = FALSE)

  write.tree(v_t,paste0(treefilepath,strsplit(file_list[j],'txt'),'nwk'))


  ### PLOTTING ###

  #d <- edg_sep(v_t)
  #deme0_edges <- d[2][[1]]
  #deme1_edges <- d[1][[1]]

  #phylotree_unp <- ladderize(v_t,right=FALSE)
  #rate_tree_colors <- c()
  ##for (i in (1:((length(phylotree_unp$node.label)+length(phylotree_unp$tip.label))-1))) { ## Ensure node labels exist
  #for (i in (1:length(phylotree_unp$edge.length))) {
  #  if (i %in% deme1_edges) {
  #    rate_tree_colors[i] <- "#BA65F1"
  #  } else {
  #    rate_tree_colors[i] <- "#17F2FD"
  #  }
  #}
  #
  #plot(phylotree_unp, show.tip.label = FALSE, show.node.label = FALSE, cex=0.5, edge.color=rate_tree_colors) 
  #


}  


### Remove in-line nodes in genealogies
### Python script for this 
### RUN THIS COMMAND IN TREEFILEPATH DIRECTORY
### IF SYSTEM FAILS TO WORK, RUN DIRECTLY THROUGH LINUX COMMAND LINE IN THE TREEFILEPATH DIRECTORY 
###
setwd(treefilepath)
output <- system(paste0("python ../phylopomp_inline_remover.py ."),intern=TRUE) 
print(output)
setwd(wd)

#### Save phylogeny as nwk file
file_list_timetree <- list.files(path=treefilepath,pattern='skip')

for (index1 in 1:length(file_list_timetree)) {
  gnn_skip1 <- read_file(paste0(treefilepath,file_list_timetree[index1]))
  gnn_skip <- system(paste0("echo ",paste0("'",gnn_skip1,"'"), " | sed 's/^(\\(.*\\));$/\\1/' "),intern=TRUE) # to remove the root and add a semicolon
  gnn_skip <- paste0(gnn_skip,';')
  
  phylotree_skip <- read.tree(text=gnn_skip)
  
  write.tree(phylotree_skip,file=paste0(treefilepath2,'phylopomp_2host_time_',index1,'.nwk'))
}

# Print session info
sessionInfo()

