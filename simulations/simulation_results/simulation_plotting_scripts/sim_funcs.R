## Script containing essential functions used in multiple analysis scripts
##

####################################################################

#install.packages('ape')
#install.packages('phytools')
#install.packages('treedater')
#install.packages('stats')

## Libraries to load
library(ape)
library(phytools)
library(treedater)
library(stats)

stopifnot(packageVersion("ape")>="5.8")
stopifnot(packageVersion("phytools")>="2.3.0")
stopifnot(packageVersion("treedater")>="1.0.2")
stopifnot(packageVersion("stats")>="4.3.3")

##### Function to log error and dump files
error_handler <- function(e, error_file) {
  # Record error timestamp
  timestamp <- Sys.time()
  pidstamp <- Sys.getpid()

  # Save error information, time, and current environment to a dump file
  dump.frames(to.file = FALSE)
  saveRDS(list(
	error = conditionMessage(e), 
	call = deparse(conditionCall(e)), 
	pid = pidstamp, 
	timestamp = timestamp, 
	last.dump = last.dump), 
  file = error_file)

 # Stop worker upon error
 stop()

}



##### Function for dating the tree using treedater in simulation
sim_date <- function(v_t, v_p, clock_type, sl, ncpus) {
  # Getting root-tip distance for all tree tips
  times <-
    c(diag(vcv(v_t)))
  names <- v_t$tip.label
  
  # time-name matrix
  sts <- setNames(times, names)  
  
  # Date tree
  t_dated <- dater(
    v_p,
    sts = sts,
    clock = clock_type,
    s = sl,
    ncpu = ncpus,
  ) 
  
  return(t_dated)
  
}



##### Function to separate tree based on deme information
##### ENSURE: DEME 0 is labeled with '_0_' in node and tip label, and DEME 1 is labeled with '_1_' in node and tip label 
edg_sep <- function(phylotree_unp) {
  phylotree_unp <- ladderize(phylotree_unp,right=FALSE)
  
  # Get node labels for demes (hosts)
  h2_n <- grep('_1_',phylotree_unp$node.label)
  h2_n <- h2_n + length(phylotree_unp$tip.label)
  t2_n <- grep('_1_',phylotree_unp$tip.label)
  
  # Isolate the relevant edges
  edg_r <- list()
  r <- c()
  for (v in 1:(length(phylotree_unp$node.label)+length(phylotree_unp$tip.label))) {
    if (v %in% h2_n) { # If deme 1
      #print(v)
      r <- c()
      r1 <- which(phylotree_unp$edge[, 2] == v) # Tip
      r2 <- which(phylotree_unp$edge[, 1] == v) # Node
      for (i in 1:(length(r2))) { # Across nodes...
        # If either tip as node has a _1_ or tip as tip has a _1_
        if ((sum(h2_n %in% phylotree_unp$edge[r2[i], 2]) + sum(t2_n %in% phylotree_unp$edge[r2[i], 2])) > 0) { 
          r <- c(r,r1,r2[i]) # Counted with both node and tip
        } else {
          r <- r1
        }# If tip has a _0_ 
      } 
    }  
    if (length(r) == 0)
      r <- NA
    edg_r[[length(edg_r) + length(r)]] <- r
  }
  
  
  # Isolate the relevant edges - ensure all relevant tips are included (sometimes not leading from a node counted above...)
  edg_r2 <- list()
  r2 <- c()
  for (v2 in 1:(length(phylotree_unp$tip.label))) {
    if (v2 %in% t2_n) {
      r2 <- which(phylotree_unp$edge[,2] == v2)
      edg_r2[[length(edg_r2)+length(r2)]] <- r2
      #print(r2)
    }
  }
  
  # function to flatten edge arrays, adapted from Stack Overflow comment (https://stackoverflow.com/a/72056420)
  flattenArray <- function(x) {
    if (is.data.frame(x)) return(list(x))
    if (!is.list(x)) return(x)
    unlist(lapply(x, flattenArray), FALSE)
  }
  
  edg_rf <- flattenArray((edg_r))
  edg_rf2 <- flattenArray((edg_r2))
  # All relevant edges isolated for deme 1
  edg_rf3 <- sort(na.omit(unique(c(edg_rf,edg_rf2))))
  
  deme1_edges <- edg_rf3
  deme0_edges <- which(!((1:(length(phylotree_unp$node.label)+length(phylotree_unp$tip.label)-1)) %in% edg_rf3))
  
  d <- cbind(list(deme0_edges),list(deme1_edges))

  return(d)
  
}



##### Getting non-numeric elements in acc_rej array
get_numbers <- function(X) {
  X[toupper(format(X,scientific=FALSE)) != tolower(format(X,scientific=FALSE))] <- NA
  return(as.double(as.character(X)))
}
