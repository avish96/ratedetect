## Script used for generating phylopomp genealogies/transmission trees
##

####################################################################
### Initialization

st <- Sys.time()

#install.packages('ape')
#install.packages('phytools')
#install.packages('stringr')
#install.packages("phylopomp",repos="https://kingaa.github.io") # install from repo for most up-to-date version

# Load libraries
library(ape)
library(phytools)
library(phylopomp)
library(stringr) 

stopifnot(packageVersion("ape")>="5.8")
stopifnot(packageVersion("phytools")>="2.3.0")
stopifnot(packageVersion("phylopomp")>="0.14.8.0")
stopifnot(packageVersion("stringr")>="1.5.1")

# Path to write transmission trees
treefilepath <- './phylopomp_2host_transmissiontrees/'
dir.create(treefilepath, recursive = TRUE)

## Number of trees desired for simulation
ntrees = 200

####################################################################

# Load any required functions
source('sim_funcs.R')

## Parameters for phylopomp simulation -- below for Two Species model
time = 5 * 365.25; t0 = 0
Beta11 = 0.3; Beta12 = 0; Beta21 = 0.3*0.01; Beta22 = 0.01 # transmission rate, or force of infection
gamma1 = 0.1; gamma2 = 0.05 # recovery rate; 1/(gamma + d) = infectious period
R0_1 = Beta11 / gamma1; R0_2 = Beta22 / gamma2 # R0 - basic reproduction number
psi1 = 0.005; psi2 = 0.025 # sampling rate; psi/gamma = odds of being sampled
omega1 = 1 / (2 * 365.25); omega2 = 1 / (0.5 * 365.25) # waning immunity rate
c1 = 1; c2 = 1; # culling probability
b1 = 1 / (7 * 365.25); b2 = 1 / (70 * 365.25) # birth rates
d1 = 1 / (7 * 365.25); d2 = 1 / (70 * 365.25) # death rates; equal to births
iota1 = 0 ; iota2 = 0 # immigration rates
phi_1 = (d1 + omega1) / gamma1; phi_2 = (d2 + omega2) / gamma2 # another sampling parameter

# Initial values for susceptible, infected, and immune/recovered populations for both species
# Calculated using S2I2R2 formulas assuming equilibrium and using final size equation we discussed earlier
N_1 = 1e4; N_2 = 1e4 # Total population size
S1_0 = N_1 / R0_1; S2_0 = N_2 # Susceptible populations
R1_0 = N_1 * ((R0_1 - 1) / (R0_1*(1 + phi_1))); R2_0 = 0  # Recovered/Immune populations
I1_0 = R1_0*phi_1; I2_0 = 0 # Infected populations

# Collect parameters
params <- list(
  time, t0 = 0,
  Beta11, Beta12, Beta21, Beta22,
  gamma1, gamma2,
  R0_1, R0_2,
  psi1, psi2, 
  omega1, omega2,
  c1, c2,
  b1, b2,
  d1, d2, 
  iota1, iota2, 
  N_1, N_2, 
  S1_0, S2_0,
  I1_0, I2_0,
  phi_1, phi_2,
  R1_0, R2_0)


## Simulate ntrees phylopomp trees 

index1 = 1
while (index1 < (ntrees+1)) {

  ## Simulate trees with above parameters
  #startt <- Sys.time()
  params |>
    with(
      phylopomp::simulate(
    "TwoSpecies",
    time = time, t0 = t0,
    Beta11 = Beta11, Beta12 = Beta12, Beta21 = Beta21, Beta22 = Beta22,
    gamma1 = gamma1, gamma2 = gamma2,
    psi1 = psi1, psi2 = psi2,
    c1 = c1, c2 = c2,
    omega1 = omega1, omega2 = omega2, # 6-month-ish waning immunity rate in chickens...
    b1 = b1, b2 = b2,
    d1 = d1, d2 = d2,
    iota1 = iota1, iota2 = iota2,
    S1_0 = S1_0, S2_0 = S2_0,
    I1_0 = I1_0, I2_0 = I2_0,
    R1_0 = R1_0, R2_0 = R2_0)
    ) |>
      freeze(seed=12*index1) -> g
  #endt <- Sys.time()
  #endt-startt
  
  # plot(g,prune=FALSE)
  

  ## Parse all transmission trees simulated
  gy <- yaml(g)
  gn <- phylopomp::newick(g,prune=TRUE,obscure = FALSE) # pruned (all extant taxa removed)
  gn_arr_1 <- strsplit(gn,";")[[1]]

  # Save those trees with at least 2 edges in each deme
  gn_arr <- c()
  for (di in 1:length(gn_arr_1)) {
    if ((str_count(gn_arr_1[di],pattern='b_0_') >= 2) && (str_count(gn_arr_1[di],pattern='b_1_') >= 2)) {
      gn_arr <- append(gn_arr,gn_arr_1[di])
    }
  }
  
  if (length(gn_arr) == 0) {
  } else {
      
      ## Extract largest simulated tree
      ind = 1
      siz1 = nchar(gn_arr[1])
      if (length(gn_arr) == 1){
        gnn <- gn_arr[ind]
      } else {
      for (i in 2:length(gn_arr)) {
        siz2 = nchar(gn_arr[i])
        if (siz2 > siz1) {
          ind = i
          siz1 <- siz2
        }
      }
      }
      gnn <- gn_arr[ind]
      
      ## Write tree 
      gnn <- paste0(gnn,';')
          
      sink(paste0(treefilepath,"phylopomp_genealogy_",index1,".txt")) 
      cat(gnn)
      sink()
      index1 = index1 + 1 
  }
}

et <- Sys.time()
print(et-st)

# Print session info
sessionInfo()
