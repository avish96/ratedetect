#!/usr/bin/bash

# Save a log file
LOGFILE="./phylopomp_run.log"

# Function to log messages
log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1" | tee -a "$LOGFILE"
}

# Record start of script
STARTTIME=$(date +%s)
log_message "Script started"

##### 

nohup R CMD BATCH phylopomp_treesim.R > phylopomp_treesim.Rout 2>&1 
nohup R CMD BATCH phylopomp_treeprep.R > phylopomp_treeprep.Rout 2>&1 

nohup R CMD BATCH master_simulation_script_phylopomp_highshape.R > master_simulation_script_phylopomp_highshape.Rout 2>&1 
nohup R CMD BATCH master_simulation_script_phylopomp_lowshape.R > master_simulation_script_phylopomp_lowshape.Rout 2>&1 
nohup R CMD BATCH master_simulation_script_phylopomp_highshape_scaled.R > master_simulation_script_phylopomp_highshape_scaled.Rout 2>&1 
nohup R CMD BATCH master_simulation_script_phylopomp_highshape_sampled_third.R > master_simulation_script_phylopomp_highshape_sampled_third.Rout 2>&1 
nohup R CMD BATCH master_simulation_script_phylopomp_highshape_sampled_twothird.R > master_simulation_script_phylopomp_highshape_sampled_twothird.Rout 2>&1 

#####

# Checking and logging if an error occurred
if [ $? -ne 0 ]; then
    log_message "ERROR: Script failed to run, check R scripts and output files"
fi

# Log script execution time
ENDTIME=$(date +%s)
DURATION=$((ENDTIME - STARTTIME))
log_message "Script ended"
log_message "Duration: $DURATION seconds"
log_message "--------------------------------------"
