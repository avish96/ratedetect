#!/usr/bin/bash

# Save a log file
LOGFILE="./bd_run.log"

# Function to log messages
log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1" | tee -a "$LOGFILE"
}

# Record start of script
STARTTIME=$(date +%s)
log_message "Script started"

######

nohup R CMD BATCH birthdeath_treeprep.R birthdeath_treeprep.Rout 2>&1

nohup R CMD BATCH master_simulation_script_birthdeath_template.R master_simulation_script_birthdeath_template.Rout 2>&1

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
