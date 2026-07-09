#!/usr/bin/bash

cd ./phylopomp_simulations/
nohup bash run_phylopomp_scripts.sh

cd ../to_etal_2016_simulations/
nohup bash run_birthdeath_scripts.sh
