# Detecting Structured Rate Variation

This repository contains all simulation and empirical workflows and coding scripts for the manuscript _A nonparametric approach for detecting substitution-rate heterogeneity in viral phylogenies. Avinash Subramanian, Adam S. Lauring, Aaron A. King, Stephen A. Smith._

Empirical analysis datasets are also housed in this repository. Simulation output files are uploaded on Zenodo, and available upon request. Upon acceptance of this manuscript, they will be publicly published.

## File Structure 
```
ratedetect
├── data
│   └── True unrooted trees with dates
│       ├── Relaxed molecular clock trees
│       └── Strict molecular clock trees
├── detect_structured_rates.R
├── empirical_analyses
│   ├── figures
│   ├── H3_analyses
│   │   ├── evenly_sampled_(no_passage_filter)
│   │   └── unevenly_sampled
│   └── worobey_etal_2014_analyses
├── scripts
│   ├── birthdeath_treeprep.R
│   ├── master_simulation_script_birthdeath_template.R
│   ├── master_simulation_script_phylopomp_template.R
│   ├── phylopomp_inline_remover.py
│   ├── phylopomp_treeprep.R
│   ├── phylopomp_treesim.R
│   ├── plot_funcs.R
│   ├── run_birthdeath_scripts.sh
│   ├── run_phylopomp_scripts.sh
│   ├── sim_funcs.R
│   └── values_for_analysis_scripts.R
├── simulations
│   ├── phylopomp_simulation_scripts
│   ├── simulation_results
│   │   ├── simulation_plotting_scripts
│   │   └── simulations_SANN
│   │       ├── figures_main
│   │       ├── figures_supp
│   │       ├── phylopomp_results
│   │       └── to_etal_results
│   └── to_etal_2016_simulation_scripts
└── treedater
```

## Method Implementation 

A general command-line implementation of our nonparametric method is provided in `detect_structured_rates.R`. All test statistics described in the paper are included, with the script outputting the corresponding p-values for detection of structured rate variation. 

The current implementation only permits a Gamma distribution for the underlying lineage rates. The script can be modified as desired to consider other test statistics 
and lineage rate distributions.

The script requires the following inputs:

- `--tp` Directory of input substitution tree (newick)
- `--td` Directory of dated tree (newick)
- `--tr` Directory of input rate tree (newick)
- `--dm1` Deme 1 label (one word)
- `--dm2` Deme 2 label (one word)
- `--nsites` Number of sites in alignment (integer)
- `--shape` Gamma shape parameter inferred for tree lineage rates (integer)
- `--scale` Gamma scale parameter inferred for tree lineage rates (integer)

For example (in directory of script, with files in same directory):

```bash
./detect_structured_rates.R --tp tree.nwk --td tree_dated.nwk --tr tree_rateogram.nwk --dm1 deme1 --dm2 deme2 --nsites 1000 --shape 0.01 --scale 100
```

## Method Citation

If you implement our method in your work, please cite:

- Subramanian A, Lauring AS, King AA, Smith SA. A nonparametric approach for detecting substitution-rate heterogeneity in viral phylogenies. 

Please include citations to the associated dependencies as well, depending on the test statistics used.

## Simulation Analyses

All simulation analysis template scripts are housed in the `scripts/` folder. 
Two sets of simulation analyses were conducted: the first with trees simulated using _phylopomp_, the second with birth-death trees simulated in To et al. (2016). 
Actual scripts used for simulation analyses are in the `simulations/` folder in the respective analysis directories.

The files to run these analyses and log execution details are `run_birthdeath_scripts.sh` and `run_phylopomp_scripts.sh`. Essential functions used during the simulation are in `sim_funcs.R`, while essential details on simulation results are extracted using the `values_for_analysis_scripts.R` file. See paper for further analysis details. 

Example tree files used for Figures 2 and S1, and the power analysis result files for simulated annealing, are uploaded in `simulations/simulation_results/simulations_SANN/` in the respective analysis directories.

### _Phylopomp_ Analyses

The simulation script for the _phylopomp_ trees is `phylopomp_treesim.R`. This script results in 200 transmission trees. Inline nodes and conversion to a phylogeny object are done using `phylopomp_treeprep.R` and `phylopomp_inline_remover.py`. The resulting _phylopomp_ trees are then used to conduct the simulation-based power analysis using the
`master_simulation_script_phylopomp_template.R` script files, with template replaced by the specific scenario being tested. 

The `values_for_analysis_scripts.R` file contains code for calculating basic statistics of outputted _phylopomp_ trees, and calculating the sampling weights used in the heterogeneous sampling through time simulations. The latter is based on the FASTA files for H3N2 HA sequences downloaded from GISAID, located in the `data/H3N2_HA_GISAID` folder. Sequences were downloaded from January 1 to December 31 of the year(s) marked in the file, with 2022 split into two files due to data download limits.

### Birth-Death Analyses

The birth-death trees used for these analyses are uploaded in `data/True unrooted trees with dates`, obtained from To et al. (2016) from their  online repository: http://www.atgc-montpellier.fr/LSD/.

The script for simulating demes on these trees is `birthdeath_treeprep.R`. The resulting birth-death trees are then used to conduct the simulation-based power analysis
using the `master_simulation_script_birthdeath_template.R` script files, with template replaced by the specific scenario being tested.

## Empirical Analyses

All empirical analysis scripts and results are housed in `empirical_analyses/`. 

The H7 data from Worobey et al. (2014) is obtained from their Dryad Digital Repository: http://dx.doi.org/10.5061/dryad.m04j9.

GISAID acknowledgement tables for H3 sequence data are provided in supplementary tables S3 - S6 of the manuscript. They are also uploaded here as EPI_SET ID files, together with the GISAID metadata files containing detailed accession code information. For H3 sequence data filtered for original passaging, they are labeled as `gisaid_epiflu_sequence_H3_HA_Avian.xls` for avian taxa and `gisaid_epiflu_sequence_H3_HA_Swine.xls` for Swine taxa. For H3 sequence data from GISIAD not filtered for passaging, they are labeled as `gisaid_epiflu_sequence_H3_nopass_HA_Avian.xls` for Avian taxa, and `gisaid_epiflu_sequence_H3_nopass_HA_Swine.xls` for Swine taxa. 

Alignment and tree-building workflows are described in `empirical_analyses_worobey_etal_2014_tree_build_commands.sh` for the Worobey et al. (2014) H7 dataset, and in 
`empirical_analyses_H3_tree_build_commands.sh` for the GISAID H3 datasets. Evenly- and unevenly-sampled H3 datasets and results are housed in their respective directories. 
The corresponding tree-dating scripts are `empirical_analyses_worobey_etal_2014_tree_dating.R`, and `empirical_analyses_H3_tree_dating.R`, 
and method application scripts (which print the p-values listed in the manuscript) are `empirical_analyses_worobey_etal_2014_method_application.R` 
and `empirical_analyses_H3_method_application.R`. See paper for analysis details.
 
## Reproducibility

To reproduce all simulation results, run `run_phylopomp_scripts.sh` and `run_birthdeath_scripts.sh` for the _phylopomp_ and birth-death analyses, respectively. To run both analyses in a single command, run `run_sims.sh`. The latter was run using `nohup` on a remote cluster for manuscript simulation results. 

To reproduce calculated tree statistics and sampling weights, run the respective code within
the `values_for_analysis_scripts.R` file.

For empirical analyses, to reproduce alignment and tree-building files, follow the step-by-step procedure in the `.sh` files. 
To reproduce tree-dating and method application results, run the respective R files from their directories.

To reproduce all simulation figures, run all R scripts within `simulations/simulation_results/simulation_plotting_scripts/`.

To reproduce all empirical figures, run all R scripts within `empirical_analyses/figures/`

All files are designed to be run from their respective working directories, with in-file file or directory paths all being relative paths.

### Package Dependencies

To reproduce analyses, the package _ete3_ v.3.1.2 is required for Python. In R, the following are required:

- _ape_ v.5.8
- _phytools_ v.2.3.0
- _castor_ v.1.8.3
- _treedater_ v.1.0.2
- _phylopomp_ v.0.14.8.0
- _matrixStats_ v.1.5.0
- _tryCatchLog_ v.1.3.1
- _future_ v.1.34.0
- _doFuture_ v.1.0.1
- _diptest_ v.0.77.1
- _stats_ v.4.3.3
- _stringr_ v.1.5.1
- _readr_ v.2.1.5
- _DescTools_ v.0.99.60

_Treedater_ was manually modified and installed to assess different optimization algorithms. 
The package modified to use simulated annealing is uploaded on this repository.

## References

1. Capella-Gutiérrez S, Silla-Martínez JM, Gabaldón T. 2009. trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses. Bioinformatics 25:1972–1973.
2. Katoh K, Standley DM. 2013. MAFFT Multiple Sequence Alignment Software Version 7: Improvements in Performance and Usability. Molecular Biology and Evolution 30:772–780.
3. Khare S, Gurry C, Freitas L, B Schultz M, Bach G, Diallo A, Akite N, Ho J, Tc Lee R, Yeo W, et al. 2021. GISAID’s Role in Pandemic Response. China CDC Weekly 3:1049–1051.
4. King AA, Lin Q, Ionides EL. 2022. Markov genealogy processes. Theoretical Population Biology 143:77–91.
5. King AA, Lin Q, Ionides EL. 2025. Exact phylodynamic likelihood via structured Markov genealogy processes. Available from: http://arxiv.org/abs/2405.17032.
6. Stamatakis A. 2014. RAxML version 8: a tool for phylogenetic analysis and post-analysis of large phylogenies. Bioinformatics 30:1312–1313.
7. To T-H, Jung M, Lycett S, Gascuel O. 2016. Fast Dating Using Least-Squares Criteria and Algorithms. Syst Biol 65:82–97.
8. Volz EM, Frost SDW. 2017. Scalable Relaxed Clock Phylogenetic Dating. Virus Evol. 3(2):vex025.
9. Worobey M, Han G-Z, Rambaut A. 2014. A synchronized global sweep of the internal genes of modern avian influenza virus. Nature 508:254–257.
