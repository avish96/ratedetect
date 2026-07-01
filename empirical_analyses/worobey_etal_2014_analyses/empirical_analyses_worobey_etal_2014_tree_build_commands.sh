### Steps and commands used to build Worobey et al. (2014) trees

## Original Worobey et al. (2014) datasets used
# H7.fnl.nex
# H7.fnl.2_rates.mcc.mod.nwk

## nexus to fasta, converting H7.fnl.nex file in Worobey et al. (2014) dataset to fasta file H7.fnl.fasta
pxs2fa -s H7.fnl.nex

## Headers in resulting fasta file changed to match names in the Worobey et al. (2014) nwk tree file, H7.fnl.2_rates.mcc.mod.nwk
# H7.fnl.fasta --> H7.fnl.headermod2.fasta

## Worobey et al. (2014) nwk tree file contains branches not in H7.fnl.nex. Pruning those branches
# H7.fnl.2_rates.mcc.mod.nwk --> H7.fnl.2_rates.mcc.mod2.nwk

## Build constrained H7 tree using the modified fasta file and constraint tree file, with prefix H7fixedT
../bin/raxmlHPC-PTHREADS-AVX2 -f e -t H7.fnl.2_rates.mcc.mod2.nwk -T 8 -m GTRGAMMA -p 12345 -s H7.fnl.headermod2.fasta -n H7fixedT 

## Tree rooted to equine clade on FigTree version 1.4.4 (https://github.com/rambaut/figtree/)
# RAxML_result.H7fixedT --> RAxML_result.H7fixedT_rooted.nwk
