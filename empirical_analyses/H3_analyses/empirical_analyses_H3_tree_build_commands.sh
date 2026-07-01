### Steps and commands used to build H3 trees


### Unevenly sampled

## GISAID datasets metadata with accession information
# gisaid_epiflu_isolates_H3_Avian.xls
# gisaid_epiflu_isolates_H3_Swine.xls

## Downloaded fasta datasets for Avian and Swine from GISAID
# gisaid_epiflu_sequence_H3_HA_Avian.fasta
# gisaid_epiflu_sequence_H3_HA_Swine.fasta

## file combined on Python
## header names also modified for easier reference
# gisaid_epiflu_sequence_H3_HA_Avian.fasta + gisaid_epiflu_sequence_H3_HA_Swine.fasta --> gisaid_epiflu_sequence_H3_HA_Avian_Swine.fasta

## file edited on Python to remove sequences below 1650 bp, duplicate sequences, and sequences with duplicate accession id 
# gisaid_epiflu_sequence_H3_HA_Avian_Swine.fasta --> gisaid_epiflu_sequence_H3_HA_Avian_Swine_lengthdupidfilt.fasta

## aligned with MAFFT
mafft gisaid_epiflu_sequence_H3_HA_Avian_Swine_lengthdupidfilt.fasta > gisaid_epiflu_sequence_H3_HA_Avian_Swine_lengthdupidfilt_aligned.fasta

## trimmed using trimAl
/Applications/trimal-trimAl/source/trimal -in gisaid_epiflu_sequence_H3_HA_Avian_Swine_lengthdupidfilt_aligned.fasta -out gisaid_epiflu_sequence_H3_HA_Avian_Swine_lengthdupidfilt_aligned_trimmed.fasta

## sequence names edited for raxml 
# gisaid_epiflu_sequence_H3_HA_Avian_Swine_lengthdupidfilt_aligned_trimmed.fasta --> gisaid_epiflu_sequence_H3_HA_Avian_Swine_lengthdupidfilt_aligned_trimmed_edit.fasta

## hanging ends removed via removal of all sites where alignment is over 40% empty
pxclsq -s gisaid_epiflu_sequence_H3_HA_Avian_Swine_lengthdupidfilt_aligned_trimmed_edit.fasta -p 0.6 > gisaid_epiflu_sequence_H3_HA_Avian_Swine_lengthdupidfilt_aligned_trimmed_edit_nohang.fasta

## tree built on raxml
../bin/raxmlHPC-PTHREADS-AVX2 -T 24 -m GTRGAMMA -p 12345 -s gisaid_epiflu_sequence_H3_HA_Avian_Swine_lengthdupidfilt_aligned_trimmed_edit_nohang.fasta -n H3AS_T 

## sequences with unusually long branches, or showing stretches of 3 or more '-' sites, manually removed from alignment
# gisaid_epiflu_sequence_H3_HA_Avian_Swine_lengthdupidfilt_aligned_trimmed_edit_nohang.fasta --> gisaid_epiflu_sequence_H3_HA_Avian_Swine_lengthdupidfilt_aligned_trimmed_edit_nohang_LBA.fasta

# tree rebuilt on raxml
../bin/raxmlHPC-PTHREADS-AVX2 -T 24 -m GTRGAMMA -p 12345 -s gisaid_epiflu_sequence_H3_HA_Avian_Swine_lengthdupidfilt_aligned_trimmed_edit_nohang_LBA.fasta -n H3AS_LBA_T 

## tree saved as newick file on FigTree version 1.4.4 (https://github.com/rambaut/figtree/)
# RAxML_bestTree.H3AS_LBA_T --> RAxML_bestTree.H3AS_LBA_T_unrooted.nwk

## name and dates extracted from alignment for use in dating in R script
# gisaid_epiflu_sequence_H3_HA_Avian_Swine_lengthdupidfilt_aligned_trimmed_edit_nohang_LBA_names.csv
# gisaid_epiflu_sequence_H3_HA_Avian_Swine_lengthdupidfilt_aligned_trimmed_edit_nohang_LBA_dates.csv



### Evenly sampled

## GISAID datasets metadata with accession information
# gisaid_epiflu_isolates_H3_nopass_HA_Avian.xls
# gisaid_epiflu_isolates_H3_nopass_HA_Swine.xls

## Downloaded fasta datasets for Avian and Swine from GISAID, without filter for passaging
# gisaid_epiflu_sequence_H3_nopass_HA_Avian.fasta
# gisaid_epiflu_sequence_H3_nopass_HA_Swine.fasta

## On Python: file combined; edited to remove sequences below 1650 bp, duplicate sequences, and sequences with duplicate accession id; and sampled with random.sample() to 100 sequences per year
## header names also modified for easier reference
# gisaid_epiflu_sequence_H3_HA_Avian.fasta + gisaid_epiflu_sequence_H3_HA_Swine.fasta --> gisaid_epiflu_sequence_H3_HA_Avian_Swine_nopass_lengthdupidfilt_sampled.fasta

## aligned with MAFFT
mafft gisaid_epiflu_sequence_H3_HA_Avian_Swine_nopass_lengthdupidfilt_sampled.fasta > gisaid_epiflu_sequence_H3_HA_Avian_Swine_nopass_lengthdupidfilt_sampled_aligned.fasta

## trimmed using trimAl
/Applications/trimal-trimAl/source/trimal -in gisaid_epiflu_sequence_H3_HA_Avian_Swine_nopass_lengthdupidfilt_sampled_aligned.fasta -out gisaid_epiflu_sequence_H3_HA_Avian_Swine_nopass_lengthdupidfilt_sampled_aligned_trimmed.fasta

## sequence names edited for raxml 
# gisaid_epiflu_sequence_H3_HA_Avian_Swine_nopass_lengthdupidfilt_sampled_aligned_trimmed.fasta --> gisaid_epiflu_sequence_H3_HA_Avian_Swine_nopass_lengthdupidfilt_sampled_aligned_trimmed_edit.fasta

## hanging ends removed via removal of all sites where alignment is over 40% empty
pxclsq -s gisaid_epiflu_sequence_H3_HA_Avian_Swine_nopass_lengthdupidfilt_sampled_aligned_trimmed_edit.fasta -p 0.6 > gisaid_epiflu_sequence_H3_HA_Avian_Swine_nopass_lengthdupidfilt_sampled_aligned_trimmed_edit_nohang.fasta

## sequences showing stretches of 3 or more '-' sites manually removed from alignment
# gisaid_epiflu_sequence_H3_HA_Avian_Swine_nopass_lengthdupidfilt_sampled_aligned_trimmed_edit_nohang.fasta --> gisaid_epiflu_sequence_H3_HA_Avian_Swine_nopass_lengthdupidfilt_sampled_aligned_trimmed_edit_nohang_LBA.fasta

# tree built on raxml
../bin/raxmlHPC-PTHREADS-AVX2 -T 24 -m GTRGAMMA -p 12345 -s gisaid_epiflu_sequence_H3_HA_Avian_Swine_nopass_lengthdupidfilt_sampled_aligned_trimmed_edit_nohang_LBA.fasta -n H3AS_nopass_samp_LBA_T 

## tree saved as newick file on FigTree version 1.4.4 (https://github.com/rambaut/figtree/)
# RAxML_bestTree.H3AS_nopass_samp_LBA_T --> RAxML_bestTree.H3AS_nopass_samp_LBA_T_unrooted.nwk

## name and dates extracted from alignment for use in dating in R script
# gisaid_epiflu_sequence_H3_HA_Avian_Swine_nopass_lengthdupidfilt_sampled_aligned_trimmed_edit_nohang_LBA_names.csv
# gisaid_epiflu_sequence_H3_HA_Avian_Swine_nopass_lengthdupidfilt_sampled_aligned_trimmed_edit_nohang_LBA_dates.csv






