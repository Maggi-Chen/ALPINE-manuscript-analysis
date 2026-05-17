#!/bin/bash

# PacBio HiFi simulation script based on work.sh
# Uses pacbio2021 error/qscore models with 99.5% accuracy

mkdir -p pacbio_simulation

# TRAC site simulations
badread simulate --reference templates/template_trac_wt.fa --quantity 6000x --length 750,10 --identity 99.5,99.9,0.5 --error_model pacbio2021 --qscore_model pacbio2021 --junk_reads 0 --random_reads 0 --chimeras 0 --glitches 0,0,0 > pacbio_simulation/simulated_trac_wt.fastq

badread simulate --reference templates/template_trac_del_small.fa --quantity 15x --length 750,10 --identity 99.5,99.9,0.5 --error_model pacbio2021 --qscore_model pacbio2021 --junk_reads 0 --random_reads 0 --chimeras 0 --glitches 0,0,0 > pacbio_simulation/simulated_trac_del_small.fastq

badread simulate --reference templates/template_trac_del_large.fa --quantity 15x --length 750,10 --identity 99.5,99.9,0.5 --error_model pacbio2021 --qscore_model pacbio2021 --junk_reads 0 --random_reads 0 --chimeras 0 --glitches 0,0,0 > pacbio_simulation/simulated_trac_del_large.fastq

badread simulate --reference templates/template_trac_ins_large.fa --quantity 15x --length 1250,10 --identity 99.5,99.9,0.5 --error_model pacbio2021 --qscore_model pacbio2021 --junk_reads 0 --random_reads 0 --chimeras 0 --glitches 0,0,0 > pacbio_simulation/simulated_trac_ins_large.fastq

badread simulate --reference templates/template_trac_ins_small.fa --quantity 15x --length 1250,10 --identity 99.5,99.9,0.5 --error_model pacbio2021 --qscore_model pacbio2021 --junk_reads 0 --random_reads 0 --chimeras 0 --glitches 0,0,0 > pacbio_simulation/simulated_trac_ins_small.fastq

badread simulate --reference templates/template_trac_snp.fa --quantity 15x --length 750,10 --identity 99.5,99.9,0.5 --error_model pacbio2021 --qscore_model pacbio2021 --junk_reads 0 --random_reads 0 --chimeras 0 --glitches 0,0,0 > pacbio_simulation/simulated_trac_snp.fastq

# TRBC site simulations
badread simulate --reference templates/template_trbc_wt.fa --quantity 6000x --length 750,10 --identity 99.5,99.9,0.5 --error_model pacbio2021 --qscore_model pacbio2021 --junk_reads 0 --random_reads 0 --chimeras 0 --glitches 0,0,0 > pacbio_simulation/simulated_trbc_wt.fastq

badread simulate --reference templates/template_trbc_del_small.fa --quantity 15x --length 750,10 --identity 99.5,99.9,0.5 --error_model pacbio2021 --qscore_model pacbio2021 --junk_reads 0 --random_reads 0 --chimeras 0 --glitches 0,0,0 > pacbio_simulation/simulated_trbc_del_small.fastq

badread simulate --reference templates/template_trbc_del_large.fa --quantity 15x --length 750,10 --identity 99.5,99.9,0.5 --error_model pacbio2021 --qscore_model pacbio2021 --junk_reads 0 --random_reads 0 --chimeras 0 --glitches 0,0,0 > pacbio_simulation/simulated_trbc_del_large.fastq

badread simulate --reference templates/template_trbc_ins_large.fa --quantity 15x --length 1250,10 --identity 99.5,99.9,0.5 --error_model pacbio2021 --qscore_model pacbio2021 --junk_reads 0 --random_reads 0 --chimeras 0 --glitches 0,0,0 > pacbio_simulation/simulated_trbc_ins_large.fastq

badread simulate --reference templates/template_trbc_ins_small.fa --quantity 15x --length 1250,10 --identity 99.5,99.9,0.5 --error_model pacbio2021 --qscore_model pacbio2021 --junk_reads 0 --random_reads 0 --chimeras 0 --glitches 0,0,0 > pacbio_simulation/simulated_trbc_ins_small.fastq

badread simulate --reference templates/template_trbc_snp.fa --quantity 15x --length 750,10 --identity 99.5,99.9,0.5 --error_model pacbio2021 --qscore_model pacbio2021 --junk_reads 0 --random_reads 0 --chimeras 0 --glitches 0,0,0 > pacbio_simulation/simulated_trbc_snp.fastq

# CD19 integrations (TRAC site)
badread simulate --reference templates/template_cd19_hdr.fa --quantity 6000x --length 2600,10 --identity 99.5,99.9,0.5 --error_model pacbio2021 --qscore_model pacbio2021 --junk_reads 0 --random_reads 0 --chimeras 0 --glitches 0,0,0 > pacbio_simulation/simulated_trac_cd19_hdr.fastq

badread simulate --reference templates/template_cd19_itr.fa --quantity 6000x --length 3100,10 --identity 99.5,99.9,0.5 --error_model pacbio2021 --qscore_model pacbio2021 --junk_reads 0 --random_reads 0 --chimeras 0 --glitches 0,0,0 > pacbio_simulation/simulated_trac_cd19_itr.fastq

# Create one-sided ITR templates (if not already created)
cat templates/template_cd19_itr_left.fa templates/template_cd19_itr_right.fa > templates/template_cd19_itr_oneside.fa
cat templates/template_cd22_itr_left.fa templates/template_cd22_itr_right.fa > templates/template_cd22_itr_oneside.fa

badread simulate --reference templates/template_cd19_itr_oneside.fa --quantity 3000x --length 2850,10 --identity 99.5,99.9,0.5 --error_model pacbio2021 --qscore_model pacbio2021 --junk_reads 0 --random_reads 0 --chimeras 0 --glitches 0,0,0 > pacbio_simulation/simulated_trac_cd19_itr_oneside.fastq

badread simulate --reference templates/template_cd19_truncated_hdr.fa --quantity 300x --length 2600,10 --identity 99.5,99.9,0.5 --error_model pacbio2021 --qscore_model pacbio2021 --junk_reads 0 --random_reads 0 --chimeras 0 --glitches 0,0,0 > pacbio_simulation/simulated_trac_cd19_truncated_hdr.fastq

# Cross-contamination: CD22 in TRAC site
badread simulate --reference templates/template_trac_cd22_hdr.fa --quantity 6000x --length 3950,10 --identity 99.5,99.9,0.5 --error_model pacbio2021 --qscore_model pacbio2021 --junk_reads 0 --random_reads 0 --chimeras 0 --glitches 0,0,0 > pacbio_simulation/simulated_trac_cd22_hdr.fastq

badread simulate --reference templates/template_trac_cd22_itr.fa --quantity 6000x --length 4500,10 --identity 99.5,99.9,0.5 --error_model pacbio2021 --qscore_model pacbio2021 --junk_reads 0 --random_reads 0 --chimeras 0 --glitches 0,0,0 > pacbio_simulation/simulated_trac_cd22_itr.fastq

# CD22 integrations (TRBC site)
badread simulate --reference templates/template_cd22_hdr.fa --quantity 6000x --length 3950,10 --identity 99.5,99.9,0.5 --error_model pacbio2021 --qscore_model pacbio2021 --junk_reads 0 --random_reads 0 --chimeras 0 --glitches 0,0,0 > pacbio_simulation/simulated_trbc_cd22_hdr.fastq

badread simulate --reference templates/template_cd22_itr.fa --quantity 6000x --length 4500,10 --identity 99.5,99.9,0.5 --error_model pacbio2021 --qscore_model pacbio2021 --junk_reads 0 --random_reads 0 --chimeras 0 --glitches 0,0,0 > pacbio_simulation/simulated_trbc_cd22_itr.fastq

badread simulate --reference templates/template_cd22_itr_oneside.fa --quantity 3000x --length 4200,10 --identity 99.5,99.9,0.5 --error_model pacbio2021 --qscore_model pacbio2021 --junk_reads 0 --random_reads 0 --chimeras 0 --glitches 0,0,0 > pacbio_simulation/simulated_trbc_cd22_itr_oneside.fastq

badread simulate --reference templates/template_cd22_truncated_hdr.fa --quantity 300x --length 2600,10 --identity 99.5,99.9,0.5 --error_model pacbio2021 --qscore_model pacbio2021 --junk_reads 0 --random_reads 0 --chimeras 0 --glitches 0,0,0 > pacbio_simulation/simulated_trbc_cd22_truncated_hdr.fastq

# Cross-contamination: CD19 in TRBC site
badread simulate --reference templates/template_trbc_cd19_hdr.fa --quantity 6000x --length 2600,10 --identity 99.5,99.9,0.5 --error_model pacbio2021 --qscore_model pacbio2021 --junk_reads 0 --random_reads 0 --chimeras 0 --glitches 0,0,0 > pacbio_simulation/simulated_trbc_cd19_hdr.fastq

badread simulate --reference templates/template_trbc_cd19_itr.fa --quantity 6000x --length 3100,10 --identity 99.5,99.9,0.5 --error_model pacbio2021 --qscore_model pacbio2021 --junk_reads 0 --random_reads 0 --chimeras 0 --glitches 0,0,0 > pacbio_simulation/simulated_trbc_cd19_itr.fastq

echo "PacBio HiFi simulation completed!"
echo "Output files saved in pacbio_simulation/ directory"
