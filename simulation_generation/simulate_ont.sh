# TRAC sites
badread simulate --reference templates/template_trac_wt.fa --quantity 6000x --length  750,10 --error_model nanopore2023 --qscore_model nanopore2023 --junk_reads  0  --random_reads 0 --chimeras 0 --glitches 0,0,0 > nanopore_simulation/simulated_trac_wt.fastq

badread simulate --reference templates/template_trac_del_small.fa  --quantity 15x --length  750,10 --error_model nanopore2023 --qscore_model nanopore2023 --junk_reads  0  --random_reads 0 --chimeras 0 --glitches 0,0,0 > nanopore_simulation/simulated_trac_del_small.fastq

badread simulate --reference templates/template_trac_del_large.fa  --quantity 15x --length  750,10 --error_model nanopore2023 --qscore_model nanopore2023 --junk_reads  0  --random_reads 0 --chimeras 0 --glitches 0,0,0 > nanopore_simulation/simulated_trac_del_large.fastq

badread simulate --reference templates/template_trac_ins_large.fa  --quantity 15x --length  1250,10 --error_model nanopore2023 --qscore_model nanopore2023 --junk_reads  0  --random_reads 0 --chimeras 0 --glitches 0,0,0 > nanopore_simulation/simulated_trac_ins_large.fastq

badread simulate --reference templates/template_trac_ins_small.fa  --quantity 15x --length  1250,10 --error_model nanopore2023 --qscore_model nanopore2023 --junk_reads  0  --random_reads 0 --chimeras 0 --glitches 0,0,0 > nanopore_simulation/simulated_trac_ins_small.fastq

badread simulate --reference templates/template_trac_snp.fa --quantity 15x --length  750,10 --error_model nanopore2023 --qscore_model nanopore2023 --junk_reads  0  --random_reads 0 --chimeras 0 --glitches 0,0,0 > nanopore_simulation/simulated_trac_snp.fastq


# TRBC WT templates
badread simulate --reference templates/template_trbc_wt.fa --quantity 6000x --length 750,10 --error_model nanopore2023 --qscore_model nanopore2023 --junk_reads 0 --random_reads 0 --chimeras 0 --glitches 0,0,0 > nanopore_simulation/simulated_trbc_wt.fastq

# TRBC Small deletions
badread simulate --reference templates/template_trbc_del_small.fa --quantity 15x --length 750,10 --error_model nanopore2023 --qscore_model nanopore2023 --junk_reads 0 --random_reads 0 --chimeras 0 --glitches 0,0,0 > nanopore_simulation/simulated_trbc_del_small.fastq

# TRBC Large deletions  
badread simulate --reference templates/template_trbc_del_large.fa --quantity 15x --length 750,10 --error_model nanopore2023 --qscore_model nanopore2023 --junk_reads 0 --random_reads 0 --chimeras 0 --glitches 0,0,0 > nanopore_simulation/simulated_trbc_del_large.fastq

# TRBC Large insertions
badread simulate --reference templates/template_trbc_ins_large.fa --quantity 15x --length 1250,10 --error_model nanopore2023 --qscore_model nanopore2023 --junk_reads 0 --random_reads 0 --chimeras 0 --glitches 0,0,0 > nanopore_simulation/simulated_trbc_ins_large.fastq

# TRBC Small insertions
badread simulate --reference templates/template_trbc_ins_small.fa --quantity 15x --length 1250,10 --error_model nanopore2023 --qscore_model nanopore2023 --junk_reads 0 --random_reads 0 --chimeras 0 --glitches 0,0,0 > nanopore_simulation/simulated_trbc_ins_small.fastq

# TRBC SNPs
badread simulate --reference templates/template_trbc_snp.fa --quantity 15x --length 750,10 --error_model nanopore2023 --qscore_model nanopore2023 --junk_reads 0 --random_reads 0 --chimeras 0 --glitches 0,0,0 > nanopore_simulation/simulated_trbc_snp.fastq


badread simulate --reference templates/template_cd19_hdr.fa --quantity 6000x --length 2600,10 --error_model nanopore2023 --qscore_model nanopore2023 --junk_reads 0 --random_reads 0 --chimeras 0  --glitches 0,0,0 > nanopore_simulation/simulated_trac_cd19_hdr.fastq

badread simulate --reference templates/template_cd19_itr.fa --quantity 6000x --length 3100,10 --error_model nanopore2023 --qscore_model nanopore2023 --junk_reads 0 --random_reads 0 --chimeras 0  --glitches 0,0,0 > nanopore_simulation/simulated_trac_cd19_itr.fastq

cat templates/template_cd19_itr_left.fa templates/template_cd19_itr_right.fa > templates/template_cd19_itr_oneside.fa
cat templates/template_cd22_itr_left.fa templates/template_cd22_itr_right.fa > templates/template_cd22_itr_oneside.fa

badread simulate --reference templates/template_cd19_itr_oneside.fa --quantity 3000x --length 2850,10 --error_model nanopore2023 --qscore_model nanopore2023 --junk_reads 0 --random_reads 0 --chimeras 0  --glitches 0,0,0 > nanopore_simulation/simulated_trac_cd19_itr_oneside.fastq

badread simulate --reference templates/template_cd19_truncated_hdr.fa --quantity 300x --length 2600,10 --error_model nanopore2023 --qscore_model nanopore2023 --junk_reads 0 --random_reads 0 --chimeras 0  --glitches 0,0,0 > nanopore_simulation/simulated_trac_cd19_truncated_hdr.fastq

badread simulate --reference templates/template_trac_cd22_hdr.fa --quantity 6000x --length 3950,10  --error_model nanopore2023 --qscore_model nanopore2023 --junk_reads 0 --random_reads 0 --chimeras 0  --glitches 0,0,0 > nanopore_simulation/simulated_trac_cd22_hdr.fastq

badread simulate --reference templates/template_trac_cd22_itr.fa --quantity 6000x --length 4500,10  --error_model nanopore2023 --qscore_model nanopore2023 --junk_reads 0 --random_reads 0 --chimeras 0  --glitches 0,0,0 > nanopore_simulation/simulated_trac_cd22_itr.fastq


# CD22 HDR insertion
badread simulate --reference templates/template_cd22_hdr.fa --quantity 6000x --length 3950,10 --error_model nanopore2023 --qscore_model nanopore2023 --junk_reads 0 --random_reads 0 --chimeras 0 --glitches 0,0,0 > nanopore_simulation/simulated_trbc_cd22_hdr.fastq

# CD22 ITR insertion  
badread simulate --reference templates/template_cd22_itr.fa --quantity 6000x --length 4500,10 --error_model nanopore2023 --qscore_model nanopore2023 --junk_reads 0 --random_reads 0 --chimeras 0 --glitches 0,0,0 > nanopore_simulation/simulated_trbc_cd22_itr.fastq

# CD22 one-sided ITR
badread simulate --reference templates/template_cd22_itr_oneside.fa --quantity 3000x --length 4200,10 --error_model nanopore2023 --qscore_model nanopore2023 --junk_reads 0 --random_reads 0 --chimeras 0 --glitches 0,0,0 > nanopore_simulation/simulated_trbc_cd22_itr_oneside.fastq

# CD22 truncated HDR
badread simulate --reference templates/template_cd22_truncated_hdr.fa --quantity 300x --length 3950,10 --error_model nanopore2023 --qscore_model nanopore2023 --junk_reads 0 --random_reads 0 --chimeras 0 --glitches 0,0,0 > nanopore_simulation/simulated_trbc_cd22_truncated_hdr.fastq

# CD22 cross-contamination (CD19 in TRBC target)
badread simulate --reference templates/template_trbc_cd19_hdr.fa --quantity 6000x --length 2600,10 --error_model nanopore2023 --qscore_model nanopore2023 --junk_reads 0 --random_reads 0 --chimeras 0 --glitches 0,0,0 > nanopore_simulation/simulated_trbc_cd19_hdr.fastq

# CD22 cross-contamination ITR (CD19 in TRBC target)
badread simulate --reference templates/template_trbc_cd19_itr.fa --quantity 6000x --length 3100,10 --error_model nanopore2023 --qscore_model nanopore2023 --junk_reads 0 --random_reads 0 --chimeras 0 --glitches 0,0,0 > nanopore_simulation/simulated_trbc_cd19_itr.fastq


