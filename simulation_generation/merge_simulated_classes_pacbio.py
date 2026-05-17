#!/usr/bin/env python3

import random
import sys
import os

def process_fastq_file(filepath, prefix, max_reads=5000):
    """
    Process a FASTQ file to extract first N reads with prefix added to read names
    Returns list of 4-line read records
    """
    reads = []
    
    try:
        with open(filepath, 'r') as f:
            read_count = 0
            while read_count < max_reads:
                # Read 4 lines for each read
                header = f.readline().strip()
                if not header:  # End of file
                    break
                
                sequence = f.readline().strip()
                plus_line = f.readline().strip()
                quality = f.readline().strip()
                
                # Extract UUID (part before space) and add prefix
                if header.startswith('@'):
                    uuid_part = header.split()[0]  # Get @UUID part
                    new_header = f"@{prefix}.{uuid_part[1:]}"  # Remove @ and add prefix
                else:
                    new_header = f"@{prefix}.{header}"
                
                # Store the 4-line read record
                reads.append([new_header, sequence, plus_line, quality])
                read_count += 1
                
    except FileNotFoundError:
        print(f"Warning: File {filepath} not found, skipping...")
        return []
    
    print(f"Processed {len(reads)} reads from {filepath} with prefix '{prefix}'")
    return reads

def main():
    # Set random seed for reproducible shuffling
    random.seed(42)
    
    # Define file mappings by target site for PacBio simulation
    trac_files = {
        # TRAC WT and variants
        'pacbio_simulation/simulated_trac_wt.fastq': 'trac_wt',
        'pacbio_simulation/simulated_trac_del_small.fastq': 'trac_del_small',
        'pacbio_simulation/simulated_trac_del_large.fastq': 'trac_del_large',
        'pacbio_simulation/simulated_trac_ins_small.fastq': 'trac_ins_small',
        'pacbio_simulation/simulated_trac_ins_large.fastq': 'trac_ins_large',
        'pacbio_simulation/simulated_trac_snp.fastq': 'trac_snp',
        
        # CD19 -> TRAC integrations
        'pacbio_simulation/simulated_trac_cd19_hdr.fastq': 'trac_cd19_hdr',
        'pacbio_simulation/simulated_trac_cd19_itr.fastq': 'trac_cd19_itr',
        'pacbio_simulation/simulated_trac_cd19_itr_oneside.fastq': 'trac_cd19_itr_oneside',
        'pacbio_simulation/simulated_trac_cd19_truncated_hdr.fastq': 'trac_cd19_truncated_hdr',
        
        # CD22 -> TRAC (cross-contamination)
        'pacbio_simulation/simulated_trac_cd22_hdr.fastq': 'trac_cd22_hdr',
        'pacbio_simulation/simulated_trac_cd22_itr.fastq': 'trac_cd22_itr'
    }
    
    trbc_files = {
        # TRBC WT and variants
        'pacbio_simulation/simulated_trbc_wt.fastq': 'trbc_wt',
        'pacbio_simulation/simulated_trbc_del_small.fastq': 'trbc_del_small',
        'pacbio_simulation/simulated_trbc_del_large.fastq': 'trbc_del_large',
        'pacbio_simulation/simulated_trbc_ins_small.fastq': 'trbc_ins_small',
        'pacbio_simulation/simulated_trbc_ins_large.fastq': 'trbc_ins_large',
        'pacbio_simulation/simulated_trbc_snp.fastq': 'trbc_snp',
        
        # CD22 -> TRBC integrations
        'pacbio_simulation/simulated_trbc_cd22_hdr.fastq': 'trbc_cd22_hdr',
        'pacbio_simulation/simulated_trbc_cd22_itr.fastq': 'trbc_cd22_itr',
        'pacbio_simulation/simulated_trbc_cd22_itr_oneside.fastq': 'trbc_cd22_itr_oneside',
        'pacbio_simulation/simulated_trbc_cd22_truncated_hdr.fastq': 'trbc_cd22_truncated_hdr',
        
        # CD19 -> TRBC (cross-contamination)
        'pacbio_simulation/simulated_trbc_cd19_hdr.fastq': 'trbc_cd19_hdr',
        'pacbio_simulation/simulated_trbc_cd19_itr.fastq': 'trbc_cd19_itr'
    }
    
    # Process TRAC files
    print("Processing TRAC site PacBio HiFi files...")
    trac_reads = []
    for filepath, prefix in trac_files.items():
        reads = process_fastq_file(filepath, prefix, max_reads=5000)
        trac_reads.extend(reads)
    
    print(f"Total TRAC reads collected: {len(trac_reads)}")
    
    # Process TRBC files
    print("\nProcessing TRBC site PacBio HiFi files...")
    trbc_reads = []
    for filepath, prefix in trbc_files.items():
        reads = process_fastq_file(filepath, prefix, max_reads=5000)
        trbc_reads.extend(reads)
    
    print(f"Total TRBC reads collected: {len(trbc_reads)}")
    
    # Shuffle reads for each target
    print("\nShuffling reads...")
    random.shuffle(trac_reads)
    random.shuffle(trbc_reads)
    
    # Write TRAC output file
    trac_output = 'pacbio_simulation/simulated_hifi_trac.fastq'
    print(f"Writing TRAC HiFi reads to {trac_output}...")
    with open(trac_output, 'w') as f:
        for read in trac_reads:
            f.write('\n'.join(read) + '\n')
    
    # Write TRBC output file
    trbc_output = 'pacbio_simulation/simulated_hifi_trbc.fastq'
    print(f"Writing TRBC HiFi reads to {trbc_output}...")
    with open(trbc_output, 'w') as f:
        for read in trbc_reads:
            f.write('\n'.join(read) + '\n')
    
    print(f"\nSuccessfully created PacBio HiFi datasets:")
    print(f"  {trac_output} with {len(trac_reads)} shuffled reads")
    print(f"  {trbc_output} with {len(trbc_reads)} shuffled reads")
    
    # Print summary statistics for TRAC
    print(f"\nTRAC HiFi summary by category:")
    trac_categories = {}
    for read in trac_reads:
        category = '_'.join(read[0].split('_')[1:3])  # Extract category from @prefix_category_uuid
        if category not in trac_categories:
            trac_categories[category] = 0
        trac_categories[category] += 1
    
    for category, count in sorted(trac_categories.items()):
        print(f"  {category}: {count} reads")
    
    # Print summary statistics for TRBC
    print(f"\nTRBC HiFi summary by category:")
    trbc_categories = {}
    for read in trbc_reads:
        category = '_'.join(read[0].split('_')[1:3])  # Extract category from @prefix_category_uuid
        if category not in trbc_categories:
            trbc_categories[category] = 0
        trbc_categories[category] += 1
    
    for category, count in sorted(trbc_categories.items()):
        print(f"  {category}: {count} reads")

if __name__ == "__main__":
    main()
