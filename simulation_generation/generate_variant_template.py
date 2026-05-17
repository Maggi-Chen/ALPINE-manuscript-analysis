#!/usr/bin/env python3

import random
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def generate_random_sequence(length):
    """Generate random DNA sequence of specified length"""
    bases = ['A', 'T', 'G', 'C']
    return ''.join(random.choice(bases) for _ in range(length))

def introduce_snp(sequence, position):
    """Introduce a SNP at the specified position"""
    bases = ['A', 'T', 'G', 'C']
    original_base = sequence[position]
    # Choose a different base
    new_bases = [b for b in bases if b != original_base]
    new_base = random.choice(new_bases)
    
    new_seq = sequence[:position] + new_base + sequence[position + 1:]
    return new_seq

def introduce_deletion(sequence, start_pos, size):
    """Introduce a deletion at the specified position"""
    end_pos = min(start_pos + size, len(sequence))
    actual_size = end_pos - start_pos
    
    # Create new sequence with deletion
    new_seq = sequence[:start_pos] + sequence[end_pos:]
    return new_seq, actual_size, start_pos, end_pos

def introduce_insertion(sequence, position, size):
    """Introduce an insertion at the specified position"""
    # Generate random sequence for insertion
    insertion_seq = generate_random_sequence(size)
    
    # Create new sequence with insertion
    new_seq = sequence[:position] + insertion_seq + sequence[position:]
    return new_seq

def get_variant_position(cleavage_pos, variant_window=20):
    """Get random position within variant window around cleavage site"""
    start_pos = max(0, cleavage_pos - variant_window)
    end_pos = cleavage_pos + variant_window
    return random.randint(start_pos, end_pos)

def get_deletion_position_and_size(sequence_length, cleavage_pos, deletion_size_range, cleavage_window=20):
    """
    Get deletion start position and size such that the deletion affects the cleavage site region.
    The deletion should either start within the cleavage window or span across it.
    Constraints: Preserve PCR primers (~30bp each end) and maintain minimum amplifiable sequence.
    """
    min_size, max_size = deletion_size_range
    
    # PCR primer constraints
    primer_buffer = 50  # Preserve ~30bp at each end for PCR primers
    min_remaining_seq = 100  # Minimum sequence length after deletion for PCR amplification
    
    # Define safe deletion region (avoid PCR primer regions)
    safe_start = primer_buffer
    safe_end = sequence_length - primer_buffer
    
    # Define the cleavage region that should be affected
    cleavage_start = max(safe_start, cleavage_pos - cleavage_window)
    cleavage_end = min(safe_end, cleavage_pos + cleavage_window)
    
    # Strategy 1: Start within cleavage window (50% of cases)
    # Strategy 2: Start before cleavage window but span into it (50% of cases)
    if random.random() < 0.5:
        # Start within cleavage window
        latest_start = max(safe_start, cleavage_end - 1)
        earliest_start = max(safe_start, cleavage_start)
        start_pos = random.randint(earliest_start, latest_start)
    else:
        # Start before cleavage window but ensure deletion spans into it
        latest_start = max(safe_start, cleavage_start - 1)
        earliest_start = safe_start
        start_pos = random.randint(earliest_start, latest_start)
        
    # Calculate maximum allowed deletion size to preserve minimum sequence
    max_del_for_min_seq = 500
    
    # Calculate maximum deletion size based on start position and primer constraints
    max_del_from_pos = safe_end - start_pos
    
    # Use the most restrictive constraint
    max_allowed_size = min(max_size, max_del_for_min_seq, max_del_from_pos)
    
    # Ensure we have a valid range
    if max_allowed_size < min_size:
        # If constraints are too tight, use a smaller deletion
        actual_size = min(min_size, max_del_from_pos, max_del_for_min_seq)
        actual_size = max(1, actual_size)  # At least 1bp deletion
    else:
        actual_size = random.randint(min_size, max_allowed_size)
    
    return start_pos, actual_size

def generate_snp_templates(reference_seq, cleavage_pos, target_name, num_templates=500):
    """Generate SNP variant templates"""
    templates = []
    
    for i in range(num_templates):
        # Get random position within ±20bp of cleavage site
        variant_pos = get_variant_position(cleavage_pos)
        
        # Introduce SNP
        variant_seq = introduce_snp(reference_seq, variant_pos)
        
        # Create SeqRecord
        record_id = f"{target_name}_SNP_{i+1}_pos{variant_pos}"
        record = SeqRecord(Seq(variant_seq), id=record_id, description=f"SNP at position {variant_pos}")
        templates.append(record)
    
    return templates

def generate_small_deletion_templates(reference_seq, cleavage_pos, target_name, num_templates=500):
    """Generate small deletion (1-50bp) variant templates"""
    templates = []
    i=0
    while len(templates) < num_templates:
    #for i in range(num_templates):
        # Get deletion position and size that affects cleavage region
        start_pos, del_size = get_deletion_position_and_size(
            len(reference_seq), cleavage_pos, (1, 50), cleavage_window=20
        )
        
        # Introduce deletion
        variant_seq, actual_size, actual_start, actual_end = introduce_deletion(reference_seq, start_pos, del_size)
        
        # Create SeqRecord
        record_id = f"{target_name}_DEL_small_{i+1}_pos{actual_start}-{actual_end}_size{actual_size}"
        record = SeqRecord(Seq(variant_seq), id=record_id,
                         description=f"Small deletion of {actual_size}bp from {actual_start} to {actual_end}")
        if min(actual_end,cleavage_pos+20)-max(actual_start,cleavage_pos-20)>0:
            templates.append(record)
            i+=1
    
    return templates

def generate_large_deletion_templates(reference_seq, cleavage_pos, target_name, num_templates=500):
    """Generate large deletion (51-200bp) variant templates"""
    templates = []
    i=0
    while len(templates) < num_templates:
    #for i in range(num_templates):
        # Get deletion position and size that affects cleavage region
        start_pos, del_size = get_deletion_position_and_size(
            len(reference_seq), cleavage_pos, (51, 200), cleavage_window=20
        )
        
        # Introduce deletion
        variant_seq, actual_size, actual_start, actual_end = introduce_deletion(reference_seq, start_pos, del_size)
        
        # Create SeqRecord
        record_id = f"{target_name}_DEL_large_{i+1}_pos{actual_start}-{actual_end}_size{actual_size}"
        record = SeqRecord(Seq(variant_seq), id=record_id,
                         description=f"Large deletion of {actual_size}bp from {actual_start} to {actual_end}")
        if min(actual_end,cleavage_pos+20)-max(actual_start,cleavage_pos-20)>0:
            templates.append(record)
            i+=1
    
    return templates

def generate_small_insertion_templates(reference_seq, cleavage_pos, target_name, num_templates=500):
    """Generate small insertion (1-50bp) variant templates"""
    templates = []
    
    for i in range(num_templates):
        # Get random position and size for insertion
        variant_pos = get_variant_position(cleavage_pos)
        ins_size = random.randint(1, 50)  # Small insertions: 1-50bp
        
        # Introduce insertion
        variant_seq = introduce_insertion(reference_seq, variant_pos, ins_size)
        
        # Create SeqRecord
        record_id = f"{target_name}_INS_small_{i+1}_pos{variant_pos}_size{ins_size}"
        record = SeqRecord(Seq(variant_seq), id=record_id, 
                         description=f"Small insertion of {ins_size}bp at position {variant_pos}")
        templates.append(record)
    
    return templates

def generate_large_insertion_templates(reference_seq, cleavage_pos, target_name, num_templates=500):
    """Generate large insertion (51-200bp) variant templates"""
    templates = []
    
    for i in range(num_templates):
        # Get random position and size for insertion
        variant_pos = get_variant_position(cleavage_pos)
        ins_size = random.randint(51, 200)  # Large insertions: 51-200bp
        
        # Introduce insertion
        variant_seq = introduce_insertion(reference_seq, variant_pos, ins_size)
        
        # Create SeqRecord
        record_id = f"{target_name}_INS_large_{i+1}_pos{variant_pos}_size{ins_size}"
        record = SeqRecord(Seq(variant_seq), id=record_id, 
                         description=f"Large insertion of {ins_size}bp at position {variant_pos}")
        templates.append(record)
    
    return templates

def save_templates_to_file(templates, output_file):
    """Save templates to FASTA file"""
    with open(output_file, 'w') as f:
        SeqIO.write(templates, f, 'fasta')
    print(f"Saved {len(templates)} templates to {output_file}")

def main():
    # Set random seed for reproducibility
    random.seed(42)
    
    # Load reference sequences
    trac_ref = open('reference/reference_trac.fa','r').read().split('\n')[1]
    trbc_ref = open('reference/reference_trbc.fa','r').read().split('\n')[1]
    
    # Define cleavage positions (middle of the reference sequences)
    trac_cleavage = len(trac_ref) // 2  # Middle of TRAC reference
    trbc_cleavage = len(trbc_ref) // 2  # Middle of TRBC reference
    
    print(f"TRAC reference length: {len(trac_ref)}, cleavage position: {trac_cleavage}")
    print(f"TRBC reference length: {len(trbc_ref)}, cleavage position: {trbc_cleavage}")
    
    # Generate variant classes for TRAC target
    print("\nGenerating TRAC variant templates...")
    
    # SNPs
    trac_snp_templates = generate_snp_templates(trac_ref, trac_cleavage, "TRAC")
    save_templates_to_file(trac_snp_templates, 'template_trac_snp.fa')
    
    # Small deletions
    trac_del_small_templates = generate_small_deletion_templates(trac_ref, trac_cleavage, "TRAC")
    save_templates_to_file(trac_del_small_templates, 'template_trac_del_small.fa')
    
    # Large deletions
    trac_del_large_templates = generate_large_deletion_templates(trac_ref, trac_cleavage, "TRAC")
    save_templates_to_file(trac_del_large_templates, 'template_trac_del_large.fa')
    
    # Small insertions
    trac_ins_small_templates = generate_small_insertion_templates(trac_ref, trac_cleavage, "TRAC")
    save_templates_to_file(trac_ins_small_templates, 'template_trac_ins_small.fa')
    
    # Large insertions
    trac_ins_large_templates = generate_large_insertion_templates(trac_ref, trac_cleavage, "TRAC")
    save_templates_to_file(trac_ins_large_templates, 'template_trac_ins_large.fa')
    
    # Generate variant classes for TRBC target
    print("\nGenerating TRBC variant templates...")
    
    # SNPs
    trbc_snp_templates = generate_snp_templates(trbc_ref, trbc_cleavage, "TRBC")
    save_templates_to_file(trbc_snp_templates, 'template_trbc_snp.fa')
    
    # Small deletions
    trbc_del_small_templates = generate_small_deletion_templates(trbc_ref, trbc_cleavage, "TRBC")
    save_templates_to_file(trbc_del_small_templates, 'template_trbc_del_small.fa')
    
    # Large deletions
    trbc_del_large_templates = generate_large_deletion_templates(trbc_ref, trbc_cleavage, "TRBC")
    save_templates_to_file(trbc_del_large_templates, 'template_trbc_del_large.fa')
    
    # Small insertions
    trbc_ins_small_templates = generate_small_insertion_templates(trbc_ref, trbc_cleavage, "TRBC")
    save_templates_to_file(trbc_ins_small_templates, 'template_trbc_ins_small.fa')
    
    # Large insertions
    trbc_ins_large_templates = generate_large_insertion_templates(trbc_ref, trbc_cleavage, "TRBC")
    save_templates_to_file(trbc_ins_large_templates, 'template_trbc_ins_large.fa')
    
    # Print summary
    print(f"\nVariant template generation complete!")
    print(f"Generated templates for each target (TRAC and TRBC):")
    print(f"  - SNPs: 500 templates each")
    print(f"  - Small deletions (1-50bp): 500 templates each") 
    print(f"  - Large deletions (51-200bp): 500 templates each")
    print(f"  - Small insertions (1-50bp): 500 templates each")
    print(f"  - Large insertions (51-200bp): 500 templates each")
    print(f"Total: 5,000 templates per target (10,000 total)")
    
    print(f"\nAll variants positioned within ±20bp of cleavage sites")
    print(f"Files saved with naming pattern: template_[target]_[variant_type].fa")

if __name__ == "__main__":
    main()
