
a=open('genomic_seq.txt','r').read().split('\n')[:-1]
allseq={}
for c in a:
    c=c.split('\t')
    allseq[c[0]]=c[1]

a=open('cd19_cDNA.fasta','r').read().split('\n')[1:-1]
cd19=''.join(a)
allseq['CD19'] = cd19

a=open('cd22_cDNA.fasta','r').read().split('\n')[1:-1]
cd22=''.join(a)
allseq['CD22'] = cd22

allseq['ITR_left'] = open('itr1_seq.txt','r').read()
allseq['ITR_right'] = open('itr2_seq.txt','r').read()
print(allseq)



# HDR: pcr - genomic - HA - kozak - transgene - polyA - HA - genomic - pcr
hdr_seq_cd19 = allseq['TRAC_PCR_left'] + allseq['TRAC_Genomic_left'] + allseq['TRAC_HA_left'] + allseq['kozak'] + allseq['CD19'] + allseq['polyA'] + allseq['TRAC_HA_right'] + allseq['TRAC_Genomic_right'] + allseq['TRAC_PCR_right']

hdr_seq_cd22 = allseq['TRBC_PCR_left'] + allseq['TRBC_Genomic_left'] + allseq['TRBC_HA_left'] + allseq['kozak'] + allseq['CD22'] + allseq['polyA'] + allseq['TRBC_HA_right'] + allseq['TRBC_Genomic_right'] + allseq['TRBC_PCR_right']


# ITR-Full: pcr - genomic - HA - ITR - ITR spacer - HA - kozak - transgene - polyA - HA - ITR spacer - ITR - HA - genomic - pcr

itr_seq_cd19  = allseq['TRAC_PCR_left'] + allseq['TRAC_Genomic_left'] + allseq['TRAC_HA_left'] + allseq['ITR_left'] + allseq['ITR_spacer'] + allseq['TRAC_HA_left'] + allseq['kozak'] + allseq['CD19'] + allseq['polyA'] + allseq['TRAC_HA_right'] + allseq['ITR_spacer'] + allseq['ITR_right'] + allseq['TRAC_HA_right'] + allseq['TRAC_Genomic_right'] + allseq['TRAC_PCR_right']

itr_seq_cd22 = allseq['TRBC_PCR_left'] + allseq['TRBC_Genomic_left'] + allseq['TRBC_HA_left'] + allseq['ITR_left'] + allseq['ITR_spacer'] + allseq['TRBC_HA_left'] + allseq['kozak'] + allseq['CD22'] + allseq['polyA'] + allseq['TRBC_HA_right'] + allseq['ITR_spacer'] + allseq['ITR_right'] + allseq['TRBC_HA_right'] + allseq['TRBC_Genomic_right'] + allseq['TRBC_PCR_right']

f=open('template_cd19_hdr.fa','w')
f.write('>template_cd19_hdr\n'+hdr_seq_cd19+'\n')
f.close()

f=open('template_cd19_itr.fa','w')
f.write('>template_cd19_itr\n'+itr_seq_cd19+'\n')
f.close()

f=open('template_cd22_hdr.fa','w')
f.write('>template_cd22_hdr\n'+hdr_seq_cd22+'\n')
f.close()

f=open('template_cd22_itr.fa','w')
f.write('>template_cd22_itr\n'+itr_seq_cd22+'\n')
f.close()

# 1-side ITR
# pcr - genomic - HA - ITR - ITR spacer - HA - kozak - transgene - polyA - HA - genomic - pcr
# or: pcr - genomic - HA - kozak - transgene - polyA - HA - ITR spacer - ITR - HA - genomic - pcr

itr_left_cd19 = allseq['TRAC_PCR_left'] + allseq['TRAC_Genomic_left'] + allseq['TRAC_HA_left'] + allseq['ITR_left'] + allseq['ITR_spacer'] + allseq['TRAC_HA_left'] + allseq['kozak'] + allseq['CD19'] + allseq['polyA'] + allseq['TRAC_HA_right'] + allseq['TRAC_Genomic_right'] + allseq['TRAC_PCR_right']
itr_right_cd19 = allseq['TRAC_PCR_left'] + allseq['TRAC_Genomic_left'] + allseq['TRAC_HA_left'] + allseq['kozak'] + allseq['CD19'] + allseq['polyA'] + allseq['TRAC_HA_right'] + allseq['ITR_spacer'] + allseq['ITR_right'] + allseq['TRAC_HA_right'] + allseq['TRAC_Genomic_right'] + allseq['TRAC_PCR_right']

itr_left_cd22 = allseq['TRBC_PCR_left'] + allseq['TRBC_Genomic_left'] + allseq['TRBC_HA_left'] + allseq['ITR_left'] + allseq['ITR_spacer'] + allseq['TRBC_HA_left'] + allseq['kozak'] + allseq['CD22'] + allseq['polyA'] + allseq['TRBC_HA_right'] + allseq['TRBC_Genomic_right'] + allseq['TRBC_PCR_right']
itr_right_cd22 = allseq['TRBC_PCR_left'] + allseq['TRBC_Genomic_left'] + allseq['TRBC_HA_left'] + allseq['kozak'] + allseq['CD22'] + allseq['polyA'] + allseq['TRBC_HA_right'] + allseq['ITR_spacer'] + allseq['ITR_right'] + allseq['TRBC_HA_right'] + allseq['TRBC_Genomic_right'] + allseq['TRBC_PCR_right']

f=open('template_cd19_itr_left.fa','w')
f.write('>template_cd19_itr_left\n'+itr_left_cd19+'\n')
f.close()
f=open('template_cd19_itr_right.fa','w')
f.write('>template_cd19_itr_right\n'+itr_right_cd19+'\n')
f.close()
f=open('template_cd22_itr_left.fa','w')
f.write('>template_cd22_itr_left\n'+itr_left_cd22+'\n')
f.close()
f=open('template_cd22_itr_right.fa','w')
f.write('>template_cd22_itr_right\n'+itr_right_cd22+'\n')
f.close()


# trac+cd22
# CD22 into TRAC (should be detected as wrong vector)
# HDR: CD22 transgene with TRAC homology arms
hdr_seq_cd22_to_trac = allseq['TRAC_PCR_left'] + allseq['TRAC_Genomic_left'] + allseq['TRAC_HA_left'] + allseq['kozak'] + allseq['CD22'] + allseq['polyA'] + allseq['TRAC_HA_right'] + allseq['TRAC_Genomic_right'] + allseq['TRAC_PCR_right']

# ITR: CD22 transgene with TRAC context and ITRs  
itr_seq_cd22_to_trac = allseq['TRAC_PCR_left'] + allseq['TRAC_Genomic_left'] + allseq['TRAC_HA_left'] + allseq['ITR_left'] + allseq['ITR_spacer'] + allseq['TRBC_HA_left'] + allseq['kozak'] + allseq['CD22'] + allseq['polyA'] + allseq['TRBC_HA_right'] + allseq['ITR_spacer'] + allseq['ITR_right'] + allseq['TRAC_HA_right'] + allseq['TRAC_Genomic_right'] + allseq['TRAC_PCR_right']

# CD19 into TRBC (should be detected as wrong vector)  
# HDR: CD19 transgene with TRBC homology arms
hdr_seq_cd19_to_trbc = allseq['TRBC_PCR_left'] + allseq['TRBC_Genomic_left'] + allseq['TRBC_HA_left'] + allseq['kozak'] + allseq['CD19'] + allseq['polyA'] + allseq['TRBC_HA_right'] + allseq['TRBC_Genomic_right'] + allseq['TRBC_PCR_right']

# ITR: CD19 transgene with TRBC context and ITRs
itr_seq_cd19_to_trbc = allseq['TRBC_PCR_left'] + allseq['TRBC_Genomic_left'] + allseq['TRBC_HA_left'] + allseq['ITR_left'] + allseq['ITR_spacer'] + allseq['TRAC_HA_left'] + allseq['kozak'] + allseq['CD19'] + allseq['polyA'] + allseq['TRAC_HA_right'] + allseq['ITR_spacer'] + allseq['ITR_right'] + allseq['TRBC_HA_right'] + allseq['TRBC_Genomic_right'] + allseq['TRBC_PCR_right']

f=open('template_trac_cd22_hdr.fa','w')
f.write('>template_trac_cd22_hdr\n'+hdr_seq_cd22_to_trac+'\n')
f.close()
f=open('template_trac_cd22_itr.fa','w')
f.write('>template_trac_cd22_itr\n'+itr_seq_cd22_to_trac+'\n')
f.close()
f=open('template_trbc_cd19_hdr.fa','w')
f.write('>template_trbc_cd19_hdr\n'+hdr_seq_cd19_to_trbc+'\n')
f.close()
f=open('template_trbc_cd19_itr.fa','w')
f.write('>template_trbc_cd19_itr\n'+itr_seq_cd19_to_trbc+'\n')
f.close()



# Reference file includes HDR + HDR + potential wrong vector insertion
wt_trac = allseq['TRAC_PCR_left'] + allseq['TRAC_Genomic_left'] + allseq['TRAC_HA_left'] +  allseq['TRAC_HA_right'] + allseq['TRAC_Genomic_right'] + allseq['TRAC_PCR_right']
wt_trbc = allseq['TRBC_PCR_left'] + allseq['TRBC_Genomic_left'] + allseq['TRBC_HA_left'] + allseq['TRBC_HA_right'] + allseq['TRBC_Genomic_right'] + allseq['TRBC_PCR_right']

f=open('reference_trac.fa','w')
f.write('>WT_TRAC\n'+wt_trac+'\n')
f.write('>Ref_HDR_CD19\n'+hdr_seq_cd19+'\n')
f.write('>Ref_ITR_CD19\n'+itr_seq_cd19+'\n')
f.write('>Ref_ITR_CD22\n'+itr_seq_cd22_to_trac+'\n')
f.close()
f=open('reference_trbc.fa','w')
f.write('>WT_TRB\n'+wt_trbc+'\n')
f.write('>Ref_HDR_CD22\n'+hdr_seq_cd22+'\n')
f.write('>Ref_ITR_CD22\n'+itr_seq_cd22+'\n')
f.write('>Ref_ITR_CD19\n'+itr_seq_cd19_to_trbc+'\n')
f.close()

# Truncated HDR: pcr - genomic - HA - kozak - transgene(random truncation) - polyA - HA - genomic - pcr
import random
def generate_hdr_truncation_templates(transgene_seq, transgene_name, construct_parts, n_templates=20):
    """Generate HDR templates with truncated transgenes only"""
    
    transgene_length = len(transgene_seq)
    templates = []
    
    # Truncation range: 10-94% for "Imperfect HDR" classification
    min_length = int(transgene_length * 0.10)  
    max_length = int(transgene_length * 0.90)  
    
    # Generate 20 random lengths
    truncation_lengths = []
    random.seed(42)  # Reproducible
    
    for i in range(n_templates * 2):  # Generate extra to ensure uniqueness
        length = random.randint(min_length, max_length)
        truncation_lengths.append(length)
    
    # Remove duplicates, sort, take first 20
    truncation_lengths = sorted(list(set(truncation_lengths)))[:n_templates]
    
    # Generate templates (mix 5' and 3' truncations)
    for i, length in enumerate(truncation_lengths):
        if i % 2 == 0:  # 5' truncation (keep beginning)
            trunc_transgene = transgene_seq[:length]
            direction = "5prime"
        else:  # 3' truncation (keep end)
            trunc_transgene = transgene_seq[-length:]
            direction = "3prime"
        
        # Assemble: FULL context + TRUNCATED transgene + FULL context
        template_seq = (construct_parts['prefix'] + 
                       trunc_transgene +  # Only transgene is truncated
                       construct_parts['suffix'])
        
        template_name = f"{transgene_name}_imperfect_HDR_{direction}_{length}bp"
        templates.append((template_name, template_seq))
    
    return templates

# Generate only HDR truncations
trac_hdr_parts={
    'prefix': allseq['TRAC_PCR_left'] + allseq['TRAC_Genomic_left'] + allseq['TRAC_HA_left'] + allseq['kozak'] ,
    'suffix': allseq['polyA'] + allseq['TRAC_HA_right'] + allseq['TRAC_Genomic_right'] + allseq['TRAC_PCR_right']
}
trbc_hdr_parts = {
    'prefix': allseq['TRBC_PCR_left'] + allseq['TRBC_Genomic_left'] + allseq['TRBC_HA_left'] + allseq['kozak'],
    'suffix':allseq['polyA'] + allseq['TRBC_HA_right'] + allseq['TRBC_Genomic_right'] + allseq['TRBC_PCR_right']
}

cd19_imperfect_hdr = generate_hdr_truncation_templates(allseq['CD19'], 'CD19', trac_hdr_parts, 20)
cd22_imperfect_hdr = generate_hdr_truncation_templates(allseq['CD22'], 'CD22', trbc_hdr_parts, 20)

# Save templates
f=open('template_cd19_truncated_hdr.fa','w')
for name, seq in cd19_imperfect_hdr:
    f.write(f">{name}\n{seq}\n")
f.close()
f=open('template_cd22_truncated_hdr.fa','w')
for name, seq in cd22_imperfect_hdr:
    f.write(f">{name}\n{seq}\n")
f.close()
print(f"\nGenerated {len(cd19_imperfect_hdr + cd22_imperfect_hdr)} imperfect HDR templates")
