a=open('/scratch/reference/hg38.fa','r').read().split('>')[1:]

chr14=a[6]
print(chr14.split('\n')[0])

chr7=a[56]
print(chr7.split('\n')[0])



trac_site=22547650
trbc_site=142791900


chr14seq=''

for c in chr14.split('\n')[1:]:
    chr14seq+=c

trac_ha_1=chr14seq[22547650-100:22547650]
trac_ha_2=chr14seq[22547650:22547650+100]

trac_pcr_1=chr14seq[22547650-325:22547650-300]
trac_pcr_2=chr14seq[22547650+300:22547650+325]

trac_genomic_1=chr14seq[22547650-300:22547650-100]
trac_genomic_2=chr14seq[22547650+100:22547650+300]


chr7seq=''
for c in chr7.split('\n')[1:]:
    chr7seq+=c
trbc_ha_1=chr7seq[142791900-100:142791900]
trbc_ha_2=chr7seq[142791900:142791900+100]

trbc_pcr_1=chr7seq[142791900-325:142791900-300]
trbc_pcr_2=chr7seq[142791900+300:142791900+325]

trbc_genomic_1=chr7seq[142791900-300:142791900-100]
trbc_genomic_2=chr7seq[142791900+100:142791900+300]

f=open('genomic_seq.txt','w')
f.write('TRAC_HA_left\t'+trac_ha_1+'\n')
f.write('TRAC_HA_right\t'+trac_ha_2+'\n')
f.write('TRAC_PCR_left\t'+trac_pcr_1+'\n')
f.write('TRAC_PCR_right\t'+trac_pcr_2+'\n')
f.write('TRAC_Genomic_left\t'+trac_genomic_1+'\n')
f.write('TRAC_Genomic_right\t'+trac_genomic_2+'\n')

f.write('TRBC_left\t'+trbc_ha_1+'\n')
f.write('TRBC_right\t'+trbc_ha_2+'\n')
f.write('TRBC_PCR_left\t'+trbc_pcr_1+'\n')
f.write('TRBC_PCR_right\t'+trbc_pcr_2+'\n')
f.write('TRBC_Genomic_left\t'+trbc_genomic_1+'\n')
f.write('TRBC_Genomic_right\t'+trbc_genomic_2+'\n')

f.close()



