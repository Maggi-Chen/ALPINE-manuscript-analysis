sampleinfo=open('PRJNA913199_SRA_Run_table.csv','r').read().split('\n')[1:-1]
pacbio={}
for c in sampleinfo:
    c=c.split(',')
    cla=c[23]
    pacbio[c[0]]=cla
    continue

    if 'HBB Unmodified' in cla:
        pacbio[c[0]]='HBB_Unmodified'
    elif 'HBB UV' in cla:
        pacbio[c[0]]='HBB_UV'
    elif 'HBB 10M' in cla:
        pacbio[c[0]]='HBB_10M'
    elif 'HBB 4M' in cla:
        pacbio[c[0]]='HBB_4M'
    elif 'HBB Untreated' in cla:
        pacbio[c[0]]='HBB_Untreated'
    elif 'HBB 200M' in cla:
        pacbio[c[0]]='HBB_200M'
    elif 'HBB' not in cla:
        continue
    else:
        print(c);quit()

ill=open('PRJNA913199_SRA_Run_table_Illumina.csv','r').read().split('\n')[1:-1]
hbb=[c for c in ill if 'HBB' in c]
illu={}
for c in hbb:
    sam=c.split(',')[0]
    cla=c.split(',')[24]
    illu[sam]=cla
    continue

    if 'HBB E7V Unmodified' in cla:
        illu[sam]='HBB_E7V_Unmodified'
    elif 'HBB E7V UV' in cla:
        illu[sam]='HBB_E7V_UV'
    elif 'HBB E7V 30X' in cla:
        illu[sam]='HBB_E7V_30X'
    elif 'HBB E7V 10X' in cla:
        illu[sam]='HBB_E7V_10X'
    elif 'HBB E7V 3X' in cla:
        illu[sam]='HBB_E7V_3X'
    elif 'HBB E7V 1X' in cla:
        illu[sam]='HBB_E7V_1X'
    elif 'HBB E7V 0.3X' in cla:
        illu[sam]='HBB_E7V_0.3X'
    elif 'HBB E7V 0.1X' in cla:
        illu[sam]='HBB_E7V_0.1X'
    elif 'HBB E7V 0.03X' in cla:
        illu[sam]='HBB_E7V_0.03X'
    elif 'HBB E7V 0.01X' in cla:
        illu[sam]='HBB_E7V_0.01X'
    elif 'HBB E7V EtOH' in cla:
        illu[sam]='HBB_E7V_EtOH'
    else:
        print(c);quit()

a=open('Merged_read_classification.txt','r').read().split('\n')[1:-1]
alpine={}

#f=open('alpine_del_freq.txt','w')
f=open('combined_all_full_group.txt','w')
f.write('\t'.join(['Sample','Group','Tool','DEL-Small','DEL-Large','DEL','INS','Indel','HDR','WT','Other'])+'\n')
for c in a:
    c=c.split('\t')
    num=c[1:]
    num=[int(m) for m in num]
    total=sum(num)
    delsmall=str(num[2]/total)
    #aldel=(num[2]+num[4])/total
    dellarge=str(num[4]/total)
    ins=str((num[3]+num[5])/total)
    hdr=str((num[6]+num[7])/total)
    wt=str((num[0]+num[1])/total)
    other=str((num[8]+num[9]+num[10])/total)

    if sum([float(x) for x in [delsmall,dellarge,ins,hdr,wt,other]])!=1:
        print(c[0],sum([float(x) for x in [delsmall,dellarge,ins,hdr,wt,other]]))

    f.write('\t'.join([c[0],pacbio[c[0]],'ALPINE',delsmall,dellarge,'0',ins,'0',hdr,wt,other])+'\n')
    #f.write(c[0]+'\t'+str(delfreq)+'\t'+str(aldel)+'\n')
    #alpine[c[0]]=delfreq
#f.close()


b=open('knock_knock_per_sample_with_HDR.csv','r').read().split('\n')[1:-1] ### This result is not used as merged table from knock-knock used "indel" class that combined insertions and small/large deletions into one group. Per-sample output files were used to count reads from sub-groups for comparison with other tools.

knock={}
for c in b:
    c=c.split(',')
    num=c[3:]
    num=[int(m) for m in num]
    total=int(c[2])
    #wt=str(num[]/total)
    indel=str((num[1]+num[2])/total)
    hdr=str((num[3]+num[4]+num[5])/total)
    wt=str(num[0]/total)
    other=str((num[6]+num[7]+num[8])/total)
    if sum([float(x) for x in [indel,hdr,wt,other]])!=1:
        print(c[0],sum([float(x) for x in [indel,hdr,wt,other]]))

    f.write('\t'.join([c[0],pacbio[c[0]],'knock-knock','0','0','0','0',indel,hdr,wt,other])+'\n')

d=open('crispresso2_per_sample_summary.csv','r').read().split('\n')[1:-1]
cris={}
for c in d:
    c=c.split(',')
    num=[c[7],c[9],c[10],c[11],c[12]]
    num=[int(x) for x in num]
    total=sum(num)
    delfreq=num[2]/total
    insfreq=num[1]/total
    wt=(num[0]+num[3])/total
    other=num[4]/total
    if sum([delfreq,insfreq,wt,other])!=1:
        print(c[0],sum([delfreq,insfreq,wt,other]))

    f.write('\t'.join([c[0],illu[c[0]],'crispresso2','0','0',str(delfreq),str(insfreq),'0','0',str(wt),str(other)])+'\n')



f.close()
