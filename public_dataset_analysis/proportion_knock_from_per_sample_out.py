sampleinfo=open('../PRJNA913199_SRA_Run_table.csv','r').read().split('\n')[1:-1]
pacbio={}
for c in sampleinfo:
    c=c.split(',')
    cla=c[23]
    pacbio[c[0]]=cla


allsample=open('filelist','r').read().split('\n')[:-1]
f=open('Merged_count_class_knock.txt','w')
for filename in allsample:
    sample=filename.split('_')[0]
    read=open(filename,'r').read().split('\n')[1:-1]
    sdel=0
    ldel=0
    ins=0
    hdr=0
    other=0
    indel=0
    wt=0

    for c in read:
        c=c.split('\t')
        if c[5]=='deletion <50 nt':
            sdel+=1; continue
        if c[5]=='deletion >=50 nt':
            ldel+=1; continue
        if c[5]=='WT':
            wt+=1; continue
        if c[5]=='insertion':
            ins+=1; continue
        if c[4] in ['incomplete HDR','HDR','donor fragment']:
            hdr+=1;continue
        if c[4] =='complex indel':
            indel+=1;continue
        print(c[4])
        other+=1

    total=sdel+ldel+ins+hdr+other+wt+indel
    sdel=str(sdel/total)
    ldel=str(ldel/total)
    ins=str(ins/total)
    indel=str(indel/total)
    hdr=str(hdr/total)
    wt=str(wt/total)
    other=str(other/total)

    f.write('\t'.join([sample,pacbio[sample],'knock-knock',sdel,ldel,'0',ins,indel,hdr,wt,other])+'\n')
f.close()
