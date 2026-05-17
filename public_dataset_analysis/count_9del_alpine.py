sampleinfo=open('../PRJNA913199_SRA_Run_table.csv','r').read().split('\n')[1:-1]
pacbio={}
for c in sampleinfo:
    c=c.split(',')
    cla=c[23]

    if 'HBB Unmodified' in cla:
        pacbio[c[0]]='Unmodified'
    elif 'HBB UV' in cla:
        pacbio[c[0]]='UV'
    elif 'HBB 10M' in cla:
        pacbio[c[0]]='10M'
    elif 'HBB 4M' in cla:
        pacbio[c[0]]='4M'
    elif 'HBB Untreated' in cla:
        pacbio[c[0]]='Untreated'
    elif 'HBB 200M' in cla:
        pacbio[c[0]]='200M'
    elif 'HBB' not in cla:
        continue
    else:
        print(c);quit()



allsample=open('filelist','r').read().split('\n')[:-1]
f=open('Merged_count_9bpdel_alpine.txt','w')
for filename in allsample:
    sample=filename.split('_')[1].split('.')[0]
    read=open(filename,'r').read().split('\n')[:-1]
    ndel=0
    for c in read:
        c=c.split('\t')
        if c[1]=='DEL-small' and c[2]=='9':
            ndel+=1
    f.write('\t'.join(['ALPINE',sample,pacbio[sample],str(ndel),str(len(read)),str(ndel/len(read)*100)])+'\n')
f.close()
