data='hifi'
out=open(data+'_trac.out','r').read().split('\n')[1:-1]

count={}

for c in out:
    true=c.split('\t')[0].split('.')[0]
    if data=='ont':
        true=true[11:]

    if 'itr_oneside' in true:
        true=true.replace('itr_oneside','itr')

    clas=c.split('\t')[4]
    tag=c.split('\t')[5]

    pair=true+'\t'+tag
#    if clas in ['']
    if pair not in count:
        count[pair]=1
    else:
        count[pair]+=1


f=open(data+'_trac.pairs','w')
for c in count:
    f.write(c+'\t'+str(count[c])+'\n')
f.close()





class_match_trac={
    "simulated_trac_cd19_hdr": "HDR",
    "simulated_trac_cd19_itr": "complex misintegration",
    "simulated_trac_cd19_truncated_hdr": "deletion in donor", 
    "simulated_trac_cd22_hdr": "complex misintegration", 
    "simulated_trac_cd22_itr": "complex misintegration",
    "simulated_trac_del_large": "deletion >=50 nt", 
    "simulated_trac_del_small": "deletion <50 nt",
    "simulated_trac_ins_large": "insertion",
    "simulated_trac_ins_small": "insertion",
    "simulated_trac_snp": "WT",
    "simulated_trac_wt": "WT",
}


class_match_trbc={
    "simulated_trbc_cd19_hdr": "complex misintegration",
    "simulated_trbc_cd19_itr": "complex misintegration",
    "simulated_trbc_cd22_truncated_hdr": "deletion in donor",
    "simulated_trbc_cd22_hdr": "HDR",
    "simulated_trbc_cd22_itr": "complex misintegration",
    "simulated_trbc_del_large": "deletion >=50 nt",
    "simulated_trbc_del_small": "deletion <50 nt",
    "simulated_trbc_ins_large": "insertion",
    "simulated_trbc_ins_small": "insertion",
    "simulated_trbc_snp": "WT",
    "simulated_trbc_wt": "WT",
}

knock_class=[]
for c in class_match_trac:
    knock_class+=[class_match_trac[c]]

recall={}
for group in class_match_trac:
    tp=0
    fn=0
    for c in count:
        if c.split('\t')[0]==group:
            if c.split('\t')[1]==class_match_trac[group]:
                tp+=count[c]
            else:
                fn+=count[c]

    recall[group]=(tp,fn)
    if tp+fn!=5000 and tp+fn!=10000:
        print(group,tp,fn)
print(recall)


precision={}
for group in knock_class:
    tp=0
    fp=0
    for c in count:
        if c.split('\t')[1]==group:
            if group==class_match_trac[c.split('\t')[0]]:
                tp+=count[c]
            else:
                fp+=count[c]
    precision[group]=(tp,fp)

print(precision)

f=open(data+'_trac.tpfp','w')
for group in class_match_trac:
    (tp,fn)=recall[group]
    (group_tp,group_fp)=precision[class_match_trac[group]]
    class_recall=tp/(tp+fn)*100
    class_precision=group_tp/(group_tp+group_fp)*100
    f1=2*class_recall*class_precision/(class_precision+class_recall)

    f.write('\t'.join([group,class_match_trac[group],str(tp),str(group_fp),str(fn),str(tp+fn),str(class_precision),str(class_recall),str(f1)])+'\n')
f.close()





#TRBC

out=open(data+'_trbc.out','r').read().split('\n')[1:-1]

count={}

for c in out:
    true=c.split('\t')[0].split('.')[0]
    if data=='ont':
        true=true[11:]

    if 'itr_oneside' in true:
        true=true.replace('itr_oneside','itr')
    clas=c.split('\t')[4]
    tag=c.split('\t')[5]
    pair=true+'\t'+tag
    if pair not in count:
        count[pair]=1
    else:
        count[pair]+=1

recall_trbc={}
precision_trbc={}

for group in class_match_trbc:
    tp=0
    fn=0
    for c in count:
        if c.split('\t')[0]==group:
            if c.split('\t')[1]==class_match_trbc[group]:
                tp+=count[c]
            else:
                fn+=count[c]

    recall_trbc[group]=(tp,fn)
    if tp+fn!=5000 and tp+fn!=10000:
        print(group,tp,fn)


for group in knock_class:
    tp=0
    fp=0
    for c in count:
        if c.split('\t')[1]==group:
            if group==class_match_trbc[c.split('\t')[0]]:
                tp+=count[c]
            else:
                fp+=count[c]
    precision_trbc[group]=(tp,fp)

f=open(data+'_trbc.tpfp','w')
for group in class_match_trbc:
    (tp,fn)=recall_trbc[group]
    (group_tp,group_fp)=precision_trbc[class_match_trbc[group]]
    class_recall=tp/(tp+fn)*100
    class_precision=group_tp/(group_tp+group_fp)*100
    f1=2*class_recall*class_precision/(class_precision+class_recall)

    f.write('\t'.join([group,class_match_trbc[group],str(tp),str(group_fp),str(fn),str(tp+fn),str(class_precision),str(class_recall),str(f1)])+'\n')
f.close()

