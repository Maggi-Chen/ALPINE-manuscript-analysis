
al=open('readname_SRR22793524.txt','r').read().split('\n')[1:-1]

kn=open('SRR22793524_knock_out','r').read().split('\n')[1:-1]


hdr=[]
alread={}
for c in al:
    cc=c.split('\t')
    name=cc[0].split('.')[1]
    alread[name]=cc[1]
    if 'HDR' in c:
        hdr+=[name]

kn_hdr=[]
knread={}
for c in kn:
    cc=c.split('\t')
    name=str(int(cc[0].split(':')[1]))
    knread[name]=cc[4]
    if 'HDR' in c:
        kn_hdr+=[name]


overlap=[c for c in hdr if c in kn_hdr]
aluniq=[c for c in hdr if c not in kn_hdr]
print(len(hdr),len(kn_hdr),len(overlap))

f=open('SRR22793524_read_both_hdr','w')
for c in overlap:
    f.write('SRR22793524.'+c+'\n')
f.close()

f=open('SRR22793524_read_alpine_hdr','w')
for c in hdr:
    if c not in kn_hdr:
        f.write('SRR22793524.'+c+'\n')
f.close()

print(overlap[0])
print(aluniq[0])


f=open('SRR22793524.sam','r')

fb=open('HDR_Both.sam','w')
fa=open('ALPINE_Unique_HDR.sam','w')

a=f.readline()
while a!='':
    if a[0]=='@':
        fa.write(a)
        fb.write(a)
        a=f.readline()
        continue

    if a.split('\t')[1] in ['256','272']:
        a=f.readline()
        continue
    name=a.split('\t')[0]
    try:
        readname=name.split('.')[1]
        
    except:
        print(a)
        quit()
    if readname in overlap:
        fb.write(a)
    if readname in aluniq:
        fa.write(a)
    a=f.readline()
f.close()
fa.close()
fb.close()
