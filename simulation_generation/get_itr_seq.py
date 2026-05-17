a=open('aav2_complete_genome.fasta','r').read().split('\n')[1:-1]
seq=''.join(a)

print(len(seq))

f=open('itr1_seq.txt','w')
f.write(seq[:145])
f.close()
f=open('itr2_seq.txt','w')
f.write(seq[4534:4679])
f.close()


