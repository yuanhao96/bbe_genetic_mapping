#! /usr/bin/python

from Bio import SeqIO
from sys import argv

info = open(argv[1],'r');

mask={};
for seq in SeqIO.parse(argv[2], 'fasta'):
    mask[seq.id] = str(seq.seq);

print('CHROM\tPOS\tREF\tALT\tSAMPLE\tMASK\n')
info.readline();
calls=info.readline();
print(calls[:-1]);


for line in info:
    result = line[:-1].split('\t');
    chrom, pos, ref, alt, sample = result[0], int(result[1]), result[2], result[3], result[4];
    if '1' in mask[chrom][pos-1-int(argv[3]):pos+int(argv[3])] :
        masking = False
    else:
        masking = True
    print(chrom+'\t'+str(pos)+'\t'+ref+'\t'+alt+'\t'+sample+'\t'+str(masking));

