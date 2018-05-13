#! /usr/bin/python

import re
from sys import argv
from Bio import SeqIO

rm={};
wm={};
combine={};

for seq in SeqIO.parse(argv[1], 'fasta'):
    rm[seq.id] = str(seq.seq);

for seq in SeqIO.parse(argv[2], 'fasta'):
    wm[seq.id] = str(seq.seq);



for seq in wm:
    combine[seq] = re.sub('[atgcnN]','1',rm[seq]);
    combine[seq] = list(re.sub('[ATGC]','0',combine[seq]));
    wm[seq] = re.sub('[atgcnN]','1',wm[seq]);

    for mask in re.finditer('1+',wm[seq]):
        combine[seq][mask.start():mask.end()]=list(mask.group());

    print('>'+seq+'\n'+''.join(combine[seq]));
