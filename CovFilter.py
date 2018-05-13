#! /usr/bin/python
from sys import argv

info = open(argv[1],'r');
fil = open(argv[2],'r');


loci = [];
fil.readline();
for line in fil:
    loc = re.search("([0-9]+)",line).group(1);
    loci.append(loc);

print('CHROM\tPOS\tREF\tALT\tSAMPLE\tMASK\tREP\tCOV\tCOVFIL')
info.readline();
calls=info.readline();
print(calls[:-1]);


for line in info:

    result = line[:-1].split('\t');
    chrom, pos, ref, alt, sample, mask, rep, cov= result[0], int(result[1]), result[2], result[3], result[4], result[5], result[6], result[7];
    covfil = False;

    if rep == 'True':
        rep_count += 1;
        if loci != []:
            if int(loci[0]) == rep_count:
                covfil = True;
                del loci[0];

        
    print(chrom+'\t'+str(pos)+'\t'+ref+'\t'+alt+'\t'+sample+'\t'+mask+'\t'+rep+'\t'+str(cov)+'\t'+str(hwp)+'\t'+pat+'\t'+str(covfil));

