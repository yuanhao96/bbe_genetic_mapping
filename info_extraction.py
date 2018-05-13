#! /usr/bin/python

from vcf import Reader
from sys import argv

path=argv[1]
reader=Reader(filename=path)

print('CHROM\tPOS\tREF\tALT\tSAMPLE\tINDEL\tHW\tMASKED\tINFORM\tREP\tA_INDEL\tA_SNP\tANNO')
calls=len(reader.samples)
output.write('calls:'+str(calls)+'\n')

for record in reader:
    info={'CHROM':record.CHROM, 'POS':record.POS, 'REF':record.REF, 'ALT':record.ALT, 'SAMPLE':[], 'INDEL':'', 'HW':'', 'MASKED':'', 'INFORM':'', 'REP':'', 'A_INDEL':'', 'A_SNP':'', 'ANNO':''}
    for sample in record.samples:
        if len(sample.data)>1:
            info['SAMPLE'].append((sample.sample,sample.data[0],sample.data[1]))
        else:
            info['SAMPLE'].append((sample.sample,sample.data[0],None))
    print(str(info['CHROM'])+'\t'+str(info['POS'])+'\t'+str(info['REF'])+'\t'+str(info['ALT'])+'\t'+str(info['SAMPLE'])+'\t'+str(info['INDEL'])+'\t'+str(info['HW'])+'\t'+str(info['MASKED'])+'\t'+str(info['INFORM'])+'\t'+str(info['REP'])+'\t'+str(info['A_INDEL'])+'\t'+str(info['A_SNP'])+'\t'+info['ANNO'])
