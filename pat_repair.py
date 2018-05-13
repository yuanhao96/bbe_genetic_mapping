#! /usr/bin/python

import re
from sys import argv

info = open(argv[1],'r');



def pat_freq(sample):
    samples = re.findall('\(.*?\)',sample);
    pool = '';
    paternal = [];
    freq=[];

    for i in samples:
        gt = str(re.findall('[0-9][/\|][0-9]',i));
        gene = re.findall('[0-9]+',gt);
        pool += ''.join(gene);
        alleles = re.findall('[A-Z]+',alt);
        alleles.append(ref);

    if len(alleles) == 2:
        whole = len(pool);
        X2 = (((whole/4)-pool.count(min(pool,key=lambda x : pool.count(x))))**2)/(whole/4) + (((whole*3/4)-pool.count(max(pool,key=lambda x : pool.count(x))))**2)/(whole*3/4);
        if X2 < 3.84:
            paternal.append(max(pool,key=lambda x : pool.count(x)));
            if int(mgene[0])-int(mgene[1]) == 0:
                for i in '01':
                    if i not in mgene:
                        paternal.append(i);
            else:
                paternal.append(max(pool,key=lambda x : pool.count(x)));
        else:
            for i in mgene:
                complement = str(1-int(i));
                paternal.append(complement);
                

    if len(alleles) == 3:
        for i in '012':
            if i not in mgene:
                paternal.append(i);
            freq.append(i+':'+str(pool.count(i)));
        if len(paternal) < 2:
            paternal.append(max(pool,key=lambda x : pool.count(x)));


    if len(alleles) == 4:
        for i in '0123':
            if i not in mgene:
                paternal.append(i);
            freq.append(i+':'+str(pool.count(i)));

    return ('/'.join(paternal),freq);

print('CHROM\tPOS\tREF\tALT\tSAMPLE\tMASK\tREP\tCOV\tCOLFIL\tSEG\tPAT')
info.readline();
calls=info.readline();
print(calls[:-1]);


for line in info:

    result = line[:-1].split('\t')
    chrom, pos, ref, alt, sample, mask, rep, cov, covfil, mdl = result[0], result[1], result[2], result[3], result[4], result[5], result[6], result[7], result[8], result[9];
    samples = re.findall("\(.*?\)",sample);

    maternal = re.search("\('maternal'.+?\)",sample);
    mgt = str(re.findall('[0-9][/\|][0-9]',maternal.group()));
    mgene = re.findall('[0-9]+',mgt);
    mat = '/'.join(mgene);
    pat = 'NA'

    pat = pat_freq(sample);

    print(chrom+'\t'+pos+'\t'+ref+'\t'+alt+'\t'+sample+'\t'+mask+'\t'+rep+'\t'+cov+'\t'+covfil+'\t'+seg+'\t'+pat[0]);
