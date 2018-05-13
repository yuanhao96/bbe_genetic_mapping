#! /usr/bin/python

from sys import argv
import re

info = open(argv[1],'r');
print('CHROM\tPOS\tREF\tALT\tSAMPLE\tMASK\tREP\tCOV')
info.readline();
calls=info.readline();
print(calls[:-1]);


for line in info:
    #info extraction
    result = line[:-1].split('\t');
    chrom, pos, ref, alt, sample, mask = result[0], int(result[1]), result[2], result[3], result[4], result[5];
    samples = re.findall("\('(?:maternal|off05|off06|off07|off_L1_83|off_L8_113|xL3_68|xL7_80).*?\)",sample);
    
    #variable claiming
    rep = True;
    alleles = re.findall('[A-Z]+',alt);
    alleles.append(ref);
    gtpool = '';
    cov = 'NA';



    #rep_test
    for i in samples:

        if not re.search('([0-9][/\|][0-9])',i):
            rep = False;
        else:
            gt = str(re.search('([0-9][/\|][0-9])',i).group(1));
            gtpool += ''.join(gt);


    if len(alleles) == 2:
        for a in '01':
            if a not in gtpool:
                rep = False;
    elif len(alleles) == 3:
        for a in '012':
            if a not in gtpool:
                rep = False;
    elif len(alleles) == 4:
        for a in '012':
            if a not in gtpool:
                rep = False;

    #cov_calculate
    if rep:
        cov = 0;
        for i in samples:
            dp = int(re.search('([0-9]+)\)',i).group(1));
            cov += dp;
        
    print(chrom+'\t'+str(pos)+'\t'+ref+'\t'+alt+'\t'+sample+'\t'+mask+'\t'+str(rep)+'\t'+str(cov))

