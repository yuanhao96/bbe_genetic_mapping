#! /usr/bin/python
from sys import argv
from scipy import stats
import re

info = open(argv[1],'r');

print('CHROM\tPOS\tREF\tALT\tSAMPLE\tMASK\tREP\tCOV\tCOVFIL\tSEG')
info.readline();
calls=info.readline();
print(calls[:-1]);


for line in info:

    result = line[:-1].split('\t')
    chrom, pos, ref, alt, sample, mask, rep, cov, covfil = result[0], result[1], result[2], result[3], result[4], result[5], result[6], result[7], result[8];
    samples = re.findall("\(.*?\)",sample);
    mdl = 'NA';

    if rep == 'True' and mask == 'True':
        gtpool = '';
        alleles = re.findall('[A-Z]+',alt);
        alleles.append(ref);

        for i in samples:
            if re.search('([0-9][/\|][0-9])',i):
                gt = re.search('([0-9][/\|][0-9])',i).group(1);
                gt = re.findall('[0-9]',gt);
                gtpool += ''.join(gt);
                if len(''.join(gt))==2:
                    gno.append(''.join(gt))

        whole = len(gtpool);
        a = gtpool.count('0')/whole;
        b = gtpool.count('1')/whole;
        c = gtpool.count('2')/whole;
        d = gtpool.count('3')/whole;
        
        aa = gno.count('00')
        ab = gno.count('01')+gno.count('10')
        ac = gno.count('02')+gno.count('20')
        ad = gno.count('03')+gno.count('30')
        bb = gno.count('11')
        bc = gno.count('12')+gno.count('21')
        bd = gno.count('13')+gno.count('31')
        cc = gno.count('22')
        cd = gno.count('23')+gno.count('32')
        dd = gno.count('33')
        

        if len(alleles) == 4:
            if min(a*a*51,a*b*51,a*c*51,a*d*51,b*b*51,b*c*51,b*d*51,c*c*51,c*d*51,d*d*51) < 5:
                alleles = '333';
            else:
                obs = [aa,ab,ac,ad,bb,bc,bd,cc,cd,dd];
                exp = [a*a*51,a*b*51,a*c*51,a*d*51,b*b*51,b*c*51,b*d*51,c*c*51,c*d*51,d*d*51];
                x2 = stats.chisquare(obs, f_exp = exp);
                mdl = x2[1];


        if len(alleles) == 3:
            if min(a*a*51,a*b*51,a*c*51,b*b*51,b*c*51,c*c*51) < 5:
                alleles = '22';
            else:
                obs = [aa,ab,ac,bb,bc,cc];
                exp = [a*a*51,a*b*51,a*c*51,b*b*51,b*c*51,c*c*51];
                x2 = stats.chisquare(obs, f_exp = exp);
                mdl = x2[1];

        if len(alleles) == 2:
            obs = [aa,ab,bb];
            exp = [a*a*51,a*b*51,b*b*51];
            x2 = stats.chisquare(obs, f_exp = exp);
            mdl = x2[1]

    print(chrom+'\t'+pos+'\t'+ref+'\t'+alt+'\t'+sample+'\t'+mask+'\t'+rep+'\t'+covfil+'\t'+str(mdl))
