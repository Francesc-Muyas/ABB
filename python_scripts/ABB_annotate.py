#!/usr/bin/env python

# python extract_info_VCF.py -infile /users/so/fmuyas/PhD/Part_1/data/200_extra_independent.txt -out /users/so/fmuyas/PhD/Part_1/data/200_extra_independent_INFO.txt
#import re
import timeit
import argparse
#import sys
import os
import warnings
import pandas as pd
import numpy
from Bio import bgzf


start = timeit.default_timer()


parser = argparse.ArgumentParser()
parser.add_argument('-i','--input_file', type=str, help='Input vcf file', required= True)
parser.add_argument('-outfile', '--outfile', type=str, help='Out vcf file', required= True)
parser.add_argument('-abb_filter', '--abb_filter', type=float, default= 0.7, help='Cutoff for ABB filter', required = False)
parser.add_argument('-abb_file', '--abb_file', type=str, help='ABB score file', required = True)

args = parser.parse_args()

infile=args.input_file
outfile=args.outfile
cutoff=args.abb_filter
ABB_FILE=args.abb_file

of1= (outfile)
OF1=open(of1,'w')


def ABB_FILTER(ABB,cutoff,FILTER):
    if (ABB != '.' and float(ABB) >= cutoff):
        return("ABB_biased")
    if (ABB == '.' or float(ABB) < cutoff):
        return(FILTER)
        

tmp = open(infile,'r')
magic_number = tmp.read(2)
print 'MESSAGE:: BGZIP file  : %s'%(magic_number=='\x1f\x8b')
tmp.close()

## ABB file
#A=pd.read_csv("/users/so/fmuyas/programs/ABB/rscripts/test2/test.txt",dtype=object,delimiter="\t")
A=pd.read_csv(ABB_FILE,dtype=object,delimiter="\t")

with open(infile) if magic_number!='\x1f\x8b' else bgzf.open(infile) as VCF:
    for variant in VCF: #to extract the information of each variant/row
        variant = variant.rstrip('\n')
        if variant.startswith('##'): #to select the name of each column
            OF1.write(variant+'\n')
            
        if variant.startswith('#CHROM'): #to select the name of each column
            name_columns= variant.split('\t')
            name_columns= name_columns[0:len(name_columns)]#there is something at the end of the line and then we remove it
            ##INFO=<ID=ABB_SCORE,Number=1,Type=Float,Description="Allele balance bias (ABB) value">
            OF1.write('##INFO=<ID=ABB,Number=1,Type=Float,Description="Allele balance bias (ABB) value">\n')
            OF1.write('##FILTER=<ID=ABB_biased,Description="Prone to systematic error based on ABB score">\n')
            OF1.write(variant+'\n')
            #print name_columns[len(name_columns)-1]
          
        if not variant.startswith('#'): #to select only the variants and not the other rows (## and #CHROM)
            columns= variant.split('\t')
            #print columns
            #samples= name_columns[9:len(name_columns)-1]
    #        print columns
            CHROM = str(columns[0])
            POS = str(columns[1])
            REF = columns[3]
            ALT = columns[4]
            FILTER = columns[6]
            INFO = columns[7]
            form = columns[8].split(':')
            
            C = A.ABB.values[(A.CHROM.values == CHROM) & (A.POS.values == POS)]
            #print C, len(C)
            
            if (len(C) > 0 and str(C[0]) != "nan"):
                ABB = C[0]
            else:
                ABB = "."
                
            
            FILT = ABB_FILTER(ABB,0.7,FILTER)
            
            #print CHROM, POS, ABB, FILTER   
            
            columns[6] = FILT
            columns[7] = INFO + ";ABB=" + str(ABB)        
            
            OF1.write('\t'.join(columns)+'\n')
            
        
#VCF.close()
OF1.close()


stop = timeit.default_timer()
print stop - start
