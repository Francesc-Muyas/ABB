#from multiprocessing import Pool
from multiprocessing import Pool, cpu_count
import timeit
import argparse
import os
import warnings
import pandas as pd
import numpy
import psutil
from Bio import bgzf

start = timeit.default_timer()

CORES=psutil.cpu_count(logical=False)

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
        
def ABB_ANNOTATE(variant):
    columns= variant.split('\t')

    CHROM = str(columns[0])
    POS = str(columns[1])
    REF = columns[3]
    ALT = columns[4]
    FILTER = columns[6]
    INFO = columns[7]
    form = columns[8].split(':')
            
    # Y represents the file with ABB scores opened as Y
    C = Y.ABB.values[(Y.CHROM.values == CHROM) & (Y.POS.values == POS)]
            
    if (len(C) > 0 and str(C[0]) != "nan"):
        ABB = C[0]
    else:
        ABB = "."
                
            
    FILT = ABB_FILTER(ABB,0.7,FILTER)
            
    #print CHROM, POS, ABB, FILTER   
            
    columns[6] = FILT
    columns[7] = INFO + ";ABB=" + str(ABB)        
               
    LINE = '\t'.join(columns)+'\n'
    return(LINE)


tmp = open(infile,'r')
magic_number = tmp.read(2)
print 'MESSAGE:: BGZIP file  : %s'%(magic_number=='\x1f\x8b')
tmp.close()

## ABB file
Y=pd.read_csv(ABB_FILE,dtype=object,delimiter="\t")

# list of tasks
tasks=[]
count_line=0
with open(infile) if magic_number!='\x1f\x8b' else bgzf.open(infile) as VCF:
    for variant in VCF:
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
            
        if not variant.startswith('#'): #to select only the variants and not the other rows (## and #CHROM)
            if (count_line < 10000 and variant is not None):
                tasks.append(variant)
                count_line = count_line + 1
                
            else:
                tasks.append(variant)

                if __name__ == '__main__':
                    p = Pool(CORES)
                    A = (p.map(ABB_ANNOTATE, (tasks)))
                    for i in A:
                        OF1.write(i)
                
                ## Re-start the counting. This step is mainly done to avoid a big usage of memory when vcf input file is super big.
                tasks=[]
                count_line=0
    
    if len(tasks) > 0:    
        if __name__ == '__main__':
            p = Pool(CORES)
            A = (p.map(ABB_ANNOTATE, (tasks)))
            for i in A:
                OF1.write(i)


#VCF.close()
OF1.close()


stop = timeit.default_timer()
print stop - start
