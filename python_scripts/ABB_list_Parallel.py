#!/usr/bin/env python

from multiprocessing import Pool, cpu_count
import re
import timeit
import argparse
import rpy2.robjects as robjects
import psutil
import os
import warnings
from Bio import bgzf

start = timeit.default_timer()

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_file', type=str, help='Input vcf file', required= True)
parser.add_argument('-outfile', '--outfile', type=str, help='Out ABB list', required= True)
#parser.add_argument('-abb', type=float, default= 0.7, help='Cutoff for ABB filter', required = False)
#parser.add_argument('-gender', type=str, default= "NA", help='Geneder of the samples: M for Males ; F for Females', required = False)

args = parser.parse_args()

infile=args.input_file
outfile=args.outfile
#gender=args.gender

CORES=psutil.cpu_count(logical=False)

# R functions
DIR = os.path.dirname(__file__)
filename = os.path.join(DIR,"../rscripts/ABB_function.r")
source = os.path.join(DIR,"../source/source.Rdata")

filename = os.path.realpath(filename)
source = os.path.realpath(source)

#print DIR
print source
print filename

r=robjects.r
r.source(filename)
#ro.r("""source('filename.R')""")
#multiprocessing.cpu_count()
def LOADING_SOURCE(DIR):
    FUNC = robjects.r['load']
    FUNC(DIR)

def DIPLOID_func(ALT_list,DP_list,ALPHA):
    alt_list = robjects.IntVector(ALT_list)
    dp_list = robjects.IntVector(DP_list)
    FUNC = robjects.r['DIPLOID_func']
    return FUNC(alt_list,dp_list,ALPHA)

def HAPLOID_func(ALT_list,DP_list,ALPHA):
    alt_list = robjects.IntVector(ALT_list)
    dp_list = robjects.IntVector(DP_list)
    FUNC = robjects.r['HAPLOID_func']
    return FUNC(alt_list,dp_list,ALPHA)

def extract_AD_info(GT,AD):
    gt = GT.split('/')
    DP = 0
    ad = AD.split(',')
    if (len(ad) < 2):
        print ad,AD
    REF_COUNT = ad[0]
    ad_2 = ad[1:]
    alt = []
    
    for i in ad_2:
        alt.append(int(i))
    
    for f in ad:
        DP = DP + int(f)
    
    ALT_COUNT = max(alt)
    
    if (len(ad) >= 2):
        return (str(DP), str(REF_COUNT), str(ALT_COUNT))

def ABB(variant):
        columns= variant.split('\t')

        CHROM = columns[0]
        POS = columns[1]
        REF = columns[3]
        ALT = columns[4]
        FILTER = columns[6]
        INFO = columns[7]
        form = columns[8].split(':')
        DP_list = []
        ALT_list = []
        
        if ("GT" in form):
            GTi=[i for i, x in enumerate(form) if x == "GT"]
        else:
            print "GT not present in vcf"
            #break
        
        if ("AD" in form):
            ADi=[i for i, x in enumerate(form) if x == "AD"]
        else:
            print "AD not present in vcf"
            return ("NO_AD")
                      
        for individual in range(9,len(columns)):#the first sample is in the 10th position of the list. It allows to move across all the samples
            #SAMPLE = name_columns[individual]
            GT = columns[individual].split(':')[GTi[0]]
            AD =columns[individual].split(':')[ADi[0]]
            
            
            if not (bool(re.search('\.',GT))): # checking for values "." in the genotype
                                
                DP, REF_COUNT, ALT_COUNT = extract_AD_info(GT,AD)
                #print GT, AD, REF_COUNT, ALT_COUNT, DP
                #print type(ALT_COUNT)
                if (int(DP) > 10): #and type(ALT_COUNT) == int and type(DP) == int):
                    DP_list.append(int(DP))
                    ALT_list.append(int(ALT_COUNT))
                
        if (len(DP_list) > 60):
            if not(bool(re.search('X',str(CHROM))) or bool(re.search('23',str(CHROM))) or bool(re.search('Y',str(CHROM))) or bool(re.search('24',str(CHROM)))):
                #print CHROM, POS, len(ALT_list), len(DP_list)
                #print DP_list
                abb = DIPLOID_func(ALT_list,DP_list,0.01)
                #print abb
                ABB = abb[3]

            elif ((bool(re.search('X',str(CHROM))) or bool(re.search('23',str(CHROM))))):
                abb  = DIPLOID_func(ALT_list,DP_list,0.01)
                ABB = abb[3]
                
            elif (bool(re.search('Y',str(CHROM))) or bool(re.search('24',str(CHROM)))):
                abb = HAPLOID_func(ALT_list,DP_list,0.01)
                ABB = abb[3]
        
            else:
                ABB = "."
        
            line = (str(CHROM), str(POS), str(ABB))
        
            LINE = '\t'.join(line)+'\n'
            return(LINE)



############################
## RUNNING PROGRAM
############################

## Open input and output files
VCF=open(infile,'r')
#VCF=open("/users/so/fmuyas/programs/ABB/test/RVAS_test.vcf",'r')


information= []
each_sample_position=[]

of1= (outfile)
OF1=open(of1,'w')
#OF1=open('/users/so/fmuyas/programs/ABB/python_scripts/parallel/result_200.vcf','w')


####### LOADING THE R SOURCE
LOADING_SOURCE(source)


# list of tasks
tasks=[]
count_line=0

####OPEN VCF FILE
tmp = open(infile,'r')
magic_number = tmp.read(2)
print 'MESSAGE:: BGZIP file  : %s'%(magic_number=='\x1f\x8b')
tmp.close()

with open(infile) if magic_number!='\x1f\x8b' else bgzf.open(infile) as VCF:
    for variant in VCF:
        variant = variant.rstrip('\n')
            
        if variant.startswith('#CHROM'): #to select the name of each column
            name_columns= variant.split('\t')
            name_columns= name_columns[0:len(name_columns)]#there is something at the end of the line and then we remove it
            ##INFO=<ID=ABB_SCORE,Number=1,Type=Float,Description="Allele balance bias (ABB) value">
            OF1.write('##INFO: ABB score list'+'\n')
            header = "CHROM\tPOS\tABB\n"
            OF1.write(header)
            
        if not variant.startswith('#'): #to select only the variants and not the other rows (## and #CHROM)
            if (count_line < 10000 and variant is not None):
                tasks.append(variant)
                count_line = count_line + 1
                
            else:
                tasks.append(variant)

                if __name__ == '__main__':
                    p = Pool(CORES)
                    A = (p.map(ABB, (tasks)))
                    #for i in A:
                    #    if i is not None:
                    #        OF1.write(i)
                    for i in A:
                        if (i == "NO_AD"):
                            print ("AD not present in the VCF")
                            break
                        elif (i is not None):
                            OF1.write(i)

                
                ## Re-start the counting. This step is mainly done to avoid a big usage of memory when vcf input file is super big.
                tasks=[]
                count_line=0
    
    if len(tasks) > 0:    
        if __name__ == '__main__':
            p = Pool(CORES)
            A = (p.map(ABB, (tasks)))
            #for i in A:
            #    if i is not None:
            #        OF1.write(i)
            for i in A:
                if (i == "NO_AD"):
                    print ("AD not present in the VCF")
                    break
                elif (i is not None):
                    OF1.write(i)


#VCF.close()
OF1.close()


stop = timeit.default_timer()
print stop - start
