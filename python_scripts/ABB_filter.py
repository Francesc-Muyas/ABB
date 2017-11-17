#!/usr/bin/env python

import re
import timeit
import argparse
import rpy2.robjects as robjects
import os
import warnings
from Bio import bgzf

start = timeit.default_timer()


parser = argparse.ArgumentParser()
parser.add_argument('-infile', type=str, help='Input vcf file', required= True)
parser.add_argument('-outfile', type=str, help='Out vcf file', required= True)
parser.add_argument('-abb', type=float, default= 0.7, help='Cutoff for ABB filter', required = False)
#parser.add_argument('-gender', type=str, default= "NA", help='Geneder of the samples: M for Males ; F for Females', required = False)

args = parser.parse_args()

infile=args.infile
outfile=args.outfile
#gender=args.gender
cutoff=args.abb

#VCF=open(infile,'r')
#VCF=open("/users/so/fmuyas/ediva_toy/project/pileup_part/test/200_additional_CLL_samples.enriched.vcf",'r')


information= []
each_sample_position=[]

of1= (outfile)
OF1=open(of1,'w')


# R functions
DIR = os.path.dirname(__file__)
filename = os.path.join(DIR,"../rscripts/ABB_function.r")
source = os.path.join(DIR,"../source/source.Rdata")

filename = os.path.realpath(filename)
source = os.path.realpath(source)

print DIR
print source
print filename

r=robjects.r
r.source(filename)
#r.source('/users/so/fmuyas/programs/ABB/rscripts/ABB_function.r')

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

def ABB_FILTER(ABB,cutoff,FILTER):
    if (ABB != '.' and float(ABB) >= cutoff):
        return("ABB_biased")
    if (ABB == '.' or float(ABB) < cutoff):
        return(FILTER)
        
def extract_AD_info(GT,AD):
    gt = GT.split('/')
    DP = 0
    ad = AD.split(',')
    REF_COUNT = ad[0]
    ad_2 = ad[1:]
    alt = []
    
    for i in ad_2:
        alt.append(int(i))
        
    
    for f in ad:
        DP = DP + int(f)
    
    ALT_COUNT = max(alt)
    
    return (str(DP), str(REF_COUNT), str(ALT_COUNT))

        
    #print DP, GT, AD



####### LOADING THE R SOURCE
LOADING_SOURCE(source)

tmp = open(infile,'r')
magic_number = tmp.read(2)
print 'MESSAGE:: BGZIP file  : %s'%(magic_number=='\x1f\x8b')
tmp.close()

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
                break
            
            if ("AD" in form):
                ADi=[i for i, x in enumerate(form) if x == "AD"]
            else:
                print "GT not present in vcf"
                break
            
            if ("GQ" in form):
                GQi=[i for i, x in enumerate(form) if x == "GQ"]
            else:
                print "GQ not present in vcf"
                break
                  
            for individual in range(9,len(columns)):#the first sample is in the 10th position of the list. It allows to move across all the samples
                SAMPLE = name_columns[individual]
                GT = columns[individual].split(':')[GTi[0]]
                AD =columns[individual].split(':')[ADi[0]]
                
                
                if not (bool(re.search('\.',GT))): # checking for values "." in the genotype
                                    
                    DP, REF_COUNT, ALT_COUNT = extract_AD_info(GT,AD)
                    if (int(DP) > 10): #and type(ALT_COUNT) == int and type(DP) == int):
                        DP_list.append(int(DP))
                        ALT_list.append(int(ALT_COUNT))
                    
                    #DP_list.append(DP)
                    #ALT_list.append(ALT_COUNT)
                    #print outline
            
    
                #OF1.write(outline)
            if (len(DP_list) > 60):
                if not(bool(re.search('X',str(CHROM))) or bool(re.search('23',str(CHROM))) or bool(re.search('Y',str(CHROM))) or bool(re.search('24',str(CHROM)))):
                    print CHROM, POS, len(DP_list), len(ALT_list), form, len(columns), DP_list[-1],ALT_list[-1]
                    abb  = DIPLOID_func(ALT_list,DP_list,0.01)
                    ABB = abb[3]
                
                elif ((bool(re.search('X',str(CHROM))) or bool(re.search('23',str(CHROM))))):
                    abb  = DIPLOID_func(ALT_list,DP_list,0.05)
                    ABB = abb[3]

                elif (bool(re.search('Y',str(CHROM))) or bool(re.search('24',str(CHROM)))):
                    abb = HAPLOID_func(ALT_list,DP_list,0.01)
                    ABB = abb[3]
                
                else:
                    ABB = '.'
                    
            else:
                ABB = "."
            
                
            FILT = ABB_FILTER(ABB,cutoff,FILTER)
                  
            columns[6] = FILT
            columns[7] = INFO + ";ABB=" + str(ABB)        
            
            OF1.write('\t'.join(columns)+'\n')
            
            

#VCF.close()
OF1.close()


stop = timeit.default_timer()
print stop - start
