#!/usr/bin/env python

# python extract_info_VCF.py -infile /users/so/fmuyas/PhD/Part_1/data/200_extra_independent.txt -out /users/so/fmuyas/PhD/Part_1/data/200_extra_independent_INFO.txt
import re
import timeit
import argparse
from Bio import bgzf

start = timeit.default_timer()

parser = argparse.ArgumentParser()
parser.add_argument('-infile', type=str, help='Input vcf file', required= True)
parser.add_argument('-outfile', type=str, help='Out vcf file', required= True)

args = parser.parse_args()

infile=args.infile
outfile=args.outfile


tmp = open(infile,'r')
magic_number = tmp.read(2)
print 'MESSAGE:: BGZIP file  : %s'%(magic_number=='\x1f\x8b')
tmp.close()


information= []
each_sample_position=[]

of1= (outfile)
OF1=open(of1,'w')

column_names=("CHROM", "POS", "REF", "ALT", "SAMPLE", "GT", "AD", "DP", "REF_COUNT", "ALT_COUNT")
column_names="\t".join(column_names)
column_names=column_names+"\n"
OF1.write(column_names)


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
####OPEN VCF FILE
with open(infile) if magic_number!='\x1f\x8b' else bgzf.open(infile) as VCF:
    for variant in VCF: #to extract the information of each variant/row
        variant = variant.rstrip('\n')
        if variant.startswith('#CHROM'): #to select the name of each column
            name_columns= variant.split('\t')
            name_columns= name_columns[0:len(name_columns)]#there is something at the end of the line and then we remove it
            #print len(name_columns)
            #print name_columns[len(name_columns)-1]
          
        if not variant.startswith('#'): #to select only the variants and not the other rows (## and #CHROM)
            columns= variant.split('\t')
            #samples= name_columns[9:len(name_columns)-1]
    #        print columns
            CHROM = columns[0]
            POS = columns[1]
            REF = columns[3]
            ALT = columns[4]
            FILTER = columns[6]
            form = columns[8].split(':')
            if ("GT" in form):
                GTi=[i for i, x in enumerate(form) if x == "GT"]
            else:
                print "GT not present in vcf"
                GTi="NO"
            
            if ("AD" in form):
                ADi=[i for i, x in enumerate(form) if x == "AD"]
            else:
                print "GT not present in vcf"
                ADi="NO"
            
            for individual in range(9,len(columns)):#the first sample is in the 10th position of the list. It allows to move across all the samples
                SAMPLE = name_columns[individual]           
                
                if (GTi != "NO"): 
                    GT = columns[individual].split(':')[GTi[0]]
                else:
                    GT = "0/0"
    
                if (ADi != "NO"): 
                    AD = columns[individual].split(':')[ADi[0]]
                else:
                    AD = "."
                                
                if (AD != "."): 
                    DP, REF_COUNT, ALT_COUNT = extract_AD_info(GT,AD)
                else:
                    DP = "."
                    REF_COUNT = "."
                    ALT_COUNT = "."
                
                
                OUTLINE = [CHROM, POS, REF, ALT,  SAMPLE, GT, AD, DP, REF_COUNT, ALT_COUNT]
                outline = '\t'.join(OUTLINE)+'\n'
                #print outline
            
    
                OF1.write(outline)
    #    

#VCF.close()
OF1.close()


stop = timeit.default_timer()

print "Time: " + str(stop - start)

