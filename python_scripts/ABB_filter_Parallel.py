#!/usr/bin/env python

from multiprocessing import Pool, cpu_count
import re
import timeit
import argparse
import psutil
import os
import warnings
from Bio import bgzf
import numpy
from scipy.stats import beta
#from scipy.stats import binom_test
from scipy.stats import binom
import scipy

start = timeit.default_timer()

parser = argparse.ArgumentParser()
parser.add_argument('-infile', type=str, help='Input vcf file', required= True)
parser.add_argument('-outfile', type=str, help='Out vcf file', required= True)
parser.add_argument('-abb_filter', type=float, default= 0.7, help='Cutoff for ABB filter', required = False)

args = parser.parse_args()

infile=args.infile
outfile=args.outfile

cutoff=args.abb_filter

CORES=psutil.cpu_count(logical=False)


information= []
each_sample_position=[]

of1= (outfile)
OF1=open(of1,'w')

## Source
DIR = os.path.dirname(__file__)
ppv_file = os.path.join(DIR,"../source/PPV.txt")

ppv_file = os.path.realpath(ppv_file)

def pBEINF(x, mu, sigma, nu, tau):
    a = mu * (1 - sigma**2)/(sigma**2)
    b = a * (1 - mu)/mu
    
    if (x > 0 and x < 1):
        cdf = nu + beta.cdf(x,a,b)
    elif (x == 0):
        cdf = nu
    elif (x == 1):
        cdf = 1 + nu + tau
        
    cdf = cdf/(1 + nu + tau)
    return(cdf)

def binom_test(A,N,P):
    AF = float(A)/N
    
    if (AF == P):
        pval = 1
    elif (AF < P):
        pval = binom.cdf(A,N,P)*2
    elif (AF > P):
        pval = (1-binom.cdf(A-1,N,P))*2
    
    return(pval)
    
def lm_model(SCORE1, SCORE2, SCORE3, OBS_l, PPV_l):
    # Predicting value basde on scores
    VAL = -4.6631-9.9920*SCORE1+28.6594*SCORE2+1.0727*SCORE3-4.9271*SCORE2*SCORE3
    
    # Min and max for PPV and OBS
    min_OBS = min(OBS_l)
    max_OBS = max(OBS_l)
        
    # Found the PPV for this val
    if (VAL <= min_OBS):
        PPV = min(PPV_l)
    elif (VAL >= max_OBS):
        PPV = max(PPV_l)
    else:
        # Values are sorted, so when we find a value >=  VAL, it is the corresponding index for the right PPV
        for PPVi in range(0,len(OBS_l)):
            if OBS_l[PPVi] >= VAL:
                break
        
        PPV = PPV_l[PPVi]
    
    return([VAL, PPV])

def DIPLOID_func(ALT,DP):
    ## Zero-one inflated parameters
    ###### PARAMETERS 0/0 GT
    MU_0 = 0.03270537
    SIGMA_0 = 0.1452151
    NU_0 = 1.689303e+01
    TAU_0 = 1.011730e-14
    
    ###### PARAMETERS 1/1 GT  
    MU_1 = 0.97263928
    SIGMA_1 =  0.1364325
    NU_1 = 2.719154e-15
    TAU_1 = 4.626247e+00
    
    # ALPHA
    ALPHA = 0.05
    
    # Allele balance
    AB = float(ALT)/DP
    
    AB2 = float(ALT-1)/DP
    
    ## Hom ref
    if (AB == 0):
        P1 = 1
    elif (AB != 0 and AB != 1):
        P1 = 1 - pBEINF(AB2, MU_0, SIGMA_0, NU_0, TAU_0)
    else:
        P1 = -1
    
    ## Het
    if (AB != 0 and AB != 1):
        #P2 = scipy.stats.binom_test(ALT, n=DP, p=0.5, alternative='two-sided')
        P2 = binom_test(ALT, DP, 0.5) # This test is only exact when P = 0.5 (This case). If not it is better to used binom_test from scipy.stats (it is too slow)
    else:
        P2 = -1

    ## Hom alt
    if (AB == 1):
        P3 = 1
    elif (AB != 1 and AB != 0):
        P3 = pBEINF(AB, MU_1, SIGMA_1, NU_1, TAU_1)
    else:
        P3 = -1

    ## Choose higher p as correct genotype
    gt = ['0/0','0/1','1/1']
    dist = [0, 0.5, 1]
    Ps = [P1, P2, P3]
    
    Pi = [x for x in range(0,len(Ps)) if Ps[x] == max(Ps)][0]
        
    ## Getting values
    GT = gt[Pi]
    P = Ps[Pi]
    DIST = abs(dist[Pi] - AB)
    
    Log = -numpy.log10(P)
    
    ## Checking if devianceis significant 
    if (Pi == 1):
        ALPHA = ALPHA/2
    
    if (P < ALPHA):
        SIG = 1
    else:
        SIG = 0
    
    return [GT, P, DIST, Log, SIG]

def HAPLOID_func(ALT,DP):
    ## Zero-one inflated parameters
    ###### PARAMETERS 0/0 GT
    MU_0 = 0.03270537
    SIGMA_0 = 0.1452151
    NU_0 = 1.689303e+01
    TAU_0 = 1.011730e-14
    
    ###### PARAMETERS 1/1 GT  
    MU_1 = 0.97263928
    SIGMA_1 =  0.1364325
    NU_1 = 2.719154e-15
    TAU_1 = 4.626247e+00
    
    # ALPHA
    ALPHA = 0.05
    
    # Allele balance
    AB = float(ALT)/DP
    
    AB2 = float(ALT-1)/DP
    
    ## Hom ref
    if (AB == 0):
        P1 = 1
    elif (AB != 0 and AB != 1):
        P1 = 1 - pBEINF(AB2, MU_0, SIGMA_0, NU_0, TAU_0)
    else:
        P1 = -1
    
    ## Het (No possible in haplody)
    P2 = -1

    ## Hom alt
    if (AB == 1):
        P3 = 1
    elif (AB != 1 and AB != 0):
        P3 = pBEINF(AB, MU_1, SIGMA_1, NU_1, TAU_1)
    else:
        P3 = -1

    ## Choose higher p as correct genotype
    gt = ['0/0','0/1','1/1']
    dist = [0, 0.5, 1]
    Ps = [P1, P2, P3]
    
    Pi = [x for x in range(0,len(Ps)) if Ps[x] == max(Ps)][0]
    
    ## Getting values
    GT = gt[Pi]
    P = Ps[Pi]
    DIST = abs(dist[Pi] - AB)
    
    Log = -numpy.log10(P)
    
    ## Checking if devianceis significant 
    if (Pi == 1):
        ALPHA = ALPHA/2
    
    if (P < ALPHA):
        SIG = 1
    else:
        SIG = 0
    
    return [GT, P, DIST, Log, SIG]

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


def ABB_FILTER(ABB,cutoff,FILTER):
    if (ABB != '.' and float(ABB) >= cutoff):
        return("ABB_biased")
    if (ABB == '.' or float(ABB) < cutoff):
        return(FILTER)
        

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
    
    if ("AD" in form):
        ADi=[i for i, x in enumerate(form) if x == "AD"]
    else:
        print "AD not present in vcf"
        return ('NO_AD')
    
    
    SAMPLE_SIZE = 1000
    COLUMN_i = range(9,len(columns))
    if len(COLUMN_i) > int(round(SAMPLE_SIZE*1.25)):
        numpy.random.seed(1234)
        COLUMN_i = numpy.random.choice(COLUMN_i, int(round(SAMPLE_SIZE*1.25)), replace=False)
    
    for individual in COLUMN_i:#the first sample is in the 10th position of the list. It allows to move across all the samples
        SAMPLE = name_columns[individual]
        GT = columns[individual].split(':')[GTi[0]]
        AD =columns[individual].split(':')[ADi[0]]
        
        
        if not (bool(re.search('\.',GT))): # checking for values "." in the genotype
                            
            DP, REF_COUNT, ALT_COUNT = extract_AD_info(GT,AD)
            if (int(DP) > 10): #and type(ALT_COUNT) == int and type(DP) == int):
                DP_list.append(int(DP))
                ALT_list.append(int(ALT_COUNT))
    
    # We need at least 60 informative samples
    if (len(DP_list) >= 60):

        if (len(DP_list) > SAMPLE_SIZE):
            DP_list = numpy.array(DP_list)
            ALT_list = numpy.array(ALT_list)
            
            numpy.random.seed(1234)
            Sample_indexes = numpy.random.choice(len(DP_list), SAMPLE_SIZE, replace = False)
            DP_list = DP_list[Sample_indexes]
            ALT_list = ALT_list[Sample_indexes]
            
        score1 = []
        score2 = []
        score3 = []
        
        # Deciding the plody
        if not (bool(re.search('Y',str(CHROM))) or bool(re.search('24',str(CHROM)))):
            FUNC = DIPLOID_func
        else:
            FUNC = HAPLOID_func
        
        for IN in range(0,len(DP_list)):
            alt = ALT_list[IN]
            dp = DP_list[IN]

            GT, P, DIST, Log, SIG =  FUNC(alt, dp)
                
            score1.append(DIST)
            score2.append(SIG)
            score3.append(Log)
    
        ## Getting scores                
        SCORE1 = numpy.mean(score1)
        SCORE2 = sum(score2)/float(len(score2))
        SCORE3 = numpy.mean(score3)
        
        ## Getting the ABB value
        VAL, ABB = lm_model(SCORE1, SCORE2, SCORE3, OBS_l, PPV_l)   
    
    else:
        ABB = "."
    
    FILT = ABB_FILTER(ABB,cutoff,FILTER)
                  
    columns[6] = FILT
    columns[7] = INFO + ";ABB=" + str(ABB)        
            
    LINE  = ('\t'.join(columns)+'\n')
    
    return(LINE)


tmp = open(infile,'r')
magic_number = tmp.read(2)
print 'MESSAGE:: BGZIP file  : %s'%(magic_number=='\x1f\x8b')
tmp.close()


OBS_l = []
Resp_l = []
PPV_l = []
with open(ppv_file) as PPV:
    for line in PPV:
        if not line.startswith('OBS'):
            line = line.rstrip('\n')
            obs, resp, ppv = line.split('\t')
            
            ## Append values to list
            OBS_l.append(float(obs))
            Resp_l.append(float(resp))
            PPV_l.append(float(ppv))


# list of tasks
tasks=[]
count_line=0      

with open(infile) if magic_number!='\x1f\x8b' else bgzf.open(infile) as VCF:
    for variant in VCF: #to extract the information of each variant/row
        variant = variant.rstrip('\n')
        if variant.startswith('##'): #to select the name of each column
            OF1.write(variant+'\n')
        
        if variant.startswith('#CHROM'): #to select the name of each column
            name_columns= variant.split('\t')
            name_columns= name_columns[0:len(name_columns)]#there is something at the end of the line and then we remove it
            
            ## Printing the header of the file
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
                    A = (p.map(ABB, (tasks)))

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

            for i in A:
                if (i == "NO_AD"):
                    print ("AD not present in the VCF")
                    break
                elif (i is not None):
                    OF1.write(i)


OF1.close()


stop = timeit.default_timer()
print 'Time:'
print stop - start
            