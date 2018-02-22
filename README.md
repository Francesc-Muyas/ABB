# ABB tool
This tool is able to detect false positive calls. It is based on a new strategy to identify systematic sequencing or alignment errors leading to false positive variant calls based on the recurrence of the allele balance bias (paper in preparation). The main applications of this tool are four:
* 'ABB_list' obtains a new callability score list based on a new subset of samples. This list label positions of the genome with values between 0 and 1, which represent the precision of being a systematic error.
* 'ABB_annotation' Annotates and filters variant from input vcfs with an existing or new ABB list under a specific threshold of ABB.
* 'ABB_filter' obtains ABB from vcfs on-fly and filter variants based on a specific ABB threshold.
* 'ABB_association' detects genes/regions and variants which ABB could explain the significant association found in case-control studies.

This tools is though to be run in DNA data, taking as input VCF files and other parameters to remove variants prone to systematic errors (highly enriched with false positive calls). Statistical analysis and results are based on the paper : paper Muyas et al.



## Get ABB tool code 
You will need to run `git clone  ` to get ABB tool. 

```
git clone https://github.com/Francesc-Muyas/ABB
```

## Install dependencies
Firstly, you will need to install some software

- [Python 2.7](https://www.python.org/download/releases/2.7/), with next python modules:
    * numpy
    * argparse
    * timeit
    * os
    * warnings
    * pandas
    * psutil
    * biopython
    * multiprocessing
    * re
    * [rpy2](https://pypi.python.org/packages/3d/9b/b76b3665936204e14174dcac4814d8c91c833e9c3164664d5e89d777dac5/rpy2-2.7.0.tar.gz). If you have problems installing this module, it can be done like:

    ```
    pip install https://pypi.python.org/packages/3d/9b/b76b3665936204e14174dcac4814d8c91c833e9c3164664d5e89d777dac5/rpy2-2.7.0.tar.gz
    ```

- R version 3.3.2, with next packages and respective dependencies:   
    * data.table
    * ggplot2
    * splines
    * gamlss.data
    * gamlss.dist
    * gamlss
    * grid
    * gridExtra
    * argparse
    * tools

- [argparse.bash](https://github.com/nhoffman/argparse-bash)

```
cd PATH/TO/ABB/source
wget https://raw.githubusercontent.com/nhoffman/argparse-bash/master/argparse.bash
chmod +x argparse.bash
```

- Download ABB list already obtained for the human exome (Hg19):

```
cd PATH/TO/ABB/source
wget https://public_docs.crg.es/sossowski/publication_data/ABB/ABB_SCORE.txt
```

- [shc](https://github.com/neurobin/shc)
Install shell script compiler. It can be downloaded and installed or you can just download a compiled binary package like next:

```
cd PATH/TO/ABB/source
wget https://github.com/neurobin/shc/releases/download/3.9.6/shc-3.9.6-bin-amd64-i386-arm64-armhf-ppc64el.tar.gz
tar -xzf shc-3.9.6-bin-amd64-i386-arm64-armhf-ppc64el.tar.gz
```

## Install/prepare the tool
This tool is based on shell, python and R scripts. Once all dependencies are installed, the tool can easily be set up with next command:
```
cd PATH/TO/ABB
bash set_up.sh -t PATH/TO/ABB -r PATH/TO/Rscript -p PATH/TO/python -a PATH/TO/argparse.bash
```

On the other hand, if you want to use it as a binary executable, you only need shc, previously downloaded, and run next command:
```
cd PATH/TO/ABB
mkdir -p bin
PATH/TO/shc -f shell/ABB_tool.sh -o bin/ABB
```

## Usage
VCF files must by in VCF [v4.1](http://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41) format. Additionally, for ABB_list, ABB_filter and ABB_association functions, vcf files must have the [Allelic Depth (AD) and Genotype (GT)](https://gatkforums.broadinstitute.org/gatk/discussion/1268/what-is-a-vcf-and-how-should-i-interpret-it) for each sample, needed to compute the allele balance and get the ABB score on-fly. The tool, as explained above has 4 main applications.

You can see all arguments just running: `bin/ABB –h `(values between [] represent default parameters).

```
usage: ABB [-h] -T {ABB_filter,ABB_list,ABB_association,ABB_annotation} -i
           INPUT_FILE [-o OUT_FOLDER] [-outfile OUTFILE] [-cases CASE_IDS]
           [-controls CONTROL_IDS] [-genes GENES_FILE] [-abb_file ABB_FILE]
           [-abb_filter ABB_FILTER]

optional arguments:
  -h, --help            show this help message and exit
  -T {ABB_filter,ABB_list,ABB_association,ABB_annotation}, --tool {ABB_filter,ABB_list,ABB_association,ABB_annotation}
                        ABB tool that you want to run
  -i INPUT_FILE, --input_file INPUT_FILE
                        VCF file to be analyzed (AD is required in the column
                        format)
  -o OUT_FOLDER, --out_folder OUT_FOLDER
                        Out folder (Full path)
  -outfile OUTFILE, --outfile OUTFILE
                        Vcf or tabulated out file name for
                        ABB_filter/ABB_annotation and ABB_list tool,
                        respectively. You have to specify the out folder with
                        the --out_folder parameter
  -cases CASE_IDS, --case_ids CASE_IDS
                        File with cases IDs (For ABB_association tool)
  -controls CONTROL_IDS, --control_ids CONTROL_IDS
                        File with control IDs (For ABB_association tool)
  -genes GENES_FILE, --genes_file GENES_FILE
                        File with gene annotation for each site (For
                        ABB_association tool)
  -abb_file ABB_FILE, --abb_file ABB_FILE
                        ABB score file (For ABB_association tool)
  -abb_filter ABB_FILTER, --abb_filter ABB_FILTER
                        ABB score cutoff [0.7]
```

* ABB_list

It obtains a new callability score list based on a new subset of samples. This list labels positions of the genome with values between 0 and 1, which represent the precision of being a systematic error. The input vcf must be a multi-sample calling vcf.
The output of this script can be used as ABB_FILE for other applications like ABB_annotation or ABB_association. To have reliable ABB values (able to detect recurrent biased allele balance positions) we recommend at least 80 samples in the multi-sample vcf file.

```
bin/ABB -T ABB_list \
    -i multi_sample.vcf \
    -o out_folder \
    -outfile ABB.list 
```

* ABB_filter

It obtains ABB from VCF on-fly and filter variants based on a specific ABB threshold. The input VCF must be a multi-sample calling vcf. It computes ABB_score on-fly and flags each variant site as ABB_biased if the obtained ABB score is greater than the threshold specified in –abb_filter. By default, this value is 0.7, representing low confidence sites. If you want to be stricter you can use any –abb_filter > 0.7 and < 1. For example, you can use -abb_filter 0.9 as hard filter, as it is used to flag very low confidence sites in the paper. You will find the flag ABB_biased in the FILTER column of the output vcf file. To have reliable ABB values (able to detect recurrent biased allele balance positions) we recommend at least 80 samples in the multi-sample vcf file.

```
bin/ABB -T ABB_filter \
    -i multi_sample.vcf \
    -o out_folder \
    -outfile filtered.vcf \
    -abb_filter ABB_threshold
```

* ABB_annotation

It annotates and filters variant from input VCFs with an existing or new ABB list under a specific threshold of ABB. Be sure that the reference genome version used for creating the ABB list is the same as the one used for the variant calling (VCF). 
This script flags each variant site as ABB_biased if the obtained ABB score from the list is greater than the threshold specified in –abb_filter. By default, this value is 0.7, representing low confidence sites. If you want to be stricter you can use whatever –abb_filter greater than 0.7. You will find the flag ABB_biased in the FILTER column of the output VCF file. In this application, the input VCF file is not needed to be multi-sample calling VCF file, it also works with a single-sample VCF file.
By default, this script uses the ABB list previously downloaded in source/ABB_SCORE.txt, however, you can use the output of the ABB_list script.

```
bin/ABB -T ABB_annotation \
    -i input.vcf \
    -o out_folder \
    -outfile filtered.vcf \
    -abb_file source/ABB_SCORE.txt
```

* ABB_association

It detects genes/regions and variants which ABB could explain the significant association found in case-control studies. Using this command line, you can run this application:

```
bin/ABB -T ABB_association \
    -i multi_sample.vcf \
    -o out_folder \
    -cases Cases.txt \
    -controls Controls.txt \
    -genes Genes.txt \
    -abb_file source/ABB_SCORE.txt
```

The input VCF file must be a multi-sample calling VCF file with the variants that previously you found associated in a Rare Variant Assocition Study (RVAS) test. Additionally, two files specifying the cases (`-cases `) or controls (`-controls `) ids are needed. These files must be one column with the same ids as the ones found in the vcf:

Cases.txt:
```
Case1_id
Case2_id
Case3_id
…
```

Controls.txt:
```
Control1_id
Control2_id
Control3_id
...
```

On the other hand, you need an easy annotation file for each variant where you specify the gene/region which the variant belongs (`-genes `). It is not necessary to be gene annotation, you can use whatever annotation you used for the association test. This tabulated file should look like this (Chrom, Positions, Gene):

Genes.txt:
```
1	114524278	OLFML3
1	185113070	TRMT1L
…
```

As previous applications, by default, the ABB list is source/ABB_SCORE.txt, but any (with same reference genome) ABB list obtained from ABB_list application can be used. 

ABB_ASSOCIATION output:

The output of this command return 4 different txt files (2 for SNVs and the other 2 for Genes) with p-values (multiple test corrected) showing the significance if that association is product or not of allele balance bias (described in main paper). 2 files represent significant SNVs/Genes and the other all SNVs/Genes with all p-values (FDR). Additionally, for each significant SNP and gene, you will find plots showing the proportion of variants not called by the variant caller, but showing significant amount of non-reference counts (MISSED) and the ones called by the variant caller (CALLED). 


SNP (ABB_ASSOCIATION.SNPs.txt and ABB_ASSOCIATION.SNPs.significant.txt) columns:

    * CHROM: Chromosome of the variant.
    * POS: Variant site coordinate.
    * ABB: ABB score for the variant site
    * Missed-Called_ratio(FDR): FDR obtained from the missed-called ratio comparison between cases and controls (Fisher test). If significant, it means that this ratio is different beween cases and controls.
    * GENE: Annotation (for ejample gene) of the variant.


GENE (ABB_ASSOCIATION.GENEs.significant.txt  and ABB_ASSOCIATION.GENEs.significant.txt) columns:

    * GENE: Tested gene.
    * Missed-Called_ratio(FDR): FDR obtained from the missed-called ratio comparison between cases and controls (Fisher test). If significant, it means that this ratio is different beween cases and controls.
    * Association_regenotyped(FDR): FDR value obtained from association chi square test between cases and control including re-genotyped variants (including the MISSED calls).
    * Association_ABB(FDR): FDR value obtained from association chi square test between cases and control including but removing significantly biased sites (significant sites in ABB_ASSOCIATION.SNPs.significant.txt).
