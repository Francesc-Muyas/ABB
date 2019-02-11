#!/bin/bash

#################################
#   TOOL FOR RUNNING ABB FILTER
#################################
#export LD_LIBRARY_PATH=/software/so/el7.2/Python-2.7.13/bin/lib:$LD_LIBRARY_PATH

export PATH=/software/so/el7.2/Python-2.7.13/bin/bin/:$PATH

echo '#################################'
echo '#   TOOL FOR RUNNING ABB ASSOCIATION TEST'
echo '#################################'
echo



####https://github.com/nhoffman/argparse-bash
#wget https://raw.githubusercontent.com/nhoffman/argparse-bash/master/argparse.bash
#chmod +x argparse.bash
#ARGPARSE="`dirname $0`"/../source/argparse.bash
ARGPARSE=PATH/TO/argparse.bash
source $ARGPARSE || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('-i','--input_vcf', required=True, help='VCF.gz file to be analyzed (AD is required in the column format)')
parser.add_argument('-o', '--out_vcf', default='$PWD', help='Out vcf file')
parser.add_argument('-abb', '--abb_filter', default=0.7, help='ABB score filter')
parser.add_argument('-t', '--threads', default=1, help='Threads to run the script')
parser.add_argument('-path', '--abb_path', required=True, help='Path to ABB tool')
EOF

echo Parameters:
echo VCF input file: "$INPUT_VCF"
echo VCF output file: "$OUT_VCF"
echo ABB score filter: "$ABB_FILTER"
echo Threads: "$THREADS"
echo ABB tool path: "$ABB_PATH"
echo

# Absolute path to ABB tool
TOOL=$(readlink -f $ABB_PATH)
# Path to python scripts
PYTHON_SCRIPTS_t="$TOOL"/python_scripts
PYTHON_SCRIPTS=$(readlink -f $PYTHON_SCRIPTS_t)
# Path to R scripts
R_SCRIPTS="$TOOL"/rscripts

#echo $TOOL
#echo $TOOLPATH
#echo $PYTHON_SCRIPTS_t
#echo $PYTHON_SCRIPTS
#echo $R_SCRIPTS
#echo

# Parameters
if [ "$THREADS" -gt 1 ];then
    out_vcf=$(dirname $OUT_VCF)
    rm -f $OUT_VCF
    rm -rf $out_vcf/temp
    rm -rf $out_vcf/out_temp
    mkdir -p $out_vcf/
    mkdir -p $out_vcf/temp
    mkdir -p $out_vcf/out_temp
       
    # Grab the header
    zcat $INPUT_VCF | head -n 10000 | grep "^#" > header
    echo done the header , split the variants now
    #grab the non header lines
    zgrep -v "^#" $INPUT_VCF > variants
       
    echo grepped the variants now splitting
    # Split the file in x lines per file to produce the amount of files desired
    N=$(wc -l variants | cut -f1 -d ' ')
    Lines_file=$(expr $N / $THREADS)
        
    split -d -l $Lines_file -a 3 variants $out_vcf/temp/vcf_
    # Reattach the header to each and clean up
    for i in $out_vcf/temp/vcf_*;do cat header $i > $i.temp && rm -f $i; echo $i ;done

    for i in $out_vcf/temp/vcf_*temp;do
        I=$(basename $i)
        echo $I
        echo python $PYTHON_SCRIPTS/ABB_filter_Parallel.py -infile $i \
            -outfile $out_vcf/out_temp/$I \
            -abb $ABB_FILTER
        python $PYTHON_SCRIPTS/ABB_filter_Parallel.py -infile $i -outfile $out_vcf/out_temp/$I -abb $ABB_FILTER &
    done
    wait
    
    # Creating the final file
    cat header > $OUT_VCF        
    for i in $out_vcf/out_temp/*temp;do
        grep -v '^#' $i >> $OUT_VCF
    done
       
    rm -rf header variants $out_vcf/temp $out_vcf/out_temp
    
    else
        # Only 1 thread used
        python $PYTHON_SCRIPTS/ABB_filter_Parallel.py -infile $INPUT_VCF -outfile $OUT_VCF -abb $ABB_FILTER

fi

end=`date +%s`
runtime=$((end-start))

rm -f argparse.bash
echo
echo "The time is: $runtime"
