#!/bin/bash

#################################
#   TOOL FOR RUNNING ABB SCORE ANNOTATION
#################################

# Python path
PYTHON=PATH/TO/python
# R path
RSCRIPT=PATH/TO/Rscript
# Absolute path to this script
SCRIPT=$(readlink -f "$0")
# Absolute path 
SCRIPTPATH=$(dirname "$SCRIPT")
# Path to python scripts
PYTHON_SCRIPTS="`dirname $0`"/../python_scripts
# Path to R scripts
R_SCRIPTS="`dirname $0`"/../rscripts

####https://github.com/nhoffman/argparse-bash
##wget https://raw.githubusercontent.com/nhoffman/argparse-bash/master/argparse.bash
##chmod +x argparse.bash
ARGPARSE=PATH/TO/argparse.bash
source $ARGPARSE || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('-i','--input_file', required=True, help='VCF file to be analyzed (AD is required in the column format)')
parser.add_argument('-outfile', '--outfile', required=True, help='Out vcf file')
parser.add_argument('-abb_file', '--abb_file', default="NO", type=str , help='ABB score file')
parser.add_argument('-abb_filter', '--abb_filter', default=0.7, help='ABB score filter')
EOF

# Parameters
if [ "$ABB_FILE" == "NO" ];then
    ABB_FILE="`dirname $0`"/../source/ABB_SCORE.txt
fi

echo Parameters:
echo VCF input file: "$INPUT_FILE"
echo VCF output file: "$OUTFILE"
echo ABB score file: "$ABB_FILE"
echo ABB score filter: "$ABB_FILTER"
echo

# VCF to temp file
OUT_FOLDER=$(dirname $OUTFILE)
OUT_folder=$(readlink -f $OUT_FOLDER)
TEMP_f=$OUT_folder/temp
mkdir -p $TEMP_f

# Checking if the file is compressed or not and grep the variant positions
gzip -t $INPUT_FILE 2>/dev/null
if [[ $? -eq 0 ]];then
    echo "Compressed file"
    zgrep -v '^#' $INPUT_FILE | awk -F'\t' -v OFS='\t' '{print $1, $2}' > $TEMP_f/POSITIONS.temp
    else
        echo "Not compressed"
        grep -v '^#' $INPUT_FILE | awk -F'\t' -v OFS='\t' '{print $1, $2}' > $TEMP_f/POSITIONS.temp
fi

# Getting ABB for each one of the desired positions
$RSCRIPT $R_SCRIPTS/ABB_annotate.r --positions_file $TEMP_f/POSITIONS.temp --abb_file $ABB_FILE --out_file $TEMP_f/ABB.temp

# Annotating and adding filter column to vcf file with ABB scores > cutoff
$PYTHON $PYTHON_SCRIPTS/ABB_annotate_Parallel.py -i $INPUT_FILE -outfile $OUTFILE --abb_file $TEMP_f/ABB.temp --abb_filter $ABB_FILTER

rm -rf $TEMP_f
