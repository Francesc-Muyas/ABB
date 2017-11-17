#!/bin/bash

#################################
#   TOOL FOR RUNNING ABB ASSOCIATION TEST
#################################

# Python path
PYTHON=/software/so/el7.2/Python-2.7.13/bin/bin/python
# R path
RSCRIPT=/usr/bin/Rscript
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
ARGPARSE="`dirname $0`"/../source/argparse.bash
source $ARGPARSE || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('-i','--input_file', required=True, help='VCF file to be analyzed (AD is required in the column format)')
parser.add_argument('-o', '--out_folder', default='$PWD', help='Out folder (Full path)')
parser.add_argument('-cases', '--case_ids', required=True, help='File with cases IDs')
parser.add_argument('-controls', '--control_ids', required=True, help='File with control IDs')
parser.add_argument('-genes', '--genes_file', required = True, help='File with gene annotation for each site')
parser.add_argument('-abb_file', '--abb_file', default="NO",type=str, help='ABB score file')
EOF

# Parameters
if [ "$ABB_FILE" == "NO" ];then
    ABB_FILE="`dirname $0`"/../source/ABB_SCORE.txt
fi

echo Parameters:
echo VCF input file: "$INPUT_FILE"
echo Cases file: "$CASE_IDS"
echo Controls file: "$CONTROL_IDS"
echo Gene annotation file: "$GENES_FILE"
echo ABB score file: "$ABB_FILE"
echo Out folder: "$OUT_FOLDER"
echo


# VCF to TSV file
OUT_folder=$(readlink -f $OUT_FOLDER)
TEMP_f=$OUT_folder/temp
mkdir -p $TEMP_f

$PYTHON $PYTHON_SCRIPTS/extract_info_VCF.py -i $INPUT_FILE -o $TEMP_f/EXPERIMENT.tsv

# TSV to ABB ASSOCIATION TEST
$RSCRIPT $R_SCRIPTS/RVAS_ABB_ASSOCIATION.r --tsv_file $TEMP_f/EXPERIMENT.tsv \
    --cases $CASE_IDS \
    --controls $CONTROL_IDS \
    --genes $GENES_FILE \
    --abb $ABB_FILE \
    --out_prefix $OUT_folder/ABB_ASSOCIATION

rm -rf $TEMP_f