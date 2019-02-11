#!/bin/bash

#################################
#   TOOL FOR RUNNING ABB ASSOCIATION TEST
#################################

echo '#################################'
echo '# ABB tool'
echo '#################################'
echo

PYTHON=PATH/TO/python
#'''
## Checking python 
#if ls $PYTHON 1> /dev/null 2>&1; then
#    echo "$PYTHON exists"
#    exit 1
#    echo
#else
#    echo "$PYTHON does not exit. Changing python to python default"
#    PYTHON=$(which python)    
#    exit 1
#fi
#'''

# Absolute path to ABB tool
ABB_PATH=PATH/TO/ABB
SCRIPT=$(readlink -f $ABB_PATH)
## Checking python 
if ls $SCRIPT 1> /dev/null 2>&1; then
    echo "$SCRIPT directory exists"
    echo
else
    echo "Error: $SCRIPT - Incorrect path to ABB"   
    exit 1
fi

####https://github.com/nhoffman/argparse-bash
##wget https://raw.githubusercontent.com/nhoffman/argparse-bash/master/argparse.bash
##chmod +x argparse.bash
ARGPARSE=PATH/TO/argparse.bash
if ls $ARGPARSE 1> /dev/null 2>&1; then
    echo "Argparser exists"
    source $ARGPARSE || exit 1
    echo
else
    echo "Argparser does not. Downloading argparser"
    wget https://raw.githubusercontent.com/nhoffman/argparse-bash/master/argparse.bash
    chmod +x argparse.bash
    ARGPARSER_temp=argparse.bash
    source $ARGPARSER_temp || exit 1
fi

argparse "$@" <<EOF || exit 1
parser.add_argument('-T','--tool', required=True, choices=['ABB_filter', 'ABB_list', 'ABB_association', 'ABB_annotation'], type=str, help='ABB tool that you want to run')
parser.add_argument('-i','--input_file', required=True, help='VCF file to be analyzed (AD is required in the column format)')
parser.add_argument('-o', '--out_folder', default='$PWD', help='Out folder (Full path)')
parser.add_argument('-outfile', '--outfile', help='Vcf or tabulated out file name for ABB_filter/ABB_annotation and ABB_list tool, respectively. You have to specify the out folder with the --out_folder parameter')
parser.add_argument('-cases', '--case_ids', required=False, help='File with cases IDs (For ABB_association tool)')
parser.add_argument('-controls', '--control_ids', required=False, help='File with control IDs (For ABB_association tool)')
parser.add_argument('-genes', '--genes_file', required=False, help='File with gene annotation for each site (For ABB_association tool)')
parser.add_argument('-abb_file', '--abb_file', default="NO", type=str, help='ABB score file (For ABB_association tool)')
parser.add_argument('-abb_filter', '--abb_filter', default=0.7, help='ABB score cutoff [0.7]')
EOF

# Parameters
if [ "$ABB_FILE" == "NO" ];then
    ABB_FILE=$ABB_PATH/source/ABB_SCORE.txt
fi

echo Parameters:
echo ABB tool to be run: "$TOOL"
echo VCF input file: "$INPUT_FILE"
echo Out folder: "$OUT_FOLDER"
echo Out file: "$OUTFILE"
echo Cases file: "$CASE_IDS"
echo Controls file: "$CONTROL_IDS"
echo Gene annotation file: "$GENES_FILE"
echo ABB score file: "$ABB_FILE"
echo ABB score filter: "$ABB_FILTER"
echo ABB tool path: "$ABB_PATH"
echo Python used: $(which $PYTHON)
echo

#### Getting paths
# Path to python scripts
PYTHON_SCRIPTS_t="$SCRIPT"/python_scripts
PYTHON_SCRIPTS=$(readlink -f $PYTHON_SCRIPTS_t)

# Path to R scripts
R_SCRIPTS="$SCRIPT"/rscripts

# Path to shell scripts
SHELL_SCRIPTS="$SCRIPT"/shell


#### Running tools
# Running ABB_filter tool
if [ "$TOOL" == "ABB_filter" ]; then
  echo "Running ABB_filter tool"
  if [ -v $OUTFILE ];then
    $PYTHON $PYTHON_SCRIPTS/ABB_filter_Parallel.py -i $INPUT_FILE -outfile $OUTFILE -abb_filter $ABB_FILTER
    
    else
        outfile=$(basename $OUTFILE)
        OUT_FILE=$OUT_FOLDER/$outfile
        $PYTHON $PYTHON_SCRIPTS/ABB_filter_Parallel.py -i $INPUT_FILE -outfile $OUT_FILE -abb_filter $ABB_FILTER
    fi
fi

# Running ABB_list tool
if [ "$TOOL" == "ABB_list" ]; then
  echo "Running ABB_list tool"
    if [ -v $OUTFILE ];then
    $PYTHON $PYTHON_SCRIPTS/ABB_list_Parallel.py -i $INPUT_FILE -outfile $OUTFILE

    else  
        outfile=$(basename $OUTFILE)
        OUT_FILE=$OUT_FOLDER/$outfile
        $PYTHON $PYTHON_SCRIPTS/ABB_list_Parallel.py -i $INPUT_FILE -outfile $OUT_FILE
    fi
fi

# Running ABB_list tool
if [ "$TOOL" == "ABB_association" ]; then
  echo "Running ABB_association tool"
  bash $SHELL_SCRIPTS/ABB_association.sh -i $INPUT_FILE -o $OUT_FOLDER -cases $CASE_IDS -controls $CONTROL_IDS -genes $GENES_FILE -abb_file $ABB_FILE
fi

# Running ABB_annotate tool
if [ "$TOOL" == "ABB_annotation" ]; then
  echo "Running ABB_annotation tool"
  if [ -v $OUTFILE ];then
    bash $SHELL_SCRIPTS/ABB_annotate.sh -i $INPUT_FILE -outfile $OUTFILE -abb_file $ABB_FILE -abb_filter $ABB_FILTER
    
    else
        outfile=$(basename $OUTFILE)
        OUT_FILE=$OUT_FOLDER/$outfile
        bash $SHELL_SCRIPTS/ABB_annotate.sh -i $INPUT_FILE -outfile $OUT_FILE -abb_file $ABB_FILE -abb_filter $ABB_FILTER
    fi
fi

rm -f ARGPARSER_temp
