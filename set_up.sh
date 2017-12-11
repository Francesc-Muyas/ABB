#!/bin/bash

#### Get options
while getopts "r:t:p:a:" opt; do
  case $opt in
    t) TOOL_PATH=$OPTARG;;
    r) RSCRIPT_PATH=$OPTARG;;
    p) PYTHON_PATH=$OPTARG;;
    a) ARGPARSE_PATH=$OPTARG;;
    \?)
      echo "Invalid option: -"$OPTARG
      echo -e "GUIDE:\n -t: PATH/TO/ABB\n -r: PATH/TO/Rscript\n -p: PATH/TO/python\n -a: PATH/TO/argparse.bash\n"
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      echo -e "GUIDE:\n -t: PATH/TO/ABB\n -r: PATH/TO/Rscript\n -p: PATH/TO/python\n -a: PATH/TO/argparse.bash\n"
      exit 1
      ;;
  esac
done

if [ -v $TOOL_PATH ];then
    echo "-t option is missed" >&2
    echo -e "GUIDE:\n -t: PATH/TO/ABB\n -r: PATH/TO/Rscript\n -p: PATH/TO/python\n -a: PATH/TO/argparse.bash\n"
    exit 1
fi

if [ -v $RSCRIPT_PATH ];then
    echo "-r option is missed" >&2
    echo -e "GUIDE:\n -t: PATH/TO/ABB\n -r: PATH/TO/Rscript\n -p: PATH/TO/python\n -a: PATH/TO/argparse.bash\n"
    exit 1
fi

if [ -v $PYTHON_PATH ];then
    echo "-p option is missed" >&2
    echo -e "GUIDE:\n -t: PATH/TO/ABB\n -r: PATH/TO/Rscript\n -p: PATH/TO/python\n -a: PATH/TO/argparse.bash\n"
    exit 1
fi

if [ -v $ARGPARSE_PATH ];then
    echo "-a option is missed" >&2
    echo -e "GUIDE:\n -t: PATH/TO/ABB\n -r: PATH/TO/Rscript\n -p: PATH/TO/python\n -a: PATH/TO/argparse.bash\n"
    exit 1
fi

# Defining paramenters
for i in shell/*sh;do
	sed -i "s|^ABB_PATH=.*|ABB_PATH=${TOOL_PATH}|" $i
	sed -i "s|^PYTHON=.*|PYTHON=${PYTHON_PATH}|" $i
	sed -i "s|^RSCRIPT=.*|RSCRIPT=${RSCRIPT_PATH}|" $i
	sed -i "s|^ARGPARSE=.*|ARGPARSE=${ARGPARSE_PATH}|"  $i
done
