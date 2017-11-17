# ABB tool, 
This tool is able to detect false positive calls. It is based on a new strategy to identify systematic sequencing or alignment errors leading to false positive variant calls based on the recurrence of the allele balance bias (paper in preparation). The main applications of this tool are four:
* 'ABB_list' obtains a new callability score list based on a new subset of samples. This list label positions of the genome with values between 0 and 1, which represent the precision of being a systematic error.
* 'ABB_annotation' Annotates and filters variant from input vcfs with an existing or new ABB list under a specific threshold of ABB.
* 'ABB_filter' obtains ABB from vcfs on-fly and filter variants based on a specific ABB threshold.
* 'ABB_association' detects genes/regions and variants which ABB could explain the significant association found in case-control studies.

This tools is though to be run in DNA data, taking as input VCF files and other parameters to remove variants prone to systematic errors (highly enriched with false positive calls). Statistical analysis and results are based on the paper : paper Muyas et al.


## Get ABB tool code 
You will need to run `git clone  ` to get ABB tool. 

Firstly, you will need to install some software
- [argparse.bash](https://github.com/nhoffman/argparse-bash)

```
cd ABB
cd source
wget https://raw.githubusercontent.com/nhoffman/argparse-bash/master/argparse.bash
chmod +x argparse.bash
```

- [shc](https://github.com/neurobin/shc)
Install shell script compiler. It can be downloaded and installed or you can just download a compiled binary package like next:

```
wget https://github.com/neurobin/shc/releases/download/3.9.6/shc-3.9.6-bin-amd64-i386-arm64-armhf-ppc64el.tar.gz
```

- [Python 2.7](https://www.python.org/download/releases/2.7/), with next python packages

+ numpy
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


