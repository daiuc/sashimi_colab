# sashimi_colab


## Installation

Note, only tested on python 3.8.15, pygenometracks 3.6

```
# create conda environment for pygenometracks

mamba create -n pygenometracks -c bioconda  pygenometracks=3.6


# copy SashimiBigWigTrack.py to your pygenometracks library path
cp ./SashimiBigWigTrack.py ./miniconda3/envs/pygenometracks/lib/python3.8/site-packages/pygenometracks/

```


use `sashimi_ingredients.py` to prepare for data.

modify `plot.sh` to plot sashimi plots

run `bash plot.sh` to plot. 


Follow [this](https://github.com/Zepeng-Mu/pyGenomeTracks/tree/master) for pyGenomeTracks installation.

Note, counts files need to be space delimited. 

Also, presect rows of counts to include just the intron cluster you want. 

You also need to index your vcf file. 

also make a list of your sample IDs (one id per line, no header) in a text file.
