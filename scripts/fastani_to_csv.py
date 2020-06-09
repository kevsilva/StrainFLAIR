#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Description:
Usage:
    fastani_to_csv.py (--inSeq <str>) (--outDir <dir>)

Options:
    -h --help
    --inSeq <file>          input file containing paths to sequences
    --outDir <dir>          Output directory
"""

import sys

try:
    import docopt
except ModuleNotFoundError:
    print("fastani_to_csv.py requires module docopt, to install: pip3 install docopt")
    sys.exit()

###############################################################################
import os
import subprocess
import pandas as pd

def main():
    
    args = docopt.docopt(__doc__)
    sequences_file = args["--inSeq"]
    output_dir = args["--outDir"]
    
    # make output folder
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # run fastani
    subprocess.check_call(f"fastani --ql {sequences_file} --rl {sequences_file} -o {output_dir}/fastani.out",shell=True)
    
    # convert to csv
    strains = []
    with open(sequences_file, 'r') as f:
        for line in f:
            strains.append(os.path.splitext(os.path.basename(line))[0])
    df = pd.DataFrame(100,columns=strains,index=strains)
    with open(f"{output_dir}/fastani.out", 'r') as f:
        for line in f:
            line = line.split()
            if line[0] != line[1]:
                i = os.path.splitext(os.path.basename(line[0]))[0]
                j = os.path.splitext(os.path.basename(line[1]))[0]
                df.loc[i,j] = line[2]
    df.to_csv(f'{output_dir}/fastani.csv')

###############################################################################
if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nStopping from keyboard")