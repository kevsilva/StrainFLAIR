#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Description:
Usage:
    read_simulator.py (--inSeq <str>) (--outDir <str>) (--number <int>) (--len <int>) (--error <float>) (--mode <str>)

Options:
    -h --help
    --inSeq <file>          input file containing paths to sequences
    --outDir <dir>          Output directory
    --number <int>          number of reads to generate
    --len <int>             read length
    --error <float>         error rate
    --mode <str>            use "vg" or "mutareads"
"""

import sys

try:
    import docopt
except ModuleNotFoundError:
    print("read_simulator.py requires module docopt, to install: pip3 install docopt")
    sys.exit()

###############################################################################
import os
import subprocess

def main():
    
    args = docopt.docopt(__doc__)
    sequences_file = args["--inSeq"]
    output_dir = args["--outDir"]
    n = int(args["--number"])
    l = int(args["--len"])
    e = float(args["--error"])
    m = args["--mode"]
    
    with open(sequences_file, 'r') as f:
        for line in f:
            
            line = line.strip()
            current_sequence = os.path.splitext(line)[0]
            current_basename = os.path.basename(current_sequence)
            
            if m == "vg":
                subprocess.check_call(f"vg construct -r {line} -m 256 > {current_sequence}.vg",stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
                subprocess.check_call(f"vg index -x {current_sequence}.xg {current_sequence}.vg",stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
                subprocess.check_call(f"vg sim -x {current_sequence}.xg -l {l} -n {n} -e {e} -a | vg view -X - > {output_dir}/vgsim_n{n}_e{e}_{current_basename}.fastq",stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
                subprocess.check_call(f"sed -i 's/@/@{current_basename}_/g' {output_dir}/vgsim_n{n}_e{e}_{current_basename}.fastq",stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
            if m == "mutareads":
                subprocess.check_call(f"mutareads {line} {output_dir}/mutareads_n{n}_e{e}_{current_basename} {n} {l} {e} 0 0",stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
                subprocess.check_call(f"perl ./softwares/fasta_to_fastq.pl {output_dir}/mutareads_n{n}_e{e}_{current_basename}.fasta > {output_dir}/mutareads_n{n}_e{e}_{current_basename}.fastq",stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
                subprocess.check_call(f"sed -i 's/@/@{current_basename}_/g' {output_dir}/mutareads_n{n}_e{e}_{current_basename}.fastq",stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)

###############################################################################
if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nStopping from keyboard")