#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Description:
Usage:
    genes_prediction.py (--inSeq <str>) (--outDir <dir>) (--outFile <str>) [--lenExtend <int>]

Options:
    -h --help
    -i --inSeq <file>          List of input sequences in fasta files
    -o --outDir <dir>          Output directory
    -f --outFile <file>        Output file name
    -e --lenExtend <int>       Extend length of the sequences [default: 0]
"""

import sys

try:
    import docopt
except ModuleNotFoundError:
    print("genes_prediction.py requires module docopt, to install: pip3 install docopt")
    sys.exit()
    
###############################################################################
import subprocess
import os
from Bio import SeqIO

def main():
    
    args = docopt.docopt(__doc__)
    fasta_list = args["--inSeq"]
    output_dir = args["--outDir"]
    output_name = args["--outFile"]
    len_extend = int(args["--lenExtend"])
        
    if not os.path.exists(output_dir):
        subprocess.check_call(f"mkdir {output_dir}",shell=True)
    
    N = 0 # initialization number of genes predicted
    with open(fasta_list) as f:
        for fasta_file in f:
            
            fasta_file = fasta_file.strip()
            fasta_basename = os.path.basename(os.path.splitext(fasta_file)[0])
            subprocess.check_call(f"prodigal -i {fasta_file} -o {output_dir}/predGenes_{fasta_basename}.txt -d {output_dir}/predGenes_{fasta_basename}.fasta",shell=True)
            p1 = subprocess.Popen(["grep",">",f"{output_dir}/predGenes_{fasta_basename}.fasta"], stdout=subprocess.PIPE)
            id_genes = p1.communicate()[0]
            id_genes = id_genes.decode().strip().split("\n")
            n = len(id_genes)
            print(f"\nNumber of genes predicted for {fasta_basename} = {n}\n")
            N += n
            subprocess.check_call(f"cat {output_dir}/predGenes_{fasta_basename}.fasta >> {output_dir}/{output_name}.fasta",shell=True)
            
            if len_extend != 0:
                record_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
                for gene in id_genes:
                    gene = gene.split("#")
                    idt = '_'.join(gene[0].split("_")[:-1])[1:]
                    start = int(gene[1])-len_extend
                    if start < 1: start = 1
                    end = int(gene[2])+len_extend
                    if end > len(record_dict[idt]): end = len(record_dict[idt])
                    with open(f"{output_dir}/{output_name}_extended{len_extend}bp.fasta", "a") as output_handle:
                        record_dict[idt].id = ' # '.join([gene[0],str(start),str(end)]+gene[3:])[1:]
                        record_dict[idt].description = ""
                        SeqIO.write(record_dict[idt][start-1:end], output_handle, "fasta")
            
    print(f"\nTotal genes predicted = {N}\n")

###############################################################################
if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nStopping from keyboard")