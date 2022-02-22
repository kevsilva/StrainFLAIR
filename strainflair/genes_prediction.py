#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess # execute vg
import os
from Bio import SeqIO # conversion to dict and to sequence
import argparse
from multiprocessing import Queue, Process, Pool, Lock, cpu_count # multiprocessing
import time # times stored in log
import logging # log.txt with times

class Timer:
    def __enter__(self):
        self.t1 = time.time()
        return self
    
    def __exit__(self, *args):
        self.t2 = time.time()
        hours, rem = divmod(self.t2-self.t1, 3600)
        minutes, seconds = divmod(rem, 60)
        self.t = "{:0>2}:{:0>2}:{:05.2f}".format(int(hours), int(minutes), seconds)

def setup_logger(name, log_path):
    l = logging.getLogger(name)
    formatter = logging.Formatter("%(asctime)s -- %(levelname)s -- %(message)s")
    fileHandler = logging.FileHandler(log_path, mode="w")
    fileHandler.setFormatter(formatter)

    l.setLevel(logging.DEBUG)
    l.addHandler(fileHandler)

def worker(q: Queue, out_dir: str, len_extend: int):
    while True:
        fasta_file = q.get()
        if not fasta_file:
            q.put(None)
            break
        predict_genes(fasta_file, out_dir, len_extend)


def predict_genes(fasta_file: str, out_dir: str, len_extend: int):
    fasta_basename = os.path.basename(os.path.splitext(fasta_file)[0])
    subprocess.check_call(f"prodigal -i {fasta_file} -o {out_dir}/predGenes_{fasta_basename}.txt -d {out_dir}/predGenes_{fasta_basename}.fasta",shell=True)
    p1 = subprocess.Popen(["grep",">",f"{out_dir}/predGenes_{fasta_basename}.fasta"], stdout=subprocess.PIPE)
    id_genes = p1.communicate()[0]
    id_genes = id_genes.decode().strip().split("\n")

    if len_extend != 0:
        record_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
        for gene in id_genes:
            gene = gene.split("#")
            idt = '_'.join(gene[0].split("_")[:-1])[1:]
            start = int(gene[1])-len_extend
            if start < 1: start = 1
            end = int(gene[2])+len_extend
            if end > len(record_dict[idt]): end = len(record_dict[idt])
            with open(f"{out_dir}/predGenes_{fasta_basename}_extended{len_extend}bp.fasta", "a") as output_handle:
                record_dict[idt].id = ' # '.join([gene[0],str(start),str(end)]+gene[3:])[1:]
                record_dict[idt].description = ""
                SeqIO.write(record_dict[idt][start-1:end], output_handle, "fasta")

def genes_prediction_main():
    
    # arguments
    parser = argparse.ArgumentParser()
    required_args = parser.add_argument_group('required arguments')
    required_args.add_argument('-s', '--in_sequences', type=str, required=True, help='Input sequence(s) in a fasta file or in a text file listing fasta file paths.')
    parser.add_argument('-o', '--out_dir', type=str, default="predicted_genes", help='Output directory. Default: "predicted_genes"')
    parser.add_argument('-l', '--len_extend', type=int, default=0, help='Number of bases to extend the length of the predicted genes at both ends. Default: 0')
    args = parser.parse_args()

    # start

    # create output directory
    if not os.path.exists(args.out_dir):
        subprocess.run(["mkdir",args.out_dir])

    # setup logger
    setup_logger("logger", f"{args.out_dir}/genes_prediction.log")
    logger = logging.getLogger("logger")

    # genes prediction
    with Timer() as _t:
        # queue initialization (for files)
        q = Queue() 

        # predict genes for each fasta in parallel
        processes = Pool(initializer=worker, initargs=(q, args.out_dir, args.len_extend))

        # fill the queue with the files to process
        if args.in_sequences.endswith(".fasta") or args.in_sequences.endswith(".fna"):
            q.put(args.in_sequences)
        else:
            with open(args.in_sequences) as f:
                for fasta_file in f:
                    q.put(fasta_file.rstrip("\n"))
        q.put(None)

        # end multiprocessing
        processes.close()
        processes.join()
    
    # concatenate
    subprocess.run(f"cat `ls -v {args.out_dir}/*.fasta | grep -v 'extended'` > {args.out_dir}/all_genes.fasta",shell=True)
    if args.len_extend != 0: subprocess.run(f"cat $(ls -v {args.out_dir}/*extended*.fasta) > {args.out_dir}/all_genes_extended.fasta",shell=True)

    logger.info(f"Multiprocessed genes prediction done in: {_t.t}")

    nb_genes = 0
    with open(f"{args.out_dir}/all_genes.fasta","r") as f:
        for line in f:
            if line.startswith(">"):
                nb_genes += 1
    logger.info(f"Total number of predicted genes: {nb_genes}")