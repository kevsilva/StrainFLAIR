#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess # execute vg
from Bio import SeqIO # conversion to dict and to sequence
import argparse # manage arguments
import os
from multiprocessing import Queue, Process, Pool, Lock, cpu_count # multiprocessing
from tempfile import TemporaryDirectory # temporary fasta and paf file for each cluster
import time # times stored in log
import logging # log file
import pickle # save d_clusters

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

def clstr2dict(clstr_file: str):
    '''
    Convert cluster file into dictionary 
    key = cluster id
    value = dict with gene_list and length of representative
    '''
    d_clusters = {}
    with open(clstr_file) as f:
        for line in f:
            line = line.strip()
            if line[0] == ">":
                cluster_name = line[1:].replace(" ","_")
                d_clusters[cluster_name] = {}
                d_clusters[cluster_name]['genes_list'] = []
            else:
                split_line = line.split()
                id_seq = split_line[2][1:-3] # delete > and ...
                d_clusters[cluster_name]['genes_list'].append(id_seq)
                # get the length of the representative sequence
                if split_line[-1] == '*':
                    d_clusters[cluster_name]['len_rep'] = int(split_line[1][:-3])
    return d_clusters

def worker(q: Queue, d_clusters: dict, d_IdToSeq: dict, out_dir: str):
    while True:
        cluster_name = q.get()
        if not cluster_name:
            q.put(None)
            break
        cluster2graph(cluster_name, d_clusters, d_IdToSeq, out_dir)

def cluster2graph(cluster_name: str, d_clusters: dict, d_IdToSeq: dict, out_dir: str):

    '''
    Graph construction
    Create a temporary fasta with all the gene sequences from the cluster specified.
    If there is only one sequence, create a linear graph with vg construct.
    If there is more than one sequence, use minimap2 for sequence alignment and seqwish to convert it into a graph.
    '''

    # creating temporary files in a temporary folder
    with TemporaryDirectory() as temp_dir:

        # create temporary merged fasta for all sequences of the cluster
        with open(f"{temp_dir}/cluster_temp.fasta", "w") as f:
            for idt in d_clusters[cluster_name]['genes_list']:
                SeqIO.write([d_IdToSeq[idt]], f, "fasta")   
    
        # build the graph
        if len(d_clusters[cluster_name]['genes_list']) == 1: # if only one sequence in the cluster, just build a linear graph with vg construct
            subprocess.run(f"vg construct -r {temp_dir}/cluster_temp.fasta -m 256 > {out_dir}/{cluster_name}.vg",stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
        else:
            subprocess.run(f"minimap2 -cx asm20 -X -t 8 {temp_dir}/cluster_temp.fasta {temp_dir}/cluster_temp.fasta | gzip > {temp_dir}/cluster_temp.paf.gz",stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
            subprocess.run(f"seqwish -s {temp_dir}/cluster_temp.fasta -p {temp_dir}/cluster_temp.paf.gz -b {temp_dir}/cluster_temp.work -g {temp_dir}/cluster_temp.gfa",shell=True)
            subprocess.run(f"vg view -Fv {temp_dir}/cluster_temp.gfa | vg mod -n -X 256 - | vg sort - > {out_dir}/{cluster_name}.vg",shell=True)
            # -n can mess up the graph, check for its integrity otherwise redo the graph without -n
            p1 = subprocess.Popen(["vg","validate",f"{out_dir}/{cluster_name}.vg"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out,err = p1.communicate()
            if err.decode() != "":
                subprocess.run(f"vg view -Fv {temp_dir}/cluster_temp.gfa | vg mod -X 256 - | vg sort - > {out_dir}/{cluster_name}.vg",shell=True)

if __name__ == "__main__":
    
    # arguments
    parser = argparse.ArgumentParser()
    required_args = parser.add_argument_group('required arguments')
    required_args.add_argument('-s', '--in_sequences', type=str, required=True, help='Input gene sequences in a fasta file.')
    required_args.add_argument('-c', '--in_clusters', type=str, required=True, help='Input clusters file (clstr ouput file from CD-HIT).')
    parser.add_argument('-o', '--out_dir', type=str, default="graphs", help='Output directory. Default: "graphs"')
    parser.add_argument('-l', '--min_length', type=float, default=0, help='Threshold on the cluster representative length. Default: 0')
    args = parser.parse_args()

    # create output directory
    if not os.path.exists(args.out_dir):
        subprocess.run(["mkdir",args.out_dir])

    # logger
    setup_logger("logger", f"{args.out_dir}/graph_construction.log")
    logger = logging.getLogger("logger")

    # start pipeline

    # sequence file into dictionary (key = sequence Id, value = sequence)
    d_IdToSeq = SeqIO.to_dict(SeqIO.parse(args.in_sequences, "fasta"))

    # cluster file into dictionary (key = cluster id, value = gene list + length of representative sequence)
    d_clusters = clstr2dict(args.in_clusters)
    pickle_out = open(args.out_dir+"/dict_clusters.pickle","wb")
    pickle.dump(d_clusters, pickle_out)
    pickle_out.close()

    with Timer() as _t:
        # queue initialization (for clusters)
        nb_clusters = 0
        q = Queue() 

        # construct graph for each cluster in parallel
        processes = Pool(initializer=worker, initargs=(q, d_clusters, d_IdToSeq, args.out_dir))

        # fill the queue with the clusters to process
        for cluster_name in d_clusters:
            if d_clusters[cluster_name]['len_rep'] >= args.min_length:
                nb_clusters += 1
                q.put(cluster_name)
        q.put(None)

        # end multiprocessing
        processes.close()
        processes.join()
    logger.info(f"Number of graphs built: {nb_clusters}")
    logger.info(f"Multiprocessed graphs building done in: {_t.t}")

