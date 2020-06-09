#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Description:
Usage:
    graphs_construction.py (--inSeq <str>) (--inClusters <str>) (--outDir <dir>) [--min_length <float>] [--verbose] [--debug]

Options:
    -h --help
    --inSeq <file>          input file containing gene sequences
    --inClusters <file>     input file containing gene clusters (CD-HIT ouput)
    --outDir <dir>          Output directory
    --min_length <float>     ignore clusters which representative sequence length is below threshold [default: 0]
    --verbose
    --debug
"""

import sys

try:
    import docopt
except ModuleNotFoundError:
    print("graphs_construction.py requires module docopt, to install: pip3 install docopt")
    sys.exit()

###############################################################################
import os
import subprocess
from Bio import SeqIO
import pickle
import time
import logging

def main():
    
    args = docopt.docopt(__doc__)
    verbose = args["--verbose"]
    debug = args["--debug"]
    sequences_file = args["--inSeq"]
    clusters_file = args["--inClusters"]
    output_dir = args["--outDir"]
    min_length = args["--min_length"]
    
    # logger
    global logger    
    setup_logger("logger", f"{output_dir}/logs.txt", verbose, debug)
    logger = logging.getLogger("logger")
    
    # set temporary folder
    cwd = os.getcwd()
    if os.path.exists(f"{cwd}/gc_tmp"):
        subprocess.check_call(f"rm -r {cwd}/gc_tmp",shell=True)
    subprocess.check_call(f"mkdir {cwd}/gc_tmp",shell=True)
    
    # make output folder
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # sequence file into dictionary (key = sequence Id, value = sequence)
    d_IdToSeq = SeqIO.to_dict(SeqIO.parse(sequences_file, "fasta"))
    
    # cluster file into dictionary (key = cluster id, value = gene_list and length of representative)
    d_clusters = {}
    with open(clusters_file) as f:
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
    
    # graph construction
    N = 0 # keep track of nodes ID
    for cluster in d_clusters:
        logger.info(f"Building graph for {cluster}...")
        with Timer() as _t:
            # check min_length argument
            if d_clusters[cluster]['len_rep'] < min_length:
                logger.info(f"{cluster} and following have been ignored, the representative sequences are too short (check -l option).")
                break # cluster file is ordered by representative sequence length so we can stop the loop
            else:
                # create fasta for previous cluster
                with open(f"{cwd}/gc_tmp/cluster_temp.fasta", "w") as f:
                    for idx in d_clusters[cluster]['genes_list']:
                        SeqIO.write([d_IdToSeq[idx]], f, "fasta")   
                # build the graph
                if len(d_clusters[cluster]['genes_list']) == 1: # if only one sequence in the cluster, just build a linear graph with vg construct
                    subprocess.check_call(f"vg construct -r {cwd}/gc_tmp/cluster_temp.fasta -m 256 | vg ids -i {N} -c - > {args.output_folder}/{cluster}.vg",shell=True)
                    subprocess.check_call(f"rm -r {cwd}/gc_tmp/*.fai",shell=True)
                else:
                    subprocess.check_call(f"minimap2 -cx asm20 -X -t 4 {cwd}/gc_tmp/cluster_temp.fasta {cwd}/gc_tmp/cluster_temp.fasta | gzip > {cwd}/gc_tmp/cluster_temp.paf.gz",shell=True)
                    subprocess.check_call(f"seqwish -s {cwd}/gc_tmp/cluster_temp.fasta -p {cwd}/gc_tmp/cluster_temp.paf.gz -b {cwd}/gc_tmp/cluster_temp.work -g {cwd}/gc_tmp/cluster_temp.gfa",shell=True)
                    subprocess.check_call(f"vg view -Fv {cwd}/gc_tmp/cluster_temp.gfa | vg mod -n -X 256 - | vg ids -i {N} -c - | vg sort - > {args.output_folder}/{cluster}.vg",shell=True)
                n = int(subprocess.check_output(['vg', 'stats', '-N', output_dir+"/"+cluster+".vg"]))
                d_clusters[cluster]['nodes_list'] = list(range(N+1,n+N+1))
                N += n
        all_time = _t.t
        logger.info(f"Graph for {cluster} done in: {all_time}")
    subprocess.check_call(f"rm -r {cwd}/gc_tmp",shell=True)
    logger.info(f"Concatenating the graphs...")
    subprocess.check_call('cat $(for i in $(ls -v '+output_dir+'); do echo '+output_dir+'/${i}; done) > '+output_dir+'/final_graph.vg',shell=True)
    subprocess.check_call(f"vg view {output_dir}/final_graph.vg > {output_dir}/final_graph.gfa",shell=True)
    
    pickle_out = open(output_dir+"/dict_clusters.pickle","wb")
    pickle.dump(d_clusters, pickle_out)
    pickle_out.close()
    
    # graph indexation
    with Timer() as _t:
        subprocess.check_call(f"vg index -x {output_dir}/final_graph.xg {output_dir}/final_graph.vg",shell=True)
    all_time = _t.t
    logger.info(f"xg indexing done in: {all_time}")
    with Timer() as _t:
        subprocess.check_call(f"vg prune {output_dir}/final_graph.vg | vg index -g {output_dir}/final_graph.gcsa -",shell=True)
    all_time = _t.t
    logger.info(f"gcsa indexing done in: {all_time}")
    with Timer() as _t:
        subprocess.check_call(f"vg snarls {output_dir}/final_graph.vg > {output_dir}/final_graph.snarls",shell=True)
    all_time = _t.t
    logger.info(f"snarls search done in: {all_time}")

###############################################################################

class Timer:
    def __enter__(self):
        self.t1 = time.time()
        return self
    
    def __exit__(self, *args):
        self.t2 = time.time()
        hours, rem = divmod(self.t2-self.t1, 3600)
        minutes, seconds = divmod(rem, 60)
        self.t = "{:0>2}:{:0>2}:{:05.2f}".format(int(hours), int(minutes), seconds)

###############################################################################

def setup_logger(name, log_path, verbose, debug):
    l = logging.getLogger(name)
    formatter = logging.Formatter("%(asctime)s -- %(levelname)s -- %(message)s")
    fileHandler = logging.FileHandler(log_path, mode="w")
    fileHandler.setFormatter(formatter)

    l.setLevel(logging.DEBUG)
    l.addHandler(fileHandler)

    if verbose or debug:
        streamHandler = logging.StreamHandler(sys.stdout)
        streamHandler.setFormatter(formatter)
        l.addHandler(streamHandler)

###############################################################################
if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nStopping from keyboard")