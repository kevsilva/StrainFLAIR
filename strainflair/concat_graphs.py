#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess # run vg
import argparse # manage arguments
import logging # log file
import time # timer
import os
import tempfile
from multiprocessing import Queue, Process, Pool, Lock, cpu_count, Manager # multiprocessing


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

def worker(q: Queue, input_dir: str):
    while True:
        l_clusters = q.get()
        if not l_clusters:
            q.put(None)
            break
        concat(l_clusters, input_dir)

def concat(l_clusters: list, input_dir: str):

    subprocess.run(['vg', 'ids', '-j', '-c']+l_clusters)
    with tempfile.NamedTemporaryFile(dir=input_dir,delete=False,suffix=".vg") as out:
        print(f"processing {out.name}")
        subprocess.run(['vg','combine'] + l_clusters, stdout=out)
    subprocess.run(['rm']+l_clusters)

def file_generator(input_dir: str):
    for entry in os.scandir(input_dir):
        if entry.name.endswith('.vg'):
            yield entry.name

def concat_graphs_main():

    # arguments
    parser = argparse.ArgumentParser()
    required_args = parser.add_argument_group('required arguments')
    required_args.add_argument('-i', '--input_dir', type=str, required=True, help='Graphs to concatenate directory.')
    parser.add_argument('-s', '--step', type=int, default=1000, help='Number of graphs to concatenate at once. Default: 1000')
    args = parser.parse_args()
        
    # logger
    setup_logger("logger", f"{args.input_dir}/concat_grahs.log")
    logger = logging.getLogger("logger")

    # start

    stage = 1
    while True:

        with Timer() as _t:
            
            print(f"start step {stage}")
            nb_files = int(subprocess.check_output(f"find {args.input_dir} -maxdepth 1 -name '*.vg' | wc -l",shell=True))
            if nb_files == 1:
                final_file = [f for f in os.listdir(args.input_dir) if f.endswith(".vg")][0]
                os.rename(f"{args.input_dir}/{final_file}", f"{args.input_dir}/all_graphs.vg")
                break

            # queue initialization
            q = Queue() 

            # for each steps in parallel
            processes = Pool(initializer=worker, initargs=(q,args.input_dir))

            print("fill the queue")
            # fill the queue with list of vg files
            ct = 0
            l_clusters = []
            for vg_file in file_generator(args.input_dir):
                l_clusters.append(f"{args.input_dir}/{vg_file}")
                ct+=1
                if ct==args.step:
                    q.put(l_clusters)
                    ct = 0
                    l_clusters = []
            if len(l_clusters)>1:
                q.put(l_clusters)
            q.put(None)

            # end multiprocessing
            processes.close()
            processes.join()
            print(f"end of step {stage}")

        logger.info(f"Step {stage} done in: {_t.t}")
        stage += 1