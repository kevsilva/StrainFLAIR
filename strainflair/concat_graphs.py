#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess
import sys
import getopt
import logging # log.txt with times
import time
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

def worker(q: Queue):
    while True:
        l_clusters = q.get()
        if not l_clusters:
            q.put(None)
            break
        concat(l_clusters)

def concat(l_clusters: list):

    subprocess.run(['vg', 'ids', '-j', '-c']+l_clusters)
    with tempfile.NamedTemporaryFile(dir=input_dir,delete=False,suffix=".vg") as out:
        print(f"processing {out.name}")
        subprocess.run(['vg','combine'] + l_clusters, stdout=out)
    subprocess.run(['rm']+l_clusters)

def file_generator(input_dir):
    for entry in os.scandir(input_dir):
        if entry.name.endswith('.vg'):
            yield entry.name

def usage():
    print(f"Usage: python3 {sys.argv[0]} -i input_dir -s step")

#if __name__ == "__main__":
def concat_graphs_main():
    # check arguments

    input_dir = None 
    step = 1000

    try:
        opts, _ = getopt.getopt(sys.argv[1:], "hi:s:")
    
    except getopt.GetoptError as err:
        # print help information and exit:
        print(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
        elif o in ("-i"):
            input_dir = a      
        elif o in ("-s"):
            step = int(a) 
        else:
            assert False, "unhandled option"
    if not input_dir or not step: 
        usage()
        exit()
        
    # logger
    setup_logger("logger", f"{input_dir}/concat_grahs_log.txt")
    logger = logging.getLogger("logger")

    # start

    stage = 1
    while True:

        with Timer() as _t:
            
            print(f"start step {stage}")
            nb_files = int(subprocess.check_output(f"find {input_dir} -maxdepth 1 -name '*.vg' | wc -l",shell=True))
            if nb_files == 1:
                final_file = [f for f in os.listdir(input_dir) if f.endswith(".vg")][0]
                os.rename(f"{input_dir}/{final_file}", f"{input_dir}/all_graphs.vg")
                break

            # queue initialization
            q = Queue() 

            # for each steps in parallel
            processes = Pool(initializer=worker, initargs=(q,))

            print("fill the queue")
            # fill the queue with list of vg files
            ct = 0
            l_clusters = []
            for vg_file in file_generator(input_dir):
                l_clusters.append(f"{input_dir}/{vg_file}")
                ct+=1
                if ct==step:
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