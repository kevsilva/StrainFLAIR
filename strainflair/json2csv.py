#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import time # to delete after test phase
import pickle # for reading the dictionnary
import numpy as np
import json
import sys # update_progress
import getopt
import os # for size files 
import re
import networkx as nx # for incompatibility reads graph/network
import itertools 
import scipy.signal # for peak detection

def update_progress(progress):
    """
    update_progress() : Displays or updates a console progress bar
    Accepts a float between 0 and 1. Any int will be converted to a float.
    A value under 0 represents a 'halt'.
    A value at 1 or bigger represents 100%
    https://stackoverflow.com/questions/3160699/python-progress-bar 
    """
    barLength = 50 # Modify this to change the length of the progress bar
    status = ""
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength*progress))
    text = "\rPercent: [{0}] {1}% {2}".format( "#"*block + "-"*(barLength-block), round(progress*100,2), status)
    sys.stderr.write(text)
    sys.stderr.flush()
    
def get_file_size(f):
    old_file_position = f.tell()
    f.seek(0, os.SEEK_END)
    size = f.tell()
    f.seek(old_file_position, os.SEEK_SET)
    return size

reverse_sign = lambda x: '-' if (x=='+') else '+'

def canonical(node_list: str):
    """
    Returns the canonical representation of a node list. 
    We define the canonical as the min between a list eg: 684619+,684620+,684618+ and its reverse eg: 684618-,684620-,684619-
    wrt to the first unsigned value (here 684619 or 684618). In this case we return 684618-,684620-,684619-
    """
    splitted_node_list = node_list.split(',')
    # if only one node, sign is always "+"
    if len(splitted_node_list) == 1:
        return splitted_node_list[0][:-1]+"+"
    else:
        if int(splitted_node_list[0][:-1]) < int(splitted_node_list[-1][:-1]): 
            return node_list
        # else
        rev_list = ','.join([val[:-1]+reverse_sign(val[-1]) for val in reversed(splitted_node_list)])
        return rev_list

def canonical_tiny(node_list: list):
    """
    Returns the canonical representation of a node list.
    For smaller lists (reads mapping paths).
    Check the sum of the difference direction of each intervalles between successive nodes
    """
    if len(node_list) == 1 or sum([t - s > 0 for s, t in zip(node_list, node_list[1:])])/(len(node_list)-1) > 0.5:
        return node_list
    return node_list[::-1]

def get_accession_number(line: str):
    accession_number = None
    line = '_'.join(line.split("_")[:-1]) # delete gene id
    line = line.split("|")
    for element in line:
        if re.match("[a-zA-Z]+_?[a-zA-Z]*[0-9]+", element):
            accession_number = element
            break
    return accession_number

def check_incomp_comp(c1:list, c2:list):
    common_nodes = list(set(c1) & set(c2))
    if common_nodes:
        ancor = common_nodes[0]
        ancor_ids_c1 = [i for i, v in enumerate(c1) if v == ancor]
        ancor_ids_c2 = [i for i, v in enumerate(c2) if v == ancor]
        for ancor_id_c1 in ancor_ids_c1:
            for ancor_id_c2 in ancor_ids_c2:
                l_left = min(len(c1[:ancor_id_c1]),len(c2[:ancor_id_c2]))
                l_right = min(len(c1[ancor_id_c1:]),len(c2[ancor_id_c2:]))
                if c1[ancor_id_c1-l_left:ancor_id_c1+l_right] == c2[ancor_id_c2-l_left:ancor_id_c2+l_right]:
                    return 1
        return -1
    # no common nodes
    return 0

def sublist(lst1:list, lst2:list):
   return set(lst1).issubset(set(lst2)) or set(lst2).issubset(set(lst1))

class Node:
    def __init__(self, sequence: str):
        self.len_sequence = len(sequence)   # we need to remind the length of the sequence of each node for statistical computations
        self.traversed_path = []            # ids (int) of the traversed paths. This corresponds to indexes in the Pangenome.paths list

class Path:
    def __init__(self):
        self.node_ids   = []                # a set of node ids (ints)
        self.cluster_id = None              # id (string) of the cluster this path belongs to 
        self.strain_ids = {}                # ids of the genes that generated this path with their counts. 
        #----
        self.unique_mapped_abundances  = [] # for each node of the path (ordered as `self.nodes`), store the coverage of mapped reads (each node of a mapped path is set to one, except the two extreme than are usually not 100% covered by the mapped sequence)
        self.total_mapped_unique_reads = 0  #  number of reads with unique mapping on this path. 
        self.total_mapped_unique_reads_normalized = 0 # number of coverage ratio reads with unique mapping on this path. Coverage ratio is the length of the read / the len of the sequence of the path
        self.multiple_mapped_abundances = []# for each node of the path (ordered as `self.nodes`), store the coverage of multiple mapped reads (normalized wrt their repartition in other paths)
        self.total_mapped_mult_reads = 0    #  number of reads with corrected multiple mapping on this path. 
        self.total_mapped_mult_reads_normalized = 0
        #----
        # self.hamming_distance   = []        # for each node of the path (ordered as `self.nodes`), store the number of substitutions when mapping reads

class Alignment:
    
    def __init__(self):
        self.mapped_nodes_id = [] # int ids of nodes
        self.mapped_nodes_cov = [] # coverage of nodes
        self.len = 0
        self.score = 0
        self.identity = 0
        self.nb_aligned = 0
        self.nb_match = 0
        self.nb_subst = 0
        self.nb_indel = 0
        self.nb_errors = 0

class Cluster:
    def __init__(self):
        self.colored_paths_id = []
        self.reads_name = [] # to delete after test
        self.unassigned_paths = []
        self.unassigned_count = 0
        self.unassigned_count_norm = 0
        self.unassigned_count_afterfilt = 0
        self.unassigned_abund = {}
        self.estimated_min_strains_before = 0
        self.estimated_min_strains_after = 0
        self.estimated_abund = 0
    
    def min_strains_inference(self):

        unassigned_paths = [canonical_tiny(p) for p in self.unassigned_paths]

        G = nx.Graph()
        list_incomp = {k: [] for k in range(len(unassigned_paths))}
        list_comp = {k: [] for k in range(len(unassigned_paths))}
        for i,read_i in enumerate(unassigned_paths):
            G.add_node(i)
            for j in range(i+1,len(unassigned_paths)):
                res = check_incomp_comp(read_i,unassigned_paths[j])
                if res==-1:
                    list_incomp[i].append(j)
                    list_incomp[j].append(i)
                    G.add_edge(i, j)
                if res==1:
                    list_comp[i].append(j)
                    list_comp[j].append(i)
        len_cliques = list(map(len,list(nx.clique.find_cliques(G))))
        self.estimated_min_strains_before = max(len_cliques) if len_cliques else 0

        unassigned_paths_filt = []
        for i in range(len(unassigned_paths)):
            comp_adjusted = [r for r in list_comp[i] if not sublist(unassigned_paths[i],unassigned_paths[r])]
            if len(comp_adjusted) >= 2:
                # threshold number comp ok
                test_abund = [1] * len(unassigned_paths[i])
                for r in list_comp[i]:
                    for n in unassigned_paths[r]:
                        if n in unassigned_paths[i]: 
                            test_abund[unassigned_paths[i].index(n)] += 1
                a = np.asarray(test_abund)
                peaks = list(scipy.signal.find_peaks(-a)[0]) # detect drops
                if not peaks: 
                    # threshold var comp ok
                    unassigned_paths_filt.append(unassigned_paths[i])
                else: # peaks detected
                    # threshold number var NOT OK
                    G.remove_node(i)
            else:
                # threshold number comp NOT OK
                G.remove_node(i)
        self.unassigned_count_afterfilt = len(unassigned_paths_filt)

        len_cliques = list(map(len,list(nx.clique.find_cliques(G))))
        self.estimated_min_strains_after = max(len_cliques) if len_cliques else 0

        estimated_abund = {}
        for p in unassigned_paths_filt:
            for n in p:
                if n not in estimated_abund: estimated_abund[n] = 0
                estimated_abund[n] += 1
        self.estimated_abund = np.mean([v for k,v in estimated_abund.items()])

class Pangenome:
    def __init__(self):
        """ 
        initializes a graph object 
        """
        self.nodes = {}                     # key: id (unsigned int), value: Node 
        self.paths = []                     # a list of ordered paths (id of a path is its ranks in this list)
        self.paths_name_to_ids = {}         # no choice: each path has a string name, eg gi|1388876906|ref|NZ_CP028116.1|_1000. Two identical paths (eg 684619+,684620+,684618+) may have distinct occurrences and thus names (as comming from distinct genes). Hence one storesfor each path name its unique path id.
        self.paths_content_to_ids = {}      # no choice: each path has a content, eg 684619+,684620+,684618+. This is the key to know its UNIQUE id (int) that is also the rank in the self.paths list
        self.species_names = set()          # store all species ids NZ_CP007592.1, NC_013654.1, ...
        self.mappingscore_freq = {}         # for each species, store the X frequence (eg. X_freq["NZ_CP007592.1"][0.1] = 12 (12 reads mapped with mappingscore between 0 and 0.1))
        self.identityscore_freq = {}        # //
        self.subst_freq = {}                # //
        self.indel_freq = {}                # //
        self.errors_freq = {}               # //
        self.clusters = []                  # a list of ordered Cluster objects (id of a cluster is its ranks in this list)

    def get_sequence_length(self, path):
        return sum([self.nodes[node_id].len_sequence for node_id in path.node_ids])

    def build_pangenome(self, gfa_file_name: str):
        """
        PARSE GFA FILE
        From the gfa file, create a dictionary of the nodes and the paths in the graph.
        S lines contain node ID and its sequence
        P lines contain path ID, list of node ID with orientation and cover of the nodes
        L lines contain links between nodes (we dont care)
        """
        print("Load the pangenome graph")
        path_id = 0
        cpt = 0     # For upgrade process
        with open(gfa_file_name, 'r') as gfa_file:
            size_file = get_file_size(gfa_file)
            while True:
                cpt+=1
                line = gfa_file.readline()

                # end of the file
                if not line:
                    break

                current_seek    = gfa_file.tell() 
                if cpt%10000 == 0: update_progress(current_seek/size_file)
      
                line = line.strip().split('\t')
                # if line S, create node
                if line[0] == 'S':
                    # S       1       ACCACGATTACGCTGGCGCTTA
                    self.nodes[int(line[1])] = Node(line[2]) 
                # if line P, create paths and add paths infos in nodes
                elif line[0] == 'P':
                    # P       gi|1388876906|ref|NZ_CP028116.1|_1000   684619+,684620+,684618+ 187M,187M,1M 
                    # path_id = line[1]

                    str_node_list = canonical(line[2])
                    strain_id = get_accession_number(line[1])
                    # If this path was already seen, we simply add this strain_id to the path.strain_ids
                    if str_node_list in self.paths_content_to_ids:
                        already_seen_path_id = self.paths_content_to_ids[str_node_list]
                        if strain_id not in self.paths[already_seen_path_id].strain_ids:
                            self.paths[already_seen_path_id].strain_ids[strain_id]=0
                        self.paths[already_seen_path_id].strain_ids[strain_id]+=1
                        self.paths_name_to_ids[line[1]] = already_seen_path_id
                        continue # nothing more to do, no incrementation of path_id, as no new path was created
                    # if first time this path is seen
                    self.species_names.add(strain_id)
                    if strain_id not in self.mappingscore_freq: self.mappingscore_freq[strain_id] = dict.fromkeys([i/100 for i in range(0,101,1)],0)
                    if strain_id not in self.identityscore_freq: self.identityscore_freq[strain_id] = dict.fromkeys([i/100 for i in range(0,101,1)],0)
                    if strain_id not in self.subst_freq: self.subst_freq[strain_id] = dict.fromkeys([i/100 for i in range(0,101,1)],0)
                    if strain_id not in self.indel_freq: self.indel_freq[strain_id] = dict.fromkeys([i/100 for i in range(0,101,1)],0)
                    if strain_id not in self.errors_freq: self.errors_freq[strain_id] = dict.fromkeys([i/100 for i in range(0,101,1)],0)
                    path = Path()
                    if strain_id not in path.strain_ids:
                        path.strain_ids[strain_id]=0
                    path.strain_ids[strain_id]+=1
                    node_list = str_node_list.split(',')
                    # seq = ""
                    for node_info in node_list:
                        node_id = int(node_info[:-1])
                        node = self.nodes[node_id]
                        # add the current path id to the traversed paths of this node
                        # we could have used a set for traversed_path, but this requires more memory
                        if path_id not in node.traversed_path:
                            node.traversed_path.append(path_id)

                        path.node_ids.append(node_id)
                        path.unique_mapped_abundances.append(0)
                        path.multiple_mapped_abundances.append(0)
                        # path.hamming_distance.append(0)
                    # assert line[1] not in self.paths_name_to_ids # todo: remove when tested that this path ids has never been seen before

                    self.paths.append(path)                                 # store this new path
                    self.paths_name_to_ids[line[1]] = path_id
                    self.paths_content_to_ids[str_node_list] = path_id
                    path_id+=1
        update_progress(1)

    def fill_clusters_info(self, pickle_file: str):
        """
        PARSE PICKLE FILE
        d_clusters: key = cluster id, value = list of colored paths names
        Add cluster id (int) in colored path object info
        Add colored path id (int) in cluster object info
        """
        d_clusters = pickle.load( open( pickle_file, "rb" ) ) # Cluster_1: ['gi|407479587|ref|NC_...58.1|_3137', 'gi|1447699251|ref|NC...95.2|_1209']
        
        # get all ids and be sure they are sorted by numerical order
        for cluster_id in range(len(d_clusters)):
            clstr = Cluster()
            for path_name in d_clusters["Cluster_"+str(cluster_id)]['genes_list']:
                # skip if the cluster in the pickle has no matching to path in graph
                if path_name not in self.paths_name_to_ids:
                    continue
                path_id = self.paths_name_to_ids[path_name]
                clstr.colored_paths_id.append(path_id)
                self.paths[path_id].cluster_id = cluster_id
            self.clusters.append(clstr)
    
    def fill_abund_dist(self, gene_match: tuple, aligned_path: Alignment, len_read: int):

        path_id,starting_node_id = gene_match
        nb_mapped_nodes = len(aligned_path.mapped_nodes_id)
        path = self.paths[path_id]
        path.total_mapped_unique_reads += 1
        path.total_mapped_unique_reads_normalized += len_read/self.get_sequence_length(path)
        # Le comptage unique "normalisé" par chemin (incrémentation de (longueur du read)/(longueur du chemin))
        for i in range(nb_mapped_nodes):
            path.unique_mapped_abundances[starting_node_id+i]+=aligned_path.mapped_nodes_cov[i]

        # update the number of errors for this path, but for unique strains
        if len(path.strain_ids) == 1:
            strain_id = list(path.strain_ids.keys())[0]
            self.mappingscore_freq[strain_id][round(aligned_path.score,2)] += 1
            self.identityscore_freq[strain_id][round(aligned_path.identity,2)] += 1
            self.subst_freq[strain_id][round(aligned_path.nb_subst,2)] += 1
            self.indel_freq[strain_id][round(aligned_path.nb_indel,2)] += 1
            self.errors_freq[strain_id][round(aligned_path.nb_errors,2)] += 1

    def get_matching_path(self, path_as_node_list):
        """
        INPUT = a node list (from alignement)
        to speed the search we get the starting node of the alignment and search only on paths crossing this start node
        to speed the seach, we use the node position in the selected path and extend the list of nodes by the same length of the alignment
        there is a match if the seed+extend list is identical to the node list of the alignment
        OUPUT = list of paths where the node list match
        """
        
        result = [] # (path_id, corresponding starting node_id, number of nodes mapped)
        start_node = path_as_node_list[0]
        paths_sel = self.nodes[start_node].traversed_path

        # for each path crossing the start node
        for path_id in paths_sel:
            colored_path = self.paths[path_id]
            all_pos_first = [i for i, v in enumerate(colored_path.node_ids) if v == start_node] # start_node from alignment may have several occurrences in the colored path

            for pos_first in all_pos_first:
                if colored_path.node_ids[pos_first:pos_first+len(path_as_node_list)] == path_as_node_list:
                    result.append((path_id, pos_first)) 

        return result
    
    def get_matching_path_old(self, path_as_node_list):
        """
        INPUT = a node list (from alignement)
        to speed the search we get the starting node of the alignment and search only on paths crossing this start node
        to speed the search, we use the node position in the selected path and extend the list of nodes by the same length of the alignment
        there is a match if the seed+extend list is identical to the node list of the alignment
        OUPUT = list of paths where the node list match
        """
        
        result = [] # (path_id, corresponding starting node_id, number of nodes mapped)
        start_node = path_as_node_list[0]
        paths_sel = self.nodes[start_node].traversed_path

        # for each path crossing the start node
        for path_id in paths_sel:
            colored_path = self.paths[path_id]
            all_pos_first = [i for i, v in enumerate(colored_path.node_ids) if v == start_node] # start_node from alignment may have several occurrences in the colored path
            # for each position of the start_node in the colored_path
            for pos_first in all_pos_first:
                current_pos_colored_path = pos_first
                found = True
                for node_id in path_as_node_list:
                    if current_pos_colored_path == len(colored_path.node_ids): 
                        found = False #  if the mapped path is longer than the colored ref path, then no mapping.
                        break 
                    if node_id != colored_path.node_ids[current_pos_colored_path]:
                        found = False
                        break
                    current_pos_colored_path+=1
                if found: 
                    result.append((path_id, pos_first))# Useless, current_pos_colored_path-id_first))     

        return result

    def print_gene_table(self, csv_file_name):
        print(f"Print gene-level results to file {csv_file_name}")
        cvs_file = open(csv_file_name, "w")
        presence_species = {} # id: species name (eg NZ_CP028116.1), value = bool present absent
        for species_name in self.species_names:
            presence_species[species_name] = False
        for species_name in presence_species:
            cvs_file.write(f"{species_name};")
        cvs_file.write(f"cluster;seq_len;nb_uniq_mapped;nb_uniq_mapped_normalized;nb_multimapped;nb_multimapped_normalized;mean_abund_uniq;mean_abund_uniq_nz;mean_abund_multiple;mean_abund_multiple_nz;ratio_covered_nodes\n")
        
        for i,path in enumerate(self.paths):
            if i%1000==0: 
                update_progress(i/len(self.paths))
            #reset species:
            presence_species = dict.fromkeys(presence_species,0)
            for strain_id,nb_occ in path.strain_ids.items():
                # gene id: gi|1388876906|ref|NZ_CP028116.1|_1
                presence_species[strain_id] = nb_occ
            #genes
            for _,nb_occ in presence_species.items():cvs_file.write(f"{nb_occ};")
            
            #nodes_abound_uniq
            #cvs_file.write("{")
            #for i, mapped_abundance in enumerate(path.unique_mapped_abundances):
            #    cvs_file.write(f"{path.node_ids[i]}:{mapped_abundance}")
            #    if i<len(path.unique_mapped_abundances)-1: cvs_file.write(',')
            #cvs_file.write("};")

            #nodes_abund_multiple
            #cvs_file.write("{")
            #for i, mapped_abundance in enumerate(path.multiple_mapped_abundances):
            #    cvs_file.write(f"{path.node_ids[i]}:{mapped_abundance}")
            #    if i<len(path.multiple_mapped_abundances)-1: cvs_file.write(',')
            #cvs_file.write("};")

            #cluster:
            cvs_file.write(f"{path.cluster_id};")

            #seq_len:
            cvs_file.write(f"{self.get_sequence_length(path)};")

            #nb unique mapped
            cvs_file.write(f"{path.total_mapped_unique_reads};")

            #nb unique mapped normalized
            cvs_file.write(f"{path.total_mapped_unique_reads_normalized};")

            #nb corrected multimapped
            cvs_file.write(f"{path.total_mapped_unique_reads+path.total_mapped_mult_reads};")

            #nb corrected multimapped normalized
            cvs_file.write(f"{path.total_mapped_unique_reads_normalized+path.total_mapped_mult_reads_normalized};")

            #mean_abund_uniq
            cvs_file.write(f"{np.mean([mapped_abundance for mapped_abundance in path.unique_mapped_abundances])};")
            
            #mean_abund_uniq_nz
            cvs_file.write(f"{np.mean([mapped_abundance for mapped_abundance in path.unique_mapped_abundances if mapped_abundance > 0])};")

            #mean_abund_multiple
            cvs_file.write(f"{np.mean([mapped_abundance for mapped_abundance in map(sum,zip(path.unique_mapped_abundances,path.multiple_mapped_abundances))])};")
            
            #mean_abund_multiple_nz
            cvs_file.write(f"{np.mean([mapped_abundance for mapped_abundance in map(sum,zip(path.unique_mapped_abundances,path.multiple_mapped_abundances)) if mapped_abundance > 0])};")

            # ratio of covered nodes
            nb_covered_nodes = 0 # for printing the ratio of covered nodes
            for i in range(len(path.unique_mapped_abundances)):
                if path.unique_mapped_abundances[i]>0 or path.multiple_mapped_abundances[i]>0:
                    nb_covered_nodes+=1
            cvs_file.write(f"{nb_covered_nodes/float(len(path.unique_mapped_abundances))}")

            #end path
            cvs_file.write("\n")
        cvs_file.close()
        update_progress(1)

    '''
    #old function
    def print_error_distribution(self, distribution_file_name):
        print(f"Print errors distribution results to file {distribution_file_name}")
        with open(distribution_file_name, "w") as distribution_file:

            for strain_id, dist in self.hamming_freq.items():
                distribution_file.write(f"Strain {strain_id}\n")
                # get the maximal value: (we need this as the `dist` dictionary is not ordered)
                max_err=-1
                for err in dist:
                    if max_err<err: max_err=err
                for nb_err in range(max_err+1):
                    if nb_err in dist: 
                        distribution_file.write(f"{nb_err}: {dist[nb_err]}\n")
                    else: 
                        distribution_file.write(f"{nb_err}: 0\n")
    
    def print_error_distribution2(self, distribution_file_name):
        print(f"Print errors distribution results to file {distribution_file_name}")
        with open(distribution_file_name, "w") as distribution_file:

            # get maximal number of errors
            max_err=-1
            for strain_id, dist in self.hamming_freq.items():
                for err in dist:
                    if max_err<err: max_err=err
            
            # print
            header = "Strain"
            for i in range(max_err+1):
                header += f";{i}"
            distribution_file.write(header+"\n")
            for strain_id, dist in self.hamming_freq.items():
                distribution_file.write(f"{strain_id}")
                for i in range(max_err+1):
                    distribution_file.write(f";{dist.get(i, 0)}")
                distribution_file.write("\n")
    '''
    def print_error_distribution(self, prefix: str, errortype: str):
        d_error = getattr(self,errortype)
        print(f"Print {errortype} distribution results to file {prefix+'_'+errortype+'.csv'}")
        with open(prefix+'_'+errortype+'.csv', "w") as distribution_file:

            header = "Strain;"+';'.join(map(str,[i/100 for i in range(0,101,1)]))
            distribution_file.write(header+"\n")
            for strain_id, dist in d_error.items():
                distribution_file.write(f"{strain_id}")
                for i in [j/100 for j in range(0,101,1)]:
                    distribution_file.write(f";{dist.get(i, 0)}")
                distribution_file.write("\n")

    def print_unassigned_info(self, unassigned_file_name: str):
        print(f"Print cluster-level unassigned reads distribution results to file {unassigned_file_name}")
        with open(unassigned_file_name, "w") as unassigned_file:

            # header
            unassigned_file.write(f"Clusters;")
            presence_species = {} # id: species name (eg NZ_CP028116.1), value = bool present absent
            for species_name in self.species_names:
                presence_species[species_name] = False
            for species_name in presence_species:
                unassigned_file.write(f"{species_name};")
            unassigned_file.write(f"nb_reads_unfilt;min_strains_unfilt;nb_reads_afterfilt;min_strains;mean_abund\n")

            for i,cluster in enumerate(self.clusters):
                unassigned_file.write(f"{'Cluster_'+str(i)};")
                for species_name in presence_species:
                    unassigned_file.write(f"{sum([self.paths[cp].strain_ids.get(species_name,0) for cp in cluster.colored_paths_id])};")
                unassigned_file.write(f"{len(cluster.unassigned_paths)};")
                unassigned_file.write(f"{cluster.estimated_min_strains_before};")
                unassigned_file.write(f"{cluster.unassigned_count_afterfilt};")
                unassigned_file.write(f"{cluster.estimated_min_strains_after};")
                unassigned_file.write(f"{cluster.estimated_abund}\n")

def metaalign2align(pangenome: Pangenome, subpath: dict, align_info: dict):
    """
    CONVERT A "META-ALIGNMENT" INTO AN ACTUAL ALIGNMENT PATH
    Alignments in the mapping file (json) are organized into "meta-nodes" (super-nodes containing the actual nodes of the pangenome graph),
    that corresponds to local alignments.
    This function is used as the final step of the BFS function to convert the best (meta)alignments into actual alignments
    (successive nodes id from the pangenome graph) and to compute node-level abundances and errors/variations.
    """
    alignment = Alignment()
    iter_nodes = [n for metanode in align_info["path"] for n in subpath[ metanode ]['path']['mapping']]
    for node in iter_nodes: # for each node get abundance
        nodeID = int(node['position']['node_id'])
        nb_match = 0
        nb_aligned = 0
        for edit in node["edit"]:
            from_len = int(edit.get("from_length", 0))
            to_len = int(edit.get("to_length", 0))
            alignment.nb_indel += abs(from_len-to_len)
            nb_aligned += min(from_len,to_len)
            nb_match += to_len if from_len == to_len and "sequence" not in edit else 0
        alignment.nb_aligned += nb_aligned
        alignment.nb_match += nb_match
        alignment.nb_subst += nb_aligned - nb_match
        abund = nb_aligned/pangenome.nodes[nodeID].len_sequence
        # if node already exists in the list, is the last item and its coverage is not 1 yet, just add the abundance to it (don't duplicate the node)
        if alignment.mapped_nodes_id and alignment.mapped_nodes_id[-1] == nodeID and alignment.mapped_nodes_cov[-1] < 1:
            alignment.mapped_nodes_cov[-1] += abund
        else:
            alignment.mapped_nodes_id.append(nodeID)
            alignment.mapped_nodes_cov.append(abund)
    alignment.len = align_info["read_len"]
    alignment.score = float(align_info["score"]/align_info["read_len"]) # normalized mapping score
    alignment.identity = float(alignment.nb_match/alignment.nb_aligned)
    alignment.nb_errors = float((alignment.nb_indel+alignment.nb_subst)/align_info["read_len"])
    alignment.nb_subst = float(alignment.nb_subst/align_info["read_len"])
    alignment.nb_indel = float(alignment.nb_indel/align_info["read_len"])
    return alignment                

def BFS(aln, pangenome: Pangenome):
    
    subpath = aln['subpath']
    starts = aln["start"] # list of starting meta-nodes (alignment tree nodes that contains pangenome graph nodes)
    
    # step 1: initialize queue
    queue = [*starts]          # contains a set of nodes not yet traversed. 
    
    current_paths = {}  # in runtime, contains all reached paths. key = last path node id, value = tuple(path, read_len, score). 
                        # This enables to quickly retreive all paths that finish on a given node and to conserve the one(s) that minimize the alignment score 
    
    for start in starts: ## there can be several roots, loop over them
        align_info = {}
        align_info["path"] = [start] # list of subpath nodes describing the path of alignement
        align_info["read_len"] = sum([int(edit.get("to_length", 0)) for node in subpath[start]['path']['mapping'] for edit in node["edit"]]) # maybe the alignment can end with a node which has successor, we need to keep track of the position in the read along the alignment
        align_info["score"] = subpath[start].get("score", 0)                                    # score field is not displayed if = 0
        current_paths[start] = [align_info]
    
    ### WALK THE WHOLE SET OF NODES. 
    ### CONSERVE FOR EACH NODE THE BEST PATH (or PATHS in case of equality) SCORE
    while queue:
        
        current_metanode = queue.pop(0) ## pop first element of the queue 
        
        ### END OF THE PROCESS
        if 'next' not in subpath[current_metanode] or all(i["read_len"]==len(aln['sequence']) for i in current_paths[current_metanode]): ## final paths obtained when there is no child or the length of the read is reached
            continue
        
        ### INSPECT ALL CHILDREN
        for child in subpath[current_metanode]['next']:
            
            ### ADD CHILD NODE IF NOT ALREADY IN THE QUEUE: 
            if child not in queue:                                              
                queue.append(child)
            
            child_info = {}
            ### keeping track of the position in the alignment
            child_info["len"] = sum([int(edit.get("to_length", 0)) for node in subpath[child]['path']['mapping'] for edit in node["edit"]]) # maybe the alignment can end with a node which has successor, we need to keep track of the position in the read along the alignment
            ### get child score (or 0 if attribute absent)
            child_info["score"] = subpath[child].get("score", 0)   
                
            ### UPDATE THE SET OF PATHS 
            # 1/ For all detected paths that end with the current_metanode (retreived thanks to `current_paths`)
            # 2/ add all new paths with their scores
            # 3/ remove updated paths from the dictionary
            # 4/ conserve only the best path(s)
            if child not in current_paths:
                current_paths[child]=[]                                                 # init all paths/scores ending with the current node
            for align_info in current_paths[current_metanode]:                              # add all new paths with their scores
                updated_info = {}
                updated_info["path"] = align_info["path"]+[child]
                updated_info["read_len"] = align_info["read_len"] + child_info["len"]
                updated_info["score"] = align_info["score"] + child_info["score"]
                if len(current_paths[child])==0:
                    current_paths[child].append(updated_info)                           # create a new path that finish with the child node
                    continue                                                            # no need to update the paths ending with child, this is the first. 
                ### conserve only the best path(s)
                best_previous_score = current_paths[child][0]["score"]
                if updated_info["score"] > best_previous_score:                         # The new path is better than all previous scores. Create a new set
                    current_paths[child]=[updated_info]                                 # create a new path that finish with the child node
                    continue
                if updated_info["score"] == best_previous_score:                        # The new path is equal to all previous scores, add the new path
                    current_paths[child].append(updated_info)                           # create a new path that finish with the child node
                    continue
                #if updated_info["score"] < best_previous_score: do nothing             # The newly created path is worse than all previous ending on child. Nothing to be done

        current_paths.pop(current_metanode)                                             # No more path end with the current node                                          
    # need an extra-step to check the best scores if paths don't end with the same node
    if len(current_paths) > 1:
        max_score = max([current_paths[end_node][0]['score'] for end_node in current_paths])
        paths_to_remove = [end_node for end_node in current_paths if current_paths[end_node][0]['score'] < max_score]
        [current_paths.pop(end_node) for end_node in paths_to_remove]
    # convert dict to list
    current_paths = list(itertools.chain.from_iterable(current_paths.values()))
    # convert meta-alignment into alignment
    final_paths = list(map(lambda x: metaalign2align(pangenome,subpath,x), current_paths))
    return final_paths

def get_all_alignments_one_read(json_file, pangenome: Pangenome, thr=0.95):
    """ 
    Singled-end: given a first read line, returns all alignemnts correponding to this read.
    Sometimes a read may occur on several successive lines, hence we concatenate the 
    corresponding alignments

    Paired-end: given a first read line, returns all alignemnts correponding to this read and its pair.
    If "name" and "paired_read_name" (mostly cases of distinct read files as opposed to interleaved)are identical, 
    both reads of the pair need to be distinguished based on their sequence.
    """
    line = json_file.readline()

    # End of file EOF
    if not line:
        return None
    
    aln = json.loads(line)
    starting_read_name = aln['name']
    starting_read_pair = aln.get('paired_read_name','')
    mapped_paths = {} # for each read pair, its sequence and list of mapped_paths
    mapped_paths["r1"] = {}
    mapped_paths["r1"]['name'] = aln['name'] # to delete after test phase
    mapped_paths["r1"]["sequence"] = aln['sequence']
    mapped_paths["r1"]["paths"] = []

    # parse only if the read has mapped
    if "subpath" in aln: 
        # get best path(s) for the alignement
        current_mapped_paths = BFS(aln, pangenome) # needs pangenome only for getting node lengths
        # store the mapped paths if the score is higher or equal to the threshold
        if  current_mapped_paths[0].score >= thr: 
            mapped_paths["r1"]["paths"] = current_mapped_paths

    while True:

        next_line_to_read_when_leaving_this_function = json_file.tell() 
        line = json_file.readline()
        if not line: 
            json_file.seek(next_line_to_read_when_leaving_this_function)
            break

        aln = json.loads(line)
        if aln['name'] != starting_read_name and aln['name'] != starting_read_pair: # next read
            json_file.seek(next_line_to_read_when_leaving_this_function)
            break

        # parse only if the read has mapped
        if "subpath" in aln: 

            # if first time seeing new sequence, paired-end case
            if aln['sequence'] != mapped_paths["r1"]["sequence"] and "r2" not in mapped_paths:
                mapped_paths["r2"] = {}
                mapped_paths["r2"]['name'] = aln['name'] # to delete after test phase
                mapped_paths["r2"]["sequence"] = aln['sequence']
                mapped_paths["r2"]["paths"] = []

            # get best path(s) for the alignement
            current_mapped_paths = BFS(aln, pangenome) # needs pangenome only for getting node lengths
            # final_paths is a list of Alignments [Alignments]
            current_best_score = current_mapped_paths[0].score

            key = "r1" if aln['sequence'] == mapped_paths["r1"]["sequence"] else "r2"
            if len(mapped_paths[key]["paths"]) == 0: # nothing already higher or equal to thr: 
                if  current_best_score >= thr: 
                    mapped_paths[key]["paths"] += current_mapped_paths
            else: # already something higher or equal to thr: 
                best_stored_score = mapped_paths[key]["paths"][0].score
                if best_stored_score < current_best_score:
                    mapped_paths[key]["paths"] = current_mapped_paths # replace the previously stored paths
                if best_stored_score == current_best_score:
                    mapped_paths[key]["paths"] += current_mapped_paths # add the current paths that have the same score
                # if best_stored_score > current_best_score: do nothing, we do not add those fund paths
    return mapped_paths

def parse_vgmpmap(json_file_name:str, pangenome: Pangenome, thr=0.95, force_se=False):

    global NB_READS
    global NB_READS_U
    global NB_READS_M
    global PAIRS
    global SINGLES
    global FILTERED
    global UNASSIGNED
    global ASSIGNED_U
    global ASSIGNED_M
    global INCONSIST
    global INCONSIST_M
    global MULTI_ALIGN
    NB_READS = 0
    NB_READS_U = 0
    NB_READS_M = 0
    PAIRS = 0
    SINGLES = 0
    FILTERED = 0
    global d_filt
    d_filt = {}
    d_filt["o104"] = 0
    d_filt["iai39"] = 0
    d_filt["third"] = 0
    UNASSIGNED = 0
    global d_unassign
    d_unassign = {}
    d_unassign["o104"] = 0
    d_unassign["iai39"] = 0
    d_unassign["third"] = 0
    ASSIGNED_U = 0
    ASSIGNED_M = 0
    INCONSIST = 0
    INCONSIST_M = 0
    MULTI_ALIGN = 0

    global SE
    global SAME_CP
    global DIFF_CP
    SE = 0
    SAME_CP = 0
    DIFF_CP = 0

    global ONE_NOALIGN
    global ONE_NOCP
    global AMBIG_CP_RESOLVED
    global AMBIG_RESOLVED
    global AMBIG_REDUCED
    ONE_NOALIGN = 0
    ONE_NOCP = 0
    AMBIG_CP_RESOLVED = 0
    AMBIG_RESOLVED = 0
    AMBIG_REDUCED = 0

    global SE_M
    SE_M = 0

    global SECOND_PASS
    SECOND_PASS = 0

    """
    PARSE MAPPING JSON FILE
    First check all alignments from the same read
    Ignore read if no alignment score > thr
    Ignore read if multiple alignment score > thr
    
    score = scoring done by vg considering bonus for matches and penalty for mismatches and gap
    identity = Portion of aligned bases that are perfect matches, or 0 if no bases are aligned.
    errors = nb of non aligned bases
    """
    
    # Optimization: we detect positions in the file of reads with unique mapping. Thus they are not tested twice
    recompute_read = set()
    
    # DO TWICE THE JOB: Once for detecting the abundance of unique mapped reads 
    # FIRST PASS/ 
    # For each read that maps uniquely: fill the abundance of 
    #  1/ each node of the mapped path (extremities are increased <= 1 for each mapped read)
    #  2/ store the abundance of each path simply in term of fully mapped reads 
    print("Parsing Alignment: first pass")
    steps = 0
    with open(json_file_name, 'r') as json_file:
        size_file = get_file_size(json_file)
        while True:
            steps += 1
            current_seek = json_file.tell() 
            if steps%1000 == 0: 
                update_progress(current_seek/size_file)

            mapped_paths = get_all_alignments_one_read(json_file, pangenome, thr)

            # end of file
            if mapped_paths == None:
                recompute_read.add(current_seek) # needed to enable the break during second pass
                break
            
            NB_READS += len(mapped_paths)
            # paired-end case
            if not force_se and len(mapped_paths) == 2:
                PAIRS += 1
                
                # 1. no align found for both reads = FILTERED
                if len(mapped_paths["r1"]["paths"]) == 0 and len(mapped_paths["r2"]["paths"]) == 0: 
                    FILTERED += 2

                    ### to delete after test
                    if "o104" in mapped_paths['r1']['name']:
                        d_filt["o104"] += 2
                    elif "iai39" in mapped_paths['r1']['name']:
                        d_filt["iai39"] += 2
                    else:
                        d_filt["third"] += 2
                    ### end

                # 1. one of the read has no align = go to single-end style
                elif len(mapped_paths["r1"]["paths"]) == 0 or len(mapped_paths["r2"]["paths"]) == 0:
                    FILTERED += 1

                    ### to delete after test
                    if "o104" in mapped_paths['r1']['name']:
                        d_filt["o104"] += 1
                    elif "iai39" in mapped_paths['r1']['name']:
                        d_filt["iai39"] += 1
                    else:
                        d_filt["third"] += 1
                    ### end

                    ONE_NOALIGN += 1
                    key_to_remove = "r1" if len(mapped_paths["r1"]["paths"]) == 0 else "r2"
                    if key_to_remove == "r1":
                        mapped_paths["r1"] = mapped_paths.pop("r2")
                    else:
                        del mapped_paths["r2"]
                
                # 1. both reads has alignments
                else:
                    
                    # For each read, parse each alignment path and get corresponding colored path
                    for key in mapped_paths:
                        mapped_paths[key]['colored_match'] = []

                        for i,aligned_path in enumerate(mapped_paths[key]["paths"]):
                            match_paths = pangenome.get_matching_path(aligned_path.mapped_nodes_id) 
                            mapped_paths[key]["colored_match"] += match_paths
                            for m in match_paths:
                                if m[0] not in mapped_paths[key]: mapped_paths[key][m[0]] = []
                                mapped_paths[key][m[0]].append(i)
                            # get_matching_path for reverse path
                            if len(aligned_path.mapped_nodes_id) > 1: # don't duplicate the result if the path has only one node
                                match_paths = pangenome.get_matching_path(aligned_path.mapped_nodes_id[::-1]) 
                                mapped_paths[key]["colored_match"] += match_paths
                                for m in match_paths:
                                    if m[0] not in mapped_paths[key]: mapped_paths[key][m[0]] = []
                                    mapped_paths[key][m[0]].append(i)

                    # 2. no cp found for both reads = UNASSIGNED
                    if len(mapped_paths["r1"]["colored_match"]) == 0 and len(mapped_paths["r2"]["colored_match"]) == 0: 
                        UNASSIGNED += 2

                        ### to delete after test
                        if "o104" in mapped_paths['r1']['name']:
                            d_unassign["o104"] += 2
                        elif "iai39" in mapped_paths['r1']['name']:
                            d_unassign["iai39"] += 2
                        else:
                            d_unassign["third"] += 2
                        ### end

                        # fill clusters with unassigned reads (but alignments with several paths possibilities for the same cluster are not treated)
                        # no shared counts for clusters as it would add noise for the incompatibility check afterwards
                        for key in mapped_paths:
                            mapped_paths[key]["clusters"] = [pangenome.paths[ pangenome.nodes[ p.mapped_nodes_id[0] ].traversed_path[0]].cluster_id for p in mapped_paths[key]["paths"]]
                        intersect_clusters = list(set(mapped_paths["r1"]["clusters"]) & set(mapped_paths["r2"]["clusters"]))
                        if len(intersect_clusters) <= 1: # if len(intersect_clusters) > 1, multiple cluster mapping, do nothing
                            if len(intersect_clusters) == 1: # update list_clstr
                                for key in mapped_paths:
                                    mapped_paths[key]["clusters"] = [x for x in mapped_paths[key]["clusters"] if x in intersect_clusters]
                            for key in mapped_paths:
                                if len(mapped_paths[key]["clusters"]) == 1: # if >1, do nothing + if duplicates, do nothing / if ==1, means there is only one path => OK
                                    traversed_cluster = mapped_paths[key]["clusters"][0]
                                    pangenome.clusters[traversed_cluster].reads_name.append(mapped_paths[key]["name"]) # to delete after test
                                    pangenome.clusters[traversed_cluster].unassigned_paths.append(mapped_paths[key]["paths"][0].mapped_nodes_id)
                                    pangenome.clusters[traversed_cluster].unassigned_count += 1
                                    pangenome.clusters[traversed_cluster].unassigned_count_norm += 1/mapped_paths[key]["paths"][0].len
                                    for i,n in enumerate(mapped_paths[key]["paths"][0].mapped_nodes_id):
                                        if n not in pangenome.clusters[traversed_cluster].unassigned_abund:
                                            pangenome.clusters[traversed_cluster].unassigned_abund[n] = 0
                                        pangenome.clusters[traversed_cluster].unassigned_abund[n] += mapped_paths[key]["paths"][0].mapped_nodes_cov[i]

                    # 2. one of the read has no cp = go to single-end style
                    # No clusters intersection here as nothing can be concluded from it
                    elif len(mapped_paths["r1"]["colored_match"]) == 0 or len(mapped_paths["r2"]["colored_match"]) == 0:
                        UNASSIGNED += 1
                        ONE_NOCP += 1
                        key_to_remove = "r1" if len(mapped_paths["r1"]["colored_match"]) == 0 else "r2"

                        ### to delete after test
                        if "o104" in mapped_paths['r1']['name']:
                            d_unassign["o104"] += 1
                        elif "iai39" in mapped_paths['r1']['name']:
                            d_unassign["iai39"] += 1
                        else:
                            d_unassign["third"] += 1
                        ### end

                        # fill clusters with unassigned reads
                        traversed_clusters = [pangenome.paths[ pangenome.nodes[ p.mapped_nodes_id[0] ].traversed_path[0]].cluster_id for p in mapped_paths[key_to_remove]["paths"]]
                        if len(traversed_clusters) == 1: # if >1, do nothing + if duplicates, do nothing
                            pangenome.clusters[traversed_clusters[0]].reads_name.append(mapped_paths[key_to_remove]["name"]) # to delete after test
                            pangenome.clusters[traversed_clusters[0]].unassigned_paths.append(mapped_paths[key_to_remove]["paths"][0].mapped_nodes_id)
                            pangenome.clusters[traversed_clusters[0]].unassigned_count += 1
                            pangenome.clusters[traversed_clusters[0]].unassigned_count_norm += 1/mapped_paths[key_to_remove]["paths"][0].len
                            for i,n in enumerate(mapped_paths[key_to_remove]["paths"][0].mapped_nodes_id):
                                if n not in pangenome.clusters[traversed_clusters[0]].unassigned_abund:
                                    pangenome.clusters[traversed_clusters[0]].unassigned_abund[n] = 0
                                pangenome.clusters[traversed_clusters[0]].unassigned_abund[n] += mapped_paths[key_to_remove]["paths"][0].mapped_nodes_cov[i]

                        # remove key from unassigned read and go to single-end style
                        if key_to_remove == "r1":
                            mapped_paths["r1"] = mapped_paths.pop("r2")
                        else:
                            del mapped_paths["r2"]
                    
                    # 2. both reads has cp
                    else:
                        
                        # colored paths matches intersection
                        list_path_id1 = [cp[0] for cp in mapped_paths['r1']['colored_match']]
                        list_path_id2 = [cp[0] for cp in mapped_paths['r2']['colored_match']]
                        colored_paths_intersect = list(set(list_path_id1) & set(list_path_id2))
                        
                        # 3. length intersection > 1 = multimapping
                        if len(colored_paths_intersect) > 1:
                            recompute_read.add(current_seek)
                            SECOND_PASS += 2
                        # 3. length intersection ==1 = no (more) ambiguity except if duplicates (= multimapping)
                        elif len(colored_paths_intersect) == 1:
                            count1 = list_path_id1.count(colored_paths_intersect[0])
                            count2 = list_path_id2.count(colored_paths_intersect[0])
                            if count1 > 1:
                                MULTI_ALIGN += 1
                            else:
                                idx_cp = [cp[0] for cp in mapped_paths['r1']["colored_match"]].index(colored_paths_intersect[0])
                                idx_align = mapped_paths['r1'][ colored_paths_intersect[0] ]
                                ### to delete after test
                                if len(idx_align) > 1:
                                    print("stop")
                                ### end
                                ASSIGNED_U += 1
                                SAME_CP += 1
                                pangenome.fill_abund_dist(mapped_paths['r1']["colored_match"][idx_cp], mapped_paths['r1']["paths"][idx_align[0]], len(mapped_paths['r1']["sequence"]))
                                if len(colored_paths_intersect) < len(mapped_paths['r1']["colored_match"]): AMBIG_CP_RESOLVED += 1
                            if count2 > 1:
                                MULTI_ALIGN += 1
                            else:
                                idx_cp = [cp[0] for cp in mapped_paths['r2']["colored_match"]].index(colored_paths_intersect[0])
                                idx_align = mapped_paths['r2'][ colored_paths_intersect[0] ]
                                ### to delete after test
                                if len(idx_align) > 1:
                                    print("stop")
                                ### end
                                ASSIGNED_U += 1
                                SAME_CP += 1
                                pangenome.fill_abund_dist(mapped_paths['r2']["colored_match"][idx_cp], mapped_paths['r2']["paths"][idx_align[0]], len(mapped_paths['r2']["sequence"]))
                                if len(colored_paths_intersect) < len(mapped_paths['r2']["colored_match"]): AMBIG_CP_RESOLVED += 1
                        # 3. no intersection = check strains intersection
                        else:
                            
                            # strains intersection
                            for key in mapped_paths:
                                mapped_paths[key]['corresp_strains'] = []
                                for i,c in enumerate(mapped_paths[key]["colored_match"]):
                                    corresp_strains = list(pangenome.paths[c[0]].strain_ids.keys())
                                    mapped_paths[key]["corresp_strains"] += corresp_strains
                                    for s in corresp_strains:
                                        if s not in mapped_paths[key]: mapped_paths[key][s] = []
                                        mapped_paths[key][s].append(i)
                            strains_intersect = list(set(mapped_paths['r1']['corresp_strains']) & set(mapped_paths['r2']['corresp_strains']))
                            
                            # reduce original list of colored path
                            if len(strains_intersect) >= 1:

                                for key in mapped_paths:
                                    update_cp = []
                                    for s in strains_intersect:
                                        update_cp += [mapped_paths[key]["colored_match"][i] for i in mapped_paths[key][s]]
                                    if len(update_cp) < len(mapped_paths[key]["colored_match"]):
                                        if len(update_cp) == 1:
                                            AMBIG_RESOLVED += 1
                                        else:
                                            AMBIG_REDUCED += 1
                                    mapped_paths[key]["colored_match"] = update_cp

                                    if len(mapped_paths[key]["colored_match"]) > 1: 
                                        recompute_read.add(current_seek)
                                        SECOND_PASS += 1
                                    else: # len(mapped_paths[key]["colored_match"]) == 1
                                        corresp_align = [mapped_paths[key]['paths'][i] for i in mapped_paths[key][ mapped_paths[key]["colored_match"][0][0] ]]
                                        if len(corresp_align) == 1:
                                            ASSIGNED_U += 1
                                            DIFF_CP += 1
                                            pangenome.fill_abund_dist(mapped_paths[key]["colored_match"][0], corresp_align[0], len(mapped_paths[key]["sequence"]))
                                        else:
                                            MULTI_ALIGN += 1
                            else:
                                INCONSIST += 2

            # single-end case
            if len(mapped_paths) == 1 or force_se:

                for key in mapped_paths:
                        
                    SINGLES += 1
                    
                    name = mapped_paths[key]['name'] # to delete after test phase
                    aligned_read = mapped_paths[key]["sequence"]
                    mapped_paths_temp = mapped_paths[key]["paths"]

                    if len(mapped_paths_temp) == 0: 
                        FILTERED += 1
                        continue # no path found

                    if len(mapped_paths_temp) > 1: 
                        recompute_read.add(current_seek)
                        SECOND_PASS += 1
                        continue # Here we deal only with reads mapping exactly one path
                    
                    # if len(mapped_paths_temp) == 1
                    aligned_path = mapped_paths_temp[0]  # for clarity
                    match_paths = pangenome.get_matching_path(aligned_path.mapped_nodes_id) 
                    if len(aligned_path.mapped_nodes_id) > 1: # don't duplicate the result if the path has only one node
                        match_paths += pangenome.get_matching_path(aligned_path.mapped_nodes_id[::-1])
                    
                    if len(match_paths) > 1: # we may have several paths corresponding to a unique alignment
                        recompute_read.add(current_seek)
                        SECOND_PASS += 1
                        continue                    
                    if len(match_paths) == 0: 
                        UNASSIGNED += 1

                        ### to delete after test
                        if "o104" in name:
                            d_unassign["o104"] += 1
                        elif "iai39" in name:
                            d_unassign["iai39"] += 1
                        else:
                            d_unassign["third"] += 1   
                        ### end
                                             
                        # add info to corresponding cluster
                        traversed_cluster = pangenome.paths[ pangenome.nodes[ aligned_path.mapped_nodes_id[0] ].traversed_path[0]].cluster_id
                        pangenome.clusters[traversed_cluster].reads_name.append(name) # to delete after test
                        pangenome.clusters[traversed_cluster].unassigned_paths.append(aligned_path.mapped_nodes_id)
                        pangenome.clusters[traversed_cluster].unassigned_count += 1
                        pangenome.clusters[traversed_cluster].unassigned_count_norm += 1/aligned_path.len
                        for i,n in enumerate(aligned_path.mapped_nodes_id):
                            if n not in pangenome.clusters[traversed_cluster].unassigned_abund:
                                pangenome.clusters[traversed_cluster].unassigned_abund[n] = 0
                            pangenome.clusters[traversed_cluster].unassigned_abund[n] += aligned_path.mapped_nodes_cov[i]
                        
                        continue
                    else:
                        ASSIGNED_U += 1
                        SE += 1
                        continue

                    # if len(match_paths) == 1
                    pangenome.fill_abund_dist(match_paths[0], aligned_path, len(aligned_read))

    update_progress(1)

    # DO TWICE THE JOB: Once for detecting the abundance of unique mapped reads 
    # Once for dealing with multimapped reads
    # SECOND PASS/ 
    # For each read that maps on several paths
    # detect the total_mapped_unique_reads of each of the mapped paths
    # This provides an abundance a,b,c eg for 3 mapped paths respectively A, B, C. 
    # For path 'A', add in each node A.multiple_mapped_abundances[node] a/(a+b+c)
    # For path 'B', add in each node B.multiple_mapped_abundances[node] b/(a+b+c)
    # For path 'C', add in each node C.multiple_mapped_abundances[node] c/(a+b+c)
    
    steps = 0
    print("Parsing Alignment: second pass")
    with open(json_file_name, 'r') as json_file:
        size_file = get_file_size(json_file)
        while True:
            steps += 1
            current_seek = json_file.tell() 
            if steps%1000==0: update_progress(current_seek/size_file)
            if current_seek not in recompute_read: 
                json_file.readline() # dont care
                continue
            mapped_paths = get_all_alignments_one_read(json_file, pangenome, thr)
            
            # end of file
            if mapped_paths == None:
                break

            # paired-end case
            if not force_se and len(mapped_paths) == 2:
                
                # 1. one of the read has no align = go to single-end style
                if len(mapped_paths["r1"]["paths"]) == 0 or len(mapped_paths["r2"]["paths"]) == 0:
                    key_to_remove = "r1" if len(mapped_paths["r1"]["paths"]) == 0 else "r2"
                    if key_to_remove == "r1":
                        mapped_paths["r1"] = mapped_paths.pop("r2")
                    else:
                        del mapped_paths["r2"]
                
                # 1. both reads has alignments
                else:

                    # For each read, parse each alignment path and get corresponding colored path
                    for key in mapped_paths:
                        mapped_paths[key]['colored_match'] = []

                        for i,aligned_path in enumerate(mapped_paths[key]["paths"]):
                            match_paths = pangenome.get_matching_path(aligned_path.mapped_nodes_id) 
                            mapped_paths[key]["colored_match"] += match_paths
                            for m in match_paths:
                                if m[0] not in mapped_paths[key]: mapped_paths[key][m[0]] = []
                                mapped_paths[key][m[0]].append(i)
                            # get_matching_path for reverse path
                            if len(aligned_path.mapped_nodes_id) > 1: # don't duplicate the result if the path has only one node
                                match_paths = pangenome.get_matching_path(aligned_path.mapped_nodes_id[::-1]) 
                                mapped_paths[key]["colored_match"] += match_paths
                                for m in match_paths:
                                    if m[0] not in mapped_paths[key]: mapped_paths[key][m[0]] = []
                                    mapped_paths[key][m[0]].append(i)

                    # 2. one of the read has no cp = go to single-end style
                    if len(mapped_paths["r1"]["colored_match"]) == 0 or len(mapped_paths["r2"]["colored_match"]) == 0:
                        key_to_remove = "r1" if len(mapped_paths["r1"]["colored_match"]) == 0 else "r2"
                        if key_to_remove == "r1":
                            mapped_paths["r1"] = mapped_paths.pop("r2")
                        else:
                            del mapped_paths["r2"]
                    # 2. both reads has cp
                    else:
                    
                        # colored paths matches intersection
                        list_path_id1 = [cp[0] for cp in mapped_paths['r1']['colored_match']]
                        list_path_id2 = [cp[0] for cp in mapped_paths['r2']['colored_match']]
                        colored_paths_intersect = list(set(list_path_id1) & set(list_path_id2))
                        
                        # 3. length intersection >= 1 ( ==1 should not be happening in second pass)
                        if len(colored_paths_intersect) >= 1:

                            for key in mapped_paths:
                                
                                # if all final cp has only one align -> OK, else too complex -> MULTI_ALIGN count
                                if all([ len(mapped_paths[key][ cp ])==1 for cp in colored_paths_intersect ]):
                                    
                                    ASSIGNED_M += 1
                                    # compute a+b+c (cf earlier comments)
                                    sum_covered_paths = sum([pangenome.paths[cp].total_mapped_unique_reads for cp in colored_paths_intersect])
                                    # fill corresponding nodes normalized abundances (a/(a+b+c) cf earlier comments
                                    for cp in mapped_paths[key]['colored_match']:
                                        if cp[0] in colored_paths_intersect:
                                            path_id = cp[0]
                                            starting_node_id = cp[1]
                                            idx_align = mapped_paths[key][ cp[0] ][0]
                                            aligned_path = mapped_paths[key]['paths'][ idx_align ]
                                            nb_mapped_nodes = len(aligned_path.mapped_nodes_id)
                                            path = pangenome.paths[path_id]
                                            if sum_covered_paths == 0: # if no unique mapped reads, equal repartition to the strains
                                                ratio = 1/len(colored_paths_intersect)
                                            else:
                                                ratio = (path.total_mapped_unique_reads)/float(sum_covered_paths)
                                            path.total_mapped_mult_reads += ratio
                                            path.total_mapped_mult_reads_normalized += ratio*len(mapped_paths[key]['sequence'])/pangenome.get_sequence_length(path)
                                            for i in range(nb_mapped_nodes):
                                                path.multiple_mapped_abundances[starting_node_id+i]+=aligned_path.mapped_nodes_cov[i]*ratio
                                else:
                                    MULTI_ALIGN += 1

                        # 3. no intersection = check strains intersection
                        else:
                            
                            # strains intersection
                            for key in mapped_paths:
                                mapped_paths[key]['corresp_strains'] = []
                                for i,c in enumerate(mapped_paths[key]["colored_match"]):
                                    corresp_strains = list(pangenome.paths[c[0]].strain_ids.keys())
                                    mapped_paths[key]["corresp_strains"] += corresp_strains
                                    for s in corresp_strains:
                                        if s not in mapped_paths[key]: mapped_paths[key][s] = []
                                        mapped_paths[key][s].append(i)
                            strains_intersect = list(set(mapped_paths['r1']['corresp_strains']) & set(mapped_paths['r2']['corresp_strains']))
                            
                            # reduce original list of colored path
                            if len(strains_intersect) >= 1:

                                for key in mapped_paths:
                                    update_cp = []
                                    for s in strains_intersect:
                                        update_cp += [mapped_paths[key]["colored_match"][i] for i in mapped_paths[key][s]]
                                    mapped_paths[key]["colored_match"] = update_cp

                                    # if len(mapped_paths[key]["colored_match"]) == 1, read has already been processed in first pass, don't process it twice
                                    if len(mapped_paths[key]["colored_match"]) > 1: 


                                        # if all final cp has only one align -> OK, else too complex -> MULTI_ALIGN count
                                        if all([ len(mapped_paths[key][ cp[0] ])==1 for cp in mapped_paths[key]["colored_match"] ]):
                                            
                                            ASSIGNED_M += 1
                                            # compute a+b+c (cf earlier comments)
                                            sum_covered_paths = sum([pangenome.paths[cp[0]].total_mapped_unique_reads for cp in mapped_paths[key]["colored_match"]])
                                            # fill corresponding nodes normalized abundances (a/(a+b+c) cf earlier comments
                                            for cp in mapped_paths[key]['colored_match']:
                                                path_id = cp[0]
                                                starting_node_id = cp[1]
                                                idx_align = mapped_paths[key][ cp[0] ][0]
                                                aligned_path = mapped_paths[key]['paths'][ idx_align ]
                                                nb_mapped_nodes = len(aligned_path.mapped_nodes_id)
                                                path = pangenome.paths[path_id]
                                                if sum_covered_paths == 0: # if no unique mapped reads, equal repartition to the strains
                                                    ratio = 1/len(mapped_paths[key]["colored_match"])
                                                else:
                                                    ratio = (path.total_mapped_unique_reads)/float(sum_covered_paths)
                                                path.total_mapped_mult_reads += ratio
                                                path.total_mapped_mult_reads_normalized += ratio*len(mapped_paths[key]['sequence'])/pangenome.get_sequence_length(path)
                                                for i in range(nb_mapped_nodes):
                                                    path.multiple_mapped_abundances[starting_node_id+i]+=aligned_path.mapped_nodes_cov[i]*ratio
                                        else:
                                            MULTI_ALIGN += 1
                            else:
                                if len(mapped_paths['r1']["colored_match"]) == 1 or len(mapped_paths['r2']["colored_match"]) == 1:
                                    INCONSIST_M += 1
                                else:
                                    INCONSIST_M += 2

            # single-end case
            if len(mapped_paths) == 1 or force_se:
                
                for key in mapped_paths:
                        
                    name = mapped_paths[key]["name"] # to delete after test phase
                    aligned_read = mapped_paths[key]["sequence"]
                    mapped_paths_temp = mapped_paths[key]["paths"]

                    # we retreive the paths corresponding to this alignments:
                    flag = False
                    for aligned_path in mapped_paths_temp:
                        found_gene_paths = pangenome.get_matching_path(aligned_path.mapped_nodes_id) 
                        if len(aligned_path.mapped_nodes_id) > 1:  # don't duplicate the result if the path has only one node
                            found_gene_paths += pangenome.get_matching_path(aligned_path.mapped_nodes_id[::-1]) 

                        if len(found_gene_paths) > 0:
                            flag=True
                        
                        # compute a+b+c (cf earlier comments)
                        sum_covered_paths = 0
                        for found_gene_path in found_gene_paths:
                            path_id = found_gene_path[0]
                            path = pangenome.paths[path_id]
                            sum_covered_paths += path.total_mapped_unique_reads # TODO: valider avec Kevin ce +1 (en cas de tout à zero)
                        
                            # update the number of mapping errors for this path
                            if len(path.strain_ids) == 1:
                                strain_id = list(path.strain_ids.keys())[0]
                                pangenome.mappingscore_freq[strain_id][round(aligned_path.score,2)] += 1
                                pangenome.identityscore_freq[strain_id][round(aligned_path.identity,2)] += 1
                                pangenome.subst_freq[strain_id][round(aligned_path.nb_subst,2)] += 1
                                pangenome.indel_freq[strain_id][round(aligned_path.nb_indel,2)] += 1
                                pangenome.errors_freq[strain_id][round(aligned_path.nb_errors,2)] += 1


                        # fill corresponding nodes normalized abundances (a/(a+b+c) cf earlier comments
                        for found_gene_path in found_gene_paths:
                            path_id = found_gene_path[0]
                            starting_node_id = found_gene_path[1]
                            nb_mapped_nodes = len(aligned_path.mapped_nodes_id)
                            path = pangenome.paths[path_id]
                            if sum_covered_paths == 0: # if no unique mapped reads, equal repartition to the strains
                                ratio = 1/len(found_gene_paths)
                            else:
                                ratio = (path.total_mapped_unique_reads)/float(sum_covered_paths)
                            path.total_mapped_mult_reads += ratio
                            path.total_mapped_mult_reads_normalized += ratio*len(aligned_read)/pangenome.get_sequence_length(path)
                            for i in range(nb_mapped_nodes):
                                path.multiple_mapped_abundances[starting_node_id+i]+=aligned_path.mapped_nodes_cov[i]*ratio
                            
                    if flag:
                        ASSIGNED_M += 1
                        SE_M += 1
                    else:
                        UNASSIGNED += 1    
                        if "o104" in name:
                            d_unassign["o104"] += 1
                        elif "iai39" in name:
                            d_unassign["iai39"] += 1
                        else:
                            d_unassign["third"] += 1
    update_progress(1)

    for i in range(len(pangenome.clusters)):
        if pangenome.clusters[i].unassigned_count != 0:
            pangenome.clusters[i].min_strains_inference()
    
def usage():
    print(f"Usage: python {sys.argv[0]} -g graph_file_name (gfa) -m mapped_file_name (json) -p dictionary_file_name (pickle) -t alignment_score_threshold -o prefix_output_files_name -f")


if __name__ == "__main__":
    '''
    graph_file = "/home/kdasilva/master_project/project_mock_v3/final_graphs/all_graphs.gfa"
    pickle_file = "/home/kdasilva/master_project/project_mock_v3/final_graphs/dict_clusters.pickle"
    mapping_file = "/home/kdasilva/master_project/project_mock_v3/mapping_mock1a_complete_allgraphs_20201119.json"
    output_file_prefix = "res"
    '''

    '''
    graph_file = "/home/kdasilva/sevens_2021/vg1.35/graphs/all_graphs.gfa"
    pickle_file = "/home/kdasilva/sevens_2021/vg1.35/graphs/dict_clusters.pickle"
    mapping_file = "/home/kdasilva/sevens_2021/vg1.35/mapping/mapping_mixpairedgenesonly_o104_iai39_bl21.json"
    output_file_prefix = "res_ecoli"
    thr = 0.95
    force_se = False
    
    '''
    graph_file = None
    pickle_file = None
    mapping_file = None
    output_file_prefix = None
    
    force_se = False
    
    try:
        opts, _ = getopt.getopt(sys.argv[1:], "hg:p:o:m:t:")
    
    except getopt.GetoptError as err:
        # print help information and exit:
        print(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
        elif o in ("-g"):
            graph_file = a
        elif o in ("-m"):
            mapping_file = a
        elif o in ("-p"):
            pickle_file = a
        elif o in ("-o"):
            output_file_prefix = a
        elif o in ("-t"):
            thr = float(a)
        elif o in ("-f"):
            force_se = True
        
        else:
            assert False, "unhandled option"
    if not graph_file or not pickle_file or not mapping_file: 
        usage()
        exit()

    output_file_csv_name = output_file_prefix+"_gene_table.csv"
    panpan = Pangenome()
    panpan.build_pangenome(graph_file)
    panpan.fill_clusters_info(pickle_file)
    parse_vgmpmap(mapping_file, panpan, thr, force_se)
    panpan.print_gene_table(output_file_csv_name)
    panpan.print_error_distribution(output_file_prefix,"mappingscore_freq")
    panpan.print_error_distribution(output_file_prefix,"identityscore_freq")
    panpan.print_error_distribution(output_file_prefix,"subst_freq")
    panpan.print_error_distribution(output_file_prefix,"indel_freq")
    panpan.print_error_distribution(output_file_prefix,"errors_freq")
    panpan.print_unassigned_info(output_file_prefix+"_unassigned.csv")

    pangenome_clusters = {}
    for i,c in enumerate(panpan.clusters):
        pangenome_clusters["cluster_"+str(i)] = {}
        pangenome_clusters["cluster_"+str(i)]["reads_name"] = panpan.clusters[i].reads_name
        pangenome_clusters["cluster_"+str(i)]["unassigned_paths"] = panpan.clusters[i].unassigned_paths
        pangenome_clusters["cluster_"+str(i)]["unassigned_abund"] = panpan.clusters[i].unassigned_abund.copy()
    pickle.dump( pangenome_clusters, open( f"allclusters.pickle", "wb" ) )

    print(f"Done, csv results are in {output_file_csv_name}, and error distribution are in {output_file_prefix+'_'+'X_freq.csv'}")
    print(f"{NB_READS} reads processed ({PAIRS} pairs, {NB_READS-PAIRS*2} singles).")
    print(f"{FILTERED} reads did not met the user-defined mapping score threshold of {thr}.")
    print(f"{UNASSIGNED} reads could not be assigned to any colored path / gene.")
    print(f"{INCONSIST+INCONSIST_M} reads have been dropped for inconsistancies of strain assignation between two reads of the same pair (paired-end only).")
    print(f"{AMBIG_CP_RESOLVED} colored path and {AMBIG_RESOLVED} strain assignation ambiguities resolved, and {AMBIG_REDUCED} strain assignation ambiguities reduced (paired-end only).")
    print(f"{MULTI_ALIGN} reads shown multiple alignments that could not be resolved.")
    print(f"In total, {FILTERED+UNASSIGNED+INCONSIST+INCONSIST_M+MULTI_ALIGN} reads were not used.")
    print(f"{ASSIGNED_U} reads have been assigned during the first step (unique mapping).")
    print(f"{ASSIGNED_M} reads have been assigned during the second step (multiple mapping).")
    print(f"In total, {ASSIGNED_U+ASSIGNED_M} ({(ASSIGNED_U+ASSIGNED_M)/NB_READS*100}%) reads were used.")
    print(d_filt)
    print(d_unassign)

   