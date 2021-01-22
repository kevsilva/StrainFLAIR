#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import pickle # for reading the dictionnary
import numpy as np
import json
import sys
import getopt
import os # for size files 
import re

# update_progress() : Displays or updates a console progress bar
## Accepts a float between 0 and 1. Any int will be converted to a float.
## A value under 0 represents a 'halt'.
## A value at 1 or bigger represents 100%
#https://stackoverflow.com/questions/3160699/python-progress-bar 
def update_progress(progress):
    barLength = 50 # Modify this to change the length of the progress bar
    status = ""
    # if isinstance(progress, int):
    #     progress = float(progress)
    # if not isinstance(progress, float):
    #     progress = 0
    #     status = "error: progress var must be float\r\n"
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
    
def file_size(f):
    old_file_position = f.tell()
    f.seek(0, os.SEEK_END)
    size = f.tell()
    f.seek(old_file_position, os.SEEK_SET)
    return size

reverse_sign = lambda x: '-' if (x=='+') else '+'

def canonical(node_list: str):
    """
    Returns the canonical representation of a node list. 
    We define the canonical as the min between a list eg: 684619+,684620+,684618+ and is reverse eg: 684618-,684620-,684619-
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

def get_accession_number(line: str):
    accession_number = None
    line = '_'.join(line.split("_")[:-1])
    line = line.split("|")
    for element in line:
        if re.match("[a-zA-Z0-9]+_[a-zA-Z0-9]+", element):
            accession_number = element
            break
    return accession_number

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
        # self.hamming_distance   = []        # for each node of the path (ordered as `self.nodes`), store the number of substitutions when mapping reads
        self.multiple_mapped_abundances = []# for each node of the path (ordered as `self.nodes`), store the coverage of multiple mapped reads (normalized wrt their repartition in other paths)
        self.total_mapped_mult_reads = 0    #  number of reads with corrected multiple mapping on this path. 
        self.total_mapped_mult_reads_normalized = 0
        

class Pangenome:
    def __init__(self):
        """ 
        initializes a graph object 
        """
        self.nodes = {}                     # key: id (unsigned int), value: node 
        self.paths = []                     # a list of ordered paths (id of a path is its ranks in this list)
        self.paths_name_to_ids = {}         # no choice: each path has a string name, eg gi|1388876906|ref|NZ_CP028116.1|_1000. Two identical paths (eg 684619+,684620+,684618+) may have distinct occurrences and thus names (as comming from distinct genes). Hence one storesfor each path name its unique path id.
        self.paths_content_to_ids = {}      # no choice: each path has a content, eg 684619+,684620+,684618+. This is the key to know its UNIQUE id (int) that is also the rank in the self.paths list
        self.species_names = set()          # store all species ids NZ_CP007592.1, NC_013654.1, ...
        self.hamming_freq = {}              # for each species, store the hamming frequence (eg. hamming_freq["NZ_CP007592.1"][3] = 12 (12 reads mapped with 3 substitutions))

    
    def get_sequence_length(self, path):
        sum_seq_len = 0
        for node_id in path.node_ids:
            sum_seq_len+=self.nodes[node_id].len_sequence
        return sum_seq_len

    def fill_pangenome(self, gfa_file_name: str):
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
            size_file = file_size(gfa_file)
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
                    if strain_id not in self.hamming_freq: 
                        self.hamming_freq[strain_id] = {}
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

    def fill_cluster_id_for_each_path(self, pickle_file):
        """
        d_clusters: key = cluster id, value = list of colored paths id
        """
        d_clusters = pickle.load( open( pickle_file, "rb" ) ) 
        # Cluster_1: ['gi|407479587|ref|NC_...58.1|_3137', 'gi|1447699251|ref|NC...95.2|_1209']

        for cluster_id, cluster in d_clusters.items():
            for path_name in cluster['genes_list']:
                # skip if the cluster in the pickle has no matching to path in graph
                if path_name not in self.paths_name_to_ids:
                    continue
                path_id = self.paths_name_to_ids[path_name]
                # assert not self.paths[path_id].cluster_id   # should not already be defined. A path belongs to a unique cluser
                self.paths[path_id].cluster_id = int(cluster_id.split("_")[-1])

    def get_matching_path2(self, path_as_node_list):
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
            id_first = [i for i, v in enumerate(colored_path.node_ids) if v == start_node] # start_node from alignment may have several occurrences in the colored path

            for id_node in id_first:
                if colored_path.node_ids[id_node:id_node+len(path_as_node_list)] == path_as_node_list:
                    result.append((path_id, id_node)) 

        return result
    
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
            all_id_first = [i for i, v in enumerate(colored_path.node_ids) if v == start_node] # start_node from alignment may have several occurrences in the colored path
            for id_first in all_id_first:
                current_pos_colored_path = id_first
                found = True
                for id_node in path_as_node_list:
                    if current_pos_colored_path == len(colored_path.node_ids): 
                        found = False #  if the mapped path is longer than the colored ref path, then no mapping.
                    
                        break 
                    if id_node != colored_path.node_ids[current_pos_colored_path]:
                        found = False
                        break
                    current_pos_colored_path+=1
                if found: 
                    result.append((path_id, id_first))# Useless, current_pos_colored_path-id_first))     

        return result

    def print_error_distribution(self, distribution_file_name):
        print(f"Print distribution results to file {distribution_file_name}")
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
                

    def print_to_csv(self, csv_file_name):
        print(f"Print results to file {csv_file_name}")
        cvs_file = open(csv_file_name, "w")
        presence_species = {} # id: species name (eg NZ_CP028116.1), value = bool present absent
        for species_name in self.species_names:
            presence_species[species_name] = False
        for species_name in presence_species:
            cvs_file.write(f"{species_name};")
        cvs_file.write(f"hamming;cluster;seq_len;nb_uniq_mapped;nb_uniq_mapped_normalized;nb_multimapped;nb_multimapped_normalized;mean_abund_uniq;mean_abund_uniq_nz;mean_abund_multiple;mean_abund_multiple_nz;ratio_covered_nodes\n")
        
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

            #hamming
            cvs_file.write("0;") # not implemented yet

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
        

            
                




class Alignment:
    
    def __init__(self):
        ## TODO: two arrays: one for node_ids, one for their coverage
        self.mapped_node_ids_cov = [] # (int ids of nodes,  of the mapped node)
    
    def addLen(self,length):
        self.len = length
    
    def addScore(self,score):
        self.score = score
    
    def addIdentity(self,identity):
        self.identity = identity
    
    def store_errors(self, nb_errors):
        self.nb_errors = nb_errors
    
    # def addVariants(self,variants):
    #     self.variants = variants

def BFS(aln, pangenome: Pangenome):
    subpath = aln['subpath']
    # step 1: initialize queue
    queue = []          # contains a set of nodes not yet traversed. 
    
    current_paths = {}  # in runtime, contains all reached paths. key = last path node id, value = tuple(path, read_len, score). 
                        # This enables to quickly retreive all paths that finish on a given node and to conserve the one(s) that minimize the alignment score 
    
    for start in range(len(subpath)): ## there can be several roots, loop over them
        try: ## for some reason aln['start'] can contain errors, "try" needed
            read_len = 0
            nb_match = 0
            nb_aligned = 0
            ## maybe the alignment can end with a node which has successor, we need to keep track of the position in the read along the alignment
            for node in subpath[start]['path']['mapping']:
                for edit in node["edit"]:
                    from_len = int(edit.get("from_length", 0))
                    to_len = int(edit.get("to_length", 0))
                    read_len += to_len
                    nb_match += to_len if from_len == to_len and "sequence" not in edit else 0
                    nb_aligned += min(from_len,to_len)
            ## end
            score = subpath[start].get("score", 0)                                    # score field is not displayed if = 0
            queue.append(start)                                                        # the queue contains only the current node 
            if start not in current_paths:                                             # init all paths/scores ending with the current node
                current_paths[start] = []                                    
            current_paths[start].append(([start],read_len,score,nb_match,nb_aligned)) # stores all paths that end with the start node.
        except IndexError: 
            continue
            
    ### WALK THE WHOLE SET OF NODES. 
    ### CONSERVE FOR EACH NODE THE BEST PATH (or PATHS in case of equality) SCORE
    while queue:

        current_node = queue.pop(0) ## pop first element of the queue 
        
        ### END OF THE PROCESS
        if 'next' not in subpath[current_node] or all(i[1]==len(aln['sequence']) for i in current_paths[current_node]): ## final paths obtained when there is no child or the length of the read is reached
            continue
        
        ### INSPECT ALL CHILDREN
        for child in subpath[current_node]['next']:
            
            ### ADD CHILD NODE IF NOT ALREADY IN THE QUEUE: 
            if child not in queue:                                              
                queue.append(child)
                
            read_len = 0
            nb_match = 0
            nb_aligned = 0
            ### keeping track of the position in the alignment
            for node in subpath[child]['path']['mapping']:
                for edit in node["edit"]:
                    from_len = int(edit.get("from_length", 0))
                    to_len = int(edit.get("to_length", 0))
                    read_len += to_len
                    nb_match += to_len if from_len == to_len and "sequence" not in edit else 0
                    nb_aligned += min(from_len,to_len)

            ### get child score (or 0 if attribute absent)
            child_score = subpath[child].get("score", 0)
                
            ### UPDATE THE SET OF PATHS 
            # 1/ For all detected paths that end with the current_node (retreived thanks to `current_paths`)
            # 2/ add all new paths with their scores
            # 3/ remove updated paths from the dictionary
            # 4/ conserve only the best path(s)
            if child not in current_paths:
                current_paths[child]=[]                                                 # init all paths/scores ending with the current node
            for (current_path, current_len, current_score, current_match, current_aligned) in current_paths[current_node]:                 # add all new paths with their scores
                updated_score = current_score + child_score
                updated_match = current_match + nb_match
                updated_aligned = current_aligned + nb_aligned                                                             # shortcut
                if len(current_paths[child])==0:
                    current_paths[child].append((current_path + [child], current_len + read_len, updated_score, updated_match, updated_aligned))        # create a new path that finish with the child node
                    continue                                                                                            # no need to update the paths ending with child, this is the first. 
                ### conserve only the best path(s)
                best_previous_score = current_paths[child][0][2]
                if updated_score > best_previous_score:                                                                 # The new path is better than all previous scores. Create a new set
                    current_paths[child]=[]                                                                             # for all paths ending on this child node
                    current_paths[child].append((current_path + [child], current_len + read_len, updated_score, updated_match, updated_aligned))        # create a new path that finish with the child node
                if updated_score == best_previous_score:                                                                # The new path is equal to all previous scores, add the new path
                    current_paths[child].append((current_path + [child], current_len + read_len, updated_score, updated_match, updated_aligned))        # create a new path that finish with the child node
                #if updated_score < best_previous_score: do nothing                                                     # The newly created path is worse than all previous ending on
                                                                                                                        # child. Nothing to be done

        current_paths.pop(current_node)      
                                                                           # No more path end with the current node
    final_paths = []
    for idx in current_paths: ## for each path among the current paths
        for tupl in current_paths[idx]: ## for each tuple we need to get back the list of nodes
            alignment = Alignment()
            for subpath_node in tupl[0]: ## parse each subpath_node
                for node in subpath[ subpath_node ]['path']['mapping']: # for each node get abundance
                    nodeID = int(node['position']['node_id'])
                    abund = 0
                    for edit in node["edit"]:
                        from_len = int(edit.get("from_length", 0))
                        to_len = int(edit.get("to_length", 0))
                        abund += min(from_len,to_len)/pangenome.nodes[nodeID].len_sequence
                    # if node already exists in the list, just add the abundance
                    current_node_list = [n[0] for n in alignment.mapped_node_ids_cov]
                    if nodeID in current_node_list:
                        node_loc = np.where(nodeID == np.array(current_node_list))[0][0]
                        alignment.mapped_node_ids_cov[node_loc] = (nodeID,alignment.mapped_node_ids_cov[node_loc][1]+abund)
                    else:
                        alignment.mapped_node_ids_cov.append((nodeID,abund))
            alignment.addLen(tupl[1])
            alignment.addScore(tupl[2])
            alignment.addIdentity(tupl[3]/tupl[4]) 
            alignment.store_errors(tupl[4]-tupl[3])                
            final_paths.append(alignment)
    return final_paths

def get_all_alignments_one_read(json_file, pangenome: Pangenome, thr=0.95):
    """ 
    given a first read line, returns all alignemnts correponding to this read.
    Sometimes a read may occur on several successive lines, hence we concatenate the 
    corresponding alignments
    """
    line = json_file.readline()

    # End of file EOF
    if not line:
        return None, None 
    
    aln = json.loads(line)
    starting_read = aln['sequence']
    mapped_paths = []
    # parse only if the read has mapped
    if "subpath" in aln: 
        # get best path(s) for the alignement
        current_mapped_paths = BFS(aln, pangenome) # needs pangenome only for getting node lengths
        # final_paths is a list of Alignments [Alignments]

        # store the mapped paths if the score is higher or equal to the threshold
        if  current_mapped_paths[0].score/len(starting_read) >= thr: 
            mapped_paths = current_mapped_paths

    while True:
        next_line_to_read_when_leaving_this_function = json_file.tell() 
        line = json_file.readline()
        if not line: 
            json_file.seek(next_line_to_read_when_leaving_this_function)
            break
        aln = json.loads(line)
        if aln['sequence'] != starting_read:
            json_file.seek(next_line_to_read_when_leaving_this_function)
            break
        # parse only if the read has mapped
        if "subpath" in aln: 
            # get best path(s) for the alignement
            current_mapped_paths = BFS(aln, pangenome) # needs pangenome only for getting node lengths
            # final_paths is a list of Alignments [Alignments]
            current_best_score = current_mapped_paths[0].score/len(aln['sequence'])
            if len(mapped_paths) == 0: # nothing already higher or equal to thr: 
                if  current_best_score >= thr: 
                    mapped_paths += current_mapped_paths
            else: # already something higher or equal to thr: 
                best_stored_score = mapped_paths[0].score/len(aln['sequence'])
                if best_stored_score < current_best_score:
                    mapped_paths = current_mapped_paths # replace the previously stored paths
                if best_stored_score == current_best_score:
                    mapped_paths += current_mapped_paths # add the current paths that have the same score
                # if best_stored_score > current_best_score: do nothing, we do not add those ofund paths
    return mapped_paths, starting_read



def parse_vgmpmap(json_file_name:str, pangenome: Pangenome, thr=0.95):
    
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
    do_not_recompute_line = set()
    
    # DO TWICE THE JOB: Once for detecting the abundance of unique mapped reads 
    # FIRST PASS/ 
    # For each read that maps uniquely: fill the abundance of 
    #  1/ each node of the mapped path (extremities are increased <= 1 for each mapped read)
    #  2/ store the abundance of each path simply in term of fully mapped reads 
    print("Parsing Alignment: first pass")
    steps = 0
    with open(json_file_name, 'r') as json_file:
        size_file = file_size(json_file)
        while True:
            steps += 1
            current_seek = json_file.tell() 
            if steps%1000 == 0: 
                update_progress(current_seek/size_file)
            mapped_paths, aligned_read = get_all_alignments_one_read(json_file, pangenome, thr)
            
            # end of file
            if mapped_paths == None:
                break

            if len(mapped_paths) == 0: 
                do_not_recompute_line.add(current_seek)
                continue # no path found

            if len(mapped_paths) > 1: 
                continue # Here we deal only with reads mapping exactly one path

            aligned_path = mapped_paths[0]  # for clarity
            aligned_path_as_nodes = [n[0] for n in aligned_path.mapped_node_ids_cov]
            found_gene_paths = pangenome.get_matching_path(aligned_path_as_nodes) 
            # don't duplicate the result if the path has only one node
            if len(aligned_path_as_nodes) > 1:
                found_gene_paths += pangenome.get_matching_path(aligned_path_as_nodes[::-1])
            # we may have several paths corresponding to a unique alignment
            if len(found_gene_paths) > 1: 
                continue
            
            do_not_recompute_line.add(current_seek) # we will not recompute those alignments during the second pass
            for found_gene_path in found_gene_paths:
                path_id = found_gene_path[0]
                starting_node_id = found_gene_path[1]
                nb_mapped_nodes = len(aligned_path.mapped_node_ids_cov)
                path = pangenome.paths[path_id]
                path.total_mapped_unique_reads += 1
                path.total_mapped_unique_reads_normalized += len(aligned_read)/pangenome.get_sequence_length(path)
                # Le comptage unique "normalisé" par chemin (incrémentation de (longueur du read)/(longueur du chemin))
                for i in range(nb_mapped_nodes):
                    path.unique_mapped_abundances[starting_node_id+i]+=aligned_path.mapped_node_ids_cov[i][1] 

                # update the number of mapping errors for this path
                for strain_id in path.strain_ids: 
                    #if strain_id not in pangenome.hamming_freq: pangenome.hamming_freq[strain_id] = {}
                    if aligned_path.nb_errors not in pangenome.hamming_freq[strain_id]:
                        pangenome.hamming_freq[strain_id][aligned_path.nb_errors] = 0 
                    pangenome.hamming_freq[strain_id][aligned_path.nb_errors]+=1
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
        size_file = file_size(json_file)
        while True:
            steps += 1
            current_seek = json_file.tell() 
            if steps%1000==0: update_progress(current_seek/size_file)
            if current_seek in do_not_recompute_line: 
                json_file.readline() # dont care
                continue
            mapped_paths, aligned_read = get_all_alignments_one_read(json_file, pangenome, thr)
            
            # end of file
            if mapped_paths == None:
                break

            # we retreive the paths corresponding to this alignments:
            for aligned_path in mapped_paths:
                aligned_path_as_nodes = [n[0] for n in aligned_path.mapped_node_ids_cov]
                found_gene_paths = pangenome.get_matching_path(aligned_path_as_nodes) 
                # don't duplicate the result if the path has only one node
                if len(aligned_path_as_nodes) > 1:
                    found_gene_paths += pangenome.get_matching_path(aligned_path_as_nodes[::-1]) 
                
                # compute a+b+c (cf earlier comments)
                sum_covered_paths = 0
                for found_gene_path in found_gene_paths:
                    path_id = found_gene_path[0]
                    path = pangenome.paths[path_id]
                    sum_covered_paths += path.total_mapped_unique_reads # TODO: valider avec Kevin ce +1 (en cas de tout à zero)
                
                    # update the number of mapping errors for this path
                    for strain_id in path.strain_ids: 
                        #if strain_id not in pangenome.hamming_freq: pangenome.hamming_freq[strain_id] = {}
                        if aligned_path.nb_errors not in pangenome.hamming_freq[strain_id]:
                            pangenome.hamming_freq[strain_id][aligned_path.nb_errors] = 0
                        pangenome.hamming_freq[strain_id][aligned_path.nb_errors]+=1


                # fill corresponding nodes normalized abundances (a/(a+b+c) cf earlier comments
                for found_gene_path in found_gene_paths:
                    path_id = found_gene_path[0]
                    starting_node_id = found_gene_path[1]
                    nb_mapped_nodes = len(aligned_path.mapped_node_ids_cov)
                    path = pangenome.paths[path_id]
                    if sum_covered_paths == 0: # if no unique mapped reads, equal repartition to the strains
                        ratio = 1/len(found_gene_paths)
                    else:
                        ratio = (path.total_mapped_unique_reads)/float(sum_covered_paths)# TODO: valider avec Kevin ce +1 (en cas de tout à zero)
                    path.total_mapped_mult_reads += ratio # TODO: valider avec Kevin
                    path.total_mapped_mult_reads_normalized += ratio*len(aligned_read)/pangenome.get_sequence_length(path)
                    for i in range(nb_mapped_nodes):
                        path.multiple_mapped_abundances[starting_node_id+i]+=aligned_path.mapped_node_ids_cov[i][1]*ratio
                    
                
                
    update_progress(1)




def usage():
    print(f"Usage: python {sys.argv[0]} -g graph_file_name (gfa) -m mapped_file_name (json) -p dictionary_file_name (pickle) -t alignment_score_threshold -o prefix_output_files_name")

    

#if __name__ == "__main__":
def json2csv_main():
    
    graph_file = "/home/kdasilva/master_project/project_mock_v3/final_graphs/all_graphs.gfa" #"final_graph.gfa"
    pickle_file = "/home/kdasilva/master_project/project_mock_v3/final_graphs/dict_clusters.pickle" #"dict_clusters.pickle"
    mapping_file = "/home/kdasilva/master_project/project_mock_v3/mapping_mock1a_complete_allgraphs_20201119.json" #"mapping_mix_o104_iai39_k12_100k.json"
    output_file_prefix = "res"
    thr = 0.95
    
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
            thr = a
        
        else:
            assert False, "unhandled option"
    if not graph_file or not pickle_file or not mapping_file: 
        usage()
        exit()


    dist_err_file_name = output_file_prefix+"_dist_err.txt"
    output_file_csv_name = output_file_prefix+".csv"
    panpan = Pangenome()
    panpan.fill_pangenome(graph_file)
    panpan.fill_cluster_id_for_each_path(pickle_file)
    parse_vgmpmap(mapping_file, panpan, thr)
    panpan.print_to_csv(output_file_csv_name)
    panpan.print_error_distribution(dist_err_file_name)

    print(f"Done, csv results are in {output_file_csv_name}, and error distribution are in {dist_err_file_name}")

            

   