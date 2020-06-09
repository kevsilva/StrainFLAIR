from collections import defaultdict
import logging
import time
import json
import numpy as np
import pandas as pd
import re
import pickle

graph_file = "/home/kdasilva/finalwork_lasso/graphs/graph_concat.gfa"
pangenome = parse_graph(graph_file)

# dico clusters (clusters -> identifiants des variants)
d_clusters = pickle.load( open( "/home/kdasilva/finalwork_lasso/graphs/dict_clusters.pickle", "rb" ) )

# dico variants (variant -> identifiants)
d_variants = {}
for key in pangenome.paths:
    nodelist = '-'.join(pangenome.paths[key].nodes)
    nodelist_r = '-'.join(pangenome.paths[key].nodes[::-1])
    if nodelist in d_variants:
        d_variants[nodelist].append(key)
    elif nodelist_r in d_variants:
        d_variants[nodelist_r].append(key)
    else:
        d_variants[nodelist] = [key]
print(f'Nombre de variants uniques : {len(d_variants)}')

# dico identifiants (identifiant -> cluster + variant)
d_idt = defaultdict(dict)
for i in d_variants:
    for j in d_variants[i]:
        d_idt[j]['variant'] = i
for i in d_clusters:
    for j in d_clusters[i]['genes_list']:
        d_idt[j]['cluster'] = i

# liste genomes
genomes = list(set([i.split("|")[-2] for i in d_idt]))
# template dataframe
df = pd.DataFrame(0,index=list(d_variants.keys()),columns=['cluster']+genomes+["abundance_uniq","abundance_adj","sum_subst"])
df['cluster'] = [ d_idt[d_variants[i][0]]['cluster'] for i in d_variants ]
for i in d_variants:
    for j in d_variants[i]:
        df.loc[i,j.split("|")[-2]] += 1

mapping_file = "/home/kdasilva/finalwork_lasso/mapping/mapping_mix_o104h4_2011c-3493_iai39_o157h7_frik2069.json"
#mapping_file = "/home/kdasilva/finalwork_lasso/mapping/mapping_mix_o104h4_2011c-3493_iai39_o157h7_frik2069_o157h7_frik2455.json"
reads, variants = parse_vgmpmap(mapping_file,thr=0.95)

# comptage des substitutions
subst = {}
for node in variants:
    for pos in variants[node]:
        if "_" not in str(pos):
            if node not in subst:
                subst[node] = 1
            else:
                subst[node] += 1

variants_subst = {}
for node in variants:
    for pos in variants[node]:
        if "_" not in str(pos):
            for nucl in variants[node][pos]:
                if len(variants[node][pos][nucl])>1:
                    if node not in variants_subst:
                        variants_subst[node] = {}
                    if pos not in variants_subst[node]:
                        variants_subst[node][pos] = {}
                    if nucl not in variants_subst[node][pos]:
                        variants_subst[node][pos][nucl] = variants[node][pos][nucl]



d_subst = {}
for i in d_variants:
    if i not in d_subst:
        d_subst[i] = 0
    list_nodes = i.split('-')
    for n in list_nodes:
        if n in subst:
            d_subst[i] += subst[n]/len(pangenome.paths[ d_variants[i][0] ].sequence) # normalisé par la taille du variant
    
    
# Nombre de reads final
print(f'Nombre de reads mappés : {len(reads)}')
# Nombre de reads par souche
d_nbreads = {}
for i in reads:
    key = '_'.join(i.split("_")[:-1])
    if key not in d_nbreads:
        d_nbreads[key] = 1
    else:
        d_nbreads[key] += 1
for i in d_nbreads:
    print(f'Nombre de reads mappés venant de {i} : {d_nbreads[i]} ({d_nbreads[i]/len(reads)*100} %)')

# get abundance for each variant
ucount, acount, stats = get_abundance(reads,True)

# finish dataframe
mat = df.copy()
mat['abundance_uniq'] = ucount.values()
mat['abundance_adj'] = acount.values()
mat['sum_subst'] = d_subst.values()
mat.to_csv(f'/home/kdasilva/finalwork_lasso/results/{mapping_file[39:-5]}.csv')

no_cpaths = pd.DataFrame(0,index=list(stats['no_cpaths_found'].keys()),columns=['nb_reads'])
no_cpaths['nb_reads'] = [len(i) for i in stats['no_cpaths_found'].values()]
no_cpaths.to_csv(f'/home/kdasilva/finalwork_lasso/results/new_paths_{mapping_file[39:-5]}.csv')

d_nbreads_nopaths = {}
for i in stats['no_cpaths_found']:
    for j in stats['no_cpaths_found'][i]:
        key = '_'.join(j[0].split("_")[:-1])
        if key not in d_nbreads_nopaths:
            d_nbreads_nopaths[key] = 1
        else:
            d_nbreads_nopaths[key] += 1
for i in d_nbreads_nopaths:
    print(f'Nombre de nouveaux chemins venant de {i} : {d_nbreads_nopaths[i]}')



d_nbreads_errors = {}
flag = False
for read in reads:
    for pair in reads[read]:
        print(read)
        
        if len(reads[read][pair]) == 1:
            align = reads[read][pair][0]
            nb_errors = len(align.variations)
            align_nodes = [n[0] for n in align.nodes]
            found = list(set(get_matching_path(align_nodes)))
            if len(found)==1:
                if "NC_000913.3" in found[0]:
                    flag = True
                    break
            variants_found = list(set([d_idt[i]['variant'] for i in found]))
        
            if len(variants_found) == 1:
                variant = variants_found[0]
                if variant not in d_nbreads_errors: d_nbreads_errors[variant] = {}
                if nb_errors not in d_nbreads_errors[variant]: d_nbreads_errors[variant][nb_errors] = 0
                d_nbreads_errors[variant][nb_errors] += 1
    if flag:
        break

d_nberrors_total = {}
d_nberrors_uniq = {}
for v in d_nbreads_errors:
    s = set([g.split("|")[3] for g in d_variants[v]])
    for e in s:
        if e not in d_nberrors_total: d_nberrors_total[e] = {}
        if e not in d_nberrors_uniq: d_nberrors_uniq[e] = {}
        for nb in d_nbreads_errors[v]:
            if nb not in d_nberrors_total[e]: d_nberrors_total[e][nb] = 0
            if nb not in d_nberrors_uniq[e]: d_nberrors_uniq[e][nb] = 0
            d_nberrors_total[e][nb] += d_nbreads_errors[v][nb]
            if len(s)==1: d_nberrors_uniq[e][nb] += d_nbreads_errors[v][nb]

strain_table_total = pd.DataFrame(0,index=list(d_nberrors_total.keys()),columns=[x for x in range(0,8)])
for i in d_nberrors_total:
    for j in d_nberrors_total[i]:
        strain_table_total.loc[i,j] = d_nberrors_total[i][j]
strain_table_total.to_csv(f'/home/kdasilva/finalwork_lasso/results/number_errors_total.csv')

strain_table_uniq = pd.DataFrame(0,index=list(d_nberrors_uniq.keys()),columns=[x for x in range(0,8)])
for i in d_nberrors_uniq:
    for j in d_nberrors_uniq[i]:
        strain_table_uniq.loc[i,j] = d_nberrors_uniq[i][j]
strain_table_uniq.to_csv(f'/home/kdasilva/finalwork_lasso/results/number_errors_uniq.csv')

######


class Edge:
    
    def __init__(self, source, target):
        self.source = source
        self.target = target

class Node:
    
    def __init__(self, ID, sequence):
        self.ID = ID
        self.sequence = sequence
        self.pathsCrossed = defaultdict(list)
        self.previous = []
        self.next = []

class Path:
    
    def __init__(self, ID):
        self.ID = ID
        self.nodes = []

class Pangenome:
    
    def __init__(self):
        self.nodes = {}
        self.edges = {}
        self.paths = {}
    
    def addNode(self, node):
        self.nodes[node.ID] = node
    
    def addEdge(self, edge):
        self.edges[(edge.source,edge.target)] = edge
    
    def addPath(self, path):
        self.paths[path.ID] = path

class Alignment:
    
    def __init__(self):
        self.nodes = []
    
    def addLen(self,length):
        self.len = length
    
    def addScore(self,score):
        self.score = score
    
    def addIdentity(self,identity):
        self.identity = identity
    
    def addVariations(self,variations):
        self.variations = variations

###############################################################################

def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join([complement[base] for base in dna[::-1]])

def parse_graph(graph_file):
    
    # PARSE GFA FILE
    # From the gfa file, create a dictionary of the nodes and the paths in the graph.
    # S lines contain node ID and its sequence
    # P lines contain path ID, list of node ID with orientation and cover of the nodes
    # L lines contain links between nodes
    # node dictionary need to contain sequence, paths crossed with orientation, previous nodes and next nodes
    # path dictionary need to contain sequence, node list, offset at the begining and the end
    
    pangenome = Pangenome()
    with open(graph_file, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            # if line S, create node
            if line[0] == 'S':
                pangenome.addNode( Node(ID=line[1],sequence=line[2]) )
            # if line P, create paths and add paths infos in nodes
            elif line[0] == 'P':
                path_id = line[1]
                path = Path(ID=path_id)
                node_list = line[2].split(',')
                seq = ""
                for idx,node_info in enumerate(node_list):
                    ori = node_info[-1]
                    node_id = node_info.rstrip(ori)
                    node = pangenome.nodes[node_id]
                    node.pathsCrossed[path.ID].append((idx,ori))
                    path.nodes.append(node.ID)
                    if ori == "+":
                        seq += node.sequence
                    else:
                        seq += reverse_complement(node.sequence)
                path.sequence = seq
                pangenome.addPath(path)
            # if line L, get parents and children of each node
            elif line[0] == "L":
                node1 = pangenome.nodes[line[1]]
                node2 = pangenome.nodes[line[3]]
                if node2.ID not in node1.next:
                    node1.next.append(node2.ID)
                if node1.ID not in node2.previous:
                    node2.previous.append(node1.ID)
                pangenome.addEdge( Edge(source=node1.ID,target=node2.ID) )
                
    return pangenome

def get_variations(mapping_node):
    for node in mapping_node['path']['mapping']:
        node_id = node['position']['node_id']
        offset = int(node['position'].get('offset',0))
        is_reverse = ('is_reverse' in node['position'])
        pos = len(pangenome.nodes[node_id].sequence)-1-offset if is_reverse else offset
        for edit in node['edit']:
            from_len = int(edit.get("from_length", 0))
            to_len = int(edit.get("to_length", 0))
            # subst or match
            if from_len == to_len: 
                # cas match
                if 'sequence' not in edit:
                    if is_reverse:
                        pos -= from_len
                    else:
                        pos += from_len
                # cas subst
                else:
                    for n in edit['sequence']:
                        yield (node_id, pos, n)
                        if is_reverse:
                            pos -= 1
                        else:
                            pos += 1
            # insert or del
            else:
                # cas deletion
                if 'sequence' not in edit:
                    for n in range(from_len):
                        yield (node_id, pos, ' ')
                        if is_reverse:
                            pos -= 1
                        else:
                            pos += 1
                # cas insertion
                else:
                    if is_reverse:
                        pos_bis = len(edit['sequence'])-1
                    else:
                        pos_bis = 0
                    for n in edit['sequence']:
                        if is_reverse: n = reverse_complement(n)
                        new_pos = str(pos-1)+'-'+str(pos)+'_'+str(pos_bis)
                        yield (node_id, new_pos, n)
                        if is_reverse:
                            pos_bis -= 1
                        else:
                            pos_bis += 1

def BFS(aln):
    
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
            variations = []
            for subpath_node in tupl[0]: ## parse each subpath_node
                for node in subpath[ subpath_node ]['path']['mapping']: # for each node get abundance
                    nodeID = node['position']['node_id']
                    abund = 0
                    for edit in node["edit"]:
                        from_len = int(edit.get("from_length", 0))
                        to_len = int(edit.get("to_length", 0))
                        abund += min(from_len,to_len)/len(pangenome.nodes[nodeID].sequence)
                    # if node already exists in the list, just add the abundance
                    current_node_list = [n[0] for n in alignment.nodes]
                    if nodeID in current_node_list:
                        node_loc = np.where(nodeID == np.array(current_node_list))[0][0]
                        alignment.nodes[node_loc] = (nodeID,alignment.nodes[node_loc][1]+abund)
                    else:
                        alignment.nodes.append((nodeID,abund))
                variations += [var for var in get_variations(subpath[subpath_node])]
            alignment.addVariations(variations)   
            alignment.addLen(tupl[1])
            alignment.addScore(tupl[2])
            alignment.addIdentity(tupl[3]/tupl[4])                 
            final_paths.append(alignment)
    return final_paths

def parse_vgmpmap(json_file, thr=0.9):
    
    # PARSE MAPPING JSON FILE
    # First check all alignments from the same read
    # Ignore read if no alignment score > thr
    # Ignore read if multiple alignment score > thr
    
    # score = scoring done by vg considering bonus for matches and penalty for mismatches and gap
    # identity = Portion of aligned bases that are perfect matches, or 0 if no bases are aligned.
    
    reads = defaultdict(dict) # clé read = id read contenant les différentes séquences (paired-end) contenant les alignements
    variants = defaultdict(dict)

    previous_read = ""
    
    with open(json_file, 'r') as aln_json:
        
        for line in aln_json:
            
            aln = json.loads(line)
            
            # parse only if the read has mapped
            if "subpath" not in aln:
                continue
            
            # get best path(s) for the alignement
            final_paths = BFS(aln)
            
            # if no path satisfying thr, skip, else add the paths
            new_score = final_paths[0].score/len(aln['sequence'])    
            if new_score >= thr:
                if aln['sequence'] not in reads[aln['name']]: 
                    reads[aln['name']][aln['sequence']] = final_paths
                else:               
                    current_score = reads[aln['name']][aln['sequence']][0].score/len(aln['sequence'])
                    # if new best score, replace
                    if new_score > current_score:
                        reads[aln['name']][aln['sequence']] = final_paths
                    # else add
                    else: 
                        if new_score == current_score:
                            reads[aln['name']][aln['sequence']] += final_paths
        
            # if all alignements checked for the current read, process variants
            if aln['name'] != previous_read and previous_read in reads:
                # but only if one alignment for the pair
                if all(len(reads[previous_read][seq]) == 1 for seq in reads[previous_read]):
                    for seq in reads[previous_read]:
                        for tupl in reads[previous_read][seq][0].variations:
                            if tupl[0] not in variants: variants[tupl[0]] = {}
                            if tupl[1] not in variants[tupl[0]]: variants[tupl[0]][tupl[1]] = {}
                            if tupl[2] not in variants[tupl[0]][tupl[1]]: variants[tupl[0]][tupl[1]][tupl[2]] = []
                            variants[tupl[0]][tupl[1]][tupl[2]].append(previous_read)

                
            previous_read = aln['name']
        
        # variants for the last line
        if aln['name'] != previous_read and previous_read in reads:
            if all(len(reads[previous_read][seq]) == 1 for seq in reads[previous_read]):
                for seq in reads[previous_read]:
                    for tupl in reads[previous_read][seq][0].variations:
                        if tupl[0] not in variants: variants[tupl[0]] = {}
                        if tupl[1] not in variants[tupl[0]]: variants[tupl[0]][tupl[1]] = {}
                        if tupl[2] not in variants[tupl[0]][tupl[1]]: variants[tupl[0]][tupl[1]][tupl[2]] = []
                        variants[tupl[0]][tupl[1]][tupl[2]].append(previous_read)
                    
    return reads,variants

def get_matching_path(node_list):
    
    # INPUT = a node list (from alignement)
    # to speed the search we get the starting node of the alignment and search only on paths crossing this start node
    # to speed the seach, we use the node position in the selected path and extend the list of nodes by the same length of the alignment
    # there is a match if the seed+extend list is identical to the node list of the alignment
    # OUPUT = list of paths where the node list match
    
    found = []
    start_node = node_list[0]
    paths_sel = pangenome.nodes[start_node].pathsCrossed
    # for each path crossing the start node
    for path in paths_sel:
        # if the node is crossed more than one time by the same path
        for idx in pangenome.nodes[start_node].pathsCrossed[path]:
            # check identity of the nodes list (not forgetting the reverse complement case)
            if pangenome.paths[path].nodes[idx[0]:idx[0]+len(node_list)] == node_list or pangenome.paths[path].nodes[idx[0]-len(node_list)+1:idx[0]+1] == node_list[::-1]:
                found.append(path)
    if len(found) == 0:
        print('no match found for 1 alignment!\n')
    
    return found

def get_abundance(reads, multimap=False):
    
    # unique count initialization
    ucount = {}
    for g in d_variants:
        ucount[g] = 0
    
    # stats initialization
    stats = {}
    stats['repartition'] = {} # TODO
    stats['no_cpaths_found'] = defaultdict(list) # key=cluster, value=(read,pair)
    
    multimapped = {}
    
    for read in reads:
        for pair in reads[read]:
            print(read)
            
            found = []
            for align in reads[read][pair]:
                align_nodes = [n[0] for n in align.nodes]
                found += get_matching_path(align_nodes)
            found = list(set(found))
            variants_found = list(set([d_idt[i]['variant'] for i in found]))
            
            if len(variants_found) == 0:
                stats['no_cpaths_found'][d_idt[list(pangenome.nodes[align_nodes[0]].pathsCrossed.keys())[0]]['cluster']].append( (read,pair) )
            elif len(variants_found) > 1:
                multimapped[(read,pair)] = variants_found
            else:
                variant = variants_found[0]
                path_seq = ''.join([pangenome.nodes[n].sequence for n in pangenome.paths[ d_variants[variant][0] ].nodes])
                ucount[variant] += len(pair)/len(path_seq) # A REVOIR EN CAS DE MAPPING IMPARFAIT
    stats['multimapped_reads'] = list(multimapped.keys())
    # ajusted count
    if multimap:
        
        acount = ucount.copy()
        
        for read,pair in multimapped.keys():
            
            current_sum = sum([ucount[c] if c in ucount else 0 for c in multimapped[(read,pair)]])

            for variant in multimapped[(read,pair)]:
                path_seq = ''.join([pangenome.nodes[n].sequence for n in pangenome.paths[ d_variants[variant][0] ].nodes])
                # if variant with no unique count, equal sharing, else ajusted count
                if current_sum == 0:
                    acount[variant] = (1/len(multimapped[(read,pair)]))*(len(pair)/len(path_seq))
                else:
                    acount[variant] += (ucount[variant]/current_sum)*(len(pair)/len(path_seq))
    
        return ucount, acount, stats
    
    else:
        return ucount, stats
