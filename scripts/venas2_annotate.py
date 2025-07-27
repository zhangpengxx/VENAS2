#!/usr/bin/env python
# coding=utf-8
# __author__ = 'zp'

import pandas as pd
import numpy as np
from datetime import datetime
import time
import argparse
import os
import random
from collections import Counter

import networkx as nx
import cdlib
from networkx.algorithms import community
from cdlib import algorithms, viz, NodeClustering

beg_time = time.time()
# parser = argparse.ArgumentParser()
# parser.add_argument('-input_name', '--input_name', dest='input_name', help="the name of input samples mutation file", required=True)
# # # parser.add_argument('-output_path', '--output_path', dest='output_path', help="the output path of the results", default= 'venas2_output')
# args = parser.parse_args()

# # #parameters
# samplefile = args.input_name
data_path = 'venas2_output/'
haplotype_df = pd.read_csv(data_path+'haplotype_final.txt',sep='\t', header=0)
# edges_df = pd.read_csv(data_path + 'edges.csv', sep=',', header=0)
#if you run venas2_optimize.py please use command below
edges_df = pd.read_csv(data_path + 'edges_optimized.csv', sep=',', header=0)
edges_df.columns = ['source', 'target', 'mutations', 'distance']

G_dir = nx.DiGraph()
G_dir.add_edges_from(edges_df[['source', 'target']].values.tolist())

G = nx.from_pandas_edgelist(edges_df,source="source",target="target")
coms = cdlib.algorithms.louvain(G, resolution=2.0,randomize=True)
n_com = len(coms.communities)

"""identify the main path"""
#find key node in each cluster
dict_deg = dict(G.degree())
keynodes_list = []
dict_clusterstokeynode = {}
dict_clusterstonodes = {}
dict_nodestoclusters = {}
dict_clusterslen = {}
clusterid = 0
for com in coms.communities:
    score=0
    keynode=0
    for node in com:
        dict_nodestoclusters[node] = clusterid 
        if score < dict_deg[node]:
            score = dict_deg[node]
            keynode = node
    # dict_keynodes[keynode] = score
    keynodes_list.append(keynode)
    dict_clusterstokeynode[clusterid] = keynode
    dict_clusterslen[clusterid]=len(com)
    dict_clusterstonodes[clusterid]= com
    clusterid += 1

#filter keyNode use threshold or cluster count
G_smaller=G.copy()
removeNode=[]
for node in G.nodes():
    if dict_deg[node] ==1:
        removeNode.append(node)

for node in removeNode:
    G_smaller.remove_node(node)

filter_keynodes = []
small_keynodes = []
threshold = len(G_smaller.nodes())/float(n_com) # use this ratio as deault
# threshold = 1100
print('threshold is ', threshold)
for keynode in keynodes_list:
    clusterid = dict_nodestoclusters[keynode]
    clusterlen = dict_clusterslen[clusterid]
    if clusterlen > threshold:
        filter_keynodes.append(keynode)
    else:
        small_keynodes.append(keynode)

print(len(filter_keynodes))

##assign the nodes in small clusters to the nearest big cluster
out_in_clusters = []
print('process the small nodes')
for cur_keynode in small_keynodes:
    # cur_keynode = small_keynodes[0]
    cur_clusterid = dict_nodestoclusters[cur_keynode]
    cur_cluster_nodes = dict_clusterstonodes[cur_clusterid]
    cur_in_clusters = []
    cur_out_clusters = []
    for each_node in cur_cluster_nodes:
        # each_node = cur_cluster_nodes[0]
        in_cluster_list = list(set(G_dir.predecessors(each_node))) #networkx.exception.NetworkXError: The node nan is not in the digraph.
        out_cluster_list = list(set(G_dir.successors(each_node)))
        in_cluster_list = [ dict_nodestoclusters[x] for x in in_cluster_list ]
        out_cluster_list = [ dict_nodestoclusters[x] for x in out_cluster_list ]
        if len(in_cluster_list) == 0:
            pass
        elif len(in_cluster_list) == 1:
            if in_cluster_list[0] == cur_clusterid:
                pass
            else:
                cur_in_clusters.append(in_cluster_list[0])
        elif len(in_cluster_list) > 1:
            cur_in_clusters += [ x for x in in_cluster_list if x != cur_clusterid ]
        ##for out clusters
        if len(out_cluster_list) == 0:
            pass
        elif len(out_cluster_list) == 1:
            if out_cluster_list[0] == cur_clusterid:
                pass
            else:
                # out_in_clusters.append(out_cluster_list[0])
                cur_out_clusters.append(out_cluster_list[0])
        elif len(out_cluster_list) > 1:
            cur_out_clusters += [ x for x in out_cluster_list if x != cur_clusterid ]
    cur_in_clusters = list(set(cur_in_clusters))
    cur_out_clusters = list(set(cur_out_clusters))
    # print(cur_clusterid, cur_keynode, cur_in_clusters, cur_out_clusters)
    ##merge the nodes in cur clusters into predecessor clusters
    if len(cur_in_clusters) == 0:
        continue
    in_clusterid = cur_in_clusters[0]
    dict_clusterstonodes[in_clusterid] += cur_cluster_nodes
    dict_clusterslen[in_clusterid] += len(cur_cluster_nodes)
    for each_node in cur_cluster_nodes:
        dict_nodestoclusters[each_node] = in_clusterid
    ##remove the corresponding clusters
    dict_clusterstokeynode.pop(cur_clusterid)
    dict_clusterslen.pop(cur_clusterid)
    dict_clusterstonodes.pop(cur_clusterid)

##process the trans nodes
print('process the trans nodes')
keynode_edges = []
# selected_nodes = []
for cur_keynode in filter_keynodes:
    # cur_keynode = filter_keynodes[7]
    cur_clusterid = dict_nodestoclusters[cur_keynode]
    cur_cluster_nodes = dict_clusterstonodes[cur_clusterid]
    cur_in_clusters = []
    cur_out_clusters = []
    for each_node in cur_cluster_nodes:
        # each_node = cur_cluster_nodes[0]
        in_cluster_list = list(set(G_dir.predecessors(each_node)))
        out_cluster_list = list(set(G_dir.successors(each_node)))
        in_cluster_list = [ dict_nodestoclusters[x] for x in in_cluster_list ]
        out_cluster_list = [ dict_nodestoclusters[x] for x in out_cluster_list ]
        if len(in_cluster_list) == 0:
            pass
        elif len(in_cluster_list) == 1:
            if in_cluster_list[0] == cur_clusterid:
                pass
            else:
                cur_in_clusters.append(in_cluster_list[0])
        elif len(in_cluster_list) > 1:
            cur_in_clusters += [ x for x in in_cluster_list if x != cur_clusterid ]
        ##for out clusters
        if len(out_cluster_list) == 0:
            pass
        elif len(out_cluster_list) == 1:
            if out_cluster_list[0] == cur_clusterid:
                pass
            else:
                cur_out_clusters.append(out_cluster_list[0])
        elif len(out_cluster_list) > 1:
            cur_out_clusters += [ x for x in out_cluster_list if x != cur_clusterid ]
    cur_in_clusters = list(set(cur_in_clusters))
    cur_out_clusters = list(set(cur_out_clusters))
    # print(cur_clusterid, cur_keynode, cur_in_clusters, cur_out_clusters)
    #for in 
    if len(cur_in_clusters) == 0:
        pass
    else:
        in_clusterid = cur_in_clusters[0]
        in_keynode = dict_clusterstokeynode[in_clusterid]
        keynode_edges.append( [in_keynode, cur_keynode] )
    #for out
    for each_out_cluster in cur_out_clusters:
        out_keynode = dict_clusterstokeynode[each_out_cluster]
        keynode_edges.append( [cur_keynode, out_keynode] )

keynode_edges_df = pd.DataFrame(keynode_edges)
keynode_edges_df = keynode_edges_df.drop_duplicates()

transnodes_list = []
for each_edge in keynode_edges_df.values:
    # each_edge = keynode_edges_df.values.tolist()[0] 
    node1, node2 = each_edge
    each_path = nx.dijkstra_path(G_smaller,node1, node2)
    if len(each_path) == 2:
        pass
    else:
        transnodes_list += each_path[1:-1]

print('keynode count ', len(set(filter_keynodes)))
print('transnode count ', len(set(transnodes_list)))

##the lineages of haps
meta_df = pd.read_csv('example/samples_metadata.txt', sep='\t', header=0)
meta_df.index = meta_df['accession_id']

##get the final nodes files
foutput = open(data_path + 'nodes_annotated.csv', 'w')
foutput.write("Id,Label,Node_value,Node_cluster,Node_cluster_len,Node_lineage,Node_lineage_raw\n")
for cur_hapidx in range(haplotype_df.shape[0]):
    cur_hap, cur_mut_count, cur_sample_names, cur_sample_count, cur_locations, cur_time_min, cur_time_max, cur_hap_pos, cur_hap_value = haplotype_df.loc[cur_hapidx]
    node = cur_hapidx
    node_hap = cur_hap
    sample_names = cur_sample_names
    if int(node) in filter_keynodes: #the real important nodes
        node_value = 1000 
    elif int(node) in transnodes_list: #the nodes link the important ones
        node_value = 100
    else:
        node_value = 0
    # if int(node) in dict_nodestoclusters.keys():
    print('process ', cur_hapidx )
    node_cluster = dict_nodestoclusters[ int(node) ]
    node_cluster_len = dict_clusterslen[ node_cluster ]
    # node_cluster = nodes_df.loc[cur_hapidx]['Node_cluster']
    # node_cluster_len = nodes_df.loc[cur_hapidx]['Node_cluster_len']
    # #node size
    # if int(node) in nodes_selected_list:
    #     node_value2 = 10
    # else:
    #     node_value2 = 5
    if cur_hapidx == 0:
        node_lineage = 'B'
        # node_lineage2 = 'B'
        node_lineage_raw = 'B'            
    else:
        # node_lineage = meta_df.loc[cur_hapidx]['lineage']
        # node_lineage2 = meta_df.loc[cur_hapidx]['lineage2']
        # node_lineage_raw = node_lineage
        sample_names = sample_names.split(';')
        if len(sample_names) >=20:
            sample_names = random.sample(sample_names, 20)
        cur_lineages = []
        for samplename in sample_names:
            cur_lineage = meta_df.loc[samplename]['pangolin_lineage']
            # cur_lineage = meta_df.loc[samplename]['lineage2']
            cur_lineages.append(cur_lineage)
        # print(cur_hapidx, cur_lineages)
        dict_pangolins = Counter(cur_lineages)
        cur_lineages = sorted(dict_pangolins.items(), key= lambda x :x[1])
        node_lineage = cur_lineages[-1][0]
        node_lineage_raw =  ';'.join([ x[0] + '-' + str(x[1]) for x in cur_lineages ])
        # node_lineage = meta_df.loc[sample_names[0]]['Pango lineage']
        # node_lineage_raw = node_lineage
    node_label = str(node)  +'-'+ node_lineage
    t = [node, node_label, node_value, node_cluster, node_cluster_len, node_lineage, node_lineage_raw]
    t= [str(x) for x in t ]
    foutput.write(','.join(t)+"\n")

foutput.close()

##if you have more than 2000 nodes please define the keynode
##if you have more than 2000 nodes please define the keynode
##if you have more than 2000 nodes please define the keynode
"""
nodes_df = pd.read_csv(data_path + 'nodes_annotated.csv', sep=',', header=0)
nodes_df_selected = nodes_df[(nodes_df['Node_value']==1000) | (nodes_df['Node_value']==100)]
nodes_df_selected.to_csv(data_path + 'nodes_mainpath.csv', sep=',', index=False, header=True)
nodes_df_selected['Node_lineage'].drop_duplicates()

# keynodes_selected = filter_keynodes + transnodes_list
nodes_selected_list = nodes_df_selected['Id'].tolist() 
edges_selected = []
for each in  edges_df.values:
    # source, target = each
    source, target, mutations, dis = each
    if source in nodes_selected_list and target in nodes_selected_list:
        edges_selected.append(each)

edges_selected = pd.DataFrame(edges_selected)
edges_selected.columns = ['source', 'target', 'mutations', 'distance']
# edges_selected.columns = ['source', 'target']
edges_selected.to_csv(data_path + 'edges_mainpath.csv', sep=',', index=False, header=True)
"""

end_time = time.time()
diff = end_time - beg_time
print(time.ctime())
print(diff, ' s')
print(diff/60, ' min')
print(diff/60/60, ' h')


