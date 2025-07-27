#!/usr/bin/env python
# coding=utf-8
# __author__ = 'zp'

import pandas as pd
import numpy as np
from datetime import datetime
from datetime import timedelta
import time
import argparse
import os
import faiss
import random
from collections import Counter
import scipy.sparse
from scipy.sparse import coo_matrix
import multiprocessing

beg_time = time.time()
# parser = argparse.ArgumentParser()
# parser.add_argument('-input_name', '--input_name', dest='input_name', help="the name of input samples mutation file", required=True)
# parser.add_argument('-cpu_count', '--cpu_count', dest='cpu_count', help="the cpu count used", default=5)
# parser.add_argument('-main_hap_rate', '--main_hap_rate', dest='main_hap_rate', help="the rate of main haps in all haps, which determines the resolution of haplotype netwwork. The range is [0,1], default value is 0.1, main hap is the center hap when cluster the haps through k-means method", default=0.1)
# parser.add_argument('-main_hap_count', '--main_hap_count', dest='main_hap_count', help="the main hap count determin the haps used in all haps, the range is [0, all haps count], default value is haps_count*0.1, main hap is the center hap when cluster the haps through k-means method", default=0.1)
# parser.add_argument('-threshold', '--threshold', dest='threshold', help="the samples count used to filter the haps with too little smaples", default=2)
# parser.add_argument('-output_path', '--output_path', dest='output_path', help="the output path of the results", default= 'venas2_output')
# parser.add_argument('-batch_mode', '--batch_mode', dest='batch_mode', help="the number of input samples", default=False)
# parser.add_argument('-initial_size', '--initial_size', dest='initial_size', help="the samples count used to construct the initial network when using batch mode", default=10000)
# parser.add_argument('-batch_size', '--batch_size', dest='batch_size', help="the samples count which added to the network each time when using batch mode", default=10000)
# args = parser.parse_args()

# #parameters
# input_name = args.input_name
# batch_mode = args.batch_mode
# initial_size = int(args.initial_size)
# batch_size = int(args.batch_size)
# can_lens = int(args.can_lens)

##load in the data
##can_lens control the number of top haps considered for each new hap in faiss
can_lens = 500
dict_nc2num = {'a':1, 't':2, 'c':3, 'g':4, '-':5, 'o':6}
dict_num2nc = {1:'a', 2:'t', 3:'c', 4:'g', 5:'-', 6:'o'}
# samplefile =  str(args.input_name)
# samplefile = 'samples'
data_path =  'example/before2020-09/'  
output_path = 'venas2_output/'
if os.path.exists(output_path):
    print('the folder exists! please check the output dir')
else:
    os.makedirs(output_path)

data = pd.read_csv(data_path +  'samples_mutations.txt', sep='\t', header=0)
data.columns = ['samplename', 'haplotype' , 'collection_date', 'location']
# data = data[data['collection_date'].notna()]
data = data.fillna('')
data['time'] = data['collection_date'].apply(lambda x: datetime.strptime(x, '%Y-%m-%d'))
# data = data.sort_values('time',ascending=True)
# data_count = data.shape[0]

def get_edge(hapidx1, hapidx2):
    # hapidx1 = parent_hapidx
    # hapidx2 = cur_hapidx
    hap1 = set(haplotype_df_updated.loc[hapidx1]['haplotype'].split(';'))
    hap2 = set(haplotype_df_updated.loc[hapidx2]['haplotype'].split(';'))
    t1 = hap1.difference(hap2)
    t2 = hap2.difference(hap1)
    # diff_count = len(t1) + len(t2)
    # diff_fn = [t1, t2]
    if t1==set() and t2!=set():
        edge_muts = ';'.join([ x for x in t2 ])
        t1_count = 0
        t2_count = len(t2)
    elif t1!=set() and t2==set():
        edge_muts = ';'.join([ 'bm:'+x for x in t1 ])
        t1_count = len(t1)
        t2_count = 0
    elif t1!=set() and t2!=set():
        edge_muts = ';'.join([ 'bm:'+x for x in t1 ]) + ';' + ';'.join([ x for x in t2 ])
        t1_count = len(t1)
        t2_count = len(t2)
    else:
        edge_muts = ''
        t1_count = t2_count = 0
    edge_diff_count = t1_count + t2_count
    return edge_muts, edge_diff_count

def convert_smaples_mutation_into_hap_df_array(data_fn):
    # data_fn = data
    haplotype_df_fn = []
    hap_pos_fn = []
    for each_haplotype, group in data_fn.groupby('haplotype'):
        # print(each_haplotype, group)
        sample_names = group['samplename'].tolist()
        sample_count = len(sample_names)
        sample_locations = ';'.join(group['location'].drop_duplicates().tolist())
        sample_times = group['time'].drop_duplicates().tolist()
        # sample_times = [ datetime.strptime(x, '%Y-%m-%d') for x in sample_times]
        sample_times_min = min(sample_times)
        sample_times_max = max(sample_times)
        if each_haplotype == '':
            muts_count = 0
            haplotype_df_fn.append([each_haplotype, muts_count, ';'.join(sample_names), sample_count, sample_locations, sample_times_min, sample_times_max, '','' ])
        else:
            each_haps = each_haplotype.split(';')
            muts_count = len(each_haps)
            hap_pos = []
            hap_value = []
            for each_cur_hap in each_haps:
                if 'del' in each_cur_hap:
                    # each_cur_hap = each_haps[0]
                    each_cur_pos = int(each_cur_hap.split('_')[0][3:])
                    each_cur_value = 5
                else:
                    each_cur_pos = int(each_cur_hap[1:-1])
                    each_cur_value = dict_nc2num[ each_cur_hap[-1:].lower() ]
                hap_pos.append(each_cur_pos)
                hap_value.append(each_cur_value)
            haplotype_df_fn.append([each_haplotype, muts_count, ';'.join(sample_names), sample_count, sample_locations, sample_times_min, sample_times_max, hap_pos, hap_value])
            hap_pos_fn += hap_pos
    haplotype_df_fn = pd.DataFrame(haplotype_df_fn)
    haplotype_df_fn.columns = ['haplotype', 'muts_count', 'sample_names', 'sample_count',  'sample_locations', 'sample_time_min', 'sample_time_max', 'hap_pos', 'hap_value']
    haplotype_df_fn = haplotype_df_fn.sort_values([ 'sample_time_min', 'muts_count','sample_count'], ascending=[True, True, False])
    # haplotype_df_fn = haplotype_df_fn.sort_values([ 'muts_count','sample_time_min','sample_count'], ascending=[True, True, False])
    # haplotype_df = pd.concat( [haplotype_df_first, haplotype_df_rest], axis=0 )
    haplotype_df_fn.index= list(range(haplotype_df_fn.shape[0]))
    hap_count_fn = haplotype_df_fn.shape[0]
    hap_pos_fn = list(sorted(set(hap_pos_fn)))
    ##convert to array
    hap_pos_count = len(hap_pos_fn)
    dict_pos_index = dict([ (x[1],x[0]) for x in enumerate(hap_pos_fn) ])
    rows = []
    cols = []
    values = []
    for each_hapidx in range(hap_count_fn):
        # each_hapidx = 41
        # print('processing ',each_hapidx)
        cur_hap_pos = haplotype_df_fn.loc[each_hapidx]['hap_pos']
        if cur_hap_pos =='':
            continue
        tmpdf = pd.DataFrame(cur_hap_pos)
        tmpdf['posidx'] = tmpdf.apply(lambda x:dict_pos_index[x[0]], axis=1)
        cur_hap_value = haplotype_df_fn.loc[each_hapidx]['hap_value']
        rows = rows + [each_hapidx]*tmpdf.shape[0]
        cols = cols + tmpdf['posidx'].tolist()
        values = values + cur_hap_value
    values = np.array(values)
    rows = np.array(rows)
    cols = np.array(cols)
    # coo = coo_matrix((values, (rows, cols)), shape=(hap_count_fn, hap_pos_count), dtype=np.int8)
    # # hap_array_df = pd.DataFrame(coo.todense())
    # # hap_array_df.columns = hap_pos_fn
    # hap_array_fn = coo.toarray()
    return haplotype_df_fn, hap_pos_fn, rows, cols, values

"""step1: get the haplotype vector form and update the latest hap"""
haplotype_df_base, hap_pos_base, rows_base, cols_base, values_base = convert_smaples_mutation_into_hap_df_array(data)
coo = coo_matrix((values_base, (rows_base, cols_base)), shape=(haplotype_df_base.shape[0], len(hap_pos_base)), dtype=np.int8)
hap_array_base = coo.toarray()
hap_array_df_base = pd.DataFrame(hap_array_base)
hap_array_df_base.columns = hap_pos_base
# np.save( output_path+'rows_base.npy', rows_base)
# np.save( output_path+'cols_base.npy', cols_base)
# np.save( output_path+'values_base.npy', values_base)
# haplotype_df_base.to_csv(output_path+'haplotype_df_base.txt',sep='\t', header=True, index=None)
# pd.DataFrame(hap_pos_base).to_csv(output_path+'hap_pos_base.txt',sep='\t', header=False, index=None)


#read in files
# data_path_coo = 'results/venas2v1_output/venas2_output_month19/before2020-09/'
# data_path_coo = output_path
# rows1 = np.load( data_path_coo+'rows_base.npy')
# cols1 = np.load( data_path_coo+'cols_base.npy')
# values1 = np.load( data_path_coo+'values_base.npy')
# haplotype_df_base = pd.read_csv(data_path_coo+'haplotype_df_base.txt',sep='\t', header=0)
# haplotype_df_base = haplotype_df_base.fillna('')
# haplotype_df_base['sample_time_min'] = haplotype_df_base['sample_time_min'].apply(lambda x: datetime.strptime(x, '%Y-%m-%d'))
# haplotype_df_base['sample_time_max'] = haplotype_df_base['sample_time_max'].apply(lambda x: datetime.strptime(x, '%Y-%m-%d'))
# # haplotype_df_base['hap_pos'] = haplotype_df_base['hap_pos'].apply(lambda row: [ int(x) for x in row[1:-1].split(',')])
# hap_pos_base = pd.read_csv(data_path_coo +  'hap_pos_base.txt', sep='\t', header=None)
# hap_pos_base = hap_pos_base[0].tolist()
# coo1 = coo_matrix((values1, (rows1, cols1)), shape=(haplotype_df_base.shape[0], len(hap_pos_base)), dtype=np.int8)
# hap_array_base = coo1.toarray()
# hap_array_df_base = pd.DataFrame(hap_array_base)
# hap_array_df_base.columns = hap_pos_base

print('the data count is ', data.shape[0], 'the haplotype count is ', haplotype_df_base.shape[0])

# orig = haplotype_df_base.copy()
# haplotype_df_base = orig.copy()
#make sure the first is ref
if '' in haplotype_df_base['haplotype'].tolist():
    idx = haplotype_df_base[haplotype_df_base['haplotype']==''].index[0]
    haplotype_df_base.loc[idx, 'sample_time_min'] = min(haplotype_df_base.loc[idx]['sample_time_min'] ,datetime.strptime('2019-12-31', '%Y-%m-%d') )
    if idx == 0:
        pass
    else:
        haps_list = haplotype_df_base.index.tolist()
        haps_list.remove(idx)
        haps_list = [idx] + haps_list
        haplotype_df_base = haplotype_df_base.loc[haps_list]
        haplotype_df_base.index = range(haplotype_df_base.shape[0])
        # hap_array_base = hap_array_base[haps_list]
        hap_array_df_base = hap_array_df_base.loc[haps_list]
        hap_array_df_base.index= range(hap_array_df_base.shape[0])
    haplotype_df_updated = haplotype_df_base.copy()
    hap_array_df_updated = hap_array_df_base.copy()
    hap_array_updated = np.array(hap_array_df_updated)
else:
    ref_df = pd.DataFrame(['', 0, 'ref', 1, 'China', datetime.strptime('2019-12-31', '%Y-%m-%d'), datetime.strptime('2019-12-31', '%Y-%m-%d'), [], [] ]).T
    ref_df.columns = ['haplotype', 'muts_count', 'sample_names', 'sample_count',  'sample_locations', 'sample_time_min', 'sample_time_max', 'hap_pos', 'hap_value']
    haplotype_df_updated = pd.concat([ref_df, haplotype_df_base ],axis=0)
    haplotype_df_updated.index = range(haplotype_df_updated.shape[0])
    tmpdf = pd.DataFrame(np.zeros((1, hap_array_df_base.shape[1])))
    tmpdf.columns = hap_array_df_base.columns 
    hap_array_df_updated = pd.concat([tmpdf, hap_array_df_base], axis=0)
    hap_array_df_updated.index= range(hap_array_df_updated.shape[0])
    hap_array_updated = np.array(hap_array_df_updated)

#process the hap pos this is for load the data
# if haplotype_df_updated.iloc[0]['haplotype']=='':
#     # idx = haplotype_df_updated[haplotype_df_updated['haplotype']==''].index[0]
#     haplotype_df_updated1 = haplotype_df_updated[haplotype_df_updated['haplotype']=='']
#     haplotype_df_updated2 = haplotype_df_updated[haplotype_df_updated['haplotype']!='']
#     haplotype_df_updated2['hap_pos'] = haplotype_df_updated2['hap_pos'].apply(lambda x: [int(t) for t in x[1:-1].split(',')])
#     haplotype_df_updated2['hap_value'] = haplotype_df_updated2['hap_value'].apply(lambda x: [int(t) for t in x[1:-1].split(',')])
#     haplotype_df_updated = pd.concat([haplotype_df_updated1, haplotype_df_updated2], axis=0)
#     haps_list =haplotype_df_updated.index.tolist()
#     haplotype_df_updated.index = range(haplotype_df_updated.shape[0])
#     hap_array_df_updated = hap_array_df_updated.loc[haps_list]
#     hap_array_df_updated.index = range(hap_array_df_updated.shape[0])
# else:
#     haplotype_df_updated['hap_pos'] = haplotype_df_updated['hap_pos'].apply(lambda x: [int(t) for t in x[1:-1].split(',')])
#     haplotype_df_updated['hap_value'] = haplotype_df_updated['hap_value'].apply(lambda x: [int(t) for t in x[1:-1].split(',')])


del haplotype_df_base
del hap_array_df_base
# """stage 2: identify the parent hap for each hap"""
print('constructing the faiss index and query...')
# beg_time = time.time()
can_lens = 500
xb = hap_array_updated.astype(np.int8)
xq = hap_array_updated.astype(np.int8)
dim= xq.shape[1]
# xb = array
# xq = array
# dim = array.shape[1]
measure = faiss.METRIC_L2
param = 'Flat'
# param = 'HNSW64'
index = faiss.index_factory(dim, param, measure)
# print(index.is_trained)
index.add(xb)
can_dis, can_ids = index.search(xq, can_lens)

# np.save( output_path+'can_dis.npy', can_dis)
# np.save( output_path+'can_ids.npy', can_ids)
# can_dis = np.load( output_path+'can_dis.npy')
# can_ids = np.load( output_path+'can_ids.npy')

"""step 3 identify the parent hap for each hap"""
def identify_edges( haplotype_df_updated, can_ids_batch, hap_ids_batch, batch_name ):
    print('start the task in the ', batch_name,'th cpu ')
    # nearest_indexs = []
    # results_logs = []
    reindex_list= []
    results_voi = []
    edges_df = pd.DataFrame(columns=[0, 1,2,3])
    hap_count_batch = can_ids_batch.shape[0]
    for cur_hap_index in range( hap_count_batch):
    # for cur_hap_index in range( 20):
        # cur_hap_index = 17
        cur_hapidx = hap_ids_batch[cur_hap_index]
        # if cur_hap_index % 1 == 0:
        #     print('processing the ', str(cur_hap_index), '/', str(cur_hapidx), '/', str(len(hap_ids)))
        cur_hap, cur_mut_count, cur_sample_names, cur_sample_count, cur_locations, cur_time_min, cur_time_max, cur_hap_pos, cur_hap_value = haplotype_df_updated.loc[cur_hapidx]
        cur_locations = cur_locations.split(';')
        cur_can_ids = can_ids_batch[cur_hap_index].tolist()
        # cur_can_ids2 = [ x for x in cur_can_ids if x!=-1 and x < cur_hapidx]
        if cur_mut_count == 0:
            continue
        if cur_mut_count == 1: #no choice but the reference hap
            parent_hapidx = 0
            edges_df.loc[edges_df.shape[0],:] = [parent_hapidx, cur_hapidx, cur_hap, cur_mut_count]
            continue
        flag_count = 0
        parent_can_list = []
        dict_parent_can = {}
        a = set(cur_hap.split(';'))
        # print('processing hap: ', cur_hap_index, cur_hapidx,'\n' ,'mutation positions: ', a, cur_hap_pos,'\n', 'values: ', cur_hap_value, '\n', [ x for x in cur_can_ids if x!=-1 ][:20], cur_time_min,cur_time_max, cur_locations)
        # a = hap_array_df_updated.loc[cur_hapidx] #first stage check the top 10 nearest haps 
        # a = a[a>0]
        for each_hap_index in range(1, can_lens):
            # each_hap_index = 10
            each_hapidx = cur_can_ids[each_hap_index]
            if each_hapidx == -1:
                break
            if flag_count >= 10:
                break
            elif each_hapidx >= cur_hapidx:
                continue
            # each_hap_pos = haplotype_df_updated.loc[each_hapidx]['hap_pos']
            flag_count += 1
            b = set(haplotype_df_updated.loc[each_hapidx]['haplotype'].split(';'))
            # b = hap_array_df_updated.loc[each_hapidx]
            # b = b[b>0]
            # diff = a.compare(b)
            # diff_count = diff.shape[0]
            if each_hapidx == 0:
                t1= a
                t2= set()
                diff_count = cur_mut_count
            else:
                t1= a.difference(b)
                t2 = b.difference(a)
                diff_count = len(t1) + len(t2)
            diff = [t1, t2]
            dict_parent_can[each_hapidx] = diff
            parent_can_list.append( [each_hapidx,  diff_count])
            # print( 'candidate parents: ', each_hapidx, diff_count, '\n', diff , haplotype_df_updated.loc[each_hapidx][['sample_time_min', 'sample_time_max', 'sample_locations', 'sample_count']].tolist())
        if len(parent_can_list) > 0:
            parent_can_df = pd.DataFrame(parent_can_list)
            parent_can_df = parent_can_df.sort_values([1], ascending=True)
            nearest_df = parent_can_df[parent_can_df[1]==parent_can_df.iloc[0][1]] # select the haps with minimum distance
            common_flag = False
            if nearest_df.shape[0] == 1: # this means there is only 1 candidate ab for abc          
                nearest_hapidx = nearest_df.iloc[0][0]
                # nearest_diff = dict_parent_can[nearest_hapidx]
                # nonzero = np.count_nonzero(nearest_diff['self'])
                # if nearest_diff.shape[0]==nonzero and nearest_diff.shape[0] <= 5:
                parent_hapidx = nearest_hapidx
                common_flag = True
            elif nearest_df.iloc[0][1] <= 3:
                common_flag = True
                select_list = nearest_df[0].tolist()
                nearest_df.index = select_list
                parent_can_df = pd.concat([ nearest_df,  haplotype_df_updated.loc[select_list][['sample_time_min', 'sample_time_max', 'sample_locations', 'sample_count']] ], axis=1)
                parent_hapidx_orig = parent_can_df.sort_values([1, 'sample_count', 'sample_time_max'], ascending=[True, False, False]).index[0]
                select_list2 = []
                for each_idx in range(parent_can_df.shape[0]):
                    each_parent_hapidx, diff_count, each_time_min, each_time_max, each_locations, sample_count = parent_can_df.iloc[each_idx]
                    each_locations = each_locations.split(';')
                    delta_time = each_time_max - cur_time_min
                    if delta_time.days >= -60 and len(set(cur_locations).intersection(set(each_locations))) >= 1:
                        select_list2.append(each_parent_hapidx)
                if len(select_list2) == 0:
                    parent_hapidx = parent_can_df.loc[select_list].sort_values([1, 'sample_count', 'sample_time_max'], ascending=[True, False, False]).index[0]
                elif len(select_list2) == 1:
                    parent_hapidx = select_list2[0]
                else: #select the minest dis else use sample count
                    # parent_can_df = parent_can_df.loc[select_list2]
                    parent_hapidx = parent_can_df.loc[select_list2].sort_values([1, 'sample_count', 'sample_time_max'], ascending=[True, False, False]).index[0]
                if parent_hapidx != parent_hapidx_orig:
                    results_voi.append([cur_hapidx, parent_hapidx, parent_hapidx_orig])
            if common_flag == False: #stage 2 continue check top 10 - 50 nearest haps 
                flag_count = 0
                parent_can_list2 = []
                # a = hap_array_df_updated.loc[cur_hapidx]
                for each_hap_index2 in range(each_hap_index, can_lens):
                    # each_hap_index2 = 1
                    each_hapidx = cur_can_ids[each_hap_index2]
                    if each_hapidx == -1:
                        break
                    if flag_count >= 30:
                        break
                    elif each_hapidx >= cur_hapidx:
                        continue
                    # each_hapidx = 212
                    # each_hap_pos = haplotype_df_updated.loc[each_hapidx]['hap_pos']
                    flag_count += 1
                    b = set(haplotype_df_updated.loc[each_hapidx]['haplotype'].split(';'))
                    if each_hapidx == 0:
                        t1= a
                        t2= set()
                        diff_count = cur_mut_count
                    else:
                        t1= a.difference(b)
                        t2 = b.difference(a)
                        diff_count = len(t1) + len(t2)
                    diff = [t1, t2]
                    dict_parent_can[each_hapidx] = diff
                    parent_can_list2.append( [each_hapidx,  diff_count])
                    # print( 'candidate parents: ', each_hapidx,diff_count, '\n',b, diff , haplotype_df_updated.loc[each_hapidx][['sample_time_min', 'sample_time_max', 'sample_locations', 'sample_count']].tolist())
                    # results_logs.append( ['candidates' ,each_hapidx, ';'.join([str(x) for x in each_hap_pos]), [ '-'.join([str(t) for t in x]) for x in zip(diff.index.tolist(), diff['self'].tolist(), diff['other'].tolist()) ] ])
                parent_can_df = pd.DataFrame(parent_can_list+parent_can_list2 )
                parent_can_df = parent_can_df.sort_values([1], ascending=True)
                nearest_df = parent_can_df[parent_can_df[1]==parent_can_df.iloc[0][1]] # select the haps with minimum distance
                if nearest_df.shape[0] == 1: # this means there is only 1 candidate 
                    parent_hapidx = nearest_df.iloc[0][0]
                else:
                    select_list = nearest_df[0].tolist()
                    nearest_df.index = select_list
                    parent_can_df = pd.concat([ nearest_df,  haplotype_df_updated.loc[select_list][['sample_time_min', 'sample_time_max', 'sample_locations', 'sample_count']] ], axis=1)
                    parent_hapidx_orig = parent_can_df.sort_values([1, 'sample_count', 'sample_time_max'], ascending=[True, False, False]).index[0]
                    select_list2 = []
                    for each_idx in range(parent_can_df.shape[0]):
                        each_parent_hapidx, diff_count, each_time_min, each_time_max, each_locations, sample_count = parent_can_df.iloc[each_idx]
                        each_locations = each_locations.split(';')
                        delta_time = each_time_max - cur_time_min
                        if delta_time.days >= -60 and len(set(cur_locations).intersection(set(each_locations))) >= 1:
                            select_list2.append(each_parent_hapidx)
                    if len(select_list2) == 0:
                        parent_hapidx = parent_can_df.loc[select_list].sort_values([1, 'sample_count', 'sample_time_max'], ascending=[True, False, False]).index[0]
                    elif len(select_list2) == 1:
                        parent_hapidx = select_list2[0]
                    else: #select the minest dis else use sample count
                        # parent_can_df = parent_can_df.loc[select_list2]
                        parent_hapidx = parent_can_df.loc[select_list2].sort_values([1, 'sample_count', 'sample_time_max'], ascending=[True, False, False]).index[0]
                    if parent_hapidx != parent_hapidx_orig:
                        results_voi.append([cur_hapidx, parent_hapidx, parent_hapidx_orig])
            # print('the final identified parent is: ', parent_hapidx)
            t1, t2 = dict_parent_can[parent_hapidx]
            if t1==set() and t2!=set():
                edge_muts = ';'.join([ 'bm:'+ x for x in t2 ])
                t1_count = 0
                t2_count = len(t2)
            elif t1!=set() and t2==set():
                edge_muts = ';'.join([ x for x in t1 ])
                t1_count = len(t1)
                t2_count = 0
            elif t1!=set() and t2!=set():
                edge_muts = ';'.join([ x for x in t1 ]) + ';' + ';'.join(['bm:'+x for x in t2 ])
                t1_count = len(t1)
                t2_count = len(t2)
            else:
                edge_muts = ''
                t1_count = t2_count = 0
            edge_diff_count = t1_count + t2_count
            # return edge_muts, edge_diff_count
            # edge_muts, edge_diff_count = get_edge(parent_hapidx, cur_hapidx)
            edges_df.loc[edges_df.shape[0], :] = [parent_hapidx, cur_hapidx, edge_muts, edge_diff_count]
        else: #this maybe caused by the can_lens is set too small just re-search for cur hap
            reindex_list.append(cur_hapidx)
    edges_df.columns = ['source', 'target', 'mutations', 'distance']
    edges_df.to_csv(output_path+'edges'+batch_name+'.csv',sep=',', index=None , header=True)
    results_voi = pd.DataFrame(results_voi)
    results_voi.to_csv(output_path+'results_voi'+batch_name+'.txt',sep='\t', index=None , header=None)
    reindex_list = pd.DataFrame(reindex_list)
    reindex_list.to_csv(output_path+'reindex_list'+batch_name+'.txt',sep='\t', index=None , header=None)

hap_ids = hap_array_df_updated.index.tolist()
hap_count = len(hap_ids)
# cpu_count = int(multiprocessing.cpu_count())
# cpu_count = int(args.cpu_count)
cpu_count = 1
batch_size = int(hap_count/cpu_count)

p = multiprocessing.Pool(cpu_count)
for i in range(cpu_count):
    # i= 0
    if i == cpu_count - 1:
        hap_ids_batch = hap_ids[ batch_size * i: ]
        can_ids_batch = can_ids[batch_size * i:]
    else:
        hap_ids_batch = hap_ids[batch_size * i : batch_size*(i+1) ]
        can_ids_batch = can_ids[batch_size * i: batch_size*(i+1)]
    batch_name = str(i)
    # print(batch_name, hap_ids_batch[:9], can_ids_batch[:9])
    p.apply_async(identify_edges, args=(haplotype_df_updated, can_ids_batch, hap_ids_batch, batch_name ))

p.close()
p.join()

print('the tasks in different cpus are over!')
edges_df = pd.DataFrame()
results_voi = pd.DataFrame()
reindex_list = pd.DataFrame()
for i in range(cpu_count):
    batch_name = str(i)
    edges_df_tmp = pd.read_csv(output_path+'edges'+batch_name+'.csv',sep=',', header=0)
    edges_df = pd.concat([edges_df, edges_df_tmp], axis=0)
    if os.path.getsize(output_path + 'results_voi'+batch_name+'.txt') == 0:
        results_voi_tmp = pd.DataFrame()
    else:
        results_voi_tmp = pd.read_csv(output_path+'results_voi'+batch_name+'.txt',sep='\t', header=None)
    results_voi = pd.concat([results_voi, results_voi_tmp], axis=0)
    if os.path.getsize(output_path + 'reindex_list'+batch_name+'.txt') == 0:
        reindex_list_tmp = pd.DataFrame()
    else:
        reindex_list_tmp = pd.read_csv(output_path+'reindex_list'+batch_name+'.txt',sep='\t', header=None)
    reindex_list = pd.concat([reindex_list, reindex_list_tmp], axis=0)

print('processing the reindex ')
# some haps can not find the colest ones from top can_len=500, so we increase the parameter
if reindex_list.shape[0]>0:
    reindex_list = reindex_list[0].tolist()
    # hap_array_batch = hap_array_updated[cur_hapidx]
    hap_array_batch = hap_array_updated[reindex_list]
    # can_lens = xb.shape[0]
    can_lens_batch = 5000
    # can_lens = int(hap_array_df_updated.shape[0]/2)
    xb = hap_array_updated.astype(np.int8)
    xq = hap_array_batch.astype(np.int8)
    if len(reindex_list)==1:
        xq = xq.reshape(1,-1)
    measure = faiss.METRIC_L2
    param = 'Flat'
    # param = 'HNSW64'
    # dim = len(hap_pos_batch)
    dim= hap_array_updated.shape[1]
    index = faiss.index_factory(dim, param, measure)
    # print(index.is_trained)
    index.add(xb)
    can_dis_batch, can_ids_batch = index.search(xq, can_lens_batch)
    # edges_df = pd.DataFrame(columns=[0, 1,2,3])
    for cur_hap_index in range(len(reindex_list)):
        # cur_hap_index = 0
        cur_hapidx = reindex_list[cur_hap_index]
        cur_hap, cur_mut_count, cur_sample_names, cur_sample_count, cur_locations, cur_time_min, cur_time_max, cur_hap_pos, cur_hap_value = haplotype_df_updated.loc[cur_hapidx]
        cur_locations = cur_locations.split(';')
        cur_can_ids = can_ids_batch[cur_hap_index].tolist()
        # cur_can_ids2 = [ x for x in cur_can_ids if x!=-1 and x < cur_hapidx]
        flag_count = 0
        parent_can_list = []
        dict_parent_can = {}
        a = set(cur_hap.split(';'))
        # print('processing hap: ', cur_hap_index, cur_hapidx,'\n' ,'mutation positions: ', a, cur_hap_pos,'\n', 'values: ', cur_hap_value, '\n', [ x for x in cur_can_ids if x!=-1 ][:20], cur_time_min,cur_time_max, cur_locations)
        # a = hap_array_df_updated.loc[cur_hapidx] #first stage check the top 10 nearest haps 
        for each_hap_index in range(1, can_lens_batch):
            # each_hap_index = 1
            each_hapidx = cur_can_ids[each_hap_index]
            if each_hapidx == -1:
                break
            if flag_count >= 10:
                break
            elif each_hapidx >= cur_hapidx:
                continue
            # each_hap_pos = haplotype_df_updated.loc[each_hapidx]['hap_pos']
            flag_count += 1
            b = set(haplotype_df_updated.loc[each_hapidx]['haplotype'].split(';'))
            if each_hapidx == 0:
                t1= a
                t2= set()
                diff_count = cur_mut_count
            else:
                t1= a.difference(b)
                t2 = b.difference(a)
                diff_count = len(t1) + len(t2)
            diff = [t1, t2]
            dict_parent_can[each_hapidx] = diff
            parent_can_list.append( [each_hapidx,  diff_count])
            # print( 'candidate parents: ', each_hapidx, diff , haplotype_df_updated.loc[each_hapidx][['sample_time_min', 'sample_time_max', 'sample_locations', 'sample_count']].tolist())
        if len(parent_can_list) > 0: 
            parent_can_df = pd.DataFrame(parent_can_list)
            parent_can_df = parent_can_df.sort_values([1], ascending=True)
            nearest_df = parent_can_df[parent_can_df[1]==parent_can_df.iloc[0][1]] # select the haps with minimum distance
            parent_hapidx = nearest_df.iloc[0][0]
            print('the final identified parent is: ', parent_hapidx)
            t1, t2 = dict_parent_can[parent_hapidx]
            if t1==set() and t2!=set():
                edge_muts = ';'.join([ 'bm:'+x for x in t2 ])
                t1_count = 0
                t2_count = len(t2)
            elif t1!=set() and t2==set():
                edge_muts = ';'.join([  x for x in t1 ])
                t1_count = len(t1)
                t2_count = 0
            elif t1!=set() and t2!=set():
                edge_muts = ';'.join([  x for x in t1 ]) + ';' + ';'.join([ 'bm:'+x for x in t2 ])
                t1_count = len(t1)
                t2_count = len(t2)
            else:
                edge_muts = ''
                t1_count = t2_count = 0
            edge_diff_count = t1_count + t2_count
            # edge_muts, edge_diff_count = get_edge(parent_hapidx, cur_hapidx)
            edges_df.loc[edges_df.shape[0], :] = [parent_hapidx, cur_hapidx, edge_muts, edge_diff_count]
        else: 
            parent_hapidx = 0
            edges_df.loc[edges_df.shape[0],:] = [parent_hapidx, cur_hapidx, cur_hap, cur_mut_count]

"""save files"""
# haplotype_df_month.to_csv(output_path+'haplotype_tmp.txt',sep='\t', header=True, index=None)
haplotype_df_updated.to_csv(output_path+'haplotype_final.txt',sep='\t', header=True, index=None)
# edges_df = pd.DataFrame(edges_list)
# edges_df.columns = ['source', 'target', 'mutations', 'distance']
edges_df.to_csv(output_path+'edges.csv',sep=',', index=None , header=True)
# results_voi = pd.DataFrame(results_voi)
results_voi.to_csv(output_path+'results_voi.txt',sep='\t', index=None , header=None)

end_time = time.time()
diff = end_time - beg_time
print(time.ctime())
print(diff, ' s')
print(diff/60, ' min')
print(diff/60/60, ' h')

