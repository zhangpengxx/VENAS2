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
import re
import faiss
import random
from scipy.sparse import coo_matrix

beg_time = time.time()
parser = argparse.ArgumentParser()
parser.add_argument('-input_name', '--input_name', dest='input_name', help="the name of input samples mutation file", required=True)
# # parser.add_argument('-output_path', '--output_path', dest='output_path', help="the output path of the results", default= 'venas2_output')
# parser.add_argument('-batch_mode', '--batch_mode', dest='batch_mode', help="the number of input samples", default=False)
# parser.add_argument('-initial_size', '--initial_size', dest='initial_size', help="the samples count used to construct the initial network when using batch mode", default=10000)
# parser.add_argument('-batch_size', '--batch_size', dest='batch_size', help="the samples count which added to the network each time when using batch mode", default=10000)
args = parser.parse_args()

# #parameters
input_name = args.input_name
# batch_mode = args.batch_mode
# initial_size = int(args.initial_size)
# batch_size = int(args.batch_size)
# can_lens = int(args.can_lens)
# input_name = '2020-10'

##load in the data
dict_nc2num = {'a':1, 't':2, 'c':3, 'g':4, '-':5, 'o':6}
dict_num2nc = {1:'a', 2:'t', 3:'c', 4:'g', 5:'-', 6:'o'}
##can_lens control the number of top haps considered for each new hap in faiss
can_lens = 500

batch_name_list = ['2020-01', '2020-02', '2020-03', '2020-04', '2020-05', '2020-06', '2020-07', '2020-08', '2020-09', '2020-10', '2020-11', '2020-12', '2021-01', '2021-02', '2021-03', '2021-04', '2021-05', '2021-06', '2021-07', '2021-08', '2021-09', '2021-10', '2021-11', '2021-12', '2022-01', '2022-02', '2022-03', '2022-04', '2022-05', '2022-06', '2022-07', '2022-08', '2022-09', '2022-10', '2022-11', '2022-12', '2023-01', '2023-02', '2023-03', '2023-04', '2023-05', '2023-06', '2023-07', '2023-08', '2023-09', '2023-10', '2023-11', '2023-12', '2024-01', '2024-02', '2024-03', '2024-04', '2024-05', '2024-06', '2024-07', '2024-08', '2024-09', '2024-10', '2024-11']

# output_path = 'results/venas2v1_output/venas2_output_month23/merged_' + str(input_name) + '/'
data_path_base = 'datasets/'
output_path = data_path_base+ 'venas2_output_month12/merged_' + str(input_name) + '/'
if os.path.exists(output_path):
    print('the folder exists! please check the output dir')
else:
    os.makedirs(output_path)

basefile = 'before2020-09'
base_path = data_path_base+ 'venas2_output_month12/'+str(basefile) + '/'

def read_haplotype(data_path_fn, samplefile):
    # data_path_fn = batch_path
    # samplefile = 'haplotype_final'
    haplotype_df = pd.read_csv( data_path_fn + samplefile + '.txt', sep='\t', header=0)
    haplotype_df = haplotype_df.fillna('')
    haplotype_df.columns = ['haplotype', 'muts_count', 'sample_names', 'sample_count','sample_locations', 'sample_time_min', 'sample_time_max', 'hap_pos', 'hap_value']
    haplotype_df['sample_time_min'] = haplotype_df['sample_time_min'].apply(lambda x:  datetime.strptime(x, '%Y-%m-%d'))
    haplotype_df['sample_time_max'] = haplotype_df['sample_time_max'].apply(lambda x:  datetime.strptime(x, '%Y-%m-%d'))
    if '' in haplotype_df['haplotype'].tolist():
    # if haplotype_df.iloc[0]['haplotype'] =='':
        haplotype_df1 = haplotype_df[haplotype_df['haplotype'] == '']
        haplotype_df2 = haplotype_df[haplotype_df['haplotype'] != '']
        haplotype_df2.loc[:,'hap_pos'] = haplotype_df2['hap_pos'].apply(lambda x: [int(t) for t in x[1:-1].split(',')])
        haplotype_df2.loc[:,'hap_value'] = haplotype_df2['hap_value'].apply(lambda x: [int(t) for t in x[1:-1].split(',')])
        haplotype_df = pd.concat([haplotype_df1, haplotype_df2], axis=0)
        # haps_list =haplotype_df.index.tolist()
        # haplotype_df.index = range(haplotype_df.shape[0])
    else:
        haplotype_df.loc[:,'hap_pos'] = haplotype_df['hap_pos'].apply(lambda x: [int(t) for t in x[1:-1].split(',')])
        haplotype_df.loc[:,'hap_value'] = haplotype_df['hap_value'].apply(lambda x: [int(t) for t in x[1:-1].split(',')])
    return haplotype_df

haplotype_df_base = read_haplotype(base_path, 'haplotype_final')
edges_df_base = pd.read_csv(base_path + 'edges.csv', sep=',', header=0)
edges_df_base['source'] = edges_df_base['source'].apply(lambda x: int(x))
edges_df_base['target'] = edges_df_base['target'].apply(lambda x: int(x))
if os.path.getsize(base_path + 'results_voi.txt') == 0:
    results_voi_base = pd.DataFrame()
else:
    results_voi_base = pd.read_csv(base_path + 'results_voi.txt', sep='\t', header=None)
    results_voi_base.columns = ['target', 'source', 'source_orig']

haplotype_df_updated = haplotype_df_base.copy()

def read_haplotype_raw(data_path_fn,samplefile):
    # data_path_fn = batch_path
    # samplefile = 'haplotype_final'
    haplotype_df = pd.read_csv( data_path_fn + samplefile + '.txt', sep='\t', header=0)
    haplotype_df = haplotype_df.fillna('')
    haplotype_df.columns = ['haplotype', 'muts_count', 'sample_names', 'sample_count','sample_locations', 'sample_time_min', 'sample_time_max', 'hap_pos', 'hap_value']
    haplotype_df['sample_time_min'] = haplotype_df['sample_time_min'].apply(lambda x:  datetime.strptime(x, '%Y-%m-%d'))
    haplotype_df['sample_time_max'] = haplotype_df['sample_time_max'].apply(lambda x:  datetime.strptime(x, '%Y-%m-%d'))
    return haplotype_df

edges_df = edges_df_base.copy()
results_voi = pd.DataFrame(columns = ['target', 'source', 'source_orig'])
# recom_list = pd.DataFrame(columns = ['idx', 'source', 'target', 'acceptors', 'bpidx', 'class'])
# recom_list_time = pd.DataFrame(columns = ['idx', 'source', 'target', 'acceptors', 'bpidx', 'class'])
# edges_df = pd.DataFrame(columns=[0, 1,2,3])
assert input_name in batch_name_list
batch_idx = batch_name_list.index(input_name)
for month in batch_name_list[8:batch_idx+1]:
# for month in batch_name_list[21:batch_idx+1]:
    # month = batch_name_list[9]
    print('processing the month ', month)
    # batch_path_raw = 'results/venas2v1_output/haps_month_main/'+str(month) + '/'
    batch_path_raw = data_path_base+ 'haps_month12/'+str(month) + '/'
    batch_path = data_path_base + 'venas2_output_month12/'+str(month) + '/'
    haplotype_df_batch = read_haplotype(batch_path, 'haplotype_final')
    edges_df_batch = pd.read_csv(batch_path + 'edges.csv', sep=',', header=0)
    edges_df_batch['source'] = edges_df_batch['source'].apply(lambda x: int(x))
    edges_df_batch['target'] = edges_df_batch['target'].apply(lambda x: int(x))
    # results_voi_batch = pd.read_csv(batch_path + 'results_voi.txt', sep='\t', header=None)
    if os.path.getsize(batch_path + 'results_voi.txt') == 0:
        results_voi_batch = pd.DataFrame()
    else:
        results_voi_batch = pd.read_csv(batch_path + 'results_voi.txt', sep='\t', header=None)
        results_voi_batch.columns = ['target', 'source', 'source_orig']

    # if os.path.getsize(batch_path + 'recom_list.txt') == 0:
    #     recom_batch = pd.DataFrame()
    # else:
    #     recom_batch = pd.read_csv(batch_path + 'recom_list.txt', sep='\t', header=None)
    #     recom_batch.columns = ['idx', 'source', 'target', 'acceptors', 'bpidx', 'class']

    # if os.path.getsize(batch_path + 'recom_list_time.txt') == 0:
    #     recom_time_batch = pd.DataFrame()
    # else:
    #     recom_time_batch = pd.read_csv(batch_path + 'recom_list_time.txt', sep='\t', header=None)
    #     recom_time_batch.columns = ['idx', 'source', 'target', 'acceptors', 'bpidx', 'class']

    # hap_index_df = pd.read_csv(batch_path + 'hap_index_df.txt',sep='\t', header=0)
    haplotype_df_batch_raw = read_haplotype_raw(batch_path_raw, 'haplotype_df')
    haplotype_df_previous = haplotype_df_updated.copy()
    #merge the existing haps into haplotype_df_previous and change it self
    print('process the haps of batch appear in previous')
    df1 = pd.DataFrame()
    df1['hapid'] = haplotype_df_previous.index
    df1['haplotype'] = haplotype_df_previous['haplotype']
    df2 = pd.DataFrame()
    df2['hapid'] = haplotype_df_batch_raw.index
    df2['haplotype'] = haplotype_df_batch_raw['haplotype']
    haplotype_df_merge = pd.merge(df1, df2, how='inner', on='haplotype')
    if haplotype_df_merge.shape[0] == 0: #there is no samples in existing haps
        pass
    else:
        for each_hap_info in haplotype_df_merge.values:
            hapid_x, haplotype, hapid_y = each_hap_info
            haplotype_df_previous.loc[hapid_x, 'sample_names'] += ';' + haplotype_df_batch_raw.loc[hapid_y, 'sample_names']
            haplotype_df_previous.loc[hapid_x, 'sample_count'] += haplotype_df_batch_raw.loc[hapid_y, 'sample_count']
            haplotype_df_previous.loc[hapid_x, 'sample_locations'] = ';'.join(set(haplotype_df_previous.loc[hapid_x, 'sample_locations'].split(';')+ haplotype_df_batch_raw.loc[hapid_y, 'sample_locations'].split(';')))
            haplotype_df_previous.loc[hapid_x, 'sample_time_min'] = min( haplotype_df_previous.loc[hapid_x, 'sample_time_min'], haplotype_df_batch_raw.loc[hapid_y, 'sample_time_min'])
            haplotype_df_previous.loc[hapid_x, 'sample_time_max'] = max( haplotype_df_previous.loc[hapid_x, 'sample_time_max'], haplotype_df_batch_raw.loc[hapid_y, 'sample_time_max'])

        # haplotype_removed_indexs = [x for x in haplotype_df_batch_raw.index.tolist() if x not in haplotype_df_merge['hapid_y'].tolist()]
    df1 = pd.DataFrame()
    df1['hapid'] = haplotype_df_previous.index
    df1['haplotype'] = haplotype_df_previous['haplotype']
    df2 = pd.DataFrame()
    # selected_hapids = edges_df_batch['target'].tolist()
    df2['hapid'] = haplotype_df_batch.index
    # df2['haplotype'] = haplotype_df_batch.loc[selected_hapids]['haplotype'].tolist()
    df2['haplotype'] = haplotype_df_batch['haplotype'].tolist()
    haplotype_df_merge = pd.merge(df1, df2, how='outer', on='haplotype')
    merged_df1 = haplotype_df_merge[haplotype_df_merge['hapid_x'].notna() & (haplotype_df_merge['hapid_y'].notna())]
    merged_df2 = haplotype_df_merge[haplotype_df_merge['hapid_x'].isna()]
    # merged_df3 = haplotype_df_merge[haplotype_df_merge['hapid_y'].isna()]
    # merged_df4 = haplotype_df_merge[haplotype_df_merge['hapid_y'].notna()]
    # merged_df2['hapid_y'] = merged_df2['hapid_y'].astype('int8')
    # hapidx_batch = merged_df2['hapid_y'].tolist()
    # merged_df4 = pd.merge( merged_df2['hapid_y'], edges_df_batch[['source','target']], how='inner', left_on='hapid_y', right_on='target' )
    # selected_hapids = merged_df4['target'].tolist()

    # first_idx = int(edges_df_batch.iloc[0]['target'])
    # if haplotype_df_batch.loc[first_idx]['haplotype'] == '':
    #     selected_hapids = list(range(first_idx+1, first_idx+edges_df_batch.shape[0]))
    # else:
    # batchnew_hapids = list(range(first_idx, first_idx+edges_df_batch.shape[0]))
    batchnew_hapids = edges_df_batch['target'].tolist()
    # batchnew_hapids = merged_df2['hapid_y'].tolist
    # overlap_list = list(set(merged_df1['hapid_y'].tolist()).intersection(set(selected_hapids)))
    selected_hapids = list(set(batchnew_hapids).difference(set(merged_df1['hapid_y'].tolist())))
    # selected_hapids = merged_df2['hapid_y'].tolist()
    # selected_hapids = edges_df_batch['target'].tolist()
    haplotype_df_new = haplotype_df_batch.loc[selected_hapids]
    haplotype_df_new_index = list(range(haplotype_df_previous.shape[0], haplotype_df_previous.shape[0] + len(selected_hapids)) )
    haplotype_df_new.index = haplotype_df_new_index
    haplotype_df_updated = pd.concat([ haplotype_df_previous, haplotype_df_new ], axis=0)
    # tmpiddf = pd.DataFrame()
    # tmpiddf['hapid_x'] = merged_df1['hapid_x'].tolist() + haplotype_df_new_index
    # tmpiddf['hapid_y'] = merged_df1['hapid_y'].tolist() + ['n'+str(x) for x in selected_hapids]
    dict1 = {}
    for idx in merged_df1.index.tolist():
        hapid1, hap, hapid2 = merged_df1.loc[idx]
        # print(hapid1, hapid2)
        dict1[int(hapid2)] = int(hapid1)

    # dict2= {}
    for idx in range(len(selected_hapids)):
        dict1[selected_hapids[idx]] = haplotype_df_new_index[idx]

    print('process the edges...')
    edges_df_batch_new = []
    for each_idx in range( edges_df_batch.shape[0]):
        # each_idx = 37
        # print('processing ', each_idx)
        source, target, mutations, dis = edges_df_batch.iloc[each_idx]
        # if source in dict1.keys() and target in dict1.keys():
        source = int(source)
        target = int(target)
        # if target in selected_hapids:
            # if source not in dict1.keys():
            #     print(each_idx,source)
        source_new = int(dict1[source])
        target_new = int(dict1[target])
        edges_df_batch_new.append([source_new, target_new, mutations, dis])

    edges_df_batch_new = pd.DataFrame(edges_df_batch_new)
    edges_df_batch_new.columns = ['source', 'target', 'mutations', 'distance']
    edges_df_batch_new = edges_df_batch_new.drop_duplicates()
    edges_df = pd.concat([edges_df, edges_df_batch_new], axis=0)
    edges_df.index = range(edges_df.shape[0]) 

    #process the results voi
    print('#process the results voi')
    if results_voi_batch.shape[0]>0:
        results_voi_batch_new = []
        for each_idx in range(results_voi_batch.shape[0]):
            target, source, source_orig = results_voi_batch.iloc[each_idx]
            target = int(target)
            source = int(source)
            source_orig = int(source_orig)
            if target in dict1.keys() and source in dict1.keys() and source_orig in dict1.keys():
                source_new = int(dict1[source])
                target_new = int(dict1[target])
                source_orig_new = int(dict1[source_orig])
                # target_new = tmpiddf[ tmpiddf['hapid_y']==target ].iloc[0]['hapid_x']
                # source_new = tmpiddf[ tmpiddf['hapid_y']==source ].iloc[0]['hapid_x']
                # source_orig_new = tmpiddf[ tmpiddf['hapid_y']==source_orig ].iloc[0]['hapid_x']
                results_voi_batch_new.append([target_new, source_new, source_orig_new])
        results_voi_batch_new = pd.DataFrame(results_voi_batch_new)
        results_voi_batch_new.columns = ['target', 'source', 'source_orig']
        results_voi = pd.concat([results_voi, results_voi_batch_new ], axis=0)

print('saving the files...')
"""save files"""
haplotype_df_updated.to_csv(output_path+'haplotype_final.txt',sep='\t', header=True, index=None)
edges_df = edges_df.drop_duplicates()
edges_df.to_csv(output_path+'edges.csv',sep=',', index=None , header=True)

results_voi = pd.concat([results_voi_base, results_voi], axis=0)
results_voi.to_csv(output_path+'results_voi.txt',sep='\t', index=None , header=True)


end_time = time.time()
diff = end_time - beg_time
print(time.ctime())
print(diff, ' s')
print(diff/60, ' min')
print(diff/60/60, ' h')
