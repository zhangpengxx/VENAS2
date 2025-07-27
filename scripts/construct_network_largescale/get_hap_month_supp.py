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
import scipy.sparse
from scipy.sparse import coo_matrix
from collections import Counter

beg_time = time.time()
parser = argparse.ArgumentParser()
parser.add_argument('-input_name', '--input_name', dest='input_name', help="the name of input samples mutation file", required=True)
# parser.add_argument('-output_path', '--output_path', dest='output_path', help="the output path of the results", default= 'venas2_output')
# parser.add_argument('-batch_mode', '--batch_mode', dest='batch_mode', help="the number of input samples", default=False)
# parser.add_argument('-initial_size', '--initial_size', dest='initial_size', help="the samples count used to construct the initial network when using batch mode", default=10000)
# parser.add_argument('-batch_size', '--batch_size', dest='batch_size', help="the samples count which added to the network each time when using batch mode", default=10000)
args = parser.parse_args()

# #parameters
input_name = args.input_name
batch_name_list = ['2020-01', '2020-02', '2020-03', '2020-04', '2020-05', '2020-06', '2020-07', '2020-08', '2020-09', '2020-10', '2020-11', '2020-12', '2021-01', '2021-02', '2021-03', '2021-04', '2021-05', '2021-06', '2021-07', '2021-08', '2021-09', '2021-10', '2021-11', '2021-12', '2022-01', '2022-02', '2022-03', '2022-04', '2022-05', '2022-06', '2022-07', '2022-08', '2022-09', '2022-10', '2022-11', '2022-12', '2023-01', '2023-02', '2023-03', '2023-04', '2023-05', '2023-06', '2023-07', '2023-08', '2023-09', '2023-10', '2023-11', '2023-12', '2024-01', '2024-02', '2024-03', '2024-04', '2024-05', '2024-06', '2024-07', '2024-08', '2024-09', '2024-10', '2024-11']

# dict_nc2num = {'a':1, 't':2, 'c':3, 'g':4, '-':5, 'o':6}
dict_nc2num = {'a':1, 't':2, 'c':3, 'g':4, '-':5, 'o':6,'A':1, 'T':2, 'C':3, 'G':4}
dict_num2nc = {1:'a', 2:'t', 3:'c', 4:'g', 5:'-', 6:'o'}
data_path = 'datasets/'+input_name +'/'
output_path = data_path
data = pd.read_csv(data_path + 'samples_mutations.txt', sep='\t', header=0)
# data.columns = ['samplename', 'mutations' , 'collection_date', 'location' ,'time']
# data = data[['samplename', 'mutations' , 'collection_date', 'location']]
data = data.fillna('')
data.columns = ['samplename', 'haplotype' , 'collection_date', 'location']
data['time'] = data['collection_date'].apply(lambda x: datetime.strptime(x, '%Y-%m-%d'))
data = data.sort_values('time', ascending=True)

data_count = data.shape[0]
if data_count > 10000:
    # threshold = 2
    threshold = data_count * 0.001
else:
    threshold = 7

data_haps = []
for each_haplotype, group in data.groupby('haplotype'):
    # if each_haplotype == 'NaN':
    # print(each_haplotype)
    sample_names = group['samplename'].tolist()
    sample_count = group.shape[0]
    # sample_count = len(sample_names)
    sample_names = ';'.join(sample_names)
    sample_locations = ';'.join(group['location'].drop_duplicates().tolist())
    sample_times = group['time'].drop_duplicates().tolist()
    # sample_times = [ datetime.strptime(x, '%Y-%m-%d') for x in sample_times]
    sample_times_min = min(sample_times)
    sample_times_max = max(sample_times)
    data_haps.append([each_haplotype, sample_names, sample_count, sample_locations, sample_times_min, sample_times_max])

data_haps = pd.DataFrame(data_haps)
data_haps.columns = ['haplotype', 'sample_names','sample_count',  'sample_locations', 'sample_time_min', 'sample_time_max']

# data_haps = data['haplotype'].drop_duplicates()
all_muts = []
# for sample_idx in range(100):
for eachidx in range(data_haps.shape[0]):
    cur_hap = data_haps.iloc[eachidx]['haplotype']
    if cur_hap == '':
        continue
    cur_hap = cur_hap.split(';')
    all_muts += cur_hap

dict_muts = Counter(all_muts)
muts_df = pd.DataFrame()
muts_df['mut'] = dict_muts.keys()
muts_df['ac'] = dict_muts.values()
muts_df_selected = muts_df[muts_df['ac']> threshold]
muts_selected = set(muts_df_selected['mut'].tolist())

hap_pos_selected = []
for eachmut in muts_selected:
    if 'del' in eachmut:
        hap_pos_selected+= list(range( int(eachmut.split('_')[0][3:]), int( eachmut.split('_')[1] )+1 ))
    else:
        hap_pos_selected.append(int(eachmut[1:-1]))

hap_pos_selected = sorted(hap_pos_selected)

haplotype_df = []
for sample_idx in range(data_haps.shape[0]):
    # sample_idx =3453
    cur_hap, cur_sample_names, cur_sample_count, cur_locations, cur_time_min, cur_time_max = data_haps.iloc[sample_idx]
    if cur_hap == '':
        hap_new = ''
        hap_count = 0
        hap_pos = []
        hap_value = []
    else:
        cur_hap = cur_hap.split(';')
        cur_hap_new = muts_selected.intersection(set(cur_hap))
        cur_hap_pos = []
        cur_hap_value = []
        for eachmut in cur_hap_new:
            if 'del' in eachmut:
                cur_hap_pos.append( int(eachmut.split('_')[0][3:]) )
                cur_hap_value.append(5)
            else:
                cur_hap_pos.append( int(eachmut[1:-1]) )
                cur_hap_value.append( dict_nc2num[ eachmut[-1:] ])
        tmpdf = pd.DataFrame()
        tmpdf['muts'] = list(cur_hap_new)
        tmpdf['pos'] = cur_hap_pos
        tmpdf['value'] = cur_hap_value
        tmpdf = tmpdf.sort_values('pos') #tmpdf maybe 0
        hap_new = ';'.join(tmpdf['muts'].tolist())
        hap_pos = tmpdf['pos'].tolist()
        hap_value = tmpdf['value'].tolist()
        # hap_pos = ';'.join( [str(x) for x in tmpdf['pos'].tolist() ])
        # hap_value = ';'.join([ str(x) for x in tmpdf['value'].tolist() ])
        hap_count =tmpdf.shape[0]
    haplotype_df.append([hap_new, hap_count, cur_sample_names, cur_sample_count, cur_locations, cur_time_min, cur_time_max, hap_pos, hap_value] )

haplotype_df = pd.DataFrame(haplotype_df)
haplotype_df.columns = ['haplotype', 'muts_count', 'sample_names', 'sample_count',  'sample_locations', 'sample_time_min', 'sample_time_max', 'hap_pos', 'hap_value']
# haplotype_df = haplotype_df.sort_values([ 'sample_time_min', 'muts_count','sample_count'], ascending=[True, True, False])
# haplotype_df = pd.concat( [haplotype_df_first, haplotype_df_rest], axis=0 )
haplotype_df.index= list(range(haplotype_df.shape[0]))
# haplotype_df2['haplotype'].drop_duplicates()
# haplotype_df2['hap_pos'].drop_duplicates()
# data_haps['haplotype'].drop_duplicates()

haplotype_df2 = []
for each_haplotype, group in haplotype_df.groupby('haplotype'):
    # if each_haplotype == 'NaN':
    if group.shape[0]>1:
        # print(each_haplotype, group)
        cur_sample_names = ';'.join(group['sample_names'].tolist())
        cur_sample_names = list(set(cur_sample_names.split(';')))
        cur_sample_count = len(cur_sample_names)
        cur_sample_names = ';'.join(cur_sample_names)
        cur_locations = ';'.join(set(';'.join(group['sample_locations'].tolist()).split(';')))
        # sample_times = group['time'].drop_duplicates().tolist()
        # sample_times = [ datetime.strptime(x, '%Y-%m-%d') for x in sample_times]
        cur_time_min = min(group['sample_time_min'])
        cur_time_max = max(group['sample_time_max'])
    else:
        cur_sample_names = group.iloc[0]['sample_names']
        cur_sample_count = group.iloc[0]['sample_count']
        cur_locations = group.iloc[0]['sample_locations']
        cur_time_min = group.iloc[0]['sample_time_min']
        cur_time_max = group.iloc[0]['sample_time_max']
    hap_count = group.iloc[0]['muts_count']
    hap_pos = group.iloc[0]['hap_pos']
    hap_value = group.iloc[0]['hap_value']
    haplotype_df2.append([each_haplotype, hap_count, cur_sample_names, cur_sample_count, cur_locations, cur_time_min, cur_time_max, hap_pos, hap_value] )

haplotype_df2 = pd.DataFrame(haplotype_df2)
haplotype_df2.columns = ['haplotype', 'muts_count', 'sample_names', 'sample_count',  'sample_locations', 'sample_time_min', 'sample_time_max', 'hap_pos', 'hap_value']
haplotype_df2 = haplotype_df2.sort_values([ 'sample_time_min', 'muts_count','sample_count'], ascending=[True, True, False])
haplotype_df2.index= list(range(haplotype_df2.shape[0]))

# # haplotype_df2[haplotype_df2['hap_pos']==haplotype_df2.loc[143]['hap_pos']]
# samples_all = []
# for eachidx in range(haplotype_df2.shape[0]):
#     # eachidx = 0
#     sample_names = haplotype_df2.loc[eachidx]['sample_names']
#     samples_all+=sample_names.split(';')
    # if haplotype_df2.loc[eachidx]['hap_pos']==[241, 3037, 14408, 23403, 26801, 28881, 28882, 28883]:
        # print(eachidx)

hap_count_fn = haplotype_df2.shape[0]
hap_pos_fn = list(sorted(set(hap_pos_selected)))
##convert to array
hap_pos_count = len(hap_pos_fn)
dict_pos_index = dict([ (x[1],x[0]) for x in enumerate(hap_pos_fn) ])
rows = []
cols = []
values = []
for each_hapidx in range(hap_count_fn):
    # each_hapidx = 1101
    print('processing ',each_hapidx)
    cur_hap_pos = haplotype_df2.loc[each_hapidx]['hap_pos']
    if cur_hap_pos ==[]:
        continue
    tmpdf = pd.DataFrame(cur_hap_pos)
    tmpdf['posidx'] = tmpdf.apply(lambda x:dict_pos_index[x[0]], axis=1)
    cur_hap_value = haplotype_df2.loc[each_hapidx]['hap_value']
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
# return haplotype_df_fn, hap_pos_fn, rows, cols, values


"""step1: get the haplotype vector form and update the latest hap"""
# haplotype_df_base, hap_pos_base, rows_base, cols_base, values_base = convert_smaples_mutation_into_hap_df_array(data)
# coo = coo_matrix((values_base, (rows_base, cols_base)), shape=(haplotype_df_base.shape[0], len(hap_pos_base)), dtype=np.int8)
# # hap_array_df = pd.DataFrame(coo.todense())
# # hap_array_df.columns = hap_pos_fn
# hap_array_base = coo.toarray()
np.save( output_path+'rows1.npy', rows)
np.save( output_path+'cols1.npy', cols)
np.save( output_path+'values1.npy', values)
haplotype_df2.to_csv(output_path+'haplotype_df.txt',sep='\t', header=True, index=None)
pd.DataFrame(hap_pos_fn).to_csv(output_path+'hap_pos.txt',sep='\t', header=False, index=None)


end_time = time.time()
diff = end_time - beg_time
print(time.ctime())
print(diff, ' s')
print(diff/60, ' min')
print(diff/60/60, ' h')

