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
from scipy.sparse import coo_matrix
import multiprocessing

beg_time = time.time()
parser = argparse.ArgumentParser()
parser.add_argument('-input_name', '--input_name', dest='input_name', help="the name of input samples mutation file", required=True)
parser.add_argument('-cpu_count', '--cpu_count', dest='cpu_count', help="the cpu count used", default=5)
# # parser.add_argument('-output_path', '--output_path', dest='output_path', help="the output path of the results", default= 'venas2_output')
# # parser.add_argument('-batch_mode', '--batch_mode', dest='batch_mode', help="the number of input samples", default=False)
# # parser.add_argument('-initial_size', '--initial_size', dest='initial_size', help="the samples count used to construct the initial network when using batch mode", default=10000)
# # parser.add_argument('-batch_size', '--batch_size', dest='batch_size', help="the samples count which added to the network each time when using batch mode", default=10000)
args = parser.parse_args()

# #parameters
input_name = args.input_name
# batch_mode = args.batch_mode
# initial_size = int(args.initial_size)
# batch_size = int(args.batch_size)
# input_name = '2020-06'

batch_name_list = ['2020-01', '2020-02', '2020-03', '2020-04', '2020-05', '2020-06', '2020-07', '2020-08', '2020-09', '2020-10', '2020-11', '2020-12', '2021-01', '2021-02', '2021-03', '2021-04', '2021-05', '2021-06', '2021-07', '2021-08', '2021-09', '2021-10', '2021-11', '2021-12', '2022-01', '2022-02', '2022-03', '2022-04', '2022-05', '2022-06', '2022-07', '2022-08', '2022-09', '2022-10', '2022-11', '2022-12', '2023-01', '2023-02', '2023-03', '2023-04', '2023-05', '2023-06', '2023-07', '2023-08', '2023-09', '2023-10', '2023-11', '2023-12', '2024-01', '2024-02', '2024-03', '2024-04', '2024-05', '2024-06', '2024-07', '2024-08', '2024-09', '2024-10', '2024-11']

assert input_name in batch_name_list
batch_idx = batch_name_list.index(input_name)
##load in the data
dict_nc2num = {'a':1, 't':2, 'c':3, 'g':4, '-':5, 'o':6}
dict_num2nc = {1:'a', 2:'t', 3:'c', 4:'g', 5:'-', 6:'o'}
##can_lens control the number of top haps considered for each new hap in faiss
can_lens = 500

# month = batch_name_list[batch_idx]
month = input_name
month1 = batch_name_list[batch_idx-1]
month2 = batch_name_list[batch_idx-2]
month3 = batch_name_list[batch_idx-3]
month4 = batch_name_list[batch_idx-4]
# month5 = batch_name_list[batch_idx-5]
# month6 = batch_name_list[batch_idx-6]
# print('processing ', month, '\n', 'previous months including: ', month6, month5, month4, month3, month2, month1)
print('processing ', month, '\n', 'previous months including: ', month4, month3, month2, month1)

# data_path =  'results/venas2v1_output/haps_month4/'
# output_path = 'results/venas2v1_output/venas2_output_month23/' + str(month) + '/'
data_path =  'datasets/'
output_path = 'datasets/' + str(month) + '/'
if os.path.exists(output_path):
    print('the folder exists! please check the output dir')
else:
    os.makedirs(output_path)

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
        edge_muts = 'tmp'
        t1_count = t2_count = 0
    edge_diff_count = t1_count + t2_count
    return edge_muts, edge_diff_count

def get_haps_info(month_fn):
    # month_fn = month
    month_path_fn = data_path +str(month_fn) + '/'
    haplotype_df = pd.read_csv( month_path_fn + 'haplotype_df.txt', sep='\t', header=0)
    haplotype_df = haplotype_df.fillna('')
    haplotype_df.columns = ['haplotype', 'muts_count', 'sample_names', 'sample_count','sample_locations', 'sample_time_min', 'sample_time_max', 'hap_pos', 'hap_value']
    haplotype_df['sample_time_min'] = haplotype_df['sample_time_min'].apply(lambda x:  datetime.strptime(x, '%Y-%m-%d'))
    haplotype_df['sample_time_max'] = haplotype_df['sample_time_max'].apply(lambda x:  datetime.strptime(x, '%Y-%m-%d'))
    hap_count_fn = haplotype_df.shape[0]
    hap_pos_df_fn = pd.read_csv( month_path_fn + 'hap_pos.txt', sep='\t', header=None)
    hap_pos_fn = hap_pos_df_fn[0].tolist()
    hap_pos_count = hap_pos_df_fn.shape[0]
    #hap array
    rows_fn = np.load( month_path_fn+'rows1.npy')
    cols_fn = np.load( month_path_fn+'cols1.npy')
    values_fn = np.load( month_path_fn+'values1.npy')
    coo_fn = coo_matrix((values_fn, (rows_fn, cols_fn)), shape=(hap_count_fn, hap_pos_count), dtype=np.int8)
    # coo_fn = coo_matrix((values_fn, (rows_fn, cols_fn)), shape=(hap_count_fn, 29903), dtype=np.int8)
    hap_array_fn = coo_fn.toarray()
    hap_array_df = pd.DataFrame(hap_array_fn)
    hap_array_df.columns = hap_pos_fn
    # hap_array_df.columns = range(29903)
    # hap_array_df = hap_array_df[hap_pos_df_fn]
    # if haplotype_df.iloc[0]['haplotype'] =='':
    if '' in haplotype_df['haplotype'].tolist():
        # idx = haplotype_df[haplotype_df['haplotype']==''].index[0]
        # haplotype_df1 = haplotype_df[haplotype_df['haplotype'].isna()]
        # haplotype_df2 = haplotype_df[haplotype_df['haplotype'].notna()]
        haplotype_df1 = haplotype_df[haplotype_df['haplotype']=='']
        haplotype_df2 = haplotype_df[haplotype_df['haplotype']!='']
        haplotype_df2['hap_pos'] = haplotype_df2['hap_pos'].apply(lambda x: [int(t) for t in x[1:-1].split(',')])
        haplotype_df2['hap_value'] = haplotype_df2['hap_value'].apply(lambda x: [int(t) for t in x[1:-1].split(',')])
        haplotype_df = pd.concat([haplotype_df1, haplotype_df2], axis=0)
        haps_list =haplotype_df.index.tolist()
        haplotype_df.index = range(haplotype_df.shape[0])
        hap_array_df = hap_array_df.loc[haps_list]
        hap_array_df.index = range(hap_array_df.shape[0])
    else:
        haplotype_df['hap_pos'] = haplotype_df['hap_pos'].apply(lambda x: [int(t) for t in x[1:-1].split(',')])
        haplotype_df['hap_value'] = haplotype_df['hap_value'].apply(lambda x: [int(t) for t in x[1:-1].split(',')])
    # hap_array_df = pd.read_csv( month_path_fn + 'hap_array_df.txt', sep='\t', header=0)
    # hap_pos = [ int(x) for x in hap_array_df.columns.tolist()]
    # hap_array_df.columns = range(29903)
    return haplotype_df, hap_array_df, hap_pos_fn

haplotype_df_month, hap_array_df_month, hap_pos_month = get_haps_info(month)

def get_haps_info2(month_fn):
    # month_fn = month
    month_path_fn = data_path +str(month_fn) + '/'
    haplotype_df = pd.read_csv( month_path_fn + 'haplotype_df.txt', sep='\t', header=0)
    haplotype_df = haplotype_df.fillna('')
    haplotype_df.columns = ['haplotype', 'muts_count', 'sample_names', 'sample_count','sample_locations', 'sample_time_min', 'sample_time_max', 'hap_pos', 'hap_value']
    haplotype_df['sample_time_min'] = haplotype_df['sample_time_min'].apply(lambda x:  datetime.strptime(x, '%Y-%m-%d'))
    haplotype_df['sample_time_max'] = haplotype_df['sample_time_max'].apply(lambda x:  datetime.strptime(x, '%Y-%m-%d'))
    hap_count_fn = haplotype_df.shape[0]
    hap_pos_df_fn = pd.read_csv( month_path_fn + 'hap_pos.txt', sep='\t', header=None)
    hap_pos_fn = hap_pos_df_fn[0].tolist()
    hap_pos_count = hap_pos_df_fn.shape[0]
    #hap array
    rows_fn = np.load( month_path_fn+'rows1.npy')
    cols_fn = np.load( month_path_fn+'cols1.npy')
    values_fn = np.load( month_path_fn+'values1.npy')
    return haplotype_df, hap_pos_fn, rows_fn, cols_fn, values_fn

haplotype_df_month1, hap_pos_month1, rows1, cols1, values1 = get_haps_info2(month1)
haplotype_df_month2, hap_pos_month2, rows2, cols2, values2 = get_haps_info2(month2)
haplotype_df_month3, hap_pos_month3, rows3, cols3, values3 = get_haps_info2(month3)
haplotype_df_month4, hap_pos_month4, rows4, cols4, values4 = get_haps_info2(month4)

tmpdf= pd.concat([haplotype_df_month1, haplotype_df_month2, haplotype_df_month3, haplotype_df_month4], axis=0)
tmpdf['haplotype'].drop_duplicates()

#change the cols
hap_pos_previous = sorted(set(hap_pos_month4 + hap_pos_month3 + hap_pos_month2 + hap_pos_month1))
hap_pos_df = pd.DataFrame(range(len(hap_pos_previous)))
hap_pos_df.index = hap_pos_previous
# hap_pos_df['pos'] = hap_pos_previous
# hap_pos_df['idx'] = range(len(hap_pos_previous))
# posmonth1 = pd.DataFrame(hap_pos_month1)
# posmonth1['pos'] = hap_pos_month1
# posmonth1['idx'] = range(len(hap_pos_month1))
hap_pos_df1=pd.DataFrame(hap_pos_month1)
cols1_new = np.array(hap_pos_df.loc[hap_pos_df1.loc[list(cols1)][0].tolist() ]).reshape(-1,)
hap_pos_df2=pd.DataFrame(hap_pos_month2)
cols2_new = np.array(hap_pos_df.loc[hap_pos_df2.loc[list(cols2)][0].tolist() ]).reshape(-1,)
hap_pos_df3=pd.DataFrame(hap_pos_month3)
cols3_new = np.array(hap_pos_df.loc[hap_pos_df3.loc[list(cols3)][0].tolist() ]).reshape(-1,)
hap_pos_df4=pd.DataFrame(hap_pos_month4)
cols4_new = np.array(hap_pos_df.loc[hap_pos_df4.loc[list(cols4)][0].tolist() ]).reshape(-1,)

coo_df1=pd.DataFrame() 
coo_df1['row'] = rows1
coo_df1['col'] = cols1_new
coo_df1['value'] = values1
coo_df2=pd.DataFrame()
coo_df2['row'] = rows2
coo_df2['col'] = cols2_new
coo_df2['value'] = values2
coo_df3=pd.DataFrame() 
coo_df3['row'] = rows3
coo_df3['col'] = cols3_new
coo_df3['value'] = values3
coo_df4=pd.DataFrame() 
coo_df4['row'] = rows4
coo_df4['col'] = cols4_new
coo_df4['value'] = values4

tmpdf = pd.concat([haplotype_df_month4[['haplotype']], haplotype_df_month3[['haplotype']], haplotype_df_month2[['haplotype']], haplotype_df_month1[['haplotype']]], axis=0)
tmpdf = tmpdf['haplotype'].drop_duplicates()
tmpdf = pd.DataFrame(tmpdf)
hapmonth1 = pd.DataFrame()
hapmonth1['haplotype'] = haplotype_df_month1['haplotype']
hapmonth1['idx'] = haplotype_df_month1.index.tolist()
hapmonth2 = pd.DataFrame()
hapmonth2['haplotype'] = haplotype_df_month2['haplotype']
hapmonth2['idx'] = haplotype_df_month2.index.tolist()
hapmonth3 = pd.DataFrame()
hapmonth3['haplotype'] = haplotype_df_month3['haplotype']
hapmonth3['idx'] = haplotype_df_month3.index.tolist()
hapmonth4 = pd.DataFrame()
hapmonth4['haplotype'] = haplotype_df_month4['haplotype']
hapmonth4['idx'] = haplotype_df_month4.index.tolist()
merged_df1=  pd.merge(tmpdf, hapmonth1 , how='left', on='haplotype')
merged_df1.columns = ['haplotype', 'month1']
merged_df2=  pd.merge(merged_df1, hapmonth2 , how='left', on='haplotype')
merged_df2.columns = ['haplotype', 'month1', 'month2']
merged_df3=  pd.merge(merged_df2, hapmonth3 , how='left', on='haplotype')
merged_df3.columns = ['haplotype', 'month1', 'month2', 'month3']
merged_df4=  pd.merge(merged_df3, hapmonth4 , how='left', on='haplotype')
merged_df4.columns = ['haplotype', 'month1', 'month2', 'month3', 'month4']
merged_df4 = merged_df4.fillna('')
del merged_df1
del merged_df2
del merged_df3

def merge_haps(hap_ids_batch, merged_df_fn, haplotype_df_month1,haplotype_df_month2,haplotype_df_month3,haplotype_df_month4, coo_df1, coo_df2, coo_df3,coo_df4, batch_name ):
    print('start the task in the ', batch_name,'th cpu ')
    # hap_ids_batch = hap_ids
    # merged_df_fn = merged_df4
    haplotype_df_fn = []
    coo_df_fn = pd.DataFrame(columns=['row', 'col', 'value'])
    # for each_idx in range(merged_df_fn.shape[0]):
    for each_idx in hap_ids_batch:
        # each_idx = 22252
        # print('processing ', each_idx, merged_df4.shape)
        hap, idx_month1, idx_month2, idx_month3, idx_month4 = merged_df_fn.iloc[each_idx]
        idx_list = [idx_month1, idx_month2, idx_month3, idx_month4]
        # idx_df= pd.DataFrame(idx_list)
        # idx_select = idx_df[idx_df[0]!=''].index.tolist()
        if hap == '':
            mut_count, hap_pos, hap_value = 0, '', ''
        else:
            if idx_month1 !='':
                mut_count, hap_pos, hap_value = haplotype_df_month1.iloc[ int(idx_month1)][['muts_count', 'hap_pos', 'hap_value']]
                tmpdf = coo_df1[coo_df1['row']==idx_month1]
                tmpdf.loc[:,'row']=each_idx
                coo_df_fn = pd.concat([coo_df_fn, tmpdf],axis=0)
            else:
                if idx_month2 != '':
                    mut_count, hap_pos, hap_value = haplotype_df_month2.iloc[ int(idx_month2)][['muts_count', 'hap_pos', 'hap_value']]
                    tmpdf = coo_df2[coo_df2['row']==idx_month2]
                    tmpdf.loc[:,'row']=each_idx
                    coo_df_fn = pd.concat([coo_df_fn, tmpdf],axis=0)
                else:
                    if idx_month3 !='':
                        mut_count, hap_pos, hap_value = haplotype_df_month3.iloc[ int(idx_month3)][['muts_count', 'hap_pos', 'hap_value']]
                        tmpdf = coo_df3[coo_df3['row']==idx_month3]
                        tmpdf.loc[:,'row']=each_idx
                        coo_df_fn = pd.concat([coo_df_fn, tmpdf],axis=0)
                    else:
                        if idx_month4 !='':
                            mut_count, hap_pos, hap_value = haplotype_df_month4.iloc[ int(idx_month4)][['muts_count', 'hap_pos', 'hap_value']]
                            tmpdf = coo_df4[coo_df4['row']==idx_month4]
                            tmpdf.loc[:,'row']=each_idx
                            coo_df_fn = pd.concat([coo_df_fn, tmpdf],axis=0)
            # hap_pos = [int(t) for t in hap_pos[1:-1].split(',')]
            # hap_value = [int(t) for t in hap_value[1:-1].split(',')]
        sample_names = []
        # sample_count = 0
        location = []
        time_min = []
        time_max = []
        for idx in idx_list:
            if idx_month1!='':
                # sample_names.append(haplotype_df_month1.iloc[ int(idx_month1)]['sample_names'].split(';'))
                sample_names1, sample_count1, location1, time_min1, time_max1 = haplotype_df_month1.iloc[ int(idx_month1)][['sample_names', 'sample_count','sample_locations', 'sample_time_min', 'sample_time_max']]
                sample_names.append(sample_names1)
                # sample_count += sample_count1
                location += location1.split(';')
                time_min.append(time_min1)
                time_max.append(time_max1) 
            if idx_month2!='':
                sample_names2, sample_count2, location2, time_min2, time_max2 = haplotype_df_month2.iloc[ int(idx_month2)][['sample_names', 'sample_count','sample_locations', 'sample_time_min', 'sample_time_max']]
                sample_names.append(sample_names2)
                # sample_count += sample_count2
                location += location2.split(';')
                time_min.append(time_min2)
                time_max.append(time_max2) 
            if idx_month3!='':
                sample_names3, sample_count3, location3, time_min3, time_max3 = haplotype_df_month3.iloc[ int(idx_month3)][['sample_names', 'sample_count','sample_locations', 'sample_time_min', 'sample_time_max']]
                sample_names.append(sample_names3)
                # sample_count += sample_count3
                location += location3.split(';')
                time_min.append(time_min3)
                time_max.append(time_max3) 
            if idx_month4!='':
                sample_names4, sample_count4, location4, time_min4, time_max4 = haplotype_df_month4.iloc[ int(idx_month4)][['sample_names', 'sample_count','sample_locations', 'sample_time_min', 'sample_time_max']]
                sample_names.append(sample_names4)
                # sample_count += sample_count4
                location += location4.split(';')
                time_min.append(time_min4)
                time_max.append(time_max4) 
        sample_names = list(set(sample_names))
        sample_count = len(sample_names)
        sample_names = ';'.join(sample_names)
        location = ';'.join(set(location))
        time_min = min(time_min)
        time_max = max(time_max)
        haplotype_df_fn.append([hap, mut_count, sample_names, sample_count, location, time_min, time_max, hap_pos, hap_value])
    haplotype_df_fn = pd.DataFrame(haplotype_df_fn)
    haplotype_df_fn.columns =['haplotype', 'muts_count', 'sample_names', 'sample_count','sample_locations', 'sample_time_min', 'sample_time_max', 'hap_pos', 'hap_value']
    haplotype_df_fn.to_csv(output_path+'haplotype_df'+str(batch_name)+'.txt',sep='\t', header=True, index=None)
    coo_df_fn.to_csv(output_path+'coo_df'+str(batch_name)+'.txt',sep='\t', header=True, index=None)
    prinnt(batch_name, ' is over!')

hap_count = merged_df4.shape[0]
hap_ids = list(range(hap_count))
# cpu_count = multiprocessing.cpu_count()
cpu_count = int(args.cpu_count)
# cpu_count = 16
batch_size = int(hap_count/cpu_count)

p = multiprocessing.Pool(cpu_count)
for i in range(cpu_count):
    # i= 0
    if i == cpu_count - 1:
        hap_ids_batch = hap_ids[ batch_size * i: ]
    else:
        hap_ids_batch = hap_ids[batch_size * i : batch_size*(i+1) ]
    batch_name = str(i)
    # print(batch_name, hap_ids_batch[:9])
    p.apply_async(merge_haps, args=(hap_ids_batch, merged_df4, haplotype_df_month1,haplotype_df_month2,haplotype_df_month3,haplotype_df_month4, coo_df1, coo_df2, coo_df3,coo_df4, batch_name))

p.close()
p.join()

print('the tasks in different cpus are over!')
haplotype_df_previous = pd.DataFrame()
coo_df_previous = pd.DataFrame()
for i in range(cpu_count):
    batch_name = str(i)
    if not os.path.exists(output_path + 'haplotype_df'+batch_name+'.txt'):
        tmpdf1 = pd.DataFrame()
    else:
        tmpdf1 = pd.read_csv(output_path+'haplotype_df'+batch_name+'.txt',sep='\t', header=0)
    if not os.path.exists(output_path + 'coo_df'+batch_name+'.txt'):
        tmpdf2 = pd.DataFrame()
    else:
        tmpdf2 = pd.read_csv(output_path+'coo_df'+batch_name+'.txt',sep='\t', header=0)
    haplotype_df_previous = pd.concat([haplotype_df_previous,  tmpdf1], axis=0)
    coo_df_previous = pd.concat([coo_df_previous,  tmpdf2], axis=0)

haplotype_df_previous.index = range(haplotype_df_previous.shape[0])
haplotype_df_previous.columns =['haplotype', 'muts_count', 'sample_names', 'sample_count','sample_locations', 'sample_time_min', 'sample_time_max', 'hap_pos', 'hap_value']
coo_df_previous.columns = ['row','col','value']
rows_previous = np.array(coo_df_previous['row'])
cols_previous = np.array(coo_df_previous['col'])
values_previous = np.array(coo_df_previous['value'])
coo_fn = coo_matrix((values_previous, (rows_previous, cols_previous)), shape=(haplotype_df_previous.shape[0], len(hap_pos_previous)), dtype=np.int8)
hap_array_fn = coo_fn.toarray()
hap_array_df_previous = pd.DataFrame(hap_array_fn)
hap_array_df_previous.columns = hap_pos_previous

#save and read in files
np.save( output_path+'rows_previous.npy', rows_previous)
np.save( output_path+'cols_previous.npy', cols_previous)
np.save( output_path+'values_previous.npy', values_previous)
haplotype_df_previous.to_csv(output_path+'haplotype_df_previous.txt',sep='\t', header=True, index=None)
pd.DataFrame(hap_pos_previous).to_csv(output_path+'hap_pos_previous.txt',sep='\t', header=False, index=None)
# data_path_coo = 'results/venas2v1_output/venas2_output_month19/before2020-09/'
# data_path_coo = output_path
# rows1 = np.load( data_path_coo+'rows_previous.npy')
# cols1 = np.load( data_path_coo+'cols_previous.npy')
# values1 = np.load( data_path_coo+'values_previous.npy')
# haplotype_df_previous = pd.read_csv(data_path_coo+'haplotype_df_previous.txt',sep='\t', header=0)
haplotype_df_previous = haplotype_df_previous.fillna('')
haplotype_df_previous['sample_time_min'] = haplotype_df_previous['sample_time_min'].apply(lambda x: datetime.strptime(x, '%Y-%m-%d'))
haplotype_df_previous['sample_time_max'] = haplotype_df_previous['sample_time_max'].apply(lambda x: datetime.strptime(x, '%Y-%m-%d'))
# haplotype_df_previous['hap_pos'] = haplotype_df_previous['hap_pos'].apply(lambda row: [ int(x) for x in row[1:-1].split(',')])
# hap_pos_previous = pd.read_csv(data_path_coo +  'hap_pos_previous.txt', sep='\t', header=None)
# hap_pos_previous = hap_pos_previous[0].tolist()
# coo1 = coo_matrix((values1, (rows1, cols1)), shape=(haplotype_df_previous.shape[0], len(hap_pos_previous)), dtype=np.int8)
# hap_array_previous = coo1.toarray()
# hap_array_df_previous = pd.DataFrame(hap_array_previous)
# hap_array_df_previous.columns = hap_pos_previous
#make sure the first is ref
if '' in haplotype_df_previous['haplotype'].tolist():
    idx = haplotype_df_previous[haplotype_df_previous['haplotype']==''].index[0]
    # haplotype_df_previous.loc[idx, 'sample_time_min'] = min(haplotype_df_previous.loc[idx]['sample_time_min'] ,datetime.strptime('2019-12-31', '%Y-%m-%d') )
    haplotype_df_previous.loc[idx, 'sample_time_min'] = '2019-12-31'
    if idx == 0:
        pass
    else:
        haps_list = haplotype_df_previous.index.tolist()
        haps_list.remove(idx)
        haps_list = [idx] + haps_list
        haplotype_df_previous = haplotype_df_previous.loc[haps_list]
        haplotype_df_previous.index = range(haplotype_df_previous.shape[0])
        # hap_array_previous = hap_array_previous[haps_list]
        hap_array_df_previous = hap_array_df_previous.loc[haps_list]
        hap_array_df_previous.index= range(hap_array_df_previous.shape[0])
else:
    ref_df = pd.DataFrame(['', 0, 'ref', 1, 'China', datetime.strptime('2019-12-31', '%Y-%m-%d'), datetime.strptime('2019-12-31', '%Y-%m-%d'), [], [] ]).T
    # ref_df = pd.DataFrame(['', 0, 'ref', 1, 'China', '2019-12-31', '2019-12-31', [], [] ]).T
    ref_df.columns = ['haplotype', 'muts_count', 'sample_names', 'sample_count',  'sample_locations', 'sample_time_min', 'sample_time_max', 'hap_pos', 'hap_value']
    haplotype_df_previous = pd.concat([ref_df, haplotype_df_previous ],axis=0)
    haplotype_df_previous.index = range(haplotype_df_previous.shape[0])
    tmpdf = pd.DataFrame(np.zeros((1, hap_array_df_previous.shape[1])))
    tmpdf.columns = hap_array_df_previous.columns 
    hap_array_df_previous = pd.concat([tmpdf, hap_array_df_previous], axis=0)
    hap_array_df_previous.index= range(hap_array_df_previous.shape[0])

if '' in haplotype_df_previous['haplotype'].tolist():
    # idx = haplotype_df_previous[haplotype_df_previous['haplotype']==''].index[0]
    haplotype_df_previous1 = haplotype_df_previous[haplotype_df_previous['haplotype']=='']
    haplotype_df_previous2 = haplotype_df_previous[haplotype_df_previous['haplotype']!='']
    haplotype_df_previous2['hap_pos'] = haplotype_df_previous2['hap_pos'].apply(lambda x: [int(t) for t in x[1:-1].split(',')])
    haplotype_df_previous2['hap_value'] = haplotype_df_previous2['hap_value'].apply(lambda x: [int(t) for t in x[1:-1].split(',')])
    haplotype_df_previous = pd.concat([haplotype_df_previous1, haplotype_df_previous2], axis=0)
    haps_list =haplotype_df_previous.index.tolist()
    haplotype_df_previous.index = range(haplotype_df_previous.shape[0])
    hap_array_df_previous = hap_array_df_previous.loc[haps_list]
    hap_array_df_previous.index = range(hap_array_df_previous.shape[0])
else:
    haplotype_df_previous['hap_pos'] = haplotype_df_previous['hap_pos'].apply(lambda x: [int(t) for t in x[1:-1].split(',')])
    haplotype_df_previous['hap_value'] = haplotype_df_previous['hap_value'].apply(lambda x: [int(t) for t in x[1:-1].split(',')])

#update the datas and convert to array
# hap_pos_updated = sorted(set(hap_pos_previous + hap_pos_month))
# hap_pos_exist = list(set(hap_pos_previous).difference(set(hap_pos_month)))
# hap_pos_overlap = list(set(hap_pos_previous).intersection(set(hap_pos_month)))
# hap_pos_new = list(set(hap_pos_month).difference(set(hap_pos_previous)))
# # sample_array = hap_array_df_previous[hap_pos_overlap]
# # sample_array_sum = np.count_nonzero(sample_array, keepdims=False, axis=1)
# # remove_list = list(np.where(sample_array_sum == 0)[0] )
# hap_pos_selected = list(sorted(set(hap_pos_overlap+hap_pos_new)))
# hap_array_df_previous = pd.concat([hap_array_df_previous, pd.DataFrame(0, index=hap_array_df_previous.index.tolist(), columns=hap_pos_new) ], axis=1)
# hap_array_df_batch = pd.concat([ pd.DataFrame(0, index=hap_array_df_month.index.tolist(), columns=hap_pos_exist), hap_array_df_month ], axis=1)
# hap_array_df_previous = hap_array_df_previous[hap_pos_selected]
# hap_array_df_batch = hap_array_df_month[hap_pos_selected]
# haplotype_df_batch = haplotype_df_month.copy()

hap_pos_updated = sorted(set(hap_pos_previous + hap_pos_month))
hap_pos_diff1 = set(hap_pos_updated).difference(set(hap_pos_previous))
hap_pos_diff2 = set(hap_pos_updated).difference(set(hap_pos_month))
hap_array_df_previous = pd.concat([hap_array_df_previous, pd.DataFrame(0, index=hap_array_df_previous.index.tolist(), columns=hap_pos_diff1) ], axis=1)
hap_array_df_previous = hap_array_df_previous[hap_pos_updated]
hap_array_df_batch = pd.concat([hap_array_df_month, pd.DataFrame(0, index=hap_array_df_month.index.tolist(), columns=hap_pos_diff2) ], axis=1)
hap_array_df_batch = hap_array_df_batch[hap_pos_updated]
haplotype_df_batch = haplotype_df_month.copy()
pd.DataFrame(hap_pos_updated).to_csv(output_path+'hap_pos_updated.txt',sep='\t', header=False, index=None)


##merge the haps into previous
haplotype_df_updated = haplotype_df_previous.copy()
hap_array_df_updated = hap_array_df_previous.copy()
# del haplotype_df_previous
# del hap_array_df_previous
# haplotype_df_updated = pd.concat([ haplotype_df_updated, haplotype_df_batch ], axis=0)
# hap_array_df_updated = pd.concat([ hap_array_df_updated, hap_array_df_batch ], axis=0)
# hap_array_updated = np.array(hap_array_df_updated, dtype='int8')
# hap_array_batch = np.array(hap_array_df_batch, dtype='int8')


df1 = pd.DataFrame()
df1['hapid'] = haplotype_df_updated.index
df1['haplotype'] = haplotype_df_updated['haplotype']
df2 = pd.DataFrame()
df2['hapid'] = haplotype_df_batch.index
df2['haplotype'] = haplotype_df_batch['haplotype']
haplotype_df_merge = pd.merge(df1, df2, how='inner', on='haplotype')
if haplotype_df_merge.shape[0] == 0: #there is no samples in existing haps
    pass
else:
    # tmpdf = haplotype_df_batch.loc[ haplotype_df_merge['hapid_y'].tolist()]
    # tmpdf['hapid'] = haplotype_df_merge['hapid_y'].tolist()
    # tmpdf.to_csv(output_path + 'haplotype_df_exists.txt',sep='\t', index=False, header=True)
    for each_hap_info in haplotype_df_merge.values:
        hapid_x, haplotype, hapid_y = each_hap_info
        haplotype_df_updated.loc[hapid_x, 'sample_names'] += ';' + haplotype_df_batch.loc[hapid_y, 'sample_names']
        haplotype_df_updated.loc[hapid_x, 'sample_count'] += haplotype_df_batch.loc[hapid_y, 'sample_count']
        haplotype_df_updated.loc[hapid_x, 'sample_locations'] = ';'.join(set(haplotype_df_updated.loc[hapid_x, 'sample_locations'].split(';')+ haplotype_df_batch.loc[hapid_y, 'sample_locations'].split(';')))
        haplotype_df_updated.loc[hapid_x, 'sample_time_min'] = min( haplotype_df_updated.loc[hapid_x, 'sample_time_min'], haplotype_df_batch.loc[hapid_y, 'sample_time_min'])
        haplotype_df_updated.loc[hapid_x, 'sample_time_max'] = max( haplotype_df_updated.loc[hapid_x, 'sample_time_max'], haplotype_df_batch.loc[hapid_y, 'sample_time_max'])
    haplotype_removed_indexs = [x for x in haplotype_df_batch.index.tolist() if x not in haplotype_df_merge['hapid_y'].tolist()]
    haplotype_df_batch = haplotype_df_batch.loc[haplotype_removed_indexs]
    hap_array_df_batch = hap_array_df_batch.loc[haplotype_removed_indexs]

hap_count_updated = haplotype_df_updated.shape[0]
hap_count_batch = haplotype_df_batch.shape[0]
haplotype_df_batch.index = range(hap_count_updated, hap_count_updated + hap_count_batch )
hap_array_df_batch.index = range(hap_count_updated, hap_count_updated + hap_count_batch )

# hap_index_df = pd.DataFrame()
# hap_index_df['origidx'] = haplotype_removed_indexs
# hap_index_df['newidx'] = haplotype_df_batch.index
# hap_index_df.to_csv(output_path + 'hap_index_df.txt',sep='\t', index=False, header=True)

haplotype_df_updated = pd.concat([ haplotype_df_updated, haplotype_df_batch ], axis=0)
hap_array_df_updated = pd.concat([ hap_array_df_updated, hap_array_df_batch ], axis=0)
hap_array_updated = np.array(hap_array_df_updated, dtype='int8')
hap_array_batch = np.array(hap_array_df_batch, dtype='int8')

hap_ids = hap_array_df_batch.index.tolist()
hap_count = len(hap_ids)
# del hap_array_df_updated
# del hap_array_df_batch


"""step2 :construct the faiss index and query to search"""
print('constructing the faiss index and query...')
xb = hap_array_updated.astype(np.int8)
xq = hap_array_batch.astype(np.int8)
measure = faiss.METRIC_L2
param = 'Flat'
# param = 'HNSW64'
# dim = len(hap_pos_batch)
dim= hap_array_batch.shape[1]
index = faiss.index_factory(dim, param, measure)
# print(xb, xq, xb.shape, xq.shape, dim,index)
# print(index.is_trained)
index.add(xb)
can_dis, can_ids = index.search(xq, can_lens)
np.save( output_path+'can_dis.npy', can_dis)
np.save( output_path+'can_ids.npy', can_ids)
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
        # results_logs.append(['processing', cur_hapidx, ';'.join([str(x) for x in cur_hap_pos]),';'.join([str(x) for x in cur_hap_value])] )
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
            if common_flag == False: #stage 2 continue check top 10 - 100 nearest haps 
                flag_count = 0
                parent_can_list2 = []
                # a = hap_array_df_updated.loc[cur_hapidx]
                for each_hap_index2 in range(each_hap_index, can_lens):
                    # each_hap_index2 = 1
                    each_hapidx = cur_can_ids[each_hap_index2]
                    if each_hapidx == -1:
                        break
                    if flag_count >= 50:
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
                edge_muts = ';'.join([ 'bm:'+x for x in t2 ])
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

# hap_ids = hap_array_df_batch.index.tolist()
# hap_count = len(hap_ids)
# cpu_count = multiprocessing.cpu_count()
cpu_count = int(args.cpu_count)
# cpu_count = 16
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
        # results_logs.append(['processing', cur_hapidx, ';'.join([str(x) for x in cur_hap_pos]),';'.join([str(x) for x in cur_hap_value])] )
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
                edge_muts = ';'.join([ x for x in t1 ])
                t1_count = len(t1)
                t2_count = 0
            elif t1!=set() and t2!=set():
                edge_muts = ';'.join([ x for x in t1 ]) + ';' + ';'.join([ 'bm:'+x for x in t2 ])
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
# edges_df['source'] = edges_df['source'].apply(lambda x: int(x))
# edges_df['target'] = edges_df['target'].apply(lambda x: int(x))

# results_voi = pd.DataFrame(results_voi)
results_voi.to_csv(output_path+'results_voi.txt',sep='\t', index=None , header=None)


end_time = time.time()
diff = end_time - beg_time
print(time.ctime())
print(diff, ' s')
print(diff/60, ' min')
print(diff/60/60, ' h')


            