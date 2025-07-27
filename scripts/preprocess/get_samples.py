#!/usr/bin/env python
# coding=utf-8
# __author__ = 'zp'

import pandas as pd
import numpy as np
import os
import re
from Bio import SeqIO
import time
from datetime import datetime
import multiprocessing
from Bio import AlignIO, SeqIO

beg_time = time.time()


meta_df2 =pd.read_csv('rawdata/genome_metadata_2024-06-12_all.txt', sep='\t', header=0)
meta_df2['time1'] = meta_df2['collection_date'].apply(lambda x: datetime.strptime(x, '%Y-%m-%d'))
meta_df2 = meta_df2.sort_values('time1',ascending=True)

meta_df2 = meta_df2[(meta_df2['host']=='Homo sapiens') | (meta_df2['host']=='Human' ) | (meta_df2['host'] == 'Humano') | (meta_df2['host'] == 'Human male') | (meta_df2['host'] == 'Human female')]
merged_df = meta_df2
merged_df2 = merged_df[merged_df['Accession ID'].notna()]
# merged_df2['Host'].drop_duplicates()
merged_df2 = merged_df2[merged_df2['Host']=='Human']
merged_df2[merged_df2['Location'].isna()]
merged_df2 = merged_df2[merged_df2['Pango lineage'].notna()]
# remember merged_df7[merged_df7['Pango lineage']!='Unassigned']
merged_df2[merged_df2['Collection date'].isna()]
merged_df2.index = range(merged_df2.shape[0])

#remove the wrong time
def identify_wrongtime(hap_ids_fn, batch_name):
    # hap_ids_fn = hap_ids[:2000]
    print('the task in ', batch_name)
    wrong_time_list = []
    # for eachidx in range(meta_df.shape[0]):
    for eachidx in hap_ids_fn:
        # eachidx = 0
        curid = meta_df.loc[eachidx]['Accession ID']
        curdate = meta_df.loc[eachidx]['Collection date']
        dcount = len(re.findall('-',curdate))
        if dcount!=2:
            wrong_time_list.append(curid)
    pd.DataFrame(wrong_time_list).to_csv('results/tmp/wrong_'+batch_name+'.txt', sep='\t', index=False, header=False)

hap_ids = meta_df.index.tolist()
hap_count = len(hap_ids)
# cpu_count = int(multiprocessing.cpu_count())
# cpu_count = int(args.cpu_count)
cpu_count = 64
batch_size = int(hap_count/cpu_count)

p = multiprocessing.Pool(cpu_count)
for i in range(cpu_count):
    # i= 0
    if i == cpu_count - 1:
        hap_ids_batch = hap_ids[ batch_size * i: ]
    else:
        hap_ids_batch = hap_ids[batch_size * i : batch_size*(i+1) ]
    batch_name = str(i)
    # print(batch_name, hap_ids_batch[:9], can_ids_batch[:9])
    p.apply_async(identify_wrongtime, args=( hap_ids_batch, batch_name ))

p.close()
p.join()


wrong_time_all = pd.read_csv('rawdata/metadata_2024_12_wrongtimelist.txt', sep='\t', header=None)
tmpdf1 = pd.DataFrame()
tmpdf1['accession_id'] = merged_df2['accession_id']
tmpdf1['idx'] = tmpdf1.index.tolist()
tmpdf2 = pd.DataFrame()
tmpdf2['accession_id'] = wrong_time_all[0]
tmpdf2['idx'] = tmpdf2.index.tolist()

merged_df3 = pd.merge(tmpdf1, tmpdf2, how='left', on='accession_id')
nowronge_list = merged_df3[merged_df3['idx_y'].isna()]['idx_x']
merged_df4 = merged_df2.loc[nowronge_list]
merged_df4.index = range(merged_df4.shape[0])


data= pd.read_csv('rawdata/vigtk_genome_qc/sample_qc/sample_qc_all.txt',sep='\t',header=0)
# tmpdf =data[['deletion_count', 'unknown_count', 'ambiguous_count']]
data['amb_sum'] = data['deletion_count']+ data['unknown_count'] + data['ambiguous_count']
# data = data[(data['seq_len']>28000) & (data['amb_sum']<1000)]
data = pd.merge( data, meta_df2[['virus_id', 'accession_id']], how='left', left_on='samplename', right_on='virus_id')
# data[data['accession_id'].isna()]

merged_df5 = pd.merge(merged_df4, data, how='left', left_on='accession_id', right_on='accession_id')
extradf = merged_df5[merged_df5['samplename'].isna()]
merged_df6 = merged_df5[merged_df5['samplename'].notna()]
merged_df6 = merged_df6[(merged_df6['seq_len']>28000) & (merged_df6['amb_sum']<1000)]

extradf2=  extradf[extradf['Sequence length']>29000]
merged_df7 = pd.concat([extradf2, merged_df6], axis=0)
merged_df7 = merged_df7[merged_df7['Pango lineage']!='Unassigned']
merged_df7.index = range(merged_df7.shape[0])
merged_df7.to_csv('results/venas2v1_output/merged_df7.txt', sep='\t', index=False, header=True)

merged_df7 = pd.read_csv('results/venas2v1_output/merged_df7.txt', sep='\t', header=0)
merged_df7 = merged_df7.fillna('')

# merged_df7[merged_df7['amb_sum']==0]
# merged_df8 = pd.merge(merged_df7, data_all, how='left', left_on='accession_id', right_on='samplename')
merged_df8 = merged_df7
sampledf1 = merged_df8[merged_df8['deletion_count']=='']
sampledf2 = merged_df8[merged_df8['deletion_count']!='']

sampledf5 = pd.merge(sampledf1, sample_info_df3, how='left', left_on=['virus_id_x', 'accession_id'], right_on=['samplename','accession_id'] )
# sampledf5[sampledf5['deletion_count_y']<50]
sampledf1[['seq_len','deletion_count','unknown_count']] = sampledf5[['seq_len_y','deletion_count_y','unknown_count_y']]
sampledf1[sampledf1['deletion_count']<50]
merged_df9 = pd.concat([sampledf1, sampledf2],axis=0)
merged_df9.index= range(merged_df9.shape[0])
merged_df9.to_csv('results/venas2v1_output/merged_df9.txt', sep='\t', index=False, header=True)

merged_df9['time1'] = merged_df9['Collection date'].apply(lambda x: datetime.strptime(x, '%Y-%m-%d'))
merged_df9['time2'] = merged_df9['Submission date'].apply(lambda x: datetime.strptime(x, '%Y-%m-%d'))

def identify_deltatime(hap_ids_fn, batch_name):
    # hap_ids_fn = hap_ids
    # tmpdf_fn = merged_df9
    print('the task in ', batch_name)
    deltas_fn = []
    for eachidx in hap_ids_fn:
    # for eachidx in range(tmpdf_fn.shape[0]):
        # eachidx = 0
        # print(eachidx)
        accid = tmpdf_fn.loc[eachidx]['accession_id']
        time1=  tmpdf_fn.loc[eachidx]['time1']
        time2=  tmpdf_fn.loc[eachidx]['time2']
        delta = time2 - time1
        print(eachidx,  delta)
        deltas_fn.append([accid, delta.days])
    pd.DataFrame(deltas_fn).to_csv('results/tmp/deltas_'+str(batch_name)+'.txt', sep='\t', index=False, header=False)

hap_ids = merged_df9.index.tolist()
hap_count = len(hap_ids)
# cpu_count = int(multiprocessing.cpu_count())
# cpu_count = int(args.cpu_count)
cpu_count = 64
batch_size = int(hap_count/cpu_count)

p = multiprocessing.Pool(cpu_count)
for i in range(cpu_count):
    # i= 0
    if i == cpu_count - 1:
        hap_ids_batch = hap_ids[ batch_size * i: ]
    else:
        hap_ids_batch = hap_ids[batch_size * i : batch_size*(i+1) ]
    batch_name = str(i)
    # print(batch_name, hap_ids_batch[:9], can_ids_batch[:9])
    p.apply_async(identify_deltatime, args=( hap_ids_batch, batch_name ))

p.close()
p.join()

deltas_all = pd.DataFrame()
for i in range(cpu_count):
    # i = 0
    print('processing ',str(i))
    batch_name = str(i)
    data_tmp = pd.read_csv('results/tmp/deltas_'+batch_name+'.txt', sep='\t', header=None)
    deltas_all = pd.concat([deltas_all, data_tmp], axis=0)

deltas_all.columns = ['accession_id', 'deltas']

merged_df10 = pd.merge(merged_df9, deltas_all, how='left', on='accession_id')
merged_df10.to_csv('results/venas2v1_output/merged_df10.txt', sep='\t', index=False, header=True)

merged_df10 = pd.read_csv('results/venas2v1_output/merged_df10.txt', sep='\t', header=0)
merged_df10 = merged_df10.fillna('')
merged_df10[(merged_df10['unknown_count']==0)]
merged_df10[(merged_df10['unknown_count']<10) & (merged_df10['deletion_count'] < 90)]
merged_df10[(merged_df10['unknown_count']<10)]

#divide innto mainn annd minonr
merged_df10['time1'] = merged_df10['Collection date'].apply(lambda x: datetime.strptime(x, '%Y-%m-%d'))
merged_df10 = merged_df10.sort_values('time1',ascending=True)
# merged_df10['amb_sum'] = merged_df10['deletion_count']+ merged_df10['unknown_count']
can_df_orig = merged_df10[merged_df10['deltas']>5* 30]
merged_df11 = merged_df10[merged_df10['deltas']<5* 30]
maindf = merged_df11[(merged_df11['unknown_count']<10) & (merged_df11['deletion_count'] < 90)]
# maindf = merged_df11[(merged_df11['unknown_count']==0)]
maindf['haplotype'].drop_duplicates()
maindf.index = range(maindf.shape[0])

data_path = 'results/venas2v10_output/'
samples_meta = maindf[['Virus name', 'Accession ID', 'virus_id_x', 'Collection date','Submission date',  'Location', 'Sequence length', 'Pango lineage', 'AA Substitutions', 'location_y']]
samples_meta.to_csv(data_path + 'samples_metadata4.txt',sep='\t', index=False, header=True)

samples_muts = maindf[['accession_id', 'haplotype', 'Collection date', 'location_y']]
samples_muts.columns = ['samplename', 'haplotype', 'collection_date', 'location']
samples_muts.to_csv(data_path + 'samples_mutations4.txt',sep='\t', index=False, header=True)


