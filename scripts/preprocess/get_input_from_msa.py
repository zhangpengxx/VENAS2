#!/usr/bin/env python
# coding=utf-8
# __author__ = 'zp'
#to prepare the samples mutation files (input) from sequence alignment 

import pandas as pd
import numpy as np
from datetime import datetime
from datetime import timedelta
import time
import argparse
import os
import re
import random
from Bio import AlignIO, SeqIO
from collections import Counter
import faiss
import scipy.sparse
from scipy.sparse import coo_matrix

beg_time = time.time()

# parser = argparse.ArgumentParser()
# parser.add_argument('-inputname', '--inputname', dest='inputname', help="the name of input samples", required=True)
# args = parser.parse_args()

##parameters
default_nc = {'a', 't', 'c', 'g'}
dict_nc2num = {'a':1, 't':2, 'c':3, 'g':4, '-':5, 'o':6} #o means others
dict_num2nc = {1:'a', 2:'t', 3:'c', 4:'g', 5:'-', 6:'o'}
selected_pos = [ x for x in range(100, 29803)]
selected_pos_count = len(selected_pos)

"""references informations"""
ref_fafile = 'example/OEAV139851.fasta'
f = open(ref_fafile, 'r')
refseq = []
for line in f.readlines():
    line = line.strip()
    if line == '>OEAV139851':
        pass
    else:
        refseq.append(line)

refseq = ''.join(refseq).lower()
refdf = pd.DataFrame(list(refseq))
refseq_len = len(refseq)
refdf = refdf.iloc[selected_pos]

# threshold_sample_count = 10
#you should have metadata
meta_df = pd.read_csv('example/samples_metadata.txt', sep='\t', header=0)
# meta_df['country'] = meta_df['location'].apply(lambda x:x.split('/')[1].strip())

##load in the data
# batch = args.batch
output_path =  'example/'
if os.path.exists(output_path):
    print('the folder exists! please check the output dir')
else:
    os.makedirs(output_path)

### build the array
ma_file =  'example/samples_1000.ma'
alignments = AlignIO.read(ma_file, "fasta")
samples_count = len(alignments)
batch_size = 10000
# batch_remainder = samples_count % batch_size
batch_times = samples_count // batch_size

dict_muta_info = {}
dict_batch_array = {}
samplename_list = []
sample_array_df_all = pd.DataFrame()
# nonatcg_samplename_list = []
for batch_idx in range(0, batch_times+1):
    # batch_idx = 0
    print('processing the array batch ', str(batch_idx),'/' ,str(batch_times))
    data_batch = alignments[batch_idx*batch_size : (batch_idx+1)*batch_size]
    cur_batch_count = len(data_batch)
    if cur_batch_count == 0:
        continue
    # sample_array_df = pd.DataFrame(0, index = range(cur_batch_count), columns=range( 29903) )
    rows = []
    cols = []
    values = []
    samplename_list_batch = []
    for sample_idx in range(cur_batch_count):
    # for sample_idx in range(10):
        # sample_idx = 0
        print('processing ',sample_idx)
        samplename = data_batch[sample_idx].id
        samplename_list_batch.append(samplename)
        samplename_list.append(samplename)
        sampleseq = data_batch[sample_idx].seq
        # print(list(sampleseq)[:50], list(sampleseq)[-50:])
        sampleseqdf = pd.DataFrame(sampleseq)
        sampleseqdf = sampleseqdf.iloc[selected_pos]
        # tmpdf = meta_df[meta_df['virus_id']==samplename]
        # print(sampleseqdf[0].tolist()[:10], sampleseqdf[0].tolist()[-10:])
        # beg_match = re.findall(r'^((\-)\2{2,})', str(sampleseq) )
        # end_match = re.findall(r'((\-)\2{2,}$)', str(sampleseq) )
        # re.match(r'^((\-)\2{2,})', str(sampleseq) )
        # re.search(r'((\-)\2{2,}$)', str(sampleseq) )
        # print(sample_idx,sampleseqdf)
        if '-' in sampleseqdf[0].tolist():  #process the seqs which have too many - in the end 
            # print(sample_idx)
            beg_match = re.findall(r'^((\-)\2{2,})', ''.join(sampleseqdf[0].tolist()) )
            end_match = re.findall(r'((\-)\2{2,}$)', ''.join(sampleseqdf[0].tolist()) )
            # middle_match = re.findall(r'((\-)\2{2,})', ''.join(sampleseqdf[0].tolist()) )
            # print(sample_idx, middle_match, [ len(x) for x in middle_match[0]])
            if len(beg_match) == 0:
                beg_count = 0
            else:
                beg_count = len(beg_match[0][0])
            if len(end_match) == 0:
                end_count = 0
            else:
                end_count = len(end_match[0][0])
        else:
            beg_count = 0
            end_count = 0
        # print(re.findall(r'^((\-)\2{2,})', str(sampleseq) ), re.findall(r'((\-)\2{2,}$)', str(sampleseq) ))
        # print(sample_idx, beg_count, end_count)
        diff = refdf.iloc[beg_count:(selected_pos_count-end_count)].compare(sampleseqdf.iloc[beg_count:(selected_pos_count-end_count)])
        # print(diff)
        if diff.shape[0] == 0:
            continue
        else:
            # if '-' in diff[0]['other'].tolist():
            #     print(sample_idx)
            # non_atcg_df = diff[0].query("other!='a' and other!='t' and other!='c' and other!='g' and other !='-' ")
            tmpdf = diff[0]['other'].replace('[^actg-]', 'o', regex=True)
            sample_vector = tmpdf.apply(lambda x:dict_nc2num[x])
            # sample_array_df.iloc[sample_idx][tmpdf.index.tolist()] = sample_vector
            # tmpdf = pd.DataFrame(cur_hap_pos)
            # tmpdf['posidx'] = tmpdf.apply(lambda x:dict_pos_index[x[0]], axis=1)
            # cur_hap_value = haplotype_df_fn.loc[each_hapidx]['hap_value']
            rows = rows + [sample_idx]*tmpdf.shape[0]
            cols = cols + sample_vector.index.tolist()
            values = values + sample_vector.values.tolist()
    values = np.array(values)
    rows = np.array(rows)
    cols = np.array(cols)
    # coo = coo_matrix((values, (rows, cols)), shape=(100, 29903), dtype=np.int8)
    coo = coo_matrix((values, (rows, cols)), shape=(cur_batch_count, 29903), dtype=np.int8)
    sample_array = coo.toarray()
    # np.save( output_path+'rows'+str(batch_idx)+'.npy', rows)
    # np.save( output_path+'cols'+str(batch_idx)+'.npy', cols)
    # np.save( output_path+'values'+str(batch_idx)+'.npy', values)
    # rows1 = np.load( output_path+'rows.npy')
    # cols1 = np.load( output_path+'cols.npy')
    # values1 = np.load( output_path+'values.npy')
    # coo1 = coo_matrix((values1, (rows1, cols1)), shape=(cur_batch_count, 29903), dtype=np.int8)
    # pd.DataFrame(samplename_list_batch).to_csv(output_path + 'samplename_list'+str(batch_idx)+'.txt', sep='\t', index=False, header=True)
    # pd.DataFrame(samplename_list_batch).to_csv(output_path2 + 'samplename_list'+str(batch_idx)+'.txt', sep='\t', index=False, header=True)
    # sample_array = np.array(sample_array_df)
    sample_array_df = pd.DataFrame(sample_array)
    dict_batch_array[batch_idx] = sample_array_df
    sample_array_df_all = pd.concat([sample_array_df_all, sample_array_df],axis=0)
    # sample_array_df.to_csv(output_path + 'samples_array_df_'+str(batch_idx)+'.txt', sep='\t', index=False, header=True)
    sample_array_sum = np.count_nonzero(sample_array, keepdims=False, axis=0)
    mut_pos_list = list(np.where(sample_array_sum > 0)[0] )
    # mut_pos_list = [ x for x in mut_pos_list if x+1 not in problem_pos]
    for each_pos in mut_pos_list:
        # each_pos = mut_pos_list[0]
        # ref_letter = refdf.loc[each_pos][0]
        # sample_array_df[each_pos].drop_duplicates()
        # sample_array_df[each_pos].value_counts()
        uniq_nc = Counter(sample_array_df[each_pos])
        # print(uniq_nc)
        for k, v in uniq_nc.items():
            if k != 0:
                # alt = dict_num2nc[k]
                alt = k
                each_ac = v
                # muta_info_df.append([each_pos, ref_letter, alt, v ])
                key = str(each_pos) + '-' +str(alt) 
                if key in dict_muta_info.keys():
                    dict_muta_info[key] = dict_muta_info[key] + v
                else:
                    dict_muta_info[key] = v

# pd.DataFrame(samplename_list).to_csv(output_path + 'samplename_list.txt', sep='\t', index=False, header=True)
# pd.DataFrame(samplename_list).to_csv(output_path2 + 'samplename_list.txt', sep='\t', index=False, header=True)
muta_info_df = []
for k, v in dict_muta_info.items():
    pos, alt = k.split('-')
    pos = int(pos)
    ref_letter = refdf.loc[pos][0]
    alt = dict_num2nc[int(alt)]
    muta_info_df.append([pos, ref_letter, alt, v])

muta_info_df = pd.DataFrame(muta_info_df)
muta_info_df.columns = ['pos', 'ref', 'alt', 'ac']
muta_info_df = muta_info_df.sort_values('pos', ascending=True )
# muta_info_df.to_csv(output_path + 'samples_muta_info_orig.txt', sep='\t', index=False, header=True)

# muta_info_df[muta_info_df['ac']>7]
# tmpdf = pd.read_csv('results/venas2v1_output/samples_haps_10v4/samples_muta_info_all.txt', sep='\t', header=0)
# tmpdf = tmpdf[tmpdf['alt']!='o']
# tmpdf = tmpdf[tmpdf['ac']>=7]
# mut_pos_list = list(sorted(tmpdf['pos'].drop_duplicates().tolist()))

# dict_pos = {}
# for idx in range(tmpdf.shape[0]):
#     # idx=  0
#     pos,ref, alt, ac= tmpdf.iloc[idx]
#     alt_num = dict_nc2num[alt]
#     key = str(pos) +'_' + str(alt_num)
#     # key = str(pos) +'-' +str(ref)+'-' + str(alt)
#     dict_pos[key] = ac

##imputation the ambiguous nucleotide
threshold_sample_count = 1
atcg_df = muta_info_df[muta_info_df['alt'] !='o']
atcg_df = atcg_df[atcg_df['ac']>= threshold_sample_count]
# atcg_df.to_csv(output_path + 'samples_muta_info_orig_selected.txt', sep='\t', index=False, header=True)
mut_pos_df = atcg_df['pos'].drop_duplicates()
mut_pos_df = pd.DataFrame(mut_pos_df)
mut_pos_list = mut_pos_df['pos'].tolist()

dict_pos = {}
for idx in range(atcg_df.shape[0]):
    # idx=  0
    pos,ref, alt, ac= atcg_df.iloc[idx]
    alt_num = dict_nc2num[alt]
    key = str(pos) +'_' + str(alt_num)
    # key = str(pos) +'-' +str(ref)+'-' + str(alt)
    dict_pos[key] = ac


sample_array_df_all.index= range(sample_array_df_all.shape[0])
sample_array_df_all = sample_array_df_all[mut_pos_list]

samples_mut = []
# for cur_hapidx in range(10):
for cur_hapidx in range(sample_array_df_all.shape[0]):
    # cur_hapidx  = 0
    print('processing ',str(cur_hapidx))
    # sample_idxs = haps_sampleids_df.loc[cur_hapidx][0]
    cur_mut_df = sample_array_df_all.loc[cur_hapidx]
    cur_mut_df = cur_mut_df[cur_mut_df>0]
    if cur_mut_df.shape[0] == 0:
        samples_mut.append([cur_hapidx, '', ''])
        continue
    #for the mut count > 0
    # cur_mut_pos = cur_mut_df.index.tolist()
    selected_pos_index = []
    for x in zip(cur_mut_df.index.tolist(), cur_mut_df.tolist()):
        key = str(x[0]) + '_' + str(x[1])
        if key in dict_pos.keys():
            # print(x)
            selected_pos_index.append(x[0])
    cur_mut_df = cur_mut_df.loc[selected_pos_index]
    samples_mut.append([cur_hapidx, ','.join([str(x) for x in cur_mut_df.index.tolist()]), ','.join([str(x) for x in cur_mut_df.tolist()]) ])

samples_mut = pd.DataFrame(samples_mut)
samples_mut.columns = ['idx', 'pos', 'value']

haps_df = []
# haps_sampleids_df = []
for each_hap, group in samples_mut.groupby(['pos', 'value']):
    # print(each_hap, group)
    sample_idxs = group['idx'].astype('str').tolist()
    haps_df.append( [each_hap[0], each_hap[1], ';'.join(sample_idxs) ] )

haps_df = pd.DataFrame(haps_df)
haps_df.columns = ['pos', 'value', 'sample_idxs']

##process the samples
print('processing the samples')
samples_info_df_notna = []
for cur_hapidx in range(haps_df.shape[0]):
    # cur_hapidx = 2
    print('processing ',cur_hapidx)
    hap_pos, hap_value, sample_idxs = haps_df.loc[cur_hapidx]
    if hap_pos == '':
        if ';' in sample_idxs:
            for sample_idx in sample_idxs.split(';'):
                samplename = samplename_list[int(sample_idx)]
                samples_info_df_notna.append([samplename, ''])
        else:
            samplename = samplename_list[int(sample_idxs)]
            samples_info_df_notna.append([samplename, ''])
        continue
    cur_mut_pos = [ int(x) for x in hap_pos.split(',')]
    cur_mut_value = [ int(x) for x in hap_value.split(',')]
    cur_mut_df = pd.DataFrame(cur_mut_value)
    cur_mut_df.index = cur_mut_pos
    haps_new = []
    for each_pos in cur_mut_pos:
        ref_letter = refdf.loc[each_pos][0]
        alt = cur_mut_df.loc[each_pos][0]
        alt = dict_num2nc[alt]
        haps_new.append([each_pos, ref_letter, alt])
    haps_new = pd.DataFrame(haps_new)
    tmpdf1 = haps_new[haps_new[2]=='-']
    tmpdf2 = haps_new[haps_new[2]!='-']
    if tmpdf1.shape[0]==0:
        tmpdf3 = tmpdf2
    else:
        haps_new2 = []
        added_pos = []
        for eachidx in range(tmpdf1.shape[0]):
            # eachidx = 0
            pos, ref, alt = tmpdf1.iloc[eachidx]
            # if alt == '-':
            if pos not in added_pos:
                i = 0
                while True:
                    pos_new = pos + i + 1
                    if pos_new in tmpdf1[0].tolist():
                        i += 1
                        added_pos.append(pos_new)
                    else:
                        break
                haps_new2.append([pos, i, '-'])
        haps_new2 = pd.DataFrame(haps_new2)
        tmpdf3 = pd.concat([tmpdf2, haps_new2],axis=0)
        tmpdf3 = tmpdf3.sort_values(0, ascending=True)
    cur_mut = []
    for eachidx in range(tmpdf3.shape[0]):
        # eachidx = 0
        each_pos, ref, alt = tmpdf3.iloc[eachidx]
        if alt != '-':
            cur_mut.append( ref.upper() + str(each_pos+1)+alt.upper())
        else:
            cur_mut.append( 'del'+str(each_pos+1) + '_' +str(each_pos+1+ref) )
    cur_mut = ';'.join(cur_mut)
    if ';' in sample_idxs:
        for sample_idx in sample_idxs.split(';'):
            samplename = samplename_list[int(sample_idx)]
            samples_info_df_notna.append([samplename, cur_mut])
    else:
        samplename = samplename_list[int(sample_idxs)]
        samples_info_df_notna.append([samplename, cur_mut])

samples_info_df_notna = pd.DataFrame(samples_info_df_notna)
samples_info_df_notna.columns = ['samplename', 'haplotype']

#get the meta info
samples_result = pd.merge(samples_info_df_notna, meta_df, how='left', left_on='samplename', right_on='virus_id')
samples_result = samples_result[['accession_id','haplotype','collection_date','country']]
samples_result.columns =['samplename','haplotype','collection_date','location']
samples_result.to_csv(output_path + 'samples_mutations.txt',sep='\t', index=False, header=True)

end_time = time.time()
diff = end_time - beg_time
print(time.ctime())
print(diff, ' s')
print(diff/60, ' min')
print(diff/60/60, ' h')


