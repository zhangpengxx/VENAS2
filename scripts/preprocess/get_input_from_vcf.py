#!/usr/bin/env python
# coding=utf-8
# __author__ = 'zp'
#to prepare the samples mutation files (input) from each vcf

import pandas as pd
import numpy as np
from datetime import datetime
from datetime import timedelta
import time
# import argparse
import os
import random

beg_time = time.time()

# @profile
# def main():
##load in the data
data_path = 'example/vcffiles/'
output_path =  'venas2_output/'
if os.path.exists(output_path):
    print('the folder exists! please check the output dir')
else:
    os.makedirs(output_path)

meta_df = pd.read_csv('example/samples_metadata.txt', sep='\t', header=0)
# samples_mutations = pd.read_csv('example/samples_mutations.txt', sep='\t', header=0)
# samples_mutations = samples_mutations.fillna('')

default_nc = {'A', 'T', 'C', 'G'}
muta_info_df = pd.DataFrame(columns=['pos', 'ref', 'alt', 'ac'])
samples_list = []
# samples_info_df = []
dict_samples_info = {}
for eachidx in range(meta_df.shape[0]):
# for eachidx in range(10):
    # eachidx = 0
    print('processing the ', str(eachidx), 'th sample')
    samplename = meta_df.loc[eachidx]['accession_id']
    muta_file = data_path + samplename + '.vcf'
    # print(muta_file)
    if os.path.exists(muta_file):
        samples_list.append(samplename)
        mutadf = pd.read_csv(muta_file, sep='\t', skiprows =3, header=0)
        if mutadf.shape[0] == 0:
            # samples_info_df.append( [sampleid, '', '' ] )
            # dict_samples_info = 
            continue
        # mutation = []
        each_pos = []
        each_muta_df = []
        for each in mutadf.values:
            # CHROM, pos, ID, ref, alt, QUAL, FILTER, INFO, FORMAT, info = mutadf.loc[0]
            CHROM, pos, ID, ref, alt, QUAL, FILTER, INFO, FORMAT, info = each
            if alt in default_nc:
                # print(pos, ref, alt)
                # if len(ref) > 1:
                # else:
                each_muta_df.append([pos, ref, alt])
                each_pos.append(pos)
        each_muta_df = pd.DataFrame(each_muta_df)
        each_muta_df.columns = ['pos', 'ref', 'alt']
        dict_samples_info[samplename] = each_muta_df
        muta_info_df = pd.merge(muta_info_df, each_muta_df, how='outer', on=['pos', 'ref', 'alt'])
        if muta_info_df.shape[0] == 0:
            pass
        else:
            muta_info_df['ac'] =  muta_info_df['ac'] + 1
        muta_info_df = muta_info_df.fillna(1)
        # each_pos_tmp = ','.join([ str(x) for x in each_pos ])
        # samples_info_df.append( [sampleid,  each_pos_tmp ] )


# samples_info_df = pd.DataFrame(samples_info_df)
# samples_info_df.columns = ['samplename', 'mutations', 'pos']
# muta_info_df.to_csv(output_path + 'samples_muta_info.txt',sep='\t', index=False, header=True)
# samples_info_df.to_csv(output_path + 'samples_info.txt',sep='\t', index=False, header=True)
 

threshold_sample_count = 1
# atcg_df = muta_info_df[muta_info_df['alt'] !='o']
# atcg_df = atcg_df[atcg_df['ac']>= threshold_sample_count]
# # atcg_df.to_csv(output_path + 'samples_muta_info_orig_selected.txt', sep='\t', index=False, header=True)
# mut_pos_df = atcg_df['pos'].drop_duplicates()
# mut_pos_df = pd.DataFrame(mut_pos_df)
# mut_pos_list = mut_pos_df['pos'].tolist()
muta_info_df  = muta_info_df[muta_info_df['ac']>= threshold_sample_count]
muta_info_df = muta_info_df.sort_values('pos', ascending=True)
dict_mutations_selected = {}
for each in muta_info_df.values:
    ac, pos, ref, alt = each
    key = str(pos) + '-'+ref + '-' + alt
    dict_mutations_selected[key] = ac
 
##correct the samples
samples_info_df = []
for eachidx in range(len(samples_list)):
# for eachidx in range(10):
    # eachidx = 1
    print('processing the ', str(eachidx), 'th sample')
    samplename = samples_list[eachidx]
    if samplename in dict_samples_info.keys():
        each_muta_df = dict_samples_info[samplename]
        each_muta = []
        for each_mut in each_muta_df.values:
            pos, ref, alt = each_mut
            key = str(pos) + '-'+ref + '-' + alt
            if key in dict_mutations_selected.keys():
                if len(ref) > 1:
                    each_muta.append(  'del' + str(pos+1) +'_' + str(pos+len(ref)-1))
                else:
                    each_muta.append( ref + str(pos) + alt)
                # each_muta.append( str(pos)+'(SNP:'+ref+'->'+alt+')')
        each_muta = ';'.join(each_muta)
        samples_info_df.append([samplename,each_muta])
    else:
        samples_info_df.append([samplename, ''])

samples_info_df = pd.DataFrame(samples_info_df)
samples_info_df.columns = ['samplename', 'mutations']

samples_result = pd.merge(samples_info_df, meta_df[['accession_id','collection_date','country']], how='left', left_on='samplename', right_on='accession_id')
samples_result = samples_result[['samplename','mutations','collection_date','country']]
samples_result.to_csv(output_path + 'samples_mutations.txt',sep='\t', index=False, header=True)


end_time = time.time()
diff = end_time - beg_time
print(time.ctime())
print(diff, ' s')
print(diff/60, ' min')
print(diff/60/60, ' h')

# # if __name__ == '__main__':
# #     main()

