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

beg_time = time.time()

data_path = 'datasets/'
# data_path = 'results/venas2v1_output/samples_haps_10v4/samples_10w_142/'
data_all = pd.read_csv(data_path + 'samples_mutations.txt', sep='\t', header=0)
data_all = data_all.fillna('')
data=  data_all
data['time'] = data['collection_date'].apply(lambda x: datetime.strptime(x, '%Y-%m-%d'))
# data = data.fillna('')
data = data.sort_values('time',ascending=True)

# time1 = datetime.strptime('2024-5-31', '%Y-%m-%d')
# time2 = datetime.strptime('2024-11-30', '%Y-%m-%d')
# batch_list = pd.date_range(start=time1, end=time2, freq='M').tolist()
#mkdir results/venas2v6_output/haps_month/before2020-09
base_time = datetime.strptime('2020-08-31', '%Y-%m-%d')
results = data[data['time']<=base_time]
results[['samplename', 'haplotype', 'collection_date', 'location']].to_csv(data_path + 'before202009/samples_mutations.txt', sep='\t', index=False, header=True)

stats = []
data_time_max = datetime.strptime('2024-5-31', '%Y-%m-%d')
base_time = datetime.strptime('2020-1-31', '%Y-%m-%d')
batch_list = pd.date_range(start=base_time, end=data_time_max, freq='M').tolist()
batch_name_list = []
batch_times = len(batch_list)
for batch_idx in range( batch_times-1):
    # batch_idx = 5
    print('processing the batch ', str(batch_idx),'/' ,str(batch_times))
    batch_name = datetime.strftime( batch_list[batch_idx+1] ,'%Y-%m' )
    beging_time = batch_list[batch_idx]
    ending_time = batch_list[batch_idx+1]
    # print(batch_idx, beging_time, ending_time, batch_name)
    # data_batch = data[ (data['time']>batch_list[batch_idx]) & (data['time']<=batch_list[batch_idx+1]) ]
    data_batch = data[ (data['time']>beging_time) & (data['time']<=ending_time) ]
    print(batch_idx, beging_time, ending_time, batch_name, data_batch.shape[0])
    stats.append([batch_idx, beging_time, ending_time, batch_name, data_batch.shape[0]])
    batch_name_list.append(batch_name)
    output_path = data_path + batch_name + '/'
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    data_batch[['samplename', 'haplotype', 'collection_date', 'location']].to_csv(output_path+'samples_mutations.txt',sep='\t', header=True, index=None)

pd.DataFrame(stats).to_csv(data_path + 'samples_stats.txt',sep='\t', header=False, index=None)


end_time = time.time()
diff = end_time - beg_time
print(time.ctime())
print(diff, ' s')
print(diff/60, ' min')
print(diff/60/60, ' h')

