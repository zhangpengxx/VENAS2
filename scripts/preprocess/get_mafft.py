
#!/usr/bin/env python
# coding=utf-8
# __author__ = 'zp'

import pandas as pd
import numpy as np
import os
import re
import time

import multiprocessing

def do_mafft(infile: str, out_file: str):
    os.system("mafft --thread 1 --quiet " + infile + " > " + out_file)

p = multiprocessing.Pool(multiprocessing.cpu_count())

samples_df = pd.read_csv('rawdata/samples.txt', sep='\t', header=None)

# batch = 9
for sampleidx in range(samples_df.shape[0]):
# for sampleidx in range(batch * 10, (batch +1)* 10):
	# samplename = non_exist_samples_df[0].iloc[sampleidx]
	samplename = samples_df.iloc[sampleidx][0]
	print('processing ', sampleidx, samplename)
	input_file = 'rawdata/vigtk_genome_fa/' + samplename + '.fasta'
	target_file =  'rawdata/vigtk_genome_ma/' + samplename + '.ma'
	p.apply_async(do_mafft, args=(input_file, target_file))

p.close()
p.join()


# hap_ids = meta_df.index.tolist()
# hap_count = len(hap_ids)
# # cpu_count = int(multiprocessing.cpu_count())
# # cpu_count = int(args.cpu_count)
# cpu_count = 64
# batch_size = int(hap_count/cpu_count)

# p = multiprocessing.Pool(cpu_count)
# for i in range(cpu_count):
#     # i= 0
#     if i == cpu_count - 1:
#         hap_ids_batch = hap_ids[ batch_size * i: ]
#     else:
#         hap_ids_batch = hap_ids[batch_size * i : batch_size*(i+1) ]
#     batch_name = str(i)
#     # print(batch_name, hap_ids_batch[:9], can_ids_batch[:9])
#     p.apply_async(identify_wrongtime, args=( hap_ids_batch, batch_name ))

# p.close()
# p.join()

"""
mafft --thread 1 --quiet vigtk_genome_fa/OEAV19990326.fasta > vigtk_genome_ma/OEAV19990326.ma
mafft --thread 1 --quiet fa2/EPI_ISL_11531673.fa  > ma/EPI_ISL_11531673.ma
"""

