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

beg_time = time.time()

meta_df =pd.read_csv('rawdata/metadata_2024.tsv', sep='\t', header=0)
# meta_df = meta_df.fillna('')
# dex(['Virus name', 'Last vaccinated', 'Passage details/history', 'Type',
#        'Accession ID', 'Collection date', 'Location',
#        'Additional location information', 'Sequence length', 'Host',
#        'Patient age', 'Gender', 'Clade', 'Pango lineage', 'Pango version',
#        'Variant', 'AA Substitutions', 'Submission date', 'Is reference?',
#        'Is complete?', 'Is high coverage?', 'Is low coverage?', 'N-Content',
#        'GC-Content']

# problem_df = pd.read_csv('/bdp-picb/student/zhangpeng/venaszp/rawdata/problematic_sites_sarsCov2.vcf', sep='\t', skiprows=88, header=0)
# problem_pos = problem_df['POS'].tolist()
# selected_pos = [ x for x in range(29903) if x+1 not in problem_pos]
selected_pos = [ x for x in range(100, 29803)]

# for batch in range(1300, 1617):
# for batch in range(10,11):
batch = 10
print('processing ',str(batch))
data_path = 'rawdata/vigtk_genome_qc/sample_ma/'
non_samples_list = []
zero_samples_list = []
sample_info_df = []
# for sampleidx in range(batch * 10, (batch +1)* 10):
for sampleidx in range(batch * 1000000 + 500000, (batch+1) * 1000000):
	# sampleidx = 0
	samplename = meta_df.iloc[sampleidx]['virus_id']
	# samplename = other_df.iloc[sampleidx]['samplename']
	# if sampleidx // 1 == 0:
	print('processing ', sampleidx, samplename)
	ma_file =  data_path + samplename + '.ma'
	if os.path.exists(ma_file):
		seqrecords_fn = list(SeqIO.parse(ma_file, "fasta"))
		sampleseq = seqrecords_fn[0].seq
		seq_df = pd.DataFrame(list(sampleseq.lower()))
		seq_df.columns = ['sample']
		# seq_new = str(sampleseq)
		# seq_orig = ''.join(seq_df['sample'])
		seq_df = seq_df.iloc[selected_pos]
		seq_new = ''.join(seq_df['sample'])
		seq_new = seq_new.lower()
		# deletion_count = seq_df[seq_df['sample']=='-'].shape[0]
		# unknown_count = seq_df[seq_df['sample']=='n'].shape[0]
		# seq_len = seq_df[(seq_df['sample']=='a')|(seq_df['sample']=='t')|(seq_df['sample']=='c')|(seq_df['sample']=='g') ].shape[0]
		# ambiguous_count = len(selected_pos) - deletion_count - unknown_count - seq_len
		deletion_count = len(re.findall('-',seq_new))
		unknown_count = len(re.findall('[n]',seq_new))
		# ambiguous_count = len(re.findall('[^atcgn-]',seq_new))
		seq_len = len(re.findall('a',seq_new)) + len(re.findall('t',seq_new)) + len(re.findall('c',seq_new)) + len(re.findall('g',seq_new))
		ambiguous_count = len(selected_pos) - deletion_count - unknown_count - seq_len
		sample_info_df.append([samplename, seq_len, deletion_count, unknown_count, ambiguous_count ])
		##write the file to merge
		# f = open(output_path + samplename + '.ma', 'w' )
		# f.write('>' + samplename + '\n')
		# f.write(seq_orig + '\n')
		# f.close()

sample_info_df = pd.DataFrame(sample_info_df)
sample_info_df.columns = ['samplename', 'seq_len', 'deletion_count', 'unknown_count', 'ambiguous_count']
sample_info_df.to_csv('rawdata/vigtk_genome_qc/sample_qc2/sample' +str(batch)+'.txt',sep='\t', index=False, header=True)
# sample_info_df.to_csv('rawdata/vigtk_genome_qc/sample_qc2/sample10_2.txt',sep='\t', index=False, header=True)

end_time = time.time()
diff = end_time - beg_time
print(time.ctime())
print(diff, ' s')
print(diff/60, ' min')
print(diff/60/60, ' h')



