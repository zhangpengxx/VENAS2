#!/bin/bash

date
script_beg_time=`date '+%s'`
echo 'script Begin....'

# for ((batch=102; batch<=102; batch++))
# do
batch=$1
merged_ma_file='results/venas2v1_output/samples_10w_ma/samples_10w_'$batch'.ma' 
# merged_ma_file='results/venas2v1_output/samples_10w_ma/samples_10w_test.ma' 
flag=0
touch $merged_ma_file
cat 'results/venas2v1_output/samples_10w_ma/samples_10w_'$batch'.txt' | while read line
# cat 'results/venas2v1_output/samples_10w_ma/test.txt' | while read line
do
	echo 'batch ' $batch ' flag ' $flag ' sample ' $line
	# cat  'rawdata/pango-designation-master20241212/recal_ma/'$line'.ma' >> $merged_ma_file 
	cat  'rawdata/vigtk_genome_qc/sample_ma/'$line'.ma' >> $merged_ma_file 
	# cat  'results/venas2v3_output/samples_10w/'$line'.txt' >> $merged_ma_file 
	flag=$[flag+1]
done
# cat rawdata/vigtk_genome_qc/sample_ma/OEAV000245.ma >> results/venas2v3_output/1000_msa.ma
echo 'all script end.....'
script_end_time=`date '+%s'`
diff=$[ script_end_time - script_beg_time ]
echo $diff
date

# done