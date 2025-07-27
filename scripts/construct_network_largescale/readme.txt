This folder is only prepared for large scale dataset, for example >1million SARS-CoV-2 sequences.

In general, we split the large samples into different month, then construct network for each month separately and merge all month results together.

1. Split the samples according to month
python3 get_hap_month.py 

2. Get the haplotype of samples in each month separately to accelerate the process
python3 get_hap_month_supp.py 

3. Construct the first network using early emerged samples. We use the samples before 2020-09 to construct the first network. This script is the same as venas2_construct.py.
python3 venas2_base.py

4. Construct the network for each month. The network only focus on the samples of in current moth and its previous four month.
python3 venas2_month.py 

5. Merge the network from each month into the final network together. 
python3 venas2_month_merge.py 







