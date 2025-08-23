# VENAS2

## Introduction

VENAS2：an updated Viral genome Evolution Network Analysis System.
Comprehensive analyses of viral genomes can provide a global picture of SARS-CoV-2 transmission and help to predict the oncoming trends of the pandemic. However, the rapid accumulation of SARS-CoV-2 genomes presents an unprecedented data size and complexity that has exceeded the capacity of existing methods in constructing evolution network through virus genotyping. The VENAS seeks to apply reliable computational algorithms to build an integrative genomic analysis system that enables researchers to trace viral mutations along the transmission routes using the daily updated SARS-CoV-2 genomes.


## Pre-requisites
VENAS2 requires python3 with following third-party library/program:
argparse
pandas
numpy
faiss
networkx
CDlib
pango 4.3.1

If you want to provide a fasta file as input file, VENAS2 also needs the MAFFT (<https://mafft.cbrc.jp/alignment/software/>) in the executable path. Then you can use the “multi_mafft.py” to perform a multi-threaded multiple sequence alignment.


## Installation
Cloning the repository via the following commands 
```
$ git clone https://github.com/zhangpengxx/VENAS2.git
```

## Basic Usage
We assume here that all the scripts are in the system path.


### Part 1: Construct the haplotype network 
We assume you have the samples_mutations.txt. You need to prepare the file following the preprocess part if you do not have this file!!!
The samples_mutations.txt contains the mutation of samples. The format is defined as following:
```
<AccessionID>\t<Haplotype>\t<Collection date>\t<Location>
```
e.g. 
EPI_ISL_784970	C241T;C1059T 2020-01-01	China
EPI_ISL_784972	C241T;C1059T,del2448_2449 2020-01-02 USA


```
python3 scripts/venas2_contruct.py
```
**Results Description:**
*	haplotype_final.txt: The file contains the haplotypes information including Accession ID, mutations, collection date, location and so on.
*	edges.csv: The edges between all haplotypes. The format is as following:
```
<Source>\t<End>\t<Mutations>\t<distance>
0	1 	del2448_2449 	1
```
The source and end mean the haplotype id. The direction of haplotype evolution is from source to end. The source-target pair are parent-child edge, which means the target is evolved from the source haplotype with corresponding mutations. The mutations and distance represent the difference and their count between source and end.


**If you have large scale dataset, for example >1million SARS-CoV-2 sequences, please construct the network following scripts in construct_network_largescale foder.**


### Part 2: Optimize the haplotype network 

remove the nodes have multiplde parents
adjust the edges with reversion mutaions which lead to the edges summary are big

```
python3 scripts/venas2_optimize.py
```
**Results Description:**
*	edges_optimized.csv: The edges


### Part 3: Annotate the haplotype network
Identify the topology communities and main path of haplotype network. This also helps to draw the network.
```
python3 scripts/venas2_annotate.py
```
**Results Description:**
*	nodes_annotated.csv: This contains the information of all haplotypes. The node id represent the haplotype id in haplotype_final.txt. The nodel value with 1000 and 0 means the main and normal nodes in the network, respectively. The node cluster and cluster len stand for the cluster and cluster size of the node. The Lineage indicates the pangolin lineage.

You can visualize the network by force-directed graph tools like Gephi (recommend) using the edges.csv and nodes_annotated.csv files.


### Preprocess the data
We have designed the process scripts for two types of input. The first is multiple sequence aligned file in ma format, and the second is variants file in vcf format. The output is the samples_mutations.txt file described in part 1. 
The reference genome used in this work is OEAV139851. 

```
python get_input_from_msa.py
python get_input_from_vcf.py
```

If you only have fasta files we recommend you process the file as following:
(1). Split the fasta samples into a folder, each sample have one fasta file and its first line is the reference sequence OEAV139851/NC_045512v2.

(2). Map each sequence to reference by Mafft.
python3 scripts/get_mafft.py

(3). sample quality control - Calculate the sequence length, deletion count, ambious count
python3 scripts/sample_qc.py

(4). Continue quality control, identfy the main samples.
#Remove samples which did not have lineage, location, date
#whose date are not year-month-day
#whose host are not human
python3 scripts/get_samples.py

(5). Merge the corresponding samples ma files into msa, if the samples too large you can use batch.
bash scripts/get_mergema.sh 


## Examples
We provide an example input file samples_mutations.txt and its corresponding metadata samples_metadata.txt for 1000 SARS-CoV-2 sequence samples. 
To help the users understand more about our preprocess part, we also provide the samples in ma and vcf format in the example foler, respectively. Tips: unzip the vcffiles.zip first.


## Contact us
If you have any questions, please feel free to contact us. 
zhangpeng@sinh.ac.cn
lingyunchao@sinh.ac.cn
rfcao@sinh.ac.cn

## About us
Bio-Med Big Data Center, CAS Key Laboratory of Computational Biology, Shanghai Institute of Nutrition and Health, Chinese Academy of Sciences

## Citation

