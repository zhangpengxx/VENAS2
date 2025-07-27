#!/usr/bin/env python
# coding=utf-8
# __author__ = 'zp'

import pandas as pd
import numpy as np
from datetime import datetime
import time
import argparse
import os
import random
from collections import Counter
import multiprocessing

beg_time = time.time()
# parser = argparse.ArgumentParser()
# parser.add_argument('-input_name', '--input_name', dest='input_name', help="the name of input samples mutation file", required=True)
# parser.add_argument('-cpu_count', '--cpu_count', dest='cpu_count', help="the cpu count used", default=5)
# args = parser.parse_args()

# # #parameters
data_path = 'venas2_output/'
haplotype_df = pd.read_csv(data_path+'haplotype_final.txt',sep='\t', header=0)
haplotype_df = haplotype_df.fillna('')
edges_df = pd.read_csv(data_path + 'edges.csv', sep=',', header=0)

## step 1: remove the circles in net 
#when merge the haps from multiple months, some edges may have two directions
target_list = sorted(set(edges_df['target'].tolist()))
can_list = range(edges_df.shape[0])

def remove_circles(can_list, batch_idx ):
    print('processing ', batch_idx)
    # edges_selected = []
    remove_list = []
    for eachidx in range(len(can_list)):
        # eachidx = 24
        # tmpdf = edges_df_main.iloc[[eachidx],:]
        eachtarget = can_list[eachidx]
        # eachtarget = target_list2[eachidx]
        # eachtarget = 101125
        # print('processing ', eachidx, eachtarget)
        tmpdf = edges_df[edges_df['target']==eachtarget]
        if tmpdf.shape[0] > 1:
            curremoves = [] #remove the two directions
            for tmpidx in tmpdf.index.tolist():
                source1, target1, mutations1, dis1 = tmpdf.loc[tmpidx]
                if source1 > target1:
                    curremoves.append(tmpidx)
            remove_list+= curremoves
            tmpdf = tmpdf.drop(curremoves)
            if tmpdf.shape[0] == 1:
                pass       
            elif tmpdf.shape[0]==2:
                # print(eachidx, tmpdf.shape)
                eachidx1, eachidx2 = tmpdf.index.tolist()
                source1, target1, mutations1, dis1 = tmpdf.iloc[0]
                source2, target2, mutations2, dis2 = tmpdf.iloc[1]
                bm_muts1 = len([ eachmut for  eachmut in mutations1.split(';') if eachmut.startswith('bm')])
                bm_muts2 = len([ eachmut for  eachmut in mutations2.split(';') if eachmut.startswith('bm')])
                if source1 == source2 and target1 == target2:
                    curremove = eachidx2
                elif bm_muts1 == 0 and bm_muts2 == 0: #both pos
                    if dis1 > dis2:
                        curremove = eachidx1
                    elif dis1<dis2:
                        curremove = eachidx2
                    else:
                        delta1 = target1 - source1
                        delta2 = target2 - source2
                        if delta1 >= delta2:
                            curremove = eachidx1
                        else:
                            curremove = eachidx2
                elif bm_muts1 == dis1 and bm_muts2 == dis2: #both neg
                    if dis1 > dis2:
                        curremove = eachidx1
                    elif dis1<dis2:
                        curremove = eachidx2
                    else:
                        delta1 = target1 - source1
                        delta2 = target2 - source2
                        if delta1 >= delta2:
                            curremove = eachidx1
                        else:
                            curremove = eachidx2
                elif bm_muts1 == 0: #1 both pos
                    curremove = eachidx2
                elif bm_muts2 == 0: #2 both neg
                    curremove = eachidx1
                elif bm_muts1 == dis1: #1 both neg
                    curremove = eachidx2
                elif bm_muts2 == dis2: #2 both neg
                    curremove = eachidx1
                else: #select the one whose pos count higher
                    pos_muts1 = dis1 - bm_muts1
                    pos_muts2 = dis2 - bm_muts2
                    if pos_muts1 >= pos_muts2:
                        curremove = eachidx2
                    else:
                        curremove = eachidx1
                    # print(eachidx, tmpdf)
                remove_list.append(curremove)
            elif tmpdf.shape[0]==3:
                # print(eachidx, tmpdf.shape)
                eachidx1, eachidx2, eachidx3 = tmpdf.index.tolist()
                source1, target1, mutations1, dis1 = tmpdf.iloc[0]
                source2, target2, mutations2, dis2 = tmpdf.iloc[1]
                source3, target3, mutations3, dis3 = tmpdf.iloc[1]
                bm_muts1 = len([ eachmut for  eachmut in mutations1.split(';') if eachmut.startswith('bm')])
                bm_muts2 = len([ eachmut for  eachmut in mutations2.split(';') if eachmut.startswith('bm')])
                bm_muts3 = len([ eachmut for  eachmut in mutations3.split(';') if eachmut.startswith('bm')])
                if bm_muts1 == 0: #1 pos
                    curremoves = [eachidx2, eachidx3 ]
                elif bm_muts2 == 0:
                    curremoves = [eachidx1, eachidx3 ]
                elif bm_muts3 == 0:
                    curremoves = [eachidx1, eachidx2 ]
                else:
                    tmpdf = tmpdf.sort_values('distance')
                    curremoves = tmpdf.index.tolist()[1:]
                remove_list += curremoves
            elif tmpdf.shape[0]>3:
                tmpdf = tmpdf.sort_values('distance')
                curremoves = tmpdf.index.tolist()[1:]
                remove_list += curremoves
            # print(eachidx, tmpdf)
    pd.DataFrame(remove_list).to_csv(data_path+'remove_list'+str(batch_idx)+'.txt', sep='t', index=False, header=False)

# target_list2 = target_list[:100]
hap_count = len(target_list)
cpu_count = int(multiprocessing.cpu_count())
cpu_count = 1
batch_size = int(hap_count/cpu_count)

p = multiprocessing.Pool(cpu_count)
for i in range(cpu_count):
    # i= 0
    if i == cpu_count - 1:
        hap_ids_batch = target_list[ batch_size * i: ]
    else:
        hap_ids_batch = target_list[batch_size * i : batch_size*(i+1) ]
    # print(hap_ids_batch)
    batch_name = str(i)
    # print(batch_name, hap_ids_batch[:9], can_ids_batch[:9])
    p.apply_async(remove_circles, args=( hap_ids_batch, batch_name ))

p.close()
p.join()

remove_df  = pd.DataFrame() 
for batch_idx in range(cpu_count):
    # batch_idx = 0
    print('processing ',batch_idx)
    curfile = data_path+'remove_list'+str(batch_idx)+'.txt'
    if os.path.getsize(curfile) == 0:
        pass
    else:
        tmpdf = pd.read_csv(curfile,sep='\t', header=None)
        remove_df = pd.concat([remove_df, tmpdf ], axis=0)
        os.remove(curfile)

if remove_df.shape[0]==0:
    remove_list = []
else:
    remove_list = remove_df[0].tolist()

# edges_df_main2 = edges_df_main.loc[edges_selected]
edges_df2 = edges_df.drop(remove_list)
edges_df2.index = range(edges_df2.shape[0])
# edges_df2.to_csv(data_path + 'edges_df2.csv', sep=',', index=False, header=True)
# edges_df2 = pd.read_csv(data_path + 'edges_df2.csv', sep=',', header=0)

##step 2: mimimize the sum of network
def identify_nodes_diff(hapidx1, hapidx2):
    # hapidx1 = 534996
    # hapidx2 = 544285
    if haplotype_df.loc[hapidx1]['haplotype']=='':
        a= set()
    else:
        a = set(haplotype_df.loc[hapidx1]['haplotype'].split(';'))
    b = set(haplotype_df.loc[hapidx2]['haplotype'].split(';'))
    t1= a.difference(b)
    t2 = b.difference(a)
    # t3 = a.intersection(b)
    diff_count = len(t1) + len(t2)
    diff_muts = list(t2) + [ 'bm:'+x for x in t1 ]
    diff_muts = ';'.join(diff_muts)
    # diff_count
    # diff_muts
    return diff_count, diff_muts

def reverse_muts(edges_df_fn, hapidxfn):
    # hapidxfn = 1197393
    source, target, mutations, dis = edges_df_fn.loc[hapidxfn]
    pos_muts = []
    # bm_muts = []
    bm_muts2 = []
    for  eachmut in mutations.split(';'):
        if eachmut.startswith('bm'):
            # bm_muts.append(eachmut)
            bm_muts2.append(eachmut[3:])
        else:
            pos_muts.append(eachmut)
    new_muts = ';'.join(bm_muts2)+ ';' + ';'.join([ 'bm:'+x for x in pos_muts])        
    output_fn = [target, source, new_muts, dis]
    return output_fn

def minimize_network1_select(can_list, batch_idx ):
    print('processing ', batch_idx)
    edges_selected = []
    # remove_list = []
    # can_list = range(10000)
    # batch_idx = '0'
    # edges_df_new = edges_df.loc[can_list].copy()
    index_list = []
    # edges_df_new = edges_df.copy()
    edges_df_new = []
    for eachidx in can_list:
        # print('processing ',eachidx)
        # source, target, mutations, dis = edges_df_new.loc[eachidx]
        source, target, mutations, dis = edges_df.loc[eachidx]
        # nodepre1, dispre1 = edges_df_main[edges_df_main['target']==source].iloc[0][['source','distance']]
        pos_muts = []
        bm_muts = []
        bm_muts2 = []
        for  eachmut in mutations.split(';'):
            if eachmut.startswith('bm'):
                bm_muts.append(eachmut)
                bm_muts2.append(eachmut[3:])
            else:
                pos_muts.append(eachmut)
        pos_count = len(pos_muts)
        bm_count = len(bm_muts)
        # bm_count = len(re.findall('bm',mutations))
        bm_rate = bm_count / dis
            #situation1 the reverse
        if bm_rate ==1:
        # if bm_rate ==1 or (bm_count >= 10 and dis>=15):
        # if bm_count > 0:
            # can_parents = []
            # can_parents_indexs = []
            # newsource = source
            # while newsource != 0:
            # pre1 = edges_df_new[edges_df_new['target']==source]
            pre1 = edges_df[edges_df['target']==source]
            # pre = edges_df_new[edges_df_new['target']==newsource]
            if pre1.shape[0]>0:
                # print(eachidx)
                # edges_selected.append(eachidx)
                sourcepre1, dispre1, mutspre1 = pre1.iloc[0][['source','distance', 'mutations']]
                diff_count1, diff_muts1 = identify_nodes_diff(sourcepre1, target )
                # pre2 = edges_df[edges_df['target']==sourcepre1]
                # sourcepre2, dispre2, mutspre2 = pre2.iloc[0][['source','distance', 'mutations']]
                # diff_count2, diff_muts2 = identify_nodes_diff(sourcepre2, target )
                first_flag = False 
                #the first order
                if diff_count1< dispre1:
                    edges_selected.append(eachidx)
                    # diff_count2, diff_muts2 = identify_nodes_diff(sourcepre1, target )
                    # edges_df_new.loc[pre1.index[0]] = [sourcepre1, target, diff_muts1, diff_count1 ]
                    # index_list.append(pre1.index[0])
                    # edges_df_new.append([sourcepre1, target, diff_muts1, diff_count1 ])
                    # #the current one
                    # cur_index= eachidx
                    # index_list.append(cur_index)
                    # each_reverse= reverse_muts(edges_df, cur_index)
                    # edges_df_new.append(each_reverse)
                    # edges_df_new.loc[eachidx] = each_reverse
                    first_flag = True
                # if first_flag == False and diff_count2< dispre2:
                #     index_list.append(pre2.index[0])
                #     edges_df_new.append([sourcepre2, target, diff_muts2, diff_count2 ])
                #     # second_flag = True
                #     #the middle one
                #     cur_index= pre1.index[0]
                #     index_list.append(cur_index)
                #     each_reverse= reverse_muts(edges_df, cur_index)
                #     edges_df_new.append(each_reverse)
                #     #the current one
                #     cur_index= eachidx
                #     index_list.append(cur_index)
                #     each_reverse= reverse_muts(edges_df, cur_index)
                #     edges_df_new.append(each_reverse)
                    # print('a')
                if bm_count>=5 and first_flag==False: #for the second order
                    # second_flag = False
                    # pre2 = edges_df_new[edges_df_new['target']==sourcepre1]
                    pre2 = edges_df[edges_df['target']==sourcepre1]
                    sourcepre2, dispre2, mutspre2 = pre2.iloc[0][['source','distance', 'mutations']]
                    diff_count2, diff_muts2 = identify_nodes_diff(sourcepre2, target )
                    if diff_count2 < dispre2:
                        # edges_df_new.loc[pre2.index[0]] = [sourcepre2, target, diff_muts2, diff_count2 ]
                        edges_selected.append(eachidx)
                        # index_list.append(pre2.index[0])
                        # edges_df_new.append([sourcepre2, target, diff_muts2, diff_count2 ])
                        # # second_flag = True
                        # #the middle one
                        # cur_index= pre1.index[0]
                        # index_list.append(cur_index)
                        # each_reverse= reverse_muts(edges_df, cur_index)
                        # edges_df_new.append(each_reverse)
                        # # edges_df_new.loc[pre1.index[0]] = each_reverse
                        # #the current one
                        # cur_index= eachidx
                        # index_list.append(cur_index)
                        # each_reverse= reverse_muts(edges_df, cur_index)
                        # edges_df_new.append(each_reverse)
                        # edges_df_new.loc[cur_index] = each_reverse
    if len(edges_selected) !=0:
        edges_selected = pd.DataFrame(edges_selected)
        # edges_df_new.columns = ['source', 'target', 'mutations', 'distance']
        # edges_df_new.index = index_list
        edges_selected.to_csv(data_path + 'edges_selected'+str(batch_idx)+'.csv', sep=',', index=False, header=False) 

# minimize_network1(hap_ids_batch, batch_name)
# target_list = sorted(set(edges_df['target'].tolist()))
edges_df = edges_df2.copy()
hapids = list(range(edges_df.shape[0]))
hap_count = len(hapids)
cpu_count = int(multiprocessing.cpu_count())
cpu_count = 1
batch_size = int(hap_count/cpu_count)

p = multiprocessing.Pool(cpu_count)
for i in range(cpu_count):
    # i= 0
    if i == cpu_count - 1:
        hap_ids_batch = hapids[ batch_size * i: ]
    else:
        hap_ids_batch = hapids[batch_size * i : batch_size*(i+1) ]
    # print(len(hap_ids_batch) )
    batch_name = str(i)
    # print(batch_name, hap_ids_batch[:9], can_ids_batch[:9])
    p.apply_async(minimize_network1_select, args=( hap_ids_batch, batch_name ))

p.close()
p.join()

##this belong to the minimize network1
## first select the edges which need to process because we need to process them at the same time

edges_df_selected  = pd.DataFrame() 
for batch_idx in range(cpu_count):
    # batch_idx = 1
    print('processing ',batch_idx)
    curfile = data_path+'edges_selected'+str(batch_idx)+'.csv'
    if os.path.exists(curfile):
        tmpdf = pd.read_csv(curfile, sep=',', header=None)
        edges_df_selected = pd.concat([edges_df_selected, tmpdf ], axis=0)
        os.remove(curfile)

beg_time = time.time()
for i in range(edges_df_selected.shape[0]):
# for eachidx in edges_df_selected[0].tolist()[:10]:
    # print('processing ',eachidx)
    # source, target, mutations, dis = edges_df_new.loc[eachidx]
    eachidx = edges_df_selected.iloc[i][0]
    print('processing ', i, '/', edges_df_selected.shape[0], eachidx)
    source, target, mutations, dis = edges_df.loc[eachidx]
    # nodepre1, dispre1 = edges_df_main[edges_df_main['target']==source].iloc[0][['source','distance']]
    pos_muts = []
    bm_muts = []
    bm_muts2 = []
    for  eachmut in mutations.split(';'):
        if eachmut.startswith('bm'):
            bm_muts.append(eachmut)
            bm_muts2.append(eachmut[3:])
        else:
            pos_muts.append(eachmut)
    pos_count = len(pos_muts)
    bm_count = len(bm_muts)
    # bm_count = len(re.findall('bm',mutations))
    bm_rate = bm_count / dis
        #situation1 the reverse
    if bm_rate ==1:
    # if bm_rate ==1 or (bm_count >= 10 and dis>=15):
    # if bm_count > 0:
        # can_parents = []
        # can_parents_indexs = []
        # newsource = source
        # while newsource != 0:
        # pre1 = edges_df_new[edges_df_new['target']==source]
        pre1 = edges_df[edges_df['target']==source]
        # pre = edges_df_new[edges_df_new['target']==newsource]
        if pre1.shape[0]>0:
            # print(eachidx)
            # edges_selected.append(eachidx)
            sourcepre1, dispre1, mutspre1 = pre1.iloc[0][['source','distance', 'mutations']]
            diff_count1, diff_muts1 = identify_nodes_diff(sourcepre1, target )
            # pre2 = edges_df[edges_df['target']==sourcepre1]
            # sourcepre2, dispre2, mutspre2 = pre2.iloc[0][['source','distance', 'mutations']]
            # diff_count2, diff_muts2 = identify_nodes_diff(sourcepre2, target )
            first_flag = False 
            #the first order
            if diff_count1< dispre1:
                # diff_count2, diff_muts2 = identify_nodes_diff(sourcepre1, target )
                edges_df.loc[pre1.index[0]] = [sourcepre1, target, diff_muts1, diff_count1 ]
                # index_list.append(pre1.index[0])
                # edges_df_new.append([sourcepre1, target, diff_muts1, diff_count1 ])
                #the current one
                cur_index= eachidx
                # index_list.append(cur_index)
                each_reverse= reverse_muts(edges_df, cur_index)
                # edges_df_new.append(each_reverse)
                edges_df.loc[eachidx] = each_reverse
                first_flag = True
            # if first_flag == False and diff_count2< dispre2:
            #     index_list.append(pre2.index[0])
            #     edges_df_new.append([sourcepre2, target, diff_muts2, diff_count2 ])
            #     # second_flag = True
            #     #the middle one
            #     cur_index= pre1.index[0]
            #     index_list.append(cur_index)
            #     each_reverse= reverse_muts(edges_df, cur_index)
            #     edges_df_new.append(each_reverse)
            #     #the current one
            #     cur_index= eachidx
            #     index_list.append(cur_index)
            #     each_reverse= reverse_muts(edges_df, cur_index)
            #     edges_df_new.append(each_reverse)
                # print('a')
            if bm_count>=5 and first_flag==False: #for the second order
                # second_flag = False
                # pre2 = edges_df_new[edges_df_new['target']==sourcepre1]
                pre2 = edges_df[edges_df['target']==sourcepre1]
                sourcepre2, dispre2, mutspre2 = pre2.iloc[0][['source','distance', 'mutations']]
                diff_count2, diff_muts2 = identify_nodes_diff(sourcepre2, target )
                if diff_count2 < dispre2:
                    edges_df.loc[pre2.index[0]] = [sourcepre2, target, diff_muts2, diff_count2 ]
                    # index_list.append(pre2.index[0])
                    # edges_df_new.append([sourcepre2, target, diff_muts2, diff_count2 ])
                    # second_flag = True
                    #the middle one
                    cur_index= pre1.index[0]
                    # index_list.append(cur_index)
                    each_reverse= reverse_muts(edges_df, cur_index)
                    # edges_df_new.append(each_reverse)
                    edges_df.loc[pre1.index[0]] = each_reverse
                    #the current one
                    cur_index= eachidx
                    # index_list.append(cur_index)
                    each_reverse= reverse_muts(edges_df, cur_index)
                    # edges_df_new.append(each_reverse)
                    edges_df.loc[cur_index] = each_reverse


# edges_df.to_csv(data_path + 'edges_df3.csv', sep=',', index=False, header=True)
# edges_df3 = pd.read_csv(data_path + 'edges_df3.csv', sep=',', header=0)

#part2  some nodes still have all bm
def minimize_network2(can_list, batch_idx ):
    print('processing ', batch_idx)
    index_list = []
    edges_df_new = []
    for eachidx in can_list:
        # print('processing ',eachidx)
        source, target, mutations, dis = edges_df.loc[eachidx]
        # nodepre1, dispre1 = edges_df_main[edges_df_main['target']==source].iloc[0][['source','distance']]
        pos_muts = []
        bm_muts = []
        bm_muts2 = []
        for  eachmut in mutations.split(';'):
            if eachmut.startswith('bm'):
                bm_muts.append(eachmut)
                bm_muts2.append(eachmut[3:])
            else:
                pos_muts.append(eachmut)
        pos_count = len(pos_muts)
        bm_count = len(bm_muts)
        # bm_count = len(re.findall('bm',mutations))
        bm_rate = bm_count / dis
        #situations where to descrease the bm count
        if (bm_count >= 10 and dis>=15):
            can_parents2 = []
            newsource = source
            loop_count = 0
            while newsource != 0 and loop_count < 20:
            # pre = edges_df[edges_df['target']==source]
            # newsource = pre.iloc[0]['source']
                pre = edges_df[edges_df['target']==newsource]
                # pre = edges_df_new[edges_df_new['target']==newsource]
                newsources = pre['source'].tolist()
                loop_count += 1
                for newsource in newsources:
                    # newsource = pre.iloc[0]['source']    
                    if newsource not in can_parents2:
                        can_parents2.append(newsource)
                        # print(newsource, target)
                        # each_parent = newsource
                        # diff_count, diff_muts = identify_nodes_diff(each_parent, target )
                        # eachbm_count = len([ x for x in diff_muts.split(';') if x.startswith('bm') ])
                        # if eachbm_count<bm_count and eachbm_count <=3:
                        #     # print( each_parent, diff_count, diff_muts)
                        #     selected_parent = each_parent
                        #     index_list.append(eachidx)
                        #     edges_df_new.append([selected_parent, target, diff_muts, diff_count ])
                        #     # edges_df_new.loc[eachidx] = [selected_parent, target, diff_muts, diff_count ]
                    else:
                        break
                # else:
                #   break
            for each_parent in can_parents2:
                # each_parent = newsource
                diff_count, diff_muts = identify_nodes_diff(each_parent, target )
                eachbm_count = len([ x for x in diff_muts.split(';') if x.startswith('bm') ])
                if eachbm_count<bm_count and eachbm_count <=3:
                    # print( each_parent, diff_count, diff_muts)
                    selected_parent = each_parent
                    index_list.append(eachidx)
                    edges_df_new.append([selected_parent, target, diff_muts, diff_count ])
                    # edges_df_new.loc[eachidx] = [selected_parent, target, diff_muts, diff_count ]
                    break
                # else:
                    # break
    if len(edges_df_new) !=0:
        edges_df_new = pd.DataFrame(edges_df_new)
        edges_df_new.columns = ['source', 'target', 'mutations', 'distance']
        edges_df_new.index = index_list
        edges_df_new.to_csv(data_path + 'edges_annotated'+str(batch_idx)+'part2.csv', sep=',', index=True, header=True)

# target_list = sorted(set(edges_df['target'].tolist()))
edges_df = edges_df.copy()
hapids = list(range(edges_df.shape[0]))
hap_count = len(hapids)
# cpu_count = int(multiprocessing.cpu_count)
cpu_count = 1
batch_size = int(hap_count/cpu_count)

p = multiprocessing.Pool(cpu_count)
for i in range(cpu_count):
    # i= 0
    if i == cpu_count - 1:
        hap_ids_batch = hapids[ batch_size * i: ]
    else:
        hap_ids_batch = hapids[batch_size * i : batch_size*(i+1) ]
    # print(len(hap_ids_batch), hap_ids_batch[:9], hap_ids_batch[-9:] )
    batch_name = str(i)
    # print(batch_name, hap_ids_batch[:9], can_ids_batch[:9])
    p.apply_async(minimize_network2, args=( hap_ids_batch, batch_name ))

p.close()
p.join()

edges_df_changed  = pd.DataFrame() 
for batch_idx in range(cpu_count):
    # batch_idx = 1
    print('processing ',batch_idx)
    curfile = data_path+'edges_annotated'+str(batch_idx)+'part2.csv'
    if os.path.exists(curfile):
        tmpdf = pd.read_csv(curfile, sep=',', index_col=0, header=0)
        edges_df_changed = pd.concat([edges_df_changed, tmpdf ], axis=0)
        os.remove(curfile)

edges_df_new2 = edges_df.copy()
edges_df_new2.loc[ edges_df_changed.index.tolist(), :] = edges_df_changed
edges_df_new2.index = range(edges_df_new2.shape[0])
# edges_df_new2.to_csv(data_path + 'edges_annotated_all.csv', sep=',', index=False, header=True)
# edges_df = pd.read_csv(data_path + 'edges_annotated_all.csv', sep=',', header=0)
# edges_df = edges_df.fillna('')

##we need to remove circles again because the preocess before may produce edges which have two directions
def remove_circles2(can_list, batch_idx ):
    print('processing ', batch_idx)
    # edges_selected = []
    remove_list = []
    for eachidx in range(len(can_list)):
        # eachidx = 24
        # tmpdf = edges_df_main.iloc[[eachidx],:]
        eachtarget = can_list[eachidx]
        # eachtarget = target_list2[eachidx]
        # eachtarget = 101125
        # print('processing ', eachidx, eachtarget)
        tmpdf = edges_df[edges_df['target']==eachtarget]
        if tmpdf.shape[0] > 1:
            curremoves = [] #remove the two directions
            for tmpidx in tmpdf.index.tolist():
                source1, target1, mutations1, dis1 = tmpdf.loc[tmpidx]
                if source1 > target1:
                    curremoves.append(tmpidx)
            remove_list+= curremoves
            tmpdf = tmpdf.drop(curremoves)
            if tmpdf.shape[0] == 1:
                pass       
            elif tmpdf.shape[0]==2:
                # print(eachidx, tmpdf.shape)
                eachidx1, eachidx2 = tmpdf.index.tolist()
                source1, target1, mutations1, dis1 = tmpdf.iloc[0]
                source2, target2, mutations2, dis2 = tmpdf.iloc[1]
                bm_muts1 = len([ eachmut for  eachmut in mutations1.split(';') if eachmut.startswith('bm')])
                bm_muts2 = len([ eachmut for  eachmut in mutations2.split(';') if eachmut.startswith('bm')])
                if source1 == source2 and target1 == target2:
                    curremove = eachidx2
                    remove_list.append(curremove)
                    break
                elif dis1 == 0:
                    curremove = eachidx2
                    remove_list.append(curremove)
                    break
                elif dis2 == 0:
                    curremove = eachidx1
                    remove_list.append(curremove)
                    break
                elif bm_muts1 == 0 and bm_muts2 == 0: #both pos
                    if dis1 > dis2:
                        curremove = eachidx1
                    elif dis1<dis2:
                        curremove = eachidx2
                    else:
                        delta1 = target1 - source1
                        delta2 = target2 - source2
                        if delta1 >= delta2:
                            curremove = eachidx1
                        else:
                            curremove = eachidx2
                    remove_list.append(curremove)
                    break
                elif bm_muts1 == dis1 and bm_muts2 == dis2: #both neg
                    if dis1 > dis2:
                        curremove = eachidx1
                    elif dis1<dis2:
                        curremove = eachidx2
                    else:
                        delta1 = target1 - source1
                        delta2 = target2 - source2
                        if delta1 >= delta2:
                            curremove = eachidx1
                        else:
                            curremove = eachidx2
                    remove_list.append(curremove)
                    break
                elif bm_muts1 == 0: #1 both pos
                    curremove = eachidx2
                    remove_list.append(curremove)
                    break
                elif bm_muts2 == 0: #2 both neg
                    curremove = eachidx1
                    remove_list.append(curremove)
                    break
                elif bm_muts1 == dis1: #1 both neg
                    curremove = eachidx2
                    remove_list.append(curremove)
                    break
                elif bm_muts2 == dis2: #2 both neg
                    curremove = eachidx1
                    remove_list.append(curremove)
                    break
                else: #select the one whose pos count higher
                    pos_muts1 = dis1 - bm_muts1
                    pos_muts2 = dis2 - bm_muts2
                    if pos_muts1 >= pos_muts2:
                        curremove = eachidx2
                    else:
                        curremove = eachidx1
                    # print(eachidx, tmpdf)
                    remove_list.append(curremove)
                    break
            # print(eachidx, tmpdf)
    pd.DataFrame(remove_list).to_csv(data_path+'remove_list'+str(batch_idx)+'.txt', sep='t', index=False, header=False)

edges_df =  edges_df_new2.copy()
target_list = list(range(edges_df.shape[0]))
# target_list2 = target_list[:100]
hap_count = len(target_list)
# cpu_count = int(multiprocessing.cpu_count)
cpu_count = 1
batch_size = int(hap_count/cpu_count)

p = multiprocessing.Pool(cpu_count)
for i in range(cpu_count):
    # i= 0
    if i == cpu_count - 1:
        hap_ids_batch = target_list[ batch_size * i: ]
    else:
        hap_ids_batch = target_list[batch_size * i : batch_size*(i+1) ]
    # print(hap_ids_batch)
    batch_name = str(i)
    # print(batch_name, hap_ids_batch[:9], can_ids_batch[:9])
    p.apply_async(remove_circles2, args=( hap_ids_batch, batch_name ))

p.close()
p.join()

remove_df  = pd.DataFrame() 
for batch_idx in range(cpu_count):
    # batch_idx = 0
    print('processing ',batch_idx)
    curfile = data_path+'remove_list'+str(batch_idx)+'.txt'
    if os.path.getsize(curfile) == 0:
        pass
    else:
        tmpdf = pd.read_csv(curfile,sep='\t', header=None)
        remove_df = pd.concat([remove_df, tmpdf ], axis=0)
        os.remove(curfile)

if remove_df.shape[0]==0:
    remove_list = []
else:
    remove_list = remove_df[0].tolist()

# edges_df_main2 = edges_df_main.loc[edges_selected]
edges_df2 = edges_df.drop(remove_list)
edges_df2.index = range(edges_df2.shape[0])
edges_df2.to_csv(data_path + 'edges_optimized.csv', sep=',', index=False, header=True)

end_time = time.time()
diff = end_time - beg_time
print(time.ctime())
print(diff, ' s')
print(diff/60, ' min')
print(diff/60/60, ' h')


