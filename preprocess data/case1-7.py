#!/usr/bin/env python
# coding: utf-8

# In[1]:


import re
from Bio import pairwise2

def alignemnt_S(seq1,seq2,total_count,GN20,target,case1,case2,case3,case4,case5,case6,case7): 
    total_count[GN20]+=1  
    tar_GCT,tar_G,bys_GCT,bys_G = 0,0,0,0
    idxes = re.finditer("\w", seq1)
    j=0
    for i in idxes:
        j+=1
        if j >= 2 and j <= 12:
            idx = i.span()[0]
            #case target
            if j == target:
                if seq2[idx] != 'A':
                    tar_GCT = 1
                    if seq2[idx] == 'G':
                        tar_G = 1    
            #case bystander 
            if j!=target:
                if seq1[idx] == 'A' and seq2[idx] != 'A':
                    bys_GCT = 1 
                    if seq2[idx] == 'G':
                        bys_G = 1 
     #1. target A->G  bystander A->A   
    if tar_G == 1 and bys_GCT == 0:        
        case1[GN20]+=1
    #2. target A->A  bystander A->A 
    if tar_GCT == 0 and bys_GCT == 0:        
        case2[GN20]+=1
#         print(seq2,tar_GCT,bys_GCT)
    #3. target A->G  any bystander A->G
    if tar_G == 1 and bys_G == 1:        
        case3[GN20]+=1
    #4. target A->A  any bystander A->G
    if tar_GCT == 0 and bys_G == 1:        
        case4[GN20]+=1
    #5. target A->G  any bystander A->!A
    if tar_G == 1 and bys_GCT == 1:        
        case5[GN20]+=1
    #6. target A->A  any bystander A->!A
    if tar_GCT == 0 and bys_GCT == 1:        
        case6[GN20]+=1       
    #7. target A->!A  bystander A->A
    if tar_GCT == 1 and bys_GCT == 0:        
        case7[GN20]+=1 
    return total_count,case1,case2,case3,case4,case5,case6,case7


# In[2]:


import pandas as pd

reference = pd.read_csv("./SORT_new_original_replace_nonreplace_pathogenic.csv", index_col=None, header=0,sep = ',',encoding="utf-8")
GN20_list,N20_list,target_list = list(reference['GN20']),list(reference['N20']),list(reference['Position'])
N20_dict,total_count,target_dict,case1,case2,case3,case4,case5,case6,case7 = {},{},{},{},{},{},{},{},{},{}
for i in zip(GN20_list,N20_list):
    N20_dict.setdefault(i[0],i[1])
    case1.setdefault(i[0],0)
    case2.setdefault(i[0],0)
    case3.setdefault(i[0],0)
    case4.setdefault(i[0],0)
    case5.setdefault(i[0],0)
    case6.setdefault(i[0],0)
    case7.setdefault(i[0],0)
    total_count.setdefault(i[0],0)
for i in zip(GN20_list,target_list):
    target_dict.setdefault(i[0],i[1])
    
barcodes = ['A1','C4','C6']
for barcode in barcodes:
    flank = open("./1_%s_N20.txt"% barcode)

    i = 0
    for line in flank:
        line = line.strip().split(',')
        GN20,N20 = line[0],line[1]

    #     alignments = pairwise2.align.globalms(N20_dict[GN20], N20,5, -4, -3, -.1,one_alignment_only = True)
        alignments = pairwise2.align.globalms(N20_dict[GN20], N20,5, -4, -5, -.1,one_alignment_only = True)
        target = target_dict[GN20]
        total_count,case1,case2,case3,case4,case5,case6,case7 = alignemnt_S(alignments[-1][0],alignments[-1][1],total_count,GN20,target,case1,case2,case3,case4,case5,case6,case7)

        if i%1000000==0:
            print(i)
        i+=1
    flank.close()

    d={
    'case1': case1,
    'case2':case2,
    'case3': case3,
    'case4':case4,
    'case5': case5,
    'case6':case6,
    'case7':case7,
    'total_count':total_count

    }
    df=pd.DataFrame(d)
    df.to_csv('N20count_%s.csv'% barcode)


# In[ ]:




