#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np


# # editing position<0

# In[2]:


df = pd.read_table('./raw_pathogenic.csv', sep = ',') 


# In[3]:


df


# In[4]:


fo = open('./raw_pathogenic.csv')
i=0

ID = []
gRNA = []
Sequence = []
Position = []
Pathogenic = []

for line in fo.readlines():  
    if (i!=0):
        seq = line.strip()
        temp = seq.split(',')   
        sequence = temp[2]
        gRNA_seq = temp[1]
        if int(temp[3])>0:           
            ID.append(temp[0])
            gRNA.append(gRNA_seq)
            Sequence.append(sequence)
            Position.append(temp[3])
            Pathogenic.append(temp[5])

        elif int(temp[3])<=0:                      
            start = sequence.index(gRNA_seq)
            new_start = start+int(temp[3])-1
            
            if sequence[new_start+22-5]=='G':
                ID.append(temp[0]+"_neg")
                gRNA.append(sequence[new_start-4:new_start+16])
                Sequence.append(sequence)
                Position.append(5)
                Pathogenic.append(temp[5])
                
            elif sequence[new_start+22-6]=='G':
                ID.append(temp[0]+"_neg")
                gRNA.append(sequence[new_start-5:new_start+15])
                Sequence.append(sequence)
                Position.append(6)
                Pathogenic.append(temp[5])
                
            elif sequence[new_start+22-7]=='G':
                ID.append(temp[0]+"_neg")
                gRNA.append(sequence[new_start-6:new_start+14])
                Sequence.append(sequence)
                Position.append(7)
                Pathogenic.append(temp[5])
                
            elif sequence[new_start+22-4]=='G':
                ID.append(temp[0]+"_neg")
                gRNA.append(sequence[new_start-3:new_start+17])
                Sequence.append(sequence)
                Position.append(4)
                Pathogenic.append(temp[5])                
                
            elif sequence[new_start+22-8]=='G' and len(sequence[new_start-7:new_start+13])>0:
                ID.append(temp[0]+"_neg")
                gRNA.append(sequence[new_start-7:new_start+13])
                Sequence.append(sequence)
                Position.append(8)
                Pathogenic.append(temp[5])                  

    i+=1
        


# In[5]:


fish_frame = pd.DataFrame({'ID':ID,
                           'gRNA':gRNA,
                           'Sequence':Sequence,
                           'Position':Position,
                          'Pathogenic':Pathogenic,
                         
                            
                          })    


# In[6]:


fish_frame


# In[7]:


fish_frame.to_csv(r'./raw_pathogenic_positive.csv', header=True, index=None, sep=',')


# In[7]:


# raw


# In[8]:


# count=0
# fifth = []
# sixth = []
# seventh = []
# fourth = []
# eighth = []

# for index, row in raw.iterrows():
#     pos = int(row['Editing position'])

#     if pos<0:
#         count+=1
#         sequence = row['Sequence context (56nt)']
#         gRNA = row['gRNA (20nt)']
#         start = sequence.index(gRNA)
#         new_start = start+pos-1
# #         print(start)
# #         print(sequence[new_start])
        
#         if sequence[new_start+22-5]=='G':
#             raw.at[index, 'gRNA (20nt)']=sequence[new_start-4:new_start+19]
#             raw.at[index, 'Editing position']= 5
# #             fifth.append(sequence[new_start-4:new_start+19])
#         elif sequence[new_start+22-6]=='G':
#             raw.at[index, 'gRNA (20nt)']=sequence[new_start-5:new_start+18]
#             raw.at[index, 'Editing position']= 6
#             #             sixth.append(sequence[new_start-5:new_start+18])
#         elif sequence[new_start+22-7]=='G':
#             raw.at[index, 'gRNA (20nt)']=sequence[new_start-6:new_start+17]
#             raw.at[index, 'Editing position']= 7
#             #             seventh.append(sequence[new_start-6:new_start+17])
#         elif sequence[new_start+22-4]=='G':
#             raw.at[index, 'gRNA (20nt)']=sequence[new_start-3:new_start+20]
#             raw.at[index, 'Editing position']= 4
#             #             fourth.append(sequence[new_start-3:new_start+20]) 
#         elif sequence[new_start+22-8]=='G':
#             if len(sequence[new_start-7:new_start+16])>0:
#                 raw.at[index, 'gRNA (20nt)']= sequence[new_start-7:new_start+16]
#                 raw.at[index, 'Editing position']= 8
#             else:
#                 raw.at[index, 'gRNA (20nt)'] = np.nan
#         else:
#              raw.at[index, 'gRNA (20nt)'] = np.nan
# #             print(new_start)
# #             print(sequence[new_start-7:new_start+16])
# #             eighth.append(sequence[new_start-7:new_start+16])            
        
# #     if row['gRNA (20nt)'][5]=='A':
        
# #     print (row['gRNA (20nt)'], row['gRNA (20nt)'][5])


# In[9]:


# for index, row in raw.iterrows():
#     pos = int(row['Editing position'])

#     if pos<0:
#         count+=1


# In[50]:


raw.dropna()


# In[25]:


fifth


# In[26]:


sixth


# In[27]:


seventh


# In[28]:


fourth


# In[29]:


eighth


# In[32]:


len(fifth)+len(sixth)+len(seventh)+len(fourth)+4


# In[51]:


raw.to_csv(r'./raw_pathogenic_positive.csv', sep=',')


# In[ ]:


for index, row in raw.iterrows():
    pos = int(row['Editing position'])
    if pos = 6:
        print

