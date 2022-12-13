#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd


# In[2]:


raw = pd.read_csv("./raw.csv", index_col=None, header=0,sep = ',',encoding="utf-8")


# In[3]:


raw


# In[4]:


raw = raw[['ID','gRNA (20nt)','Sequence context (56nt)','Editing position','Pathogenic nucleotide (gRNA orientation)','CLNSIG',]]


# In[5]:


raw


# In[6]:


raw = raw.dropna()


# In[7]:


raw['CLNSIG'].unique()


# In[8]:


raw.shape


# In[9]:


raw_pathogenic = raw[(raw['CLNSIG']!='Likely_pathogenic')]
raw_pathogenic = raw_pathogenic[(raw_pathogenic['CLNSIG']!='Likely_pathogenic,_drug_response')]


# In[11]:


raw_pathogenic


# In[12]:


raw_pathogenic.to_csv(r'./raw_pathogenic.csv', sep=',',index = False)


# In[ ]:




