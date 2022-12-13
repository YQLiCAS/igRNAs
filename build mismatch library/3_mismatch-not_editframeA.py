#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import random


# 1: 2-6,替换一个不为编辑点A的，为非T非本身，如果是3C，不改为G <br>
# 2: 2-6,替换一个不为编辑点A的，也不为T的，为T<br>
# 3：2-6，删除一个不为编辑点A的<br>
# 4：2-6，插入一个<br>
# 5: 7-9,替换一个不为编辑点A的，为非T非本身，如果是3C，不改为G<br>
# 6: 7-9,替换一个不为编辑点A的，也不为T的，为T<br>
# 7：7-9，删除一个不为编辑点A的<br>
# 8：7-9，插入一个<br>
# 9: 2-9, 改两个

# In[2]:


# df = pd.read_table('./raw_pathogenic.csv', sep = ',') 


# In[3]:


# df = pd.read_table('./raw_pathogenic_positive_6,7_2A.csv', sep = ',') 


# # 1: 2-6,替换一个不为编辑框内A的或编辑点A的，为非T非本身，如果是3C，不改为G

# In[4]:


fo = open('./raw_pathogenic_positive_6,7_2A.csv')


# In[5]:


def r(i,k):
    i.remove(k)
    return i


# In[6]:


"sgwregsrh".count('r')


# In[7]:


ID = []
gRNA = []
Sequence = []
Position = []
Pathogenic = []
Site = []
Old_sgRNA = []

i=0
for line in fo.readlines():
    
    ##2-6
    for site in range(2,7): 
        
        if (i!=0):
            seq = line.strip()
            temp = seq.split(',')   
            sequence = temp[2]
            gRNA_seq = temp[1]
            pos = temp[3]

            ##编辑框内A不编辑
            if site < 9 and site > 3 and gRNA_seq[site-1]=='A':
                continue;
                
            ##不为编辑点A     
            if int(site)!=int(pos): 
                ##C3G
                if site==3 and gRNA_seq[site-1]=='C': 
                    Base = 'A'

                else:
                ##非T
                    Bases = ['A','G','C']    
                    k = gRNA_seq[site-1] 
                    ##非本身
                    if k != 'T':
                        Bases = r(Bases,k)
                    for Base in Bases:                                    
                        ID.append(temp[0]+"_replaced_"+str(site)+'by_'+Base)
                        new_gRNA = gRNA_seq[0:site-1]+Base+gRNA_seq[site-1+1:]
                        gRNA.append(new_gRNA)
                        Old_sgRNA.append(gRNA_seq)
                        Site.append(site)
                        Position.append(temp[3])
                        Sequence.append(sequence)
                        Pathogenic.append(temp[4])

    i+=1


# In[8]:


fish_frame1 = pd.DataFrame({'ID':ID,
                            'gRNA':gRNA, 
                            'Sequence':Sequence, 
                            'Position':Position,
                            'Pathogenic':Pathogenic,
                            'Site':Site,
                            'Old_sgRNA':Old_sgRNA,})   


# In[9]:


fish_frame1


# In[10]:


fish_frame1.to_csv(r'./fish_frame1.csv', encoding='utf_8_sig',header=True, index=None, sep=',')


# # 2: 2-6,替换一个不为编辑点A的或不为编辑框内A的，也不为T的，为T

# In[11]:


fo = open('./raw_pathogenic_positive_6,7_2A.csv')


# In[12]:


ID = []
gRNA = []
Sequence = []
Position = []
Pathogenic = []
Site = []
Old_sgRNA = []

i=0
for line in fo.readlines():
    
    ##2-6
    for site in range(2,7): 
        if (i!=0):
            seq = line.strip()
            temp = seq.split(',')   
            sequence = temp[2]
            gRNA_seq = temp[1]
            pos = temp[3]
            
            ##编辑框内A不编辑
            if site < 9 and site > 3 and gRNA_seq[site-1]=='A':
                continue;
                
            ##不为编辑点A     
            if int(site)!=int(pos): 
                Base = 'T'   
                k = gRNA_seq[site-1] 
                ##非本身
                if k != 'T':                                  
                    ID.append(temp[0]+"_replaced_"+str(site)+'by_'+Base)
                    new_gRNA = gRNA_seq[0:site-1]+Base+gRNA_seq[site-1+1:]
                    gRNA.append(new_gRNA)
                    Old_sgRNA.append(gRNA_seq)
                    Site.append(site)
                    Position.append(temp[3])
                    Sequence.append(sequence)
                    Pathogenic.append(temp[4])

    i+=1


# In[13]:


len(ID)


# In[14]:


fish_frame2 = pd.DataFrame({'ID':ID,
                            'gRNA':gRNA, 
                            'Sequence':Sequence, 
                            'Position':Position,
                            'Pathogenic':Pathogenic,
                            'Site':Site,
                            'Old_sgRNA':Old_sgRNA,})   


# In[15]:


fish_frame2


# In[16]:


fish_frame2.to_csv(r'./fish_frame2.csv', encoding='utf_8_sig',header=True, index=None, sep=',')


# # 3：2-6，删除一个不为编辑点A的，或编辑框内A

# In[17]:


fo = open('./raw_pathogenic_positive_6,7_2A.csv')

ID = []
gRNA = []
Sequence = []
Position = []
Pathogenic = []
Site = []
Old_sgRNA = []

i=0
for line in fo.readlines():
    
    ##2-6
    for site in range(2,7): 
        if (i!=0):
            seq = line.strip()
            temp = seq.split(',')   
            sequence = temp[2]
            gRNA_seq = temp[1]
            pos = temp[3]
            
            ##编辑框内A不编辑
            if site < 9 and site > 3 and gRNA_seq[site-1]=='A':
                continue;
                
            ##不为编辑点A     
            if int(site)!=int(pos): 
                Base = ''                                 
                ID.append(temp[0]+"_Indel_"+str(site))
                new_gRNA = gRNA_seq[0:site-1]+Base+gRNA_seq[site-1+1:]
                gRNA.append(new_gRNA)
                Old_sgRNA.append(gRNA_seq)
                Site.append(site)
                Position.append(temp[3])
                Sequence.append(sequence)
                Pathogenic.append(temp[4])

    i+=1


# In[18]:


fish_frame3 = pd.DataFrame({'ID':ID,
                            'gRNA':gRNA, 
                            'Sequence':Sequence, 
                            'Position':Position,
                            'Pathogenic':Pathogenic,
                            'Site':Site,
                            'Old_sgRNA':Old_sgRNA,})   


# In[19]:


fish_frame3


# In[20]:


fish_frame3.to_csv(r'./fish_frame3.csv', encoding='utf_8_sig',header=True, index=None, sep=',')


# # 4: 2-6，插入一个

# In[21]:


fo = open('./raw_pathogenic_positive_6,7_2A.csv')

ID = []
gRNA = []
Sequence = []
Position = []
Pathogenic = []
Site = []
Old_sgRNA = []

i=0
for line in fo.readlines():
    
    ##2-6
    for site in range(2,7): 
        if (i!=0):
            seq = line.strip()
            temp = seq.split(',')   
            sequence = temp[2]
            gRNA_seq = temp[1]
            pos = temp[3]
            

#             ##不为编辑点A     
#             if int(site)!=int(pos): 
            Bases = ['A','T','G','C'] 
            for Base in Bases:                    
                ID.append(temp[0]+"_Insert_"+str(site)+'_by_'+Base)
                new_gRNA = gRNA_seq[0:site-1]+Base+gRNA_seq[site-1:]
                gRNA.append(new_gRNA)
                Old_sgRNA.append(gRNA_seq)
                Site.append(site)
                Position.append(temp[3])
                Sequence.append(sequence)
                Pathogenic.append(temp[4])

    i+=1


# In[22]:


fish_frame4 = pd.DataFrame({'ID':ID,
                            'gRNA':gRNA, 
                            'Sequence':Sequence, 
                            'Position':Position,
                            'Pathogenic':Pathogenic,
                            'Site':Site,
                            'Old_sgRNA':Old_sgRNA,})   


# In[23]:


fish_frame4


# In[24]:


fish_frame4.to_csv(r'./fish_frame4.csv', encoding='utf_8_sig',header=True, index=None, sep=',')


# # 5: 7-9,替换一个不为编辑点A的或编辑框内A，为非T非本身，如果是3C，不改为G

# In[25]:


fo = open('./raw_pathogenic_positive_6,7_2A.csv')

ID = []
gRNA = []
Sequence = []
Position = []
Pathogenic = []
Site = []
Old_sgRNA = []

i=0
for line in fo.readlines():
    
    ##7-9
    for site in range(7,10): 
        
        if (i!=0):
            seq = line.strip()
            temp = seq.split(',')   
            sequence = temp[2]
            gRNA_seq = temp[1]
            pos = temp[3]
            
            ##编辑框内A不编辑
            if site < 9 and site > 3 and gRNA_seq[site-1]=='A':
                continue;
                
            ##不为编辑点A     
            if int(site)!=int(pos): 
                ##C3G
                if site==3 and gRNA_seq[site-1]=='C': 
                    Base = 'A'

                else:
                ##非T
                    Bases = ['A','G','C']    
                    k = gRNA_seq[site-1] 
                    ##非本身
                    if k != 'T':
                        Bases = r(Bases,k)
                    for Base in Bases:                                    
                        ID.append(temp[0]+"_replaced_"+str(site)+'by_'+Base)
                        new_gRNA = gRNA_seq[0:site-1]+Base+gRNA_seq[site-1+1:]
                        gRNA.append(new_gRNA)
                        Old_sgRNA.append(gRNA_seq)
                        Site.append(site)
                        Position.append(temp[3])
                        Sequence.append(sequence)
                        Pathogenic.append(temp[4])

    i+=1


# In[26]:


fish_frame5 = pd.DataFrame({'ID':ID,
                            'gRNA':gRNA, 
                            'Sequence':Sequence, 
                            'Position':Position,
                            'Pathogenic':Pathogenic,
                            'Site':Site,
                            'Old_sgRNA':Old_sgRNA,})   


# In[27]:


fish_frame5


# In[28]:


fish_frame5.to_csv(r'./fish_frame5.csv', encoding='utf_8_sig',header=True, index=None, sep=',')


# # 6: 7-9,替换一个不为编辑点A的或编辑框内A，也不为T的，为T<br>

# In[29]:


fo = open('./raw_pathogenic_positive_6,7_2A.csv')

ID = []
gRNA = []
Sequence = []
Position = []
Pathogenic = []
Site = []
Old_sgRNA = []

i=0
for line in fo.readlines():
    
    ##7-9
    for site in range(7,10): 
        if (i!=0):
            seq = line.strip()
            temp = seq.split(',')   
            sequence = temp[2]
            gRNA_seq = temp[1]
            pos = temp[3]
            
            ##编辑框内A不编辑
            if site < 9 and site > 3 and gRNA_seq[site-1]=='A':
                continue;
                
            ##不为编辑点A     
            if int(site)!=int(pos): 
                Base = 'T'   
                k = gRNA_seq[site-1] 
                ##非本身
                if k != 'T':                                  
                    ID.append(temp[0]+"_replaced_"+str(site)+'by_'+Base)
                    new_gRNA = gRNA_seq[0:site-1]+Base+gRNA_seq[site-1+1:]
                    gRNA.append(new_gRNA)
                    Old_sgRNA.append(gRNA_seq)
                    Site.append(site)
                    Position.append(temp[3])
                    Sequence.append(sequence)
                    Pathogenic.append(temp[4])

    i+=1


# In[30]:


fish_frame6 = pd.DataFrame({'ID':ID,
                            'gRNA':gRNA, 
                            'Sequence':Sequence, 
                            'Position':Position,
                            'Pathogenic':Pathogenic,
                            'Site':Site,
                            'Old_sgRNA':Old_sgRNA,})   


# In[31]:


fish_frame6


# In[32]:


fish_frame6.to_csv(r'./fish_frame6.csv', encoding='utf_8_sig',header=True, index=None, sep=',')


# # 7：7-9，删除一个不为编辑点A或编辑框内A的

# In[33]:


fo = open('./raw_pathogenic_positive_6,7_2A.csv')

ID = []
gRNA = []
Sequence = []
Position = []
Pathogenic = []
Site = []
Old_sgRNA = []

i=0
for line in fo.readlines():
    
    ##7-9
    for site in range(7,10): 
        if (i!=0):
            seq = line.strip()
            temp = seq.split(',')   
            sequence = temp[2]
            gRNA_seq = temp[1]
            pos = temp[3]
            
            ##编辑框内A不编辑
            if site < 9 and site > 3 and gRNA_seq[site-1]=='A':
                continue;
                
            ##不为编辑点A     
            if int(site)!=int(pos): 
                Base = ''                                 
                ID.append(temp[0]+"_Indel_"+str(site))
                new_gRNA = gRNA_seq[0:site-1]+Base+gRNA_seq[site-1+1:]
                gRNA.append(new_gRNA)
                Old_sgRNA.append(gRNA_seq)
                Site.append(site)
                Position.append(temp[3])
                Sequence.append(sequence)
                Pathogenic.append(temp[4])

    i+=1


# In[34]:


fish_frame7 = pd.DataFrame({'ID':ID,
                            'gRNA':gRNA, 
                            'Sequence':Sequence, 
                            'Position':Position,
                            'Pathogenic':Pathogenic,
                            'Site':Site,
                            'Old_sgRNA':Old_sgRNA,})   


# In[35]:


fish_frame7


# In[36]:


fish_frame7.to_csv(r'./fish_frame7.csv', encoding='utf_8_sig',header=True, index=None, sep=',')


# # 8：7-9，插入一个

# In[37]:


fo = open('./raw_pathogenic_positive_6,7_2A.csv')

ID = []
gRNA = []
Sequence = []
Position = []
Pathogenic = []
Site = []
Old_sgRNA = []

i=0
for line in fo.readlines():
    
    ##7-9
    for site in range(7,10): 
        if (i!=0):
            seq = line.strip()
            temp = seq.split(',')   
            sequence = temp[2]
            gRNA_seq = temp[1]
            pos = temp[3]

#             ##不为编辑点A     
#             if int(site)!=int(pos): 
            Bases = ['A','T','G','C'] 
            for Base in Bases:                    
                ID.append(temp[0]+"_Insert_"+str(site)+'_by_'+Base)
                new_gRNA = gRNA_seq[0:site-1]+Base+gRNA_seq[site-1:]
                gRNA.append(new_gRNA)
                Old_sgRNA.append(gRNA_seq)
                Site.append(site)
                Position.append(temp[3])
                Sequence.append(sequence)
                Pathogenic.append(temp[4])

    i+=1


# In[38]:


fish_frame8 = pd.DataFrame({'ID':ID,
                            'gRNA':gRNA, 
                            'Sequence':Sequence, 
                            'Position':Position,
                            'Pathogenic':Pathogenic,
                            'Site':Site,
                            'Old_sgRNA':Old_sgRNA,})   


# In[39]:


fish_frame8


# In[40]:


fish_frame8.to_csv(r'./fish_frame8.csv', encoding='utf_8_sig',header=True, index=None, sep=',')


# # 9: 2-9, 改两个

# In[50]:


fo = open('./raw_pathogenic_positive_6,7_2A.csv')

ID = []
gRNA = []
Sequence = []
Position = []
Pathogenic = []
Site = []
Old_sgRNA = []

i=0
for line in fo.readlines():  

    
    if (i!=0):
        id1 = ''
        id2 = ''
        new_gRNA_1 = ''
        new_gRNA_2 = ''
        
        seq = line.strip()
        temp = seq.split(',')   
        sequence = temp[2]
        gRNA_seq = temp[1]
        pos = temp[3]

        ##不为编辑点A                             
        pos1=get_pos()
        pos2=get_pos()
        while int(pos1) == int(pos) or (pos1<9 and pos1>3 and gRNA_seq[pos1-1]=='A'):
            pos1 = get_pos()
        while int(pos2) == int(pos1) or int(pos2) == int(pos) or (pos2<9 and pos2>3 and gRNA_seq[pos2-1]=='A'):
            pos2=get_pos()
        if int(pos1)>int(pos2):
            tmp = pos1
            pos1 = pos2
            pos2 = tmp

        operate1 = get_cate()
        operate2 = get_cate()
        #1:indel,2:insert,3:edit
        if operate1 == 1:
            while int(pos1) == int(pos) or (pos1<9 and pos1>3 and gRNA_seq[pos1-1]=='A'):
                pos1 = get_pos()
            while int(pos2) == int(pos1) or int(pos2) == int(pos) or (pos2<9 and pos2>3 and gRNA_seq[pos2-1]=='A'):
                pos2=get_pos()
            new_gRNA_1 = gRNA_seq[0:pos1-1]
            id1 = temp[0]+"_indel_"+ str(pos1)
        elif operate1 == 2:
            Base = random.choice(['A','T','G','C'])
            new_gRNA_1 = gRNA_seq[0:pos1-1]+Base+gRNA_seq[pos1-1]
            id1 = temp[0]+"_insert_"+ str(pos1)+"_by_"+Base
        elif operate1 == 3:    
            while int(pos1) == int(pos) or (pos1<8 and pos1>2 and gRNA_seq[pos1]=='A'):
                pos1 = get_pos()
            while int(pos2) == int(pos1) or int(pos2) == int(pos) or (pos2<9 and pos2>3 and gRNA_seq[pos2-1]=='A'):
                pos2=get_pos()
            k = gRNA_seq[pos1-1] 
            Bases = ['A','T','G','C']
            Bases = r(Bases,k)
            Base = random.choice(Bases)
            new_gRNA_1 = gRNA_seq[0:pos1-1]+Base
            id1 = temp[0]+"_replaced_"+ str(pos1)+"_by_"+Base
            
        #1:indel,2:insert,3:edit
        if operate2 == 1:
            while int(pos1) == int(pos) or (pos1<8 and pos1>2 and gRNA_seq[pos1]=='A'):
                pos1 = get_pos()
            while int(pos2) == int(pos1) or int(pos2) == int(pos) or (pos2<9 and pos2>3 and gRNA_seq[pos2-1]=='A'):
                pos2=get_pos()          
            new_gRNA_2 = gRNA_seq[pos1:pos2-1]+gRNA_seq[pos2:]
            id2 = temp[0]+"_indel_"+ str(pos2)
        elif operate2 == 2:           
            Base = random.choice(['A','T','G','C'])
            new_gRNA_2 = gRNA_seq[pos1:pos2-1]+Base+gRNA_seq[pos2-1:]
            id2 = temp[0]+"_insert_"+ str(pos2)+"_by_"+Base
        elif operate2 == 3:  
            while int(pos1) == int(pos) or (pos1<8 and pos1>2 and gRNA_seq[pos1]=='A'):
                pos1 = get_pos()
            while int(pos2) == int(pos1) or int(pos2) == int(pos) or (pos2<9 and pos2>3 and gRNA_seq[pos2-1]=='A'):
                pos2=get_pos()            
            k = gRNA_seq[pos2-1] 
            Bases = ['A','T','G','C']
            Bases = r(Bases,k)
            Base = random.choice(Bases)
            new_gRNA_2 = gRNA_seq[pos1:pos2-1]+Base+gRNA_seq[pos2:]               
            id2 = temp[0]+"_replaced_"+ str(pos2)+"_by_"+Base                

#         print(new_gRNA_1+new_gRNA_2)

        ID.append(id1+'___'+id2)
        new_gRNA = new_gRNA_1+new_gRNA_2
        gRNA.append(new_gRNA)
        Old_sgRNA.append(gRNA_seq)
        Site.append(str(pos1)+'__'+str(pos2))
        Position.append(temp[3])
        Sequence.append(sequence)
        Pathogenic.append(temp[4])               

    i+=1


# In[51]:


fish_frame9 = pd.DataFrame({'ID':ID,
                            'gRNA':gRNA, 
                            'Sequence':Sequence, 
                            'Position':Position,
                            'Pathogenic':Pathogenic,
                            'Site':Site,
                            'Old_sgRNA':Old_sgRNA,})   


# In[52]:


fish_frame9


# In[53]:


fish_frame9.to_csv(r'./fish_frame9.csv', encoding='utf_8_sig',header=True, index=None, sep=',')


# In[54]:


fish_frame1 = pd.read_csv("./fish_frame1.csv", index_col=None, header=0,sep = ',',encoding="utf-8")
fish_frame2 = pd.read_csv("./fish_frame2.csv", index_col=None, header=0,sep = ',',encoding="utf-8")
fish_frame3 = pd.read_csv("./fish_frame3.csv", index_col=None, header=0,sep = ',',encoding="utf-8")
fish_frame4 = pd.read_csv("./fish_frame4.csv", index_col=None, header=0,sep = ',',encoding="utf-8")
fish_frame5 = pd.read_csv("./fish_frame5.csv", index_col=None, header=0,sep = ',',encoding="utf-8")
fish_frame6 = pd.read_csv("./fish_frame6.csv", index_col=None, header=0,sep = ',',encoding="utf-8")
fish_frame7 = pd.read_csv("./fish_frame7.csv", index_col=None, header=0,sep = ',',encoding="utf-8")
fish_frame8 = pd.read_csv("./fish_frame8.csv", index_col=None, header=0,sep = ',',encoding="utf-8")
fish_frame9 = pd.read_csv("./fish_frame9.csv", index_col=None, header=0,sep = ',',encoding="utf-8")


# In[55]:


raw = pd.read_csv("./raw_pathogenic_positive_6,7_2A.csv", index_col=None, header=0,sep = ',',encoding="utf-8")


# In[56]:


raw['Old_sgRNA']=raw['gRNA']
raw['Site']='original'


# In[57]:


raw.head()


# In[58]:


all_mismatch = pd.concat([raw,fish_frame1,fish_frame2,fish_frame3,fish_frame4,fish_frame5,fish_frame6,fish_frame7,fish_frame8,fish_frame9])


# In[62]:


all_mismatch_dd = all_mismatch.drop_duplicates()


# In[63]:


all_mismatch_dd.to_csv(r'./dd_all_noA_mismatch.csv', encoding='utf_8_sig',header=True, index=None, sep=',')


# In[64]:


all_mismatch_dd


# In[61]:


raw_mismatch = pd.concat([raw,mismatch])


# In[74]:


raw_mismatch.to_csv(r'./raw_mismatch.csv', encoding='utf_8_sig',header=True, index=None, sep=',')


# In[75]:


raw_mismatch


# 1: 2-6,替换一个不为编辑点A的，为非T非本身，如果是3C，不改为G <br>
# 2: 2-6,替换一个不为编辑点A的，也不为T的，为T<br>
# 3：2-6，删除一个不为编辑点A的<br>
# 4：2-6，插入一个<br>
# 5: 7-9,替换一个不为编辑点A的，为非T非本身，如果是3C，不改为G<br>
# 6: 7-9,替换一个不为编辑点A的，也不为T的，为T<br>
# 7：7-9，删除一个不为编辑点A的<br>
# 8：7-9，插入一个<br>
# 9: 2-9, 改两个

# In[48]:


def random_index(rate):
    start = 0
    index = 0
    randnum = random.randint(1, sum(rate))
    for index, scope in enumerate(rate):
        start += scope
        if randnum <= start:
            break
    return index


# In[49]:


pos = 0
num = 0
cate = 0


# # mismatch position

# In[43]:


def get_pos():
    arr = ['2-6', '7-9']
    rate = [90, 10]
    value = arr[random_index(rate)]
    if value == '2-6':
        pos = random.choice([2,3,4,5,6])
    elif value == '7-9':
        pos = random.choice([7,8,9])
    return pos


# # mismatch number

# In[44]:


def get_num():
    arr = ['1','2','3']
    rate = [93,6,1]
    num = int(arr[random_index(rate)])
    return num


# # mismatch category

# In[46]:


# def get_cate():
#     #1:indel,2:insert,3:edit_not_T,4:edit_T
#     arr = ['1','2','3','4']
#     rate = [20,20,30,30]
#     cate = int(arr[random_index(rate)])
#     return cate
def get_cate():
    #1:indel,2:insert,3:edit_not_T,4:edit_T
    arr = ['1','2','3']
    rate = [20,20,60]
    cate = int(arr[random_index(rate)])
    return cate


# In[ ]:





# In[ ]:





# In[30]:





# 

# In[31]:





# In[32]:





# In[33]:





# 

# In[ ]:




