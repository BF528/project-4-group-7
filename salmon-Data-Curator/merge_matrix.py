#!/usr/bin/env python
# coding: utf-8

# In[49]:


pwd 


# In[50]:


cd /projectnb/bf528/users/group7/project4/Data/fastq_counts/alevin_matrix


# In[54]:


pwd


# In[55]:


from scipy.io import mmread
import pandas as pd


# In[56]:


#reading in matrix generated by alevin
alevin_df1 = mmread("quants_mat1.mtx").toarray()
alevin_df2 = mmread("quants_mat2.mtx").toarray()
alevin_df3 = mmread("quants_mat3.mtx").toarray()

#columns for dataframe 
with open("quants_mat1_cols.txt") as cols1:
    r_columns1 = cols1.readlines()
    columns_1= [elem.strip("\n") for elem in r_columns1]
#print(columns_1)

with open("quants_mat2_cols.txt") as cols2:
    r_columns2 = cols2.readlines()
    columns_2= [elem.strip("\n") for elem in r_columns2]

with open("quants_mat3_cols.txt") as cols3:
    r_columns3 = cols3.readlines()
    columns_3= [elem.strip("\n") for elem in r_columns3]

#rows for dataframe 
with open("quants_mat1_rows.txt") as rows1:
    r_rows1 = rows1.readlines()
    rows_1= [elem.strip("\n") for elem in r_rows1]


with open("quants_mat2_rows.txt") as rows2:
    r_rows2 = rows2.readlines()
    rows_2= [elem.strip("\n") for elem in r_rows2]

with open("quants_mat3_rows.txt") as rows3:
    r_rows3 = rows3.readlines()
    rows_3= [elem.strip("\n") for elem in r_rows3]


# In[57]:


df1 = pd.DataFrame(alevin_df1, columns=columns_1, index=rows_1)


# In[58]:


df1


# In[59]:


df2 = pd.DataFrame(alevin_df2, columns=columns_2, index=rows_2)


# In[60]:


df2


# In[61]:


df3 = pd.DataFrame(alevin_df3, columns=columns_3, index=rows_3)


# In[62]:


df3


# In[63]:


#transpose data frame to have genes as index and barcodes as column headers 
df1_transposed = df1.T
df2_transposed = df2.T
df3_transposed = df3.T


# In[64]:


df1_transposed


# In[65]:


df2_transposed


# In[66]:


df3_transposed


# In[67]:


#concatenate all the data frames to get a umi count matrix for all libraries 
umi_count = pd.concat([df1_transposed, df2_transposed, df3_transposed], axis=1)


# In[68]:


umi_count


# In[69]:


umi_count.size


# In[70]:


umi_count.to_csv('UMI_counts_matrix.csv')


# In[ ]:





# In[ ]:





# In[ ]:




