#!/usr/bin/env python
# coding: utf-8

# In[1]:


import GEOparse
import pandas as pd
import pylab as pl
import seaborn as sns #(heat map need to do from this package for variable selection)
import matplotlib.pyplot as plt
import scipy as spy
import numpy as np
pl.rcParams['figure.figsize'] = (14, 10)
pl.rcParams['ytick.labelsize'] = 12
pl.rcParams['xtick.labelsize'] = 11
pl.rcParams['axes.labelsize'] = 23
pl.rcParams['legend.fontsize'] = 20
sns.set_style('ticks')
c1, c2, c3, c4 = sns.color_palette("Set1", 4)


# In[2]:


controls = ['GSM1337848',
            'GSM1337846',
            'GSM1337851',
            'GSM1337852',
            'GSM1337875',
            'GSM1337890',
            'GSM1337950']


# In[3]:


gse = GEOparse.get_GEO("GSE55489")


# In[4]:


gsm = GEOparse.get_GEO("GSM1338078")


# In[5]:


gsm.metadata.items()


# In[6]:


gse.gpls['GPL1261'].columns


# In[7]:


gse.gsms["GSM1337845"].columns


# In[8]:


pivoted_control_samples = gse.pivot_samples('VALUE')
print(pivoted_control_samples.shape)

#CALCULATE AVG OF EVERY ROW, HISTOGRAM


# In[9]:


pivoted_control_samples.iloc[1]


# In[10]:


pivoted_control_samples.hist("GSM1337845")
sns.despine(offset=10, trim=True)


# In[11]:


print (gse.phenotype_data[["title", "source_name_ch1", "submission_date"]]) #Works without summary, why doesn't work with summary, #will split control and treated eventually


# In[12]:


print (gse.phenotype_data) 


# FINDING UNIQUE IN SOURCE_NAME_CH1

# In[13]:


a = len(gse.phenotype_data["source_name_ch1"].unique())


# In[14]:


print(a) 


# In[15]:


print("GSM example:")
for gsm_name, gsm in gse.gsms.items():
    print("Name: ", gsm_name)
    print("Metadata:",)
for key, value in gsm.metadata.items():
    print(" - %s : %s" % (key, ", ".join(value)))
    print ("Table data:",)
    print (gsm.table.head())
    break
print()
print("GPL example:")
for gpl_name, gpl in gse.gpls.items():
    print("Name: ", gpl_name)
    print("Metadata:",)
for key, value in gpl.metadata.items():
    print(" - %s : %s" % (key, ", ".join(value)))
    print("Table data:",)
    print(gpl.table.head())
    break


# In[16]:


#FILTERING OUT UNEXPRESSED GENE SAMPLES WITH PIVOTED SAMPLES BASE
pivoted_control_samples_average = pivoted_control_samples.median(axis=1)
print("Number of probes before filtering: ", len(pivoted_control_samples_average))
expression_threshold = pivoted_control_samples_average.quantile(0.20) #60-70% OF GENES LOW EXPRESSION, ANALYZE, housekeeping genes fall away when we analyze which genes different in diff conditions 
expressed_probes = pivoted_control_samples_average[pivoted_control_samples_average >= expression_threshold].index.tolist()
print("Number of probes above threshold: ", len(expressed_probes))
samples = gse.pivot_samples("VALUE").ix[expressed_probes]
print(pivoted_control_samples_average)


# In[19]:


print(pivoted_control_samples_average)


# In[17]:


plt.title("Gene Expression Values")
plt.xlabel('Average Expression in Mice')
plt.ylabel('Number of Genes')
#LABEL X AND Y AXIS
n = plt.hist(pivoted_control_samples_average, bins = 25)


# In[20]:


print (n) #


# ----------------------------------------------------------------Standard Deviation----------------------------------------------------------------------------------------

# In[21]:


print(pivoted_control_samples_average.mean())


# In[22]:


pivoted_control_samples_std = pivoted_control_samples.std(axis=1)
print(pivoted_control_samples_std)
plt.title('Frequency of Genes vs Standard Deviation')
plt.xlabel('Standard Deviation')
plt.ylabel('Genes')
plt.hist(pivoted_control_samples_std, bins = 50) 

----------------------------------------------TABLE OF TREATMENT & NONTREATED  VS STEATOSIS OR NO STEATOSIS----------------------------------------------------------------
# In[23]:


#TREATMENT VEHICLE OR INH
treatment_table = gse.phenotype_data["source_name_ch1"]
#print(treatment_table)
treatment_tbl = gse.phenotype_data["source_name_ch1"].value_counts()
#print(treatment_tbl)
#treatment_tbl = treatment_tbl.str.split(", ")
#print(treatment_tbl[1])
vals = treatment_table.values
#print(type(vals))
#print(vals)
list1 = vals.tolist()
#out_arr = list1.split(",")
out_arr = [x.split(", ") for x in list1]
out_arr_arr = np.array(out_arr)
#print(out_arr_arr)

vtreat = print(np.count_nonzero( out_arr_arr[:, 1] == "vehicle treated"))
itreat = print(np.count_nonzero( out_arr_arr[:, 1] == "isoniazid treated"))
#numpy.count_nonzero(a, axis=None, *, keepdims=False)
"""
pseudocode: 
    - take the values in dictionary treatment_table
    - split the values based on ", " to get a part with 
    vehicle treated or isoniazid treated as separate from the 
    liver value
    - sort by samples by "vehicle treated" and "isoniazid treated"
    - create two new lists containing the samples GSM number for
    vehicle treated and inh treated
    - count all samples that were vehicle treated and all samples
    that were isoniazid treated 
        - return as tuple?
        
    OR
    
    - use treatment_tbl and find the sum of the counts of 
    the uniques with vehicle treated in them, and isoniazid 
    treated
    - return as a tuple (vehicle treated, isoniazid treated)
        """


# In[24]:


treatment_table = gse.phenotype_data["source_name_ch1"]
print(treatment_table)


# In[25]:


#STEATOSIS YES OR NO
#print(gsm.metadata.items())

# cnt_ste_yes = 0
# for gsm in gse.gsms.keys():
#     print (type(gsm))
#     print(gsm.metadata)
#     print (gsm.metadata["characteristics_ch1"][5])
# for gsm in gse.gsms[value]:
#      if (gsm.metadata["characteristics_ch1"][5] == "Yes"):
#         cnt_ste_yes += 1
# print(cnt_ste_yes)
print (gsm.metadata)

#print(type(gse.gsms))

cnt_ste_yes = 0 
gsm_yes = []
gsm_no = []
cnt_ste_no = 0
substring = "Yes"
for i in gse.gsms.keys():
    gsm = GEOparse.get_GEO(i, silent = True)
    if substring in gsm.metadata["characteristics_ch1"][5]:
        cnt_ste_yes += 1
        gsm_yes.append(i)
    else:
        cnt_ste_no += 1
        gsm_no.append(i)
print (cnt_ste_yes)
print (cnt_ste_no)
print(gsm_yes)
print (gsm_no)
        
"""
pseudocode: 
    - call object "steatosis" from gsm.metadata.items()
    dictionary
    - group

Create matrix of different strains and treatment
"""


# In[26]:


#CONDITIONAL DATA
ysubstring = "Yes"
nsubstring = "No"
inh = "Isoniazid"
veh = "Vehicle"
cnt_ste_yes_ison_yes = 0
cnt_ste_yes_veh_yes = 0
cnt_ste_no_ison_yes = 0
cnt_ste_no_veh_yes = 0
#print (gsm.metadata["characteristics_ch1"])
for i in gse.gsms.keys():
    gsm = GEOparse.get_GEO(i, silent = True)
    if  ysubstring in gsm.metadata["characteristics_ch1"][5] :
        if inh in gsm.metadata["characteristics_ch1"][4]:
            cnt_ste_yes_ison_yes += 1
        elif veh in gsm.metadata["characteristics_ch1"][4]:
            cnt_ste_yes_veh_yes += 1
    elif nsubstring in gsm.metadata["characteristics_ch1"][5] :
        if inh in gsm.metadata["characteristics_ch1"][4]:
            cnt_ste_no_ison_yes += 1
        elif veh in gsm.metadata["characteristics_ch1"][4]:
            cnt_ste_no_veh_yes += 1
print(cnt_ste_yes_ison_yes)
print(cnt_ste_yes_veh_yes)
print(cnt_ste_no_ison_yes)
print(cnt_ste_no_veh_yes)


# In[27]:


# for i in gse.gsms.keys():
#     gsm = GEOparse.get_GEO(i, silent = True)
#     print(gsm.metadata["characteristics_ch1"][4])


# In[28]:


print(gse.gsms.keys())


# In[29]:


#FILTERING OUT UNEXPRESSED GENE SAMPLES WITH METADATA BASE
print(gpl.table.head())

-------------------------------------------------------BAR GRAPH TREATMENT & STRAIN VS STEATOSIS --------------------------------------------------------------------------------------------------
# In[30]:


#FINDING UNIQUE STRAINS
strains = []

for i in gse.gsms.keys():
    gsm = GEOparse.get_GEO(i, silent = True)
    strains.append(gsm.metadata["characteristics_ch1"][3] + ", " + gsm.metadata["characteristics_ch1"][4])
    strains_set = set(strains)  
print(strains_set)


# df = pd.DataFrame(strains_set, columns=["A", "B"])


# rows = []
# for i in range(3):
#     rows.append([i, i + 1])
#     print(rows)
#     OUTPUT
#     [[0, 1], [1, 2], [2, 3]]
#     df = pd.DataFrame(rows, columns=["A", "B"])
#     print(df)

    


# In[31]:


cnt_yes = []
cnt_no = []
for i in gse.gsms.keys():
    gsm = GEOparse.get_GEO(i, silent = True)
    if ysubstring in gsm.metadata["characteristics_ch1"][5]:
        cnt_yes.append (["yes" + ", "  + gsm.metadata["characteristics_ch1"][3] + gsm.metadata["characteristics_ch1"][4]])
    elif nsubstring in gsm.metadata["characteristics_ch1"][5]:
        cnt_no.append ("no" + ", " + gsm.metadata["characteristics_ch1"][3] + gsm.metadata["characteristics_ch1"][4])
# print(cnt_yes)
# print(cnt_no)


# In[32]:


df_yes = pd.DataFrame(cnt_yes, columns=["Name"])
df_yes[['steatosis','category']] = df_yes.Name.str.split(",", expand=True)


# In[33]:


#df[['steatosis','category']]


# In[34]:


df_yes=df_yes[['steatosis','category']].groupby(['category']).agg('count')
df_yes=df_yes.rename(columns={'steatosis': "steatosis_yes"})
df_yes


# In[35]:


df_no = pd.DataFrame(cnt_no, columns=["Name"])
df_no[['steatosis','category']] = df_no.Name.str.split(",", expand=True)


# In[36]:


df_no=df_no[['steatosis','category']].groupby(['category']).agg('count')
df_no=df_no.rename(columns={'steatosis': "steatosis_no"})
df_no


# In[37]:


df_no=df_no.reset_index()
df_yes=df_yes.reset_index()


# In[38]:


print(len(df_no.category.unique()))
print(len(df_yes.category.unique()))


# In[39]:


pd.set_option('display.max_rows', 10)
strain_cnt=pd.merge(df_no, df_yes, how='outer', on='category').fillna(0)


# In[40]:


strain_cnt


# In[41]:


import matplotlib.pyplot as plt
import numpy as np


labels = strain_cnt.category
s_no = strain_cnt.steatosis_no
s_yes = strain_cnt.steatosis_yes

x = np.arange(len(labels))  # the label locations
width = 0.3  # the width of the bars

fig, ax = plt.subplots()
plt.rcParams["figure.figsize"] = (20,3)
plt.xticks(rotation=90)
rects1 = ax.bar(labels, s_no, width, label='steatosis ')
rects2 = ax.bar(labels, s_yes, width, label='strain ')
# rects1 = ax.bar(x - width/2, men_means, width, label='Men')
# rects2 = ax.bar(x + width/2, women_means, width, label='Women')

# # Add some text for labels, title and custom x-axis tick labels, etc.
# ax.set_ylabel('Scores')
# ax.set_title('Scores by group and gender')
# ax.set_xticks(x)
# ax.set_xticklabels(labels)
# ax.legend()

# ax.bar_label(rects1, padding=3)
# ax.bar_label(rects2, padding=3)

#fig.tight_layout()
#from matplotlib.pyplot import figure

plt.show()


# In[42]:


number_of_unique_values = len(strains_set)
print(number_of_unique_values)


# In[42]:


print(gse.phenotype_data["source_name_ch1"][0])


# ______________ THREE-WAY ANOVA _______________
# 
# in a loop create a temporary dataframe that has the variables from the meta data and then add expression value of one gene. then do the anova test and save the result . Then when it loops the second time, it will rewrite the temp dataframe and add the expression values of the next gene. do this for every gene and save the results for each gene.

# In[43]:


steatosis = []
trtment = []
anova_strain = []
gsm_list = []
for i in gse.gsms.keys():
    gsm = GEOparse.get_GEO(i, silent = True)
    gsm_list.append(gsm)
    steatosis.append (gsm.metadata["characteristics_ch1"][5])
    trtment.append(gsm.metadata["characteristics_ch1"][4])
    anova_strain.append(gsm.metadata["characteristics_ch1"][3])


# In[44]:


#anova_df = pd.DataFrame([gsm, steatosis, trtment, anova_strain], )
anova_dict = {"gsm": gsm_list, "steatosis": steatosis, "treatment": trtment, "strain": anova_strain}
anova_df = pd.DataFrame(anova_dict)


# In[45]:


anova_df.head()


# In[46]:


anova_df.shape


# In[47]:


pivoted_control_samples = gse.pivot_samples('VALUE')
filtered_genes = pivoted_control_samples.reset_index()[pivoted_control_samples.reset_index()["ID_REF"].isin(expressed_probes)]


# In[48]:


#for i in pivoted_control_samples.reset_index()["ID_REF"]: 
filtered_genes.iloc[1]
# pivoted_control_samples.iloc[1]


# In[49]:


print(type(expressed_probes))


# In[50]:


pivoted_control_samples.reset_index()


# In[51]:


anova_df.head()


# In[54]:


#filtered_genes.head()
filtered_genes_withoutindex= filtered_genes.set_index(['ID_REF'])


# In[55]:


filtered_genes_withoutindex


# In[85]:


filtered_genes_withoutindex.filter(like = '1416956_at', '1417232_at', '1417340_at', '1417655_a_at', '1417933_at', '1418404_at', axis = 0).mean(axis = 1)
#mean(axis=1)


# In[87]:


filtered_genes_withoutindex.columns = filtered_genes_withoutindex.columns.str.strip()
print(filtered_genes_withoutindex.mean(axis = 1))


# In[88]:


#filtered_genes_withoutindex.set_option("display.max_rows", None, "display.max_columns", None)  
#filtered_genes_withoutindex.groupby(ID_REF).agg mean(axis = 1)
filtered_genes_withoutindex.groupby(ID_REF).mean(axis = 1)


# In[89]:


prac = filtered_genes_withoutindex.iloc[135]


# In[90]:


filtered_genes_withoutindex.iloc[134]


# In[91]:


prac


# In[92]:


title=filtered_genes['ID_REF'][135]
title


# In[93]:


#for i in range(len(filtered_genes)):
 ###########   print (i)


# In[144]:


anova_df_temp=anova_df.copy(deep=True) 


# In[145]:


anova_df_temp


# In[126]:


filtered_genes_withoutindex


# In[153]:


anova_df_temp


# In[247]:


anova_df_temp2 = anova_df_temp.copy()
str(anova_df_temp2['gsm'].loc[0])
anova_df_temp2['gsm'] = [str(i)[9:-1] for i in anova_df_temp2['gsm']]
display(anova_df_temp2['gsm'])


# In[249]:


df100 = pd.DataFrame(columns = anova_df_temp2[anova_df_temp2['steatosis'].str[-3:]=='Yes']['gsm'])
filtered_genes_withoutindex2 = filtered_genes_withoutindex.copy()
list(df100.columns)
filtered_genes_withoutindex2_new = filtered_genes_withoutindex2[list(df100.columns)]


# In[258]:


filtered_genes_withoutindex2_new
genesid = ['1416956_at','1417232_at', '1417340_at', '1417655_a_at', '1417933_at', '1418404_at', '1418628_at', '1419906_at', '1421947_at', '1423101_at', '1428995_at', '1429616_at', '1431768_a_at', '1432531_at', '1434310_at', '1441280_at', '1442425_at', '1448872_at']
for i in genesid:
    print(np.mean(filtered_genes_withoutindex2_new.loc[i]))


# In[260]:


filtered_genes_withoutindex
genesid = ['1416956_at','1417232_at', '1417340_at', '1417655_a_at', '1417933_at', '1418404_at', '1418628_at', '1419906_at', '1421947_at', '1423101_at', '1428995_at', '1429616_at', '1431768_a_at', '1432531_at', '1434310_at', '1441280_at', '1442425_at', '1448872_at']
for i in genesid:
    print(np.mean(filtered_genes_withoutindex.loc[i]))


# In[155]:


anova_df_temp.loc[anova_df_temp['steatosis'] == 'steatosis: Yes']


# In[158]:


anova_df_temp2


# In[131]:


anova_df_temp.loc[anova_df_temp['steatosis'] == 'steatosis: Yes']


# In[123]:


print(gsm_list)


# In[115]:


display(anova_df_temp)


# In[96]:


anova_df['strain'].unique()


# In[97]:


for i in range(len(filtered_genes_withoutindex)): 
    title=filtered_genes['ID_REF'][i]
    anova_df_temp["values"]=np.array(filtered_genes_withoutindex.iloc[i,:])
    #print(anova_df_temp)
    anova_df_temp=anova_df_temp.drop([title], axis=1)
    #df[i(14...) + 'gene']= value
    #filtered_gene


# In[98]:


filtered_genes['ID_REF'].index


# In[99]:


filtered_genes['ID_REF'][134:137]


# In[100]:


import statsmodels.api as sm
from statsmodels.formula.api import ols
import warnings 
warnings.filterwarnings('ignore')

for i in range(len(filtered_genes_withoutindex)):
    if i in filtered_genes['ID_REF'].index: 
        title=filtered_genes['ID_REF'][i]
        anova_df_temp["values"]=np.array(filtered_genes_withoutindex.iloc[i,:])
        #print(anova_df_temp)
        model = ols('values ~ C(treatment) + C(strain) + C(steatosis) + C(treatment):C(strain) + C(steatosis):C(strain) + C(steatosis):C(treatment)' , data = anova_df_temp).fit()
        #print(sm.stats.anova_lm(model, typ=2))
        anova_df_temp=anova_df_temp.drop(["values"], axis=1)
    #df[i(14...) + 'gene']= value
    #filtered_gene
    else: 
        continue


# In[65]:


import statsmodels.api as sm
from statsmodels.formula.api import ols
values_df= pd.DataFrame(columns=['C(treatment)','C(strain)','C(steatosis)','C(treatment):C(strain)','C(steatosis):C(strain)','C(steatosis):C(treatment)', 'Residual'])
import warnings
warnings.filterwarnings('ignore')

for i in range(len(filtered_genes_withoutindex)):
    if i in filtered_genes['ID_REF'].index: 
        if i%1000==0:
            print(i)
        title=filtered_genes['ID_REF'][i]
        anova_df_temp["values"]=np.array(filtered_genes_withoutindex.iloc[i,:])
        #print(anova_df_temp)
        model = ols('values ~ C(treatment) + C(strain) + C(steatosis) + C(treatment):C(strain) + C(steatosis):C(strain) + C(steatosis):C(treatment)' , data = anova_df_temp).fit()
        #print(sm.stats.anova_lm(model, typ=2))
        anova_df_temp=anova_df_temp.drop(["values"], axis=1)
        #print(sm.stats.anova_lm(model, typ=2)[0,0]) PR(>F) 
        #print(type(sm.stats.anova_lm(model, typ = 2)[[1]][["PR(>F)"]]))
        data=sm.stats.anova_lm(model, typ = 2)[["PR(>F)"]].T.rename(index={"PR(>F)": title})
        values_df=values_df.append(data)
    #df[i(14...) + 'gene']= value
    #filtered_gene
    else: 
        continue


# In[66]:


# DROP RESIDUAL
values_df
values_df=values_df.drop(["Residual"], axis=1)
values_df


# In[67]:


new_cut_off = 0.01/len(values_df)
new_cut_off


# In[75]:


#values_df_hyp_sig_corr


# In[68]:


values_df_sig = values_df[(values_df<new_cut_off).sum(axis=1)>=5]
values_df_sig


# In[153]:


filtered_genes_withoutindex.loc[values_df_sig.index.values.tolist()].head(18)


# In[152]:


filtered_genes_withoutindex.loc[values_df_sig.index.values.tolist()]


# In[70]:


filtered_genes_withoutindex.loc[values_df_sig.index.values.tolist()].T.reset_index()


# In[71]:



#create classification model, all the genes that passed multiple hypothesis testing. remove rows from filtered_genes_withoutindex or create new matrix with just those rows
#transpose the dataframe, and add steatosis, strain, treatment
#split model into training and test (used to find the best model), and validation (10% -- only at the end when we decide the best model) before
#cross validation
#start coding logistic regression model

#what is human ortholog for mice genes


# In[ ]:


#need to find average expression of 18 genes in a dataframe with expression of only gsms with steatosis yes


# In[72]:


anova_df_temp


# In[86]:





# In[ ]:


print(gsm_list)


# In[73]:


modeldf = pd.concat([anova_df, filtered_genes_withoutindex.loc[values_df_sig.index.values.tolist()].T.reset_index()], axis=1)
modeldf=modeldf.drop(columns='gsm')
modeldf


# In[74]:


modeldf=modeldf.replace({'steatosis: Yes': 1})
modeldf=modeldf.replace({'steatosis: No': 0})
modeldf=modeldf.replace({'treatment: Isoniazid': 1})
modeldf=modeldf.replace({'treatment: Vehicle': 0})
modeldf


# In[76]:


one_hot = pd.get_dummies(modeldf['strain'])
one_hot


# In[77]:


# Import PCA 
from sklearn.decomposition import PCA
# Import StandardScaler
from sklearn.preprocessing import StandardScaler


# In[78]:


filtered_genes_withoutindex


# In[79]:


anova_df_temp


# In[82]:


from sklearn.utils import shuffle
from sklearn.model_selection import train_test_split

modeldf = shuffle(modeldf)

X=modeldf.loc[:, modeldf.columns!='steatosis']
X=X.loc[:, X.columns!='name']
X=X.loc[:, X.columns!='strain']
#X=pd.concat([X, one_hot], axis=1)
y=modeldf.loc[:, modeldf.columns=='steatosis']

X_hyp, X_val, y_hyp, y_val = train_test_split(X, y, test_size=0.10, random_state=1)
#X_train, X_test, y_train, y_test = train_test_split(X_hyp, y_hyp, test_size=0.2, random_state=1)
# X_test, X_val, y_test, y_val = train_test_split(X_test, y_test, test_size=0.5, random_state=1)


# In[302]:


X


# In[272]:


X_hyp


# In[84]:


X_hyp.shape


# In[83]:


X_val.shape


# In[84]:


y_val


# In[85]:


#X_train.shape


# In[95]:


import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report, confusion_matrix, accuracy_score, f1_score, roc_auc_score


clf=LogisticRegression(random_state = 1).fit(X_hyp,y_hyp)
y_pred=clf.predict(X_val)
print(clf.coef_)
print(clf)
# import seaborn as sns
#sns.regplot(np.array(X_train), np.array(y_train), logistic=True)


# In[96]:


#APPLY BEST HYPERPARAMETERS ON TRAINING SET


# In[97]:


print(accuracy_score(y_val, y_pred))
print(f1_score(y_val, y_pred))
print(roc_auc_score(y_val, y_pred))


# In[98]:


# LOGISTIC REGRESSION BEFORE GRID SEARCH CONFUSION MATRIX
# Let's call the confusion matrix method
from sklearn import metrics
cnf_matrix = metrics.confusion_matrix(y_val, y_pred)
# import required modules
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
get_ipython().run_line_magic('matplotlib', 'inline')

class_names=[0,1] # name  of classes
fig, ax = plt.subplots()
tick_marks = np.arange(len(class_names))
plt.xticks(tick_marks, class_names)
plt.yticks(tick_marks, class_names)
# create heatmap
sns.heatmap(pd.DataFrame(cnf_matrix), annot=True, cmap="YlGnBu" ,fmt='g')
ax.xaxis.set_label_position("top")
plt.tight_layout()
plt.title('Logistic Regression Confusion Matrix', y=1.1)
plt.ylabel('Actual')
plt.xlabel('Predicted')


# In[304]:


print(accuracy_score(y_val, y_pred))
print(f1_score(y_val, y_pred))
print(roc_auc_score(y_val, y_pred))


# In[ ]:





# In[305]:


#TRYING NO PARA WITH AVG LOGISTIC REGRESSION
acc_scr = []
f1_scr = []
rocauc_scr = []
for i in range(100):
    X_copy = X.copy()
    y_copy = y.copy()
    X_copy, y_copy = shuffle(X_copy, y_copy, random_state=1)
    X_hyp, X_val, y_hyp, y_val = train_test_split(X_copy, y, test_size=0.1, random_state=1)
    clf=LogisticRegression().fit(X_hyp,y_hyp)
    y_pred=clf.predict(X_val)
    acc_scr.append(accuracy_score(y_val, y_pred))
    f1_scr.append(f1_score(y_val, y_pred))
    rocauc_scr.append(roc_auc_score(y_val, y_pred))


# In[306]:


import statistics as stats
print(stats.mean(acc_scr))
print(stats.mean(f1_scr))
print(stats.mean(rocauc_scr))


# In[307]:


acc_scr


# In[426]:


a = np.array([1, 5, 6, 7, -8])
b = np.array([5, 6, -7, 10, 0])

a_sh, b_sh = shuffle(a, b, random_state=0)
print(a_sh)
print(b_sh)


# In[308]:


# Let's call the confusion matrix method
from sklearn import metrics
cnf_matrix = metrics.confusion_matrix(y_val, y_pred)
# import required modules
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
get_ipython().run_line_magic('matplotlib', 'inline')

class_names=[0,1] # name  of classes
fig, ax = plt.subplots()
tick_marks = np.arange(len(class_names))
plt.xticks(tick_marks, class_names)
plt.yticks(tick_marks, class_names)
# create heatmap
sns.heatmap(pd.DataFrame(cnf_matrix), annot=True, cmap="YlGnBu" ,fmt='g')
ax.xaxis.set_label_position("top")
plt.tight_layout()
plt.title('Confusion matrix', y=1.1)
plt.ylabel('Actual label')
plt.xlabel('Predicted label')


# In[ ]:


LogisticRegression(
    fit_intercept=True,
    intercept_scaling=1,
    class_weight=None,
    random_state=None,
    max_iter=100,
    multi_class='auto',
    verbose=0,
    warm_start=False,
    n_jobs=None,
    l1_ratio=None,)


# In[102]:


{'class_weight' : {'1': 0.45, '0': 0.55} }


# In[309]:


#from sklearn import svm, datasets
from sklearn.model_selection import GridSearchCV
#iris = datasets.load_iris()
parameters = {'penalty':['l1', 'l2', 'elasticnet', 'none'], 'solver':['newton-cg', 'lbfgs', 'liblinear', 'sag', 'saga'], 'dual':[True, False], 'tol' = [0.0001, 0.01, 0.00001], 'C':[0.01, 0.1, 1, 5, 9, 10],
             'fit_intercept':[True, False], 'intercept_scaling': [0.01, 0.1, 1, 10, 100], 'class_weight': [dict, 'balanced']
             'max_iter': [80, 100, 200, 500, 1000], 'multi_class': ['auto', 'ovr', 'multinomial'], 'warm_start': [True, False]}
clf_gs = GridSearchCV(LogisticRegression(), parameters, cv = 5, scoring = 'f1', random_state = 0, return_train_score=True)
clf_gs.fit(X_hyp,y_hyp)
# clf.fit(iris.data, iris.target)
# GridSearchCV(estimator=SVC(),
#              param_grid={'C': [1, 10], 'kernel': ('linear', 'rbf')})
# sorted(clf_gs.cv_results_.keys())
# ['mean_fit_time', 'mean_score_time', 'mean_test_score',...
#  'param_C', 'param_kernel', 'params',...
#  'rank_test_score', 'split0_test_score',...
# #  'split2_test_score', ...
# #  'std_fit_time', 'std_score_time', 'std_test_score']
pd.DataFrame(clf_gs.cv_results_)


# In[106]:


from sklearn.model_selection import GridSearchCV
#iris = datasets.load_iris()
parameters = {'penalty':['l1', 'l2', 'elasticnet'], 'solver':['newton-cg', 'lbfgs', 'liblinear', 'sag', 'saga'], 'dual':[True, False], 'tol' : [0.0001, 0.01, 0.00001], 'C':[0.01, 0.1, 1, 5, 9, 10],
             'intercept_scaling': [0.01, 0.1, 1, 10, 100], 'class_weight': [dict, 'balanced'],
             'max_iter': [80, 100, 200, 500, 1000]}
lg_gs = GridSearchCV(LogisticRegression(), parameters, cv = 5, scoring = 'f1', return_train_score=True)
lg_gs.fit(X_hyp,y_hyp)


# In[109]:


print(lg_gs.best_params_)
print(lg_gs.best_score_)


# In[110]:


y_pred_lg_gs = lg_gs.predict(X_val)
print(accuracy_score(y_val, y_pred_lg_gs))
print(f1_score(y_val, y_pred_lg_gs))
print(roc_auc_score(y_val, y_pred_lg_gs))


# In[111]:


# GRID SEARCH CV LOGISTIC REGRESSION
from sklearn import metrics
cnf_matrix_lg_gs = metrics.confusion_matrix(y_val, y_pred_lg_gs)
# import required modules
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
get_ipython().run_line_magic('matplotlib', 'inline')

class_names=[0,1] # name  of classes
fig, ax = plt.subplots()
tick_marks = np.arange(len(class_names))
plt.xticks(tick_marks, class_names)
plt.yticks(tick_marks, class_names)
# create heatmap
sns.heatmap(pd.DataFrame(cnf_matrix_lg_gs), annot=True, cmap="YlGnBu" ,fmt='g')
ax.xaxis.set_label_position("top")
plt.tight_layout()
plt.title('Logistic Regression Confusion Matrix After GridSearch', y=1.1)
plt.ylabel('Actual')
plt.xlabel('Predicted')


# In[112]:


def plot_search_results(grid):
    """
    Params: 
        grid: A trained GridSearchCV object.
    """
    ## Results from grid search
    results = grid.cv_results_
    means_test = results['mean_test_score']
    stds_test = results['std_test_score']
    means_train = results['mean_train_score']
    stds_train = results['std_train_score']

    ## Getting indexes of values per hyper-parameter
    masks=[]
    masks_names= list(grid.best_params_.keys())
    for p_k, p_v in grid.best_params_.items():
        masks.append(list(results['param_'+p_k].data==p_v))

    params=grid.param_grid

    ## Ploting results
    fig, ax = plt.subplots(1,len(params),sharex='none', sharey='all',figsize=(20,5))
    fig.suptitle('Score per parameter')
    fig.text(0.04, 0.5, 'MEAN SCORE', va='center', rotation='vertical')
    pram_preformace_in_best = {}
    for i, p in enumerate(masks_names):
        m = np.stack(masks[:i] + masks[i+1:])
        pram_preformace_in_best
        best_parms_mask = m.all(axis=0)
        best_index = np.where(best_parms_mask)[0]
        x = np.array(params[p])
        y_1 = np.array(means_test[best_index])
        e_1 = np.array(stds_test[best_index])
        y_2 = np.array(means_train[best_index])
        e_2 = np.array(stds_train[best_index])
        ax[i].errorbar(x, y_1, e_1, linestyle='--', marker='o', label='test')
        ax[i].errorbar(x, y_2, e_2, linestyle='-', marker='^',label='train' )
        ax[i].set_xlabel(p.upper())

    plt.legend()
    plt.show()


# In[113]:


plot_search_results(lg_gs)


# In[ ]:


#OPTUNA FOR LOGISTIC REGRESSION
import optuna


# In[116]:


lg_gs_df = pd.DataFrame(lg_gs.cv_results_)
lg_gs_df.columns
import seaborn as sns
sns.set_theme(style="ticks", palette="pastel")

# Load the example tips dataset
# tips = sns.load_dataset("tips")
plt.title('Logistic Regression')
plt.ylim([0, 0.8])
# Draw a nested boxplot to show bills by day and time
sns.boxplot(  y="mean_test_score",
            data=lg_gs_df)
sns.despine(offset=10, trim=True)


# In[155]:


#clf_gs.cv_results_.


# In[310]:


lg_gs = lg_gs.best_estimators_


# In[377]:


lg_gs.best_estimators_


# In[320]:


acc_scr = []
f1_scr = []
rocauc_scr = []
for i in range(100):
    X = shuffle(X)
    X_hyp, X_val, y_hyp, y_val = train_test_split(X, y, test_size=0.1, random_state=1)
    y_pred_gs = clf_gs.predict(X_val)
    acc_scr.append(accuracy_score(y_val, y_pred_gs))
    f1_scr.append(f1_score(y_val, y_pred_gs))
    rocauc_scr.append(roc_auc_score(y_val, y_pred_gs))
print(stats.mean(acc_scr))
print(stats.mean(f1_scr))
print(stats.mean(rocauc_scr))


# In[ ]:





# In[ ]:


#CREATE CONFUSION MATRIX


# In[317]:


# Let's call the confusion matrix method
from sklearn import metrics
cnf_matrix = metrics.confusion_matrix(y_val, y_pred_gs)
# import required modules
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
get_ipython().run_line_magic('matplotlib', 'inline')

class_names=[0,1] # name  of classes
fig, ax = plt.subplots()
tick_marks = np.arange(len(class_names))
plt.xticks(tick_marks, class_names)
plt.yticks(tick_marks, class_names)
# create heatmap
sns.heatmap(pd.DataFrame(cnf_matrix), annot=True, cmap="YlGnBu" ,fmt='g')
ax.xaxis.set_label_position("top")
plt.tight_layout()
plt.title('Confusion matrix', y=1.1)
plt.ylabel('Actual label')
plt.xlabel('Predicted label')


# In[143]:





# In[318]:


clf2=LogisticRegression(penalty = 'l2', solver= 'newton-cg', C= 5).fit(X_hyp,y_hyp)
y_pred2=clf2.predict(X_val)


# In[319]:


accuracy_score(y_val, y_pred2)


# In[151]:


RandomForestClassifier


# In[117]:


#RANDOM FOREST
from sklearn.ensemble import RandomForestClassifier
rf = RandomForestClassifier(n_estimators=100, random_state = 1) #, criterion='gini', max_depth=None, min_samples_split=2, min_samples_leaf=1, min_weight_fraction_leaf=0.0, max_features='auto', max_leaf_nodes=None, min_impurity_decrease=0.0, min_impurity_split=None, bootstrap=True, oob_score=False, n_jobs=None, random_state=None, verbose=0, warm_start=False, class_weight=None, ccp_alpha=0.0, max_samples=None)
rf.fit(X_hyp, y_hyp)


# In[118]:


# X = shuffle(X)
# X_hyp, X_val, y_hyp, y_val = train_test_split(X, y, test_size=0.1, random_state=1)
ypred_rf = rf.predict(X_val)
print(accuracy_score(y_val, ypred_rf))
print(f1_score(y_val, ypred_rf))
print(roc_auc_score(y_val, ypred_rf))
# print(stats.mean(acc_scr))
# print(stats.mean(f1_scr))
# print(stats.mean(rocauc_scr))


# In[119]:


#CONFUSION MATRIX RANDOM FOREST
from sklearn import metrics
cnf_matrix_rf = metrics.confusion_matrix(y_val, ypred_rf)
# import required modules
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
get_ipython().run_line_magic('matplotlib', 'inline')

class_names=[0,1] # name  of classes
fig, ax = plt.subplots()
tick_marks = np.arange(len(class_names))
plt.xticks(tick_marks, class_names)
plt.yticks(tick_marks, class_names)
# create heatmap
sns.heatmap(pd.DataFrame(cnf_matrix_rf), annot=True, cmap="Greens" ,fmt='g')
ax.xaxis.set_label_position("top")
plt.tight_layout()
plt.title('Random Forest Confusion matrix', y=1.1)
plt.ylabel('Actual')
plt.xlabel('Predicted')


# In[120]:


#RANDOM FOREST GS
#from sklearn import svm, datasets
from sklearn.model_selection import GridSearchCV
#iris = datasets.load_iris()
#BOOTSTRAP? 
parameters = {'n_estimators': [5, 20, 50, 70, 100, 200, 300],"max_depth": [3, 6, 9], "max_features": [3, 6, 9], "criterion": ["gini", "entropy"], 'min_samples_leaf': [1, 3, 5, 7, 9, 20]} 
# svc = svm.SVC()
rf_gs = GridSearchCV(RandomForestClassifier(), parameters, cv = 5, return_train_score= True, scoring = 'f1')
rf_gs.fit(X_hyp,y_hyp)
# clf.fit(iris.data, iris.target)
# GridSearchCV(estimator=SVC(),
#              param_grid={'C': [1, 10], 'kernel': ('linear', 'rbf')})
# sorted(clf_gs.cv_results_.keys())
# ['mean_fit_time', 'mean_score_time', 'mean_test_score',...
#  'param_C', 'param_kernel', 'params',...
#  'rank_test_score', 'split0_test_score',...
# #  'split2_test_score', ...
# #  'std_fit_time', 'std_score_time', 'std_test_score']
pd.DataFrame(rf_gs.cv_results_).head()


# In[121]:


#SCORING
print(rf_gs.best_params_)
print(rf_gs.best_score_)
print(rf_gs.best_estimator_)

y_pred_gs_rf = rf_gs.predict(X_val)
print(accuracy_score(y_val, y_pred_gs_rf))
print(f1_score(y_val, y_pred_gs_rf))
print(roc_auc_score(y_val, y_pred_gs_rf))


# In[122]:


#CONFUSION MATRIX RANDOM FOREST AFTER GRID SEARCH
from sklearn import metrics
cnf_matrix_rf_gs = metrics.confusion_matrix(y_val, y_pred_gs_rf)
# import required modules
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
get_ipython().run_line_magic('matplotlib', 'inline')

class_names=[0,1] # name  of classes
fig, ax = plt.subplots()
tick_marks = np.arange(len(class_names))
plt.xticks(tick_marks, class_names)
plt.yticks(tick_marks, class_names)
# create heatmap
sns.heatmap(pd.DataFrame(cnf_matrix_rf_gs), annot=True, cmap="Greens" ,fmt='g')
ax.xaxis.set_label_position("top")
plt.tight_layout()
plt.title('Random Forest Confusion Matrix After Grid Search', y=1.1)
plt.ylabel('Actual')
plt.xlabel('Predicted')


# In[123]:


#BOX PLOT
rf_gs_df = pd.DataFrame(rf_gs.cv_results_)
rf_gs_df.columns
import seaborn as sns
sns.set_theme(style="ticks", palette="pastel")

# Load the example tips dataset
# tips = sns.load_dataset("tips")
plt.title('Random Forest')
plt.ylim([0, 0.8])
# Draw a nested boxplot to show bills by day and time
sns.boxplot(y="mean_test_score",
            data=rf_gs_df)
sns.despine(offset=10, trim=True)


# In[ ]:





# In[421]:


#ALREADY EXECUTED
def plot_search_results(grid):
    """
    Params: 
        grid: A trained GridSearchCV object.
    """
    ## Results from grid search
    results = grid.cv_results_
    means_test = results['mean_test_score']
    stds_test = results['std_test_score']
    means_train = results['mean_train_score']
    stds_train = results['std_train_score']

    ## Getting indexes of values per hyper-parameter
    masks=[]
    masks_names= list(grid.best_params_.keys())
    for p_k, p_v in grid.best_params_.items():
        masks.append(list(results['param_'+p_k].data==p_v))

    params=grid.param_grid

    ## Ploting results
    fig, ax = plt.subplots(1,len(params),sharex='none', sharey='all',figsize=(20,5))
    fig.suptitle('Score per parameter')
    fig.text(0.04, 0.5, 'MEAN SCORE', va='center', rotation='vertical')
    pram_preformace_in_best = {}
    for i, p in enumerate(masks_names):
        m = np.stack(masks[:i] + masks[i+1:])
        pram_preformace_in_best
        best_parms_mask = m.all(axis=0)
        best_index = np.where(best_parms_mask)[0]
        x = np.array(params[p])
        y_1 = np.array(means_test[best_index])
        e_1 = np.array(stds_test[best_index])
        y_2 = np.array(means_train[best_index])
        e_2 = np.array(stds_train[best_index])
        ax[i].errorbar(x, y_1, e_1, linestyle='--', marker='o', label='test')
        ax[i].errorbar(x, y_2, e_2, linestyle='-', marker='^',label='train' )
        ax[i].set_xlabel(p.upper())

    plt.legend()
    plt.show()


# In[124]:


plot_search_results(rf_gs)


# In[422]:


plot_search_results(clf_gs)


# In[125]:


pd.DataFrame(rf_gs.cv_results_)


# In[190]:


pd.DataFrame(clf_gs_rf.cv_results_).columns


# In[126]:


from sklearn.svm import SVC
svclassifier = SVC(kernel='linear', gamma='scale')
svclassifier.fit(X_hyp, y_hyp)


# In[127]:


y_pred_svm = svclassifier.predict(X_val)
print(accuracy_score(y_val, y_pred_svm))
print(f1_score(y_val, y_pred_svm))
print(roc_auc_score(y_val, y_pred_svm))


# In[129]:


#CONFUSION MATRIX RANDOM FOREST AFTER GRID SEARCH
from sklearn import metrics
cnf_matrix_svm = metrics.confusion_matrix(y_val, y_pred_svm)
# import required modules
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
get_ipython().run_line_magic('matplotlib', 'inline')

class_names=[0,1] # name  of classes
fig, ax = plt.subplots()
tick_marks = np.arange(len(class_names))
plt.xticks(tick_marks, class_names)
plt.yticks(tick_marks, class_names)
# create heatmap
sns.heatmap(pd.DataFrame(cnf_matrix_svm), annot=True, cmap="Purples" ,fmt='g')
ax.xaxis.set_label_position("top")
plt.tight_layout()
plt.title('SVM Confusion Matrix', y=1.1)
plt.ylabel('Actual')
plt.xlabel('Predicted')


# In[128]:


svclassifier2 = SVC(kernel = 'rbf')
svclassifier2.fit(X_hyp, y_hyp)
y_pred_svm2 = svclassifier2.predict(X_val)
print(accuracy_score(y_val, y_pred_svm2))
print(f1_score(y_val, y_pred_svm2))
print(roc_auc_score(y_val, y_pred_svm2))


# In[131]:


#SVM
from sklearn.svm import SVC
# Let's try different kernels with GridSearch
param_dist = {'kernel' : ['linear', 'poly', 'rbf', 'sigmoid'], 
             'C': [0.01, 0.1, 0.5, 1, 2, 3, 4, 5, 7], 'gamma': ['scale','auto'], 
             'degree': [2, 3, 4]}
# Intiate GridSearch
svm_gs = GridSearchCV(svclassifier, param_dist, cv=5, scoring = 'f1', return_train_score = True)
svm_gs.fit(X_hyp,y_hyp)
# clf.fit(iris.data, iris.target)
# GridSearchCV(estimator=SVC(),
#              param_grid={'C': [1, 10], 'kernel': ('linear', 'rbf')})
# sorted(clf_gs.cv_results_.keys())
# ['mean_fit_time', 'mean_score_time', 'mean_test_score',...
#  'param_C', 'param_kernel', 'params',...
#  'rank_test_score', 'split0_test_score',...
# #  'split2_test_score', ...
# #  'std_fit_time', 'std_score_time', 'std_test_score']
pd.DataFrame(svm_gs.cv_results_).head()


# In[132]:


plot_search_results(svm_gs)


# In[133]:


y_pred_svm_gs = svm_gs.predict(X_val)


# In[134]:


print(svm_gs.best_params_)
print(svm_gs.best_score_)
print(svm_gs.best_estimator_)
# y_pred_svm_gs = svm_gs.predict(X_val)
print(accuracy_score(y_val, y_pred_svm_gs))
print(f1_score(y_val, y_pred_svm_gs))
print(roc_auc_score(y_val, y_pred_svm_gs))


# In[154]:


y_pred_svm_gs


# In[162]:


checkin = y_val.to_numpy
y_val[0:10]


# In[165]:


y_val[10:20]


# In[166]:


y_val[20:25]


# In[135]:


#CONFUSION MATRIX SVM AFTER GRID SEARCH
from sklearn import metrics
cnf_matrix_svm_gs = metrics.confusion_matrix(y_val, y_pred_svm_gs)
# import required modules
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
get_ipython().run_line_magic('matplotlib', 'inline')

class_names=[0,1] # name  of classes
fig, ax = plt.subplots()
tick_marks = np.arange(len(class_names))
plt.xticks(tick_marks, class_names)
plt.yticks(tick_marks, class_names)
# create heatmap
sns.heatmap(pd.DataFrame(cnf_matrix_svm_gs), annot=True, cmap="Purples" ,fmt='g')
ax.xaxis.set_label_position("top")
plt.tight_layout()
plt.title('SVM Confusion Matrix After Grid Search', y=1.1)
plt.ylabel('Actual')
plt.xlabel('Predicted')


# In[136]:


#BOX PLOT SVM AFTER GRID SEARCH
svm_gs_df = pd.DataFrame(svm_gs.cv_results_)
svm_gs_df.columns
import seaborn as sns
sns.set_theme(style="ticks", palette="pastel")

# Load the example tips dataset
# tips = sns.load_dataset("tips")
plt.title('SVM')
plt.ylim([0, 0.8])
# Draw a nested boxplot to show bills by day and time
sns.boxplot(y="mean_test_score",
            data=svm_gs_df)
sns.despine(offset=10, trim=True)


# In[443]:


y_val


# In[444]:


y_pred_svm_gs


# In[137]:


#GRADIENT BOOSTING CLASSIFIER
from sklearn.ensemble import GradientBoostingClassifier
gb_clf = GradientBoostingClassifier(n_estimators=250, learning_rate=0.05, max_depth=2.5, random_state = 1)
gb_clf.fit(X_hyp, y_hyp)
y_pred_gb = gb_clf.predict(X_val)
print(accuracy_score(y_val, y_pred_gb))
print(f1_score(y_val, y_pred_gb))
print(roc_auc_score(y_val, y_pred_gb))


# In[138]:


#CONFUSION MATRIX GRADIENT BOOSTING CLASSIFIER
from sklearn import metrics
cnf_matrix_gb_clf = metrics.confusion_matrix(y_val, y_pred_gb)
# import required modules
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
get_ipython().run_line_magic('matplotlib', 'inline')

class_names=[0,1] # name  of classes
fig, ax = plt.subplots()
tick_marks = np.arange(len(class_names))
plt.xticks(tick_marks, class_names)
plt.yticks(tick_marks, class_names)
# create heatmap
sns.heatmap(pd.DataFrame(cnf_matrix_gb_clf), annot=True, cmap="Oranges" ,fmt='g')
ax.xaxis.set_label_position("top")
plt.tight_layout()
plt.title('Gradient Boosting Classifier Confusion Matrix', y=1.1)
plt.ylabel('Actual')
plt.xlabel('Predicted')


# In[146]:


#GRADIENT BOOSTING GRID SEARCH
from sklearn.svm import SVC
# Let's try different kernels with GridSearch
param_dist = {'n_estimators' : [100, 150, 200, 250, 300, 350, 400], 'max_depth': [2, 3, 4, 5, 6, 7, 8], 
              'subsample': [0.6, 0.7, 0.8, 0.9, 1.0, 1.1], 'learning_rate': [0.001, 0.01, 0.1, 0.3, 0.5]}
# Intiate GridSearch
gb_gs = GridSearchCV(GradientBoostingClassifier(), param_dist, cv=5, scoring = 'f1', return_train_score = True)
gb_gs.fit(X_hyp,y_hyp)
# clf.fit(iris.data, iris.target)
# GridSearchCV(estimator=SVC(),
#              param_grid={'C': [1, 10], 'kernel': ('linear', 'rbf')})
# sorted(clf_gs.cv_results_.keys())
# ['mean_fit_time', 'mean_score_time', 'mean_test_score',...
#  'param_C', 'param_kernel', 'params',...
#  'rank_test_score', 'split0_test_score',...
# #  'split2_test_score', ...
# #  'std_fit_time', 'std_score_time', 'std_test_score']
pd.DataFrame(gb_gs.cv_results_).head()


# In[141]:


y_pred_gb_gs = gb_gs.predict(X_val)
print(accuracy_score(y_val, y_pred_gb_gs))
print(f1_score(y_val, y_pred_gb_gs))
print(roc_auc_score(y_val, y_pred_gb_gs))


# In[143]:


#CONFUSION MATRIX GRADIENT BOOSTING CLASSIFIER GRID SEARCH
from sklearn import metrics
cnf_matrix_gb_gs = metrics.confusion_matrix(y_val, y_pred_gb_gs)
# import required modules
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
get_ipython().run_line_magic('matplotlib', 'inline')

class_names=[0,1] # name  of classes
fig, ax = plt.subplots()
tick_marks = np.arange(len(class_names))
plt.xticks(tick_marks, class_names)
plt.yticks(tick_marks, class_names)
# create heatmap
sns.heatmap(pd.DataFrame(cnf_matrix_gb_gs), annot=True, cmap="Oranges" ,fmt='g')
ax.xaxis.set_label_position("top")
plt.tight_layout()
plt.title('Gradient Boosting Classifier After Grid Search Confusion Matrix', y=1.1)
plt.ylabel('Actual')
plt.xlabel('Predicted')


# In[147]:


plot_search_results(gb_gs)


# In[148]:


#BOX PLOT AFTER GRID SEARCH FOR GRADIENT BOOSTER CLASSIFIER
gb_gs_df = pd.DataFrame(gb_gs.cv_results_)
gb_gs_df.columns
import seaborn as sns
sns.set_theme(style="ticks", palette="pastel")

# Load the example tips dataset
# tips = sns.load_dataset("tips")
plt.title('Gradient Boosting Classifier')
plt.ylim([0, 0.8])
# Draw a nested boxplot to show bills by day and time
sns.boxplot(y="mean_test_score",
            data=gb_gs_df)
sns.despine(offset=10, trim=True)


# In[149]:


y_pred_gb_gs = gb_gs.predict(X_val)
print(accuracy_score(y_val, y_pred_gb_gs))
print(f1_score(y_val, y_pred_gb_gs))
print(roc_auc_score(y_val, y_pred_gb_gs))


# In[142]:


print(gb_gs.best_score_)
print(gb_gs.best_estimator_)
print(gb_gs.best_params_)


# In[448]:


gb_clf2 = GradientBoostingClassifier(n_estimators=200, subsample = 0.6, learning_rate=0.1, max_depth=3)
gb_clf2.fit(X_hyp, y_hyp)
y_pred_gb2 = gb_clf2.predict(X_val)
print(accuracy_score(y_val, y_pred_gb2))
print(f1_score(y_val, y_pred_gb2))
print(roc_auc_score(y_val, y_pred_gb2))


# In[456]:


#KNN CLASSIFIER
from sklearn.neighbors import KNeighborsClassifier
param_dist = {'n_neighbors': [3, 5, 7, 9]}
knn_gs = GridSearchCV(KNeighborsClassifier(), param_dist, cv = 5, scoring = 'accuracy', return_train_score = True)
knn_gs.fit(X_hyp, y_hyp)
pd.DataFrame(knn_gs.cv_results_)


# In[457]:


y_pred_knn_gs = knn_gs.predict(X_val)
print(accuracy_score(y_val, y_pred_knn_gs))
print(f1_score(y_val, y_pred_knn_gs))
print(roc_auc_score(y_val, y_pred_knn_gs))


# In[463]:


from sklearn.neighbors import KNeighborsClassifier
knn = KNeighborsClassifier()
knn.fit(X_hyp, y_hyp)
y_pred_knn = knn.predict(X_val)
print(accuracy_score(y_val, y_pred_knn))
print(f1_score(y_val, y_pred_knn))
print(roc_auc_score(y_val, y_pred_knn))


# In[196]:


anova_df_temp
values = 1415670_at
model = ols('values ~ C(treatment) + C(strain) + C(treatment):C(strain)' , data = anova_df_temp).fit()
    #anova_df_temp=anova_df_temp.drop([title], axis=1)
print(sm.stats.anova_lm(model, typ=2))


# In[ ]:





# In[197]:


import statsmodels.api as sm
from statsmodels.formula.api import ols

for i in range(2): 
    title=filtered_genes['ID_REF'][i]
    anova_df_temp[values]=np.array(filtered_genes_withoutindex.iloc[i,:])
    #print(anova_df_temp)
    model = ols('values ~ C(steatosis) + C(trtment) + C(strain) + C(steatosis):C(trtment) + C(steatosis):C(strain) + C(strain):C(trtment)' , data = anova_df_temp).fit()
    anova_df_temp=anova_df_temp.drop([values], axis=1)
    #print(sm.stats.anova_lm(model, typ=2))


# In[101]:


filtered_genes_withoutindex


# In[279]:


anova_df


# In[ ]:




