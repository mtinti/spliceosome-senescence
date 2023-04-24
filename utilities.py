import warnings
warnings.filterwarnings("ignore")
#define helping function
import os
from tqdm import tqdm_notebook
from sklearn.decomposition import PCA
from sklearn.manifold import MDS
from adjustText import adjust_text
from matplotlib.lines import Line2D
from Bio import SeqIO
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
import matplotlib.pyplot as plt
import seaborn as sns
import missingno as msno
import matplotlib
plt.style.use('ggplot')

#get only the gene id from
#the new TryTripDB format
def clean_id(temp_id):
    temp_id = temp_id.split(':')[0]
    if temp_id.count('.')>2:
        temp_id = '.'.join(temp_id.split('.')[0:3])
    return temp_id

#helper function to print out
#the protein removed at each threshold
def print_result(start_df_shape, shape_before, df, what):
    removed = shape_before[0]- df.shape[0]
    removed_from_beginning = start_df_shape[0]-df.shape[0]
    if removed > 0:
        print ('removed ',removed, what )  
        print ('tot ', removed_from_beginning, ' entries removed' )
        print ('---------------')
    else:
        print (what)
        print ('nothing removed')
        print ('---------------')

#remove rubbish entires from a
#maxquant output
def clean_df(df, id_by_site=True, rev_database=True, contaminant=True, unique_pep_threshold=2):  
    before,start = df.shape,df.shape
    print('starting from:', before)
    if id_by_site:
        #remove Only identified by site
        before,start = df.shape,df.shape
        col = 'Only identified by site'
        df = df[df[col] != '+'] 
        print_result(start, before, df, col)
    
    if rev_database:
        #remove hits from reverse database
        before = df.shape
        col = 'Reverse'
        df = df[df[col] != '+']
        print_result(start, before, df, col)
     

    if contaminant:
        #remove contaminants (mainly keratine and bsa)
        before = df.shape
        col = 'Potential contaminant'
        df = df[df[col] != '+']
        print_result(start, before, df, col)
    
    ##remove protein groups with less thatn 2 unique peptides
    before = df.shape
    col = 'Peptide counts (unique)'
    df['unique_int'] = [int(n.split(';')[0]) for n in df[col]]
    df = df[df['unique_int'] >= unique_pep_threshold]
    print_result(start, before, df, col)
    return df  

#format legend of hist plots 
#with lines instead of boxes
def hist_legend(ax, title = False):
    handles, labels = ax.get_legend_handles_labels()
    new_handles = [Line2D([], [], c=h.get_edgecolor()) for h in handles]
    ax.legend(handles=new_handles, labels=labels, 
    title=title,loc='center left', bbox_to_anchor=(1, 0.5))  

#rename some of maxquant output 
#columns
def mod_df(df, desc_dict=None):
    df['Gene_id'] = [clean_id(n.split(':')[0].split(';')[0])
                     for n in df['Protein IDs']]
    if desc_dict:
        
        df['Gene_desc'] = [desc_dict[n].split('=')[1].strip() 
                           for n in df['Gene_id']]
        df['Gene_id_all'] = [';'.join([ clean_id(a) for a in n.split(';')]) 
                             for n in df['Protein IDs' ]  ]  
        df['Gene_desc_all'] = ['; '.join([ desc_dict[a].split('=')[1].strip() 
                                              if a in desc_dict else 'none' for a in n.split(';')]) 
                               for n in df['Gene_id_all' ]  ]   

#create a dictionary id -> description
#from a trytripdb fasta file
def make_desc_dict(path_to_file):
    desc_dict = {}
    
    with open(path_to_file, "r") as handle:
        a=0
        for record in SeqIO.parse(handle, "fasta"):
            a+=1
            temp_id = clean_id(record.id).strip()
            temp_desc = record.description.split('|')[4].strip()
            desc_dict[temp_id]=temp_desc
    return desc_dict



                           
#make pca plot from pandas df
def make_pca(in_df, palette, ax, top=500):
    cols = in_df.columns
    pca = PCA(n_components=2)
    
    sorted_mean = in_df.mean(axis=1).sort_values()
    select = sorted_mean.tail(top)
    #print(top)
    in_df = in_df.loc[select.index.values]
    pca.fit(in_df)
    temp_df = pd.DataFrame()
    temp_df['pc_1']=pca.components_[0]
    temp_df['pc_2']=pca.components_[1]
    temp_df.index = cols
    print(pca.explained_variance_ratio_)
    temp_df['color']=palette
    #fig,ax=plt.subplots(figsize=(12,6))
    temp_df.plot(kind='scatter',x='pc_1',y='pc_2',s=30, c=temp_df['color'], ax=ax)
    #print(temp_df.index.values)
       
    texts = [ax.text(temp_df.iloc[i]['pc_1'], 
                       temp_df.iloc[i]['pc_2'],
                       cols[i])
                       for i in range(temp_df.shape[0])]
    
    
    adjust_text(texts, arrowprops=dict(arrowstyle='->', color='red'),ax=ax)
    ax.set_title('PCA', size=14)
    ax.set_xlabel('PC1_{:.3f}'.format(pca.explained_variance_ratio_[0]),size=12)
    ax.set_ylabel('PC2_{:.3f}'.format(pca.explained_variance_ratio_[1]),size=12)
    ax.yaxis.label.set_size(12)
    ax.xaxis.label.set_size(12)
    

#make mds plot from pandas df  
def make_mds(in_df, palette, ax, top=500):
    cols = in_df.columns
    pca = MDS(n_components=2,metric=True)
    
    sorted_mean = in_df.mean(axis=1).sort_values()
    select = sorted_mean.tail(top)
    #print(top)
    in_df = in_df.loc[select.index.values]
    temp_df = pd.DataFrame(pca.fit_transform(in_df.T),
                                 index=cols,columns =['pc_1','pc_2'] )
    
    temp_df['color']=palette
    
    temp_df.plot(kind='scatter',x='pc_1',y='pc_2',s=50, c=temp_df['color'], ax=ax)
    #print(temp_df.index.values)
       
    texts = [ax.text(temp_df.iloc[i]['pc_1'], 
                       temp_df.iloc[i]['pc_2'],
                       cols[i])
                       for i in range(temp_df.shape[0])]
    
    adjust_text(texts, arrowprops=dict(arrowstyle='->', color='red'),ax=ax)
    ax.set_title('MDS',size=14)
    ax.set_xlabel('DIM_1',size=12)
    ax.set_ylabel('DIM_2',size=12)
    ax.yaxis.label.set_size(12)
    ax.xaxis.label.set_size(12)

#format legend of hist plots 
#with lines instead of boxes
def hist_legend(ax):
    handles, labels = ax.get_legend_handles_labels()
    new_handles = [Line2D([], [], c=h.get_edgecolor()) for h in handles]
    ax.legend(handles=new_handles, labels=labels)

#get a random distribution of numbers 
#around the minimum value 
#of a columns (greather than zero)
#with small std
def get_random(in_col):
    mean_random = in_col[in_col>0].min()
    std_random = mean_random*0.25
    random_values = np.random.normal(mean_random, 
                                     scale=std_random, 
                                     size=in_col.shape[0])
    return  random_values

#add a small random value to each element
#of a cloumn, optionally plots the distribution
#of the random values
def impute(in_col, ax=False):
    random_values = get_random(in_col)
    if ax:
        np.log10(pd.Series(random_values)).plot(kind='hist',histtype='step', 
                          density=True,ax=ax,label=in_col.name)  
    
    fake_col = in_col.copy()
    fake_col = fake_col+random_values
    index = in_col[in_col==0].index.values 
    in_col.loc[index] = fake_col.loc[index] 
    return in_col    

#replace missing values with zeros
def replace_nan(col):
    col = col.replace('NaN', np.nan)
    col = col.fillna(0)
    return col

#normalization of dataframe
#to account for uneven loading
def norm_loading(df):
    col_sum = df.median(axis=0)
    print(col_sum)
    target = np.mean(col_sum)
    print(target)
    norm_facs = target / col_sum
    print(norm_facs)
    data_norm = df.multiply(norm_facs, axis=1)
    return data_norm

#essentially a scatter plot with the option 
#of annootating group of genes 
def make_vulcano(df, ax, x='-Log10PValue', 
                 y='Log2FC',
                 fc_col = 'Log2FC',
                 fc_limit=False,
                 
                 pval_col = 'PValue',
                 pval_limit=False,
                 
                 annot_index=pd.Series(), 
                 annot_names=pd.Series(),
                 title='Volcano',
                 legend_title='',
                 
                 label_for_selection = None,
                 label_for_all = None,
                 add_text = True,
                 do_adjust_text=True,
                 rolling_mean = False,
                 alpha_main=0.05,
                point_size_selection=1,point_size_all=1):
    

    if fc_limit and pval_limit:
        upper = df[df[fc_col]>fc_limit].copy()
        lower = df[df[fc_col]<-fc_limit].copy()
        
        upper = upper[upper[pval_col]<pval_limit]
        lower = lower[lower[pval_col]<pval_limit]
         
    elif pval_limit:
        upper = df[df[pval_col]<pval_limit].copy()
        lower = df[df[pval_col]<pval_limit].copy()
        
    elif fc_limit:
        upper = df[df[fc_col]>fc_limit].copy()
        lower = df[df[fc_col]<-fc_limit].copy()

    else: 
        print('no selection')    
            
    
    to_remove = []
    if 'upper' in locals() and upper.shape[0]>0:
        #print(upper.head())
        upper.plot(
        kind='scatter',x=x,y=y, ax=ax, 
        c='r', label='Up'.format(fc_limit), alpha=0.5, zorder=5)
        to_remove.append(upper)
        
    if 'lower' in locals() and lower.shape[0]>0:     
        lower.plot(
        kind='scatter',x=x,y=y, ax=ax, 
        c='g', label='Down'.format(fc_limit), alpha=0.5, zorder=5)       
        to_remove.append(lower)


    if len(annot_index) > 0:
        df.loc[annot_index].plot(kind='scatter', x=x, y=y, c='r', 
                                 s=point_size_selection, ax=ax, label=label_for_selection, alpha=1, zorder=10)
        to_remove.append(df.loc[annot_index])
                

 

    
    if len(to_remove)>0:
        to_remove=pd.concat(to_remove)
        idx = df.index.difference(to_remove.index)
        df.loc[idx].plot(kind='scatter', x=x, y=y, ax=ax, 
                         alpha=alpha_main,c='b', zorder=0, label=label_for_all,
                        s=point_size_all)
    else:
        df.plot(kind='scatter', x=x, y=y, ax=ax, 
                alpha=alpha_main,c='b', zorder=0,label=label_for_all,s=point_size_all)
    
    if rolling_mean:
        df = df.sort_values(x,ascending=False)
        df['rolling_mean'] = df[y].rolling(100).mean()
        #print(df.head())
        temp = df[['rolling_mean',x]]
        temp=temp.dropna()
        temp.plot(ax=ax, x=x, y='rolling_mean', label = 'rolling mean', c='r',alpha=0.3)
    
        ax.set_xlim(df[x].min()-df[x].min()*0.01,
                 df[x].max()+df[x].min()*0.01)


    if add_text:
        texts = [ax.text(df.loc[i][x], df.loc[i][y],name)
                               for i,name in zip(annot_index,annot_names)]
        #print(texts)
        if do_adjust_text:
            #print('adjusting text')
            adjust_text(texts, arrowprops=dict(arrowstyle='->',
                                               color='red'),
                        ax=ax)

    
    ax.legend(loc='upper center', bbox_to_anchor=(0.8, 0.8), title=legend_title)
    
    ax.set_title(title)
    ax.yaxis.label.set_size(12)
    ax.xaxis.label.set_size(12)
    

#helper function to visualize the correlation between experiments
def plot_correlation(df, figname='corr_prot'):
    #function to annotate the axes with
    #the pearson correlation coefficent
    def corrfunc(x, y, **kws):
        corr = np.corrcoef(x, y)
        r = corr[0][1]
        ax = plt.gca()
        ax.annotate("p = {:.2f}".format(r),
                    xy=(.1, .9), xycoords=ax.transAxes)
    
    #prepare the seaborn grid and plot
    g = sns.PairGrid(df.dropna(), palette=["red"], height=1.8, aspect=1.5)
    g.map_upper(plt.scatter, s=5)
    g.map_diag(sns.distplot, kde=False)
    g.map_lower(sns.kdeplot, cmap="Blues_d")
    g.map_lower(corrfunc)
    sns.set(font_scale=1.1)