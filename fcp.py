import os
import numpy as np
import pandas as pd
import datatable as dt
import matplotlib.pyplot as plt
import seaborn as sn
os.chdir('/mnt/e/Results/sequencing/Mouse_immune_system_development/202107_v7/DEV_latest/figs_mouse_DEV_v7_noSkin/Orginal_WSY/')

v7 = dt.fread('./mouse_DEV_v7_noSkin_mean.csv') # datatable读取大文件更快
del v7[:,'expressed_genes_avg_umi']
v7 = v7.to_pandas() # 读取成功后转换为pandas格式，因为用惯了pandas的数据处理方法
v7.index = v7.iloc[:,0]
v7.drop(columns='C0',inplace = True) #注意这里datatable读取的index需要在pandas里重设
Ann = pd.read_table('./mouse_DEV_v7_noSkin_ann.txt')
Ann = Ann.loc[:,'Populations']

#-------------------set Ann specific color-------------------------------------
Anncolor = { 'BM':'#008b8b',
 'B_a':'#1e90ff',
 'B_b':'#0000ff',
 'B_c':'#00ffff',
 'B_d':'#6495ed',
 'Basophil':'#00cdcd',
 'DC':'#daa520',
 'Doublets':'#ff8c00',
 'ILC':'#ffd700',
 'Mast':'#8fbc8f',
 'Mf':'#9932cc',
 'Mf_MHCII':'#8b008b',
 'Mo':'#e9967a',
 'NK':'#ffb6c1',
 'NK_a':'#8b5f65',
 'Neu':'#ff1493',
 'Neu_a':'#ff3030',
 'Neu_b':'#ff7f50',
 'Neu_c':'#d2691e',
 'Neu_d':'#00ced1',
 'Neu_e':'#2f4f4f',
 'Neu_f':'#483d8b',
 'Plasma':'#8b3e2f',
 'Progenitor':'#cd5c5c',
 'Progenitor_b':'#bc8f8f',
 'T':'#add8e6',
 'T_Cd24a':'#ff4040',
 'T_Cd24a_Cd8a':'#c71585',
 'T_Cd8a':'#66cdaa',
 'T_Cd8a_Nkg7':'#db7093',
 'mDC':'#556b2f',
 'pDC':'#cd3333'}
# keys = list(set(list(Ann)))
# keynum = len(keys)
# cmap = plt.get_cmap("tab20")#tab20c,gnuplot
# colors = [cmap(i) for i in np.linspace(0, 1, keynum)]
# for i in range(keynum):
#     Anncolor[keys[i]] = colors[i]
Typepool = {'B':['B_c','B_a','B_b','B_d'],
            'T_ILC':['ILC','NK','NK_a','T','T_Cd24a','T_Cd24a_Cd8a','T_Cd8a','T_Cd8a_Nkg7'],
            'Neu':['Neu','Neu_a','Neu_b','Neu_c','Neu_d','Neu_e','Neu_f'],
            'Mof':['Mo','Mf_MHCII','Mf',],
            'Progenitor':['Progenitor','Progenitor_b'],
            'DC':['mDC','DC','pDC'],
            'T':['T','T_Cd24a','T_Cd24a_Cd8a','T_Cd8a','T_Cd8a_Nkg7'],
           'other':['Basophil','BM','Plasma','Mast']}
#----------------------------------------------------------------------------    
def getmckey(self): # 根据Ann.text对不同元素分组建立字典索引，方便后面画图进行按组分颜色
    Anndic = {}
    keys = list(set(list(Ann))) # cell type list
    keynum = len(keys) # cell type num
    for i in range(keynum):
        Anndic[keys[i]] = list(map(lambda x:x+1,self[self == keys[i]].index.tolist())) # map apply to mc (specific cell type) index
    return Anndic    
Anndic = getmckey(Ann) # mc name 和mc index 的字典映射

def GeneTranslate(self,cap = True):
    self = self.replace(" ","") # delete space
    if ',' in self:
        self = self.split(',')
#     self = list(set(self))   # 去重，但set会打乱顺序
    if cap == True:
        self = [i.capitalize() for i in self]
#     for i in range(0,len(self)):
#         self[i] = self[i].capitalize() # 转换为小鼠基因格式,!!!注意某些小鼠基因不光只是首字母大写，比如Tcrg-C1，H2-DM，经过此函数转换后会出现不匹配
    return self

def FindList(self,mcV):
    f = []
    mismatch = 0
    for i in range(0,len(self)):
        if self[i] in list(mcV.index):
            f.append(self[i])
        else:
            print('*{}* Not Found'.format(self[i]))
            mismatch+=1
    print('Total',mismatch,'genes mismatched')
    print('{} gene(s) found'.format(len(f)))
    return f

def CreateMatrix(self,mcV):
    GeneMatrix = pd.DataFrame(np.zeros(shape = (len(self),len(mcV.columns))))
    GeneMatrix.index = [i for i in self]
    GeneMatrix.columns = [i for i in range(1,len(mcV.columns)+1)]
    for i in range(len(self)):
        GeneMatrix.iloc[i,:] = mcV.loc[self[i],:]
    return GeneMatrix

def drawmarker(self,label =None,mcname = None):  # label adds vlines reference vertical line
    GL = self.index # genelist
    mcNum = len(self.columns) # metacell numbers
    gene_num = len(GL) # control plot height and legend color point size and its text fontsize

    height = gene_num * 1.2
    width = mcNum*0.3
    if width>50:
        width = 50
    default = 'darkcyan'

    plt.figure(facecolor = 'w', figsize = (width,height))
    gs = plt.GridSpec(gene_num, 3, hspace = 0.2, wspace = 0.0,width_ratios = [0.5,16,2])
    if mcname==None:
        mcname = list(Anndic.keys())# if mcname == None: # if not input mcname, mcname should be all mc.

#-----------------------------plot gene in y axis-----------------------------------------------
    for i in range(gene_num):
        ax0 = plt.subplot(gs[i,0])
        ax0.set_ylabel(GL[i], fontsize = 50, rotation = 'horizontal', va = 'center', ha = 'right') # family:Arial
        plt.yticks([])
        plt.xticks([])
        ax0.spines['left'].set_visible(False)
        ax0.spines['top'].set_visible(False)
        ax0.spines['bottom'].set_visible(False)
#-----------------------------plot bar groupby metacell---------------------------------------------        
    for i in range(gene_num):
        ax1 = plt.subplot(gs[i,1])
        
        plt.xlim(-1,mcNum+1)         
        ymin = min(self.iloc[i,:])

        ymax = max(self.iloc[i,:])
        ax1.set_ylim(ymin,ymax)  
#         ax1.set_ylim(0,ymax)  
#         ax1.hlines(1,0,mcNum,color = 'grey',linestyle='dashed')
        for t in range(len(mcname)):
            absolute_index = [s for s in filter(lambda x:x in Anndic[mcname[t]],self)]
            relative_index = getlabel(absolute_index,self)
            ax1.bar(relative_index,self.loc[GL[i],absolute_index],linewidth = 0,width =0.8,color =Anncolor[mcname[t]])
            ax1.yaxis.label.set_fontsize(45)
        if label != None:
            ax1.vlines(label,ymin,ymax,color = 'r',linestyle = 'dashed')
            
        if mcNum>400:
            size = 4
        elif mcNum>200:
            size = 8
        elif mcNum>100:
            size = 10
        elif mcNum<=100:
            size = 15
        if i !=gene_num-1: #非尾列不用画ticks
            plt.xticks([])
        else:
            plt.xticks(np.arange(mcNum),self.columns,rotation='vertical',fontsize = size) # label mc index
#-----------------------------plot colorbar---------------------------------------------  
    ax2 = plt.subplot(gs[:,2])
    La=0
    for mc in mcname:
        ax2.scatter(1,La+1,color=Anncolor[mc],s = 15*gene_num,label = mc,edgecolors = 'black')
        ax2.text(2,La+0.95,s= mc,ha='center',fontsize=gene_num*2)
        La+=1
        plt.xticks([])
        plt.yticks([])
        plt.xlim(0.8,2)
        plt.ylim(-1,len(mcname)+2)
        ax2.spines['right'].set_visible(False)
        ax2.spines['top'].set_visible(False)
        ax2.spines['bottom'].set_visible(False)
        

# get absolute metacell number for inputing label
def getlabel(self,mc_mt):
    final = [] # absoluate num list
    ref = [i for i in mc_mt.columns]
    for i in self:
        final.append(ref.index(i))
    return final

def Findmclist(Ann = Ann,keyword = []):
    Ann_unique = set(list(Ann))
    mc_cols = []
    for i in keyword:
        mclist = Ann[Ann == i].index.tolist()
        mc_cols.extend(mclist)
    print('Total {} metacells has been found in all'.format(len(mc_cols)))
    return mc_cols

def fcp(mcV,glist,mcname = None,output = None,label = None,wcsv = False,plot = False):
    '''
    fcp aims to find matrix composed of defined gene list and defined metacell list in one version of metacell_mean_csv
    mcV -> metacell version, e.g. v7
    mcname -> metacell name keyword list
    label could be int list-> create vines in specific metacell, e.g. metacell 175 of 1-751 labeled 175
    wcsv -> whether to save csv.    output name -> csv name'''
    
#     glist = GeneTranslate(glist) 
    glist = FindList(glist,mcV) 
    if glist == None:
        print('Gene matched is none, recheck, please!')
        return None
    mc_mt = CreateMatrix(glist,mcV) # 基于输入基因list 查找创建相应矩阵
    if mcname!=None:
        mclist = Findmclist(Ann = Ann,keyword = mcname)
        mc_mt = mc_mt.iloc[:,mclist]
    if wcsv == True:
        mc_mt.to_csv('./{}.csv'.format(output)) 
    if plot == True:
        drawmarker(mc_mt,label = getlabel(label,mc_mt),mcname = mcname)
    return mc_mt

