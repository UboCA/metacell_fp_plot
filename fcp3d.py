from fcp import *

path__v7_general = '/mnt/e/Results/sequencing/Mouse_immune_system_development/202107_v7/DEV_latest/figs_mouse_DEV_v7_noSkin/General/'
path_v7_wsy = '/mnt/e/Results/sequencing/Mouse_immune_system_development/202107_v7/DEV_latest/figs_mouse_DEV_v7_noSkin/Orginal_WSY/'
path_v7_1 = '/mnt/e/Results/sequencing/Mouse_immune_system_development/202110_new_v7/figs_mouse_DEV_v7_noSkin/'
os.chdir(path__v7_general)
TimPt = pd.read_csv('./v7_TimPt_Mt.csv',header = 0,index_col = 0)
TisPt = pd.read_csv('./v7_TisPt_Mt.csv',header = 0,index_col = 0)
TimCt = pd.read_csv('./v7_TimCt_Mt.csv',header = 0,index_col = 0)
TisCt = pd.read_csv('./v7_TisCt_Mt.csv',header = 0,index_col = 0)

TimPt.index.name = None
TisPt.index.name = None

os.chdir(path_v7_wsy)
Ann = pd.read_table('./mouse_DEV_v7_noSkin_ann.txt',index_col=0) # first annotation
Ann = Ann.loc[:,'Populations']
AnnList = set(list(Ann))

Ann_1 = pd.read_table(path_v7_1+'mouse_DEV_v7_noSkin_ann.txt')  # new annotation
Ann_1 = Ann_1.loc[:,'Populations']
AnnList1 = set(list(Ann_1))


# Anndic
Anncolor = { 'BM':'#008b8b', 'B_a':'#1e90ff', 'B_b':'#0000ff', 'B_c':'#00ffff', 'B_d':'#6495ed', 
            'Basophil':'#00cdcd', 'DC':'#daa520', 'Doublets':'#ff8c00', 'ILC':'#ffd700', 
            'Mast':'#8fbc8f', 'Mf':'#9932cc', 'Mf_MHCII':'#8b008b', 'Mo':'#e9967a', 'NK':'#ffb6c1', 'NK_a':'#8b5f65',
            'Neu':'#ff1493', 'Neu_a':'#ff3030', 'Neu_b':'#ff7f50', 'Neu_c':'#d2691e', 'Neu_d':'#00ced1', 'Neu_e':'#2f4f4f', 'Neu_f':'#483d8b',
            'Plasma':'#8b3e2f', 'Progenitor':'#cd5c5c', 'Progenitor_b':'#bc8f8f',
            'T':'#add8e6', 'T_Cd24a':'#ff4040', 'T_Cd24a_Cd8a':'#c71585', 'T_Cd8a':'#66cdaa', 'T_Cd8a_Nkg7':'#db7093',
            'mDC':'#556b2f', 'pDC':'#cd3333'}
Typepool = {'B':['B_c','B_a','B_b','B_d'],
            'T_ILC':['ILC','NK','NK_a','T','T_Cd24a','T_Cd24a_Cd8a','T_Cd8a','T_Cd8a_Nkg7'],
            'T':['T','T_Cd24a','T_Cd24a_Cd8a','T_Cd8a','T_Cd8a_Nkg7'],
            'Neu':['Neu','Neu_a','Neu_b','Neu_c','Neu_d','Neu_e','Neu_f'],
            'Mof':['Mo','Mf_MHCII','Mf',],
            'Progenitor':['Progenitor','Progenitor_b'],
            'DC':['mDC','DC','pDC'],
           'other':['Basophil','BM','Plasma','Mast']}

def tocolor(self): # 从fcp返回的矩阵读取col的绝对index，进而在Anndict获得注释名称和画图颜色
    collect = []
    for i in self.columns.to_list():
        for j in list(Anndic.keys()):
            if i in Anndic[j]:
                collect.append(j)
    colorlist = []
    for i in collect:
        colorlist.append(Anncolor[i])
    return colorlist

# plot
def fch3D(mc,tis = False,tim = False,center = 10,vmin = 1,vmax= 10,cmap_name='BuGn', ncolors = 10,pad = 0.03,aspect = 50,shrink = 1):
    sn.set(font_scale=2,style="whitegrid")
    abid = mc.columns.to_list()
    width = len(abid) # mc number
    reid = np.arange(width)
    gn = len(mc.index) # gene number
    height = 1+gn # tis*13+tim*5+
    plt.figure(facecolor = 'w', figsize = (0.2*width,height))
    gs = plt.GridSpec(height, 1, hspace = 0.2, wspace = 0.0)



    ax1 = plt.subplot(gs[0,0])
    ax1.bar(reid,0.5,color=tocolor(mc),width=1)
    plt.yticks([])
    plt.xticks([])
    ax1.spines['left'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_visible(False);
    plt.margins(0,2,tight=True)  # 去除轴外沿后所占的空白，削一些y轴的颜色条

    ax2 = plt.subplot(gs[1:height,0])
    sn.heatmap(mc,vmin = vmin,vmax = vmax,ax=ax2,
                cmap=plt.cm.get_cmap(cmap_name,ncolors), # YlOrRd
                cbar = True,
                linecolor='white',linewidths=0.02,
                cbar_kws ={'label':'Footprint value',
                     'orientation':'horizontal',
                     'ticks': np.array([1,3,5,10,15,25]),
                     'format':'%.1f',
                      'pad':pad,
                        'fraction':0.1,
                           'aspect':aspect,
                          })
    #                 cbar_ax=fig.add_axes([1, 1, 50, 20]), # 暂时调不出来
    plt.yticks(rotation=0,size = 30)
    plt.xticks([])# size = 8 
# color 色条排序问题，暂时无法解决    
def fcc3D(mc,tis = False,tim = False,z_score = None,cmap='Reds',row_cluster = True,col_cluster = True,standard_score = 0,width=10,height=10):
    sn.set(font_scale=2,style="whitegrid")
    abid = mc.columns.to_list()
    width = len(abid) # mc number
    reid = np.arange(width)
    gn = len(mc.index) # gene number
    height = 1+gn # tis*13+tim*5+
    plt.figure(facecolor = 'w', figsize = (0.2*width,height))
    gs = plt.GridSpec(height, 1, hspace = 0.2, wspace = 0.0)

    ax1 = plt.subplot(gs[0,0])
    ax1.bar(reid,0.5,color=tocolor(mc),width=1)
    plt.yticks([])
    plt.xticks([])
    ax1.spines['left'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_visible(False);
    plt.margins(0,2,tight=True)  # 去除轴外沿后所占的空白，削一些y轴的颜色条
    sn.clustermap(mc,
    pivot_kws=None,
    method='average',
    metric='euclidean',
    z_score=z_score,
    cmap='Reds',
    standard_scale=0,
    figsize=(width,height),
    cbar_kws=None,
    row_cluster=True,
    col_cluster=True,
    row_linkage=None,
    col_linkage=None,
    row_colors=None,
    col_colors=None,
    mask=None,
    dendrogram_ratio=0.2,
    colors_ratio=0.03,
    cbar_pos=(0.02, 0.8, 0.05, 0.18),
    tree_kws=None,
    )
