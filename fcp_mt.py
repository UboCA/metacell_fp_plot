# load
import os
import seaborn as sn
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import datatable as dt


os.chdir('/mnt/e/Results/sequencing/Mouse_immune_system_development/202107_v7/DEV_latest/figs_mouse_DEV_v7_noSkin/Orginal_WSY/')

Ann = pd.read_table('./mouse_DEV_v7_noSkin_ann.txt')
Ann = Ann.loc[:,'Populations']
AnnList = set(list(Ann))

os.chdir('/mnt/e/Results/sequencing/Mouse_immune_system_development/202107_v7/DEV_latest/figs_mouse_DEV_v7_noSkin/General/')
TimPt = pd.read_csv('./v7_TimPt_Mt.csv',header = 0)
TisPt = pd.read_csv('./v7_TisPt_Mt.csv',header = 0)


# 输入Ann [pd.DataFrame] 得到所有的mcname映射mcindex 的字典，方便后面画图进行按组分颜色
def getmckey(self): 
    Anndic = {}
    keys = list(set(list(self)))
    keynum = len(keys)
    for i in range(keynum):
        Anndic[keys[i]] = list(map(lambda x:x+1,self[self == keys[i]].index.tolist()))
    return Anndic   
