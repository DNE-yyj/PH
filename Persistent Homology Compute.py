# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 12:58:28 2023

@author: 12739
"""

import tkinter as tk
from tkinter.filedialog import askdirectory
import os
import numpy as np
from gtda.homology import VietorisRipsPersistence
from pymatgen.core import Structure
import pandas as pd
from pymatgen.transformations.advanced_transformations import CubicSupercellTransformation
from gtda.diagrams import PairwiseDistance
import time



def computer_one_life(diagram):
    Bar0Death = []; Bar1Birth = []; Bar1Death = []; Bar2Birth = []; Bar2Death = []
    for line in diagram[0]:
        birth, death,dim = line[0],line[1],line[2]
        # Birth
        if dim == 1:
            Bar1Birth.append(birth)
        elif dim == 2:
            Bar2Birth.append(birth)
        # Death
        if death == float('inf'): continue
        if dim == 0:
            Bar0Death.append(death)
        elif dim == 1:
            Bar1Death.append(death)
        elif dim == 2:
            Bar2Death.append(death)
    Bar0Death = np.asarray(Bar0Death); Bar1Birth = np.asarray(Bar1Birth); Bar1Death = np.asarray(Bar1Death); Bar2Birth = np.asarray(Bar2Birth); Bar2Death = np.asarray(Bar2Death);
    Feature_2 = []
    # Betti0
    if len(Bar0Death) > 0:
        Feature_2.append(len(Bar0Death))#0D生死对数0
        Feature_2.append(np.sum(Bar0Death, axis=0))#总生存时间1
        Feature_2.append(np.max(Bar0Death, axis=0)-np.min(Bar0Death, axis=0))#死亡时间之间的最大差距2
        Feature_2.append(np.mean(Bar0Death, axis=0))#平均生存时间3
        Feature_2.append(np.std(Bar0Death, axis=0))#死亡时间标准差4
        Feature_2.append(np.max(Bar0Death, axis=0))#0D最大寿命以及最大死亡时间5
        Feature_2.append(np.min(Bar0Death, axis=0))#最小死亡时间6
        Feature_2.append(np.sum(Bar0Death, axis=0)/len(Bar0Death)) #死亡时间平均间隔7
        i=0
        for death in Bar0Death:
            if "{:.2f}".format(death) != "{:.2f}".format(Bar0Death[0]):
                i=i+1
            else:    
                continue
        Feature_2.append(i)#生死对类型数量8        
    else:
        Feature_2.extend([0.]*5)
    # Betti1
    if len(Bar1Death) > 0:
        
        Feature_2.append(len(Bar1Death))#1D生死对数9
        Feature_2.append(np.sum(Bar1Death - Bar1Birth, axis=0))#总生存时间10
        Feature_2.append(np.mean(Bar1Death - Bar1Birth, axis=0))#平均生存时间11
        Feature_2.append(np.std(Bar1Death - Bar1Birth, axis=0))#12
        Feature_2.append(np.max(Bar1Death - Bar1Birth, axis=0))#最长生存时间13
        Feature_2.append(np.min(Bar1Death - Bar1Birth, axis=0))
        Feature_2.append(np.mean(Bar1Birth, axis=0))#平均出生时间15
        Feature_2.append(np.std(Bar1Birth, axis=0))#16
        Feature_2.append(np.max(Bar1Birth, axis=0))#最晚出生时间17
        Feature_2.append(np.max(Bar1Birth, axis=0)-np.min(Bar1Birth, axis=0))#出生时间之间的最大差距18
        Feature_2.append(np.mean(Bar1Death, axis=0))#平均死亡时间19
        Feature_2.append(np.std(Bar1Death, axis=0))#出生时间标准差20
        Feature_2.append(np.max(Bar1Death, axis=0))#最晚死亡时间21
        Feature_2.append(np.min(Bar1Death, axis=0))#22
        Feature_2.append(np.sum(Bar1Death, axis=0)/len(Bar1Death)) #死亡时间平均间隔23
        Feature_2.append(np.sum(Bar1Birth, axis=0)/len(Bar1Birth)) #出生时间平均间隔24

       
    else:
        Feature_2.extend([0.]*15)
    # Betti2
    if len(Bar2Death) > 0:
        Feature_2.append(len(Bar2Death))#2D生死对数25
        Feature_2.append(np.max(Bar2Birth, axis=0)-np.min(Bar2Birth, axis=0))#出生时间之间的最大差距26
        Feature_2.append(np.sum(Bar2Death - Bar2Birth, axis=0))#总生存时间27
        Feature_2.append(np.mean(Bar2Death - Bar2Birth, axis=0))#28
        Feature_2.append(np.std(Bar2Death - Bar2Birth, axis=0))#29
        Feature_2.append(np.max(Bar2Death - Bar2Birth, axis=0))#30
        Feature_2.append(np.min(Bar2Death - Bar2Birth, axis=0))#31
        Feature_2.append(np.mean(Bar2Birth, axis=0))#32
        Feature_2.append(np.std(Bar2Birth, axis=0))#33
        Feature_2.append(np.max(Bar2Birth, axis=0))#34
        Feature_2.append(np.min(Bar2Birth, axis=0))#35(1D无)
        Feature_2.append(np.mean(Bar2Death, axis=0))#36
        Feature_2.append(np.std(Bar2Death, axis=0))#37
        Feature_2.append(np.max(Bar2Death, axis=0))#38
        Feature_2.append(np.min(Bar2Death, axis=0))#39
        Feature_2.append(np.sum(Bar2Death, axis=0)/len(Bar2Death)) #死亡时间平均间隔40
        Feature_2.append(np.sum(Bar2Birth, axis=0)/len(Bar2Birth)) #出生时间平均间隔41
        
    else:
        Feature_2.extend([0.]*15)
    return Feature_2




def supercell(filename):
    supercell_structure = Structure.from_file(filename)    
    supercell_structure = CubicSupercellTransformation(min_length=24).apply_transformation(supercell_structure)
    xyz=supercell_structure.cart_coords   
    VR = VietorisRipsPersistence(homology_dimensions=[0,1,2],collapse_edges=True)  # Parameter explained in the text
    diagram = VR.fit_transform([xyz])
    feature=computer_one_life(diagram)
  
    return feature


def barcode_feature():
    filename = os.listdir(path.get())
    barcodes=[]
    nos=[]
    j=0
    text.insert('insert','\n')
    text.insert('insert','Computing......')
    for i in filename:       
        filename = os.path.join(path.get(), i).replace("//", "\\")     
        diag_fea=supercell(filename)
        j=j+1
        string_print='\n'+' NO '+str(j)+'个 '+i+" barcode compute finished"
        text.insert('insert',string_print)
        barcodes.append(diag_fea)
        nos.append(i)
    
    
    columns_name=['0D_the number of birth-death pair',
                  '0D_total death time',
                  '0D_the maximum gap between death times',
                  '0D_mean death time',
                  '0D_std death time',
                  '0D_the maximum death time',
                  '0D_the minimum deah time',
                  '0D_the average death times',
                  '0D_the number of birth-death pair types',
                  '1D_the number of birth-death pair',
                  '1D_total survival time',
                  '1D_mean survival time',
                  '1D_std survival time',
                  '1D_the maximum survival time',
                  '1D_the minimum survival time',
                  '1D_mean birth time',
                  '1D_std birth time',
                  '1D_the maximum birth time',
                  '1D_the maximum gap between birth times',
                  '1D_mean death time',
                  '1D_std death time',
                  '1D_the maximum death time',
                  '1D_the minimum death time',
                  '1D_the average death times',
                  '1D_the average birth times',
                  '2D_the number of birth-death pair',
                  '2D_the maximum gap between birth times',
                  '2D_total survival time',
                  '2D_mean survival time',
                  '2D_std survival time',
                  '2D_the maximum survival time',
                  '2D_the minimum survival time',
                  '2D_mean birth time',
                  '2D_std birth time',
                  '2D_the maximum birth time',
                  '2D_the minimum birth time',
                  '2D_mean death time',
                  '2D_std death time',
                  '2D_the maximum death time',
                  '2D_the minimum death time',
                  '2D_the average death times',
                  '2D_the average birth times',
                  
                  ]
    feature_all = np.asarray(barcodes)
    ID_feature_df = pd.DataFrame(data=feature_all,
                                 columns=columns_name,
                                 index=nos)
    ID_feature_df.index.name = 'ID'
    outpath=path1.get().replace("//", "\\")+'barcodes_fature.csv'
    ID_feature_df.to_csv(outpath)
    tk.messagebox.showinfo('提示','Computer Finised')


def computer_distance():
    filename = os.listdir(path.get())
    points=[]
    nos=[]
    text.insert('insert','\n')
    text.insert('insert','Computing......')
    for i in filename:   
        filename = os.path.join(path.get(), i).replace("//", "\\")     
        supercell_structure = Structure.from_file(filename)
        #将得到的cif文件结构转换为立方超级细胞晶格
        supercell_structure = CubicSupercellTransformation(min_length=24).apply_transformation(supercell_structure)
        xyz=supercell_structure.frac_coords
        points.append(xyz)
        nos.append(i)
                
    homology_dimensions = [0, 1, 2]

    persistence = VietorisRipsPersistence(
            metric="euclidean",
            homology_dimensions=homology_dimensions,
            n_jobs=6,
            collapse_edges=True,
        )
    diagrams = persistence.fit_transform(points)
    PD=PairwiseDistance()
    dis_landscape=PD.fit_transform(diagrams)


    ID_feature_df = pd.DataFrame(data=dis_landscape,
                                 columns=nos,
                                 index=nos)
    ID_feature_df.index.name = 'no'
    outpath=path1.get().replace("//", "\\")+'distance.csv'
    ID_feature_df.to_csv(outpath)
    tk.messagebox.showinfo('提示','Computer Finised')



########选取输入文件路径###########
def selectPath():
    path_ = askdirectory() #使用askdirectory()方法返回文件夹的路径
    if path_ == "":
        path.get() #当打开文件路径选择框后点击"取消" 输入框会清空路径，所以使用get()方法再获取一次路径
    else:
        path_ = path_.replace("/", "\\")  # 实际在代码中执行的路径为“\“ 所以替换一下
        path.set(path_)

def openPath():
    dir = os.path.dirname(path.get()+"\\")
    os.system('start ' + dir)

###########获取输出文件路径############
def selectPath1():
    path_ = askdirectory() #使用askdirectory()方法返回文件夹的路径
    if path_ == "":
        path1.get() #当打开文件路径选择框后点击"取消" 输入框会清空路径，所以使用get()方法再获取一次路径
    else:
        path_ = path_.replace("/", "\\")  # 实际在代码中执行的路径为“\“ 所以替换一下
        path1.set(path_)


def openPath1():
    dir = os.path.dirname(path1.get()+"\\")
    os.system('start ' + dir)


def ternimal_print(msg, info_type):
    text.insert('end', "\n%s [%s] %s" % (time.strftime('%Y-%m-%d %H:%M:%S'), info_type.upper(), msg))
    text.update()

##########可视化界面###############
root = tk.Tk()
root.title("Persistent Homology Compute APP ")
fr1=tk.Frame(root) # 不设置边线宽，无法显示
fr1.pack(padx=10,pady=5)
fr2=tk.Frame(root) # 不设置边线宽，无法显示
fr2.pack(padx=10,pady=5)
fr3=tk.Frame(root) # 不设置边线宽，无法显示
fr3.pack(padx=(50,0),pady=5)




path = tk.StringVar()
label1=tk.Label(fr1, text="输入文件路径:",font=('微软雅黑', 11, 'bold')).grid(row=0, column=0)
inputpath=tk.Entry(fr1, textvariable=path,state="readonly").grid(row=0, column=1,ipadx=200)
button1=tk.Button(fr1, text="路径选择", command=selectPath).grid(row=0, column=2)
button2=tk.Button(fr1, text="打开文件位置", command=openPath).grid(row=0, column=3)


path1 = tk.StringVar()
label2=tk.Label(fr2, text="输出文件路径:",font=('微软雅黑', 11, 'bold')).grid(row=1, column=0)
outputpath=tk.Entry(fr2, textvariable=path1,state="readonly").grid(row=1, column=1,ipadx=200)
button3=tk.Button(fr2, text="路径选择", command=selectPath1).grid(row=1, column=2)
button4=tk.Button(fr2, text="打开文件位置", command=openPath1).grid(row=1, column=3)

button5 = tk.Button(fr3, text='computer barcode feature',command=barcode_feature)
button5.pack(side='left',ipadx=30,padx=5)
button6 = tk.Button(fr3, text='computer distance',command=computer_distance)
button6.pack(side='left',ipadx=30,padx=5)




fm_t = tk.Frame(root)     # text
fm_t.pack(fill='both')
l2 = tk.Label(fm_t, text='实时终端控制台', font=('微软雅黑', 15, 'bold'), width=200, justify='left', anchor='w')   # justify控制对其方向，anchor控制位置 共同使文本靠左
l2.pack()

s2 = tk.Scrollbar(fm_t)      # 设置垂直滚动条
b2 = tk.Scrollbar(fm_t, orient='horizontal')    # 水平滚动条
s2.pack(side='right', fill='y')     # 靠右，充满Y轴
b2.pack(side='bottom', fill='x')    # 靠下，充满x轴

#实时输出进程
text = tk.Text(fm_t, font=('Consolas', 9), undo=True, autoseparators=False,
       wrap='none', xscrollcommand=b2.set, yscrollcommand=s2.set)  # , state=DISABLED, wrap='none'表示不自动换行
text.pack(fill='both', expand='yes')
text.insert('end', 'Successfully connected to window')
s2.config(command=text.yview)  # Text随着滚动条移动被控制移动
b2.config(command=text.xview)




root.mainloop()

