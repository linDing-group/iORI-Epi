
import os
import pandas as pd
import numpy as np
import sys

def plusFea(file1,file2,dst_file):

    df = pd.read_csv(file1,sep=',')
    cl=df.columns[:].tolist()
    df2 = pd.read_csv(file2,sep=',')
    cl2 = df2.columns[1:].tolist()

    # df3 = pd.read_csv(file3,sep=',')
    # cl3 = df3.columns[1:].tolist()

    #all_list = [df[cl],df2[cl2],df3[cl3]]
    all_list = [df[cl],df2[cl2]]
    all_df=pd.concat(all_list,axis=1)

    all_df.to_csv(dst_file, index=None)


if __name__ == '__main__':

    os.system('Rscript script/ChIAPET-fea1.R')

    posfile = 'data/ORI/MCF-7/pos.bed'
    negfile = 'data/ORI/MCF-7/neg.bed'

    dicORI = {}
    dicORI['class'] = []
    pos = open(posfile).readlines()
    neg = open(negfile).readlines()

    for j in range(len(pos)):
        dicORI['class'].append('1')

    for j in range(len(neg)):
        dicORI['class'].append('0')

    HMpath = 'data/ChIA-PET/MCF7'
    filename = os.listdir(HMpath)
    for each in filename:
        eachname = each.split('_')[0]
        dicORI[eachname] = []
        valPos = os.popen('bedtools coverage -a '+posfile+' -b '+HMpath+'/'+each).readlines()
        valNeg = os.popen('bedtools coverage -a '+negfile+' -b '+HMpath+'/'+each).readlines()

        for i in valPos:
            i = i.split('\t')[3]
            dicORI[eachname].append(i)
        for i in valNeg:
            i = i.split('\t')[3]
            dicORI[eachname].append(i)
    #print(dicORI)

    dataframe = pd.DataFrame(dicORI)
    dataframe.to_csv(r"ChIA-PET2_Fea.csv",sep=',',index=False)
    plusFea('ChIA-PET1_Fea.csv','ChIA-PET2_Fea.csv','ChIA-PET_Fea.csv')
    os.system('Rscript script/RunRF.R')



     