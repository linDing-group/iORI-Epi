import pandas as pd
import numpy as np
import sys

def plusFea(file1,file2,file3,dst_file):

    df = pd.read_csv(file1,sep=',')
    cl=df.columns[:].tolist()
    df2 = pd.read_csv(file2,sep=',')
    cl2 = df2.columns[1:].tolist()

    df3 = pd.read_csv(file3,sep=',')
    cl3 = df3.columns[1:].tolist()

    all_list = [df[cl],df2[cl2],df3[cl3]]
    #all_list = [df[cl],df2[cl2]]
    all_df=pd.concat(all_list,axis=1)

    all_df.to_csv(dst_file, index=None)

# file1 = sys.argv[1]
# file2 = sys.argv[2]
# file3 = sys.argv[3]
# dst_file = sys.argv[4]

file1 = 'test_Epimarks_Fea.csv'
file2 = 'test_3D_Fea.csv'
file3 = 'test_Motif_Fea.csv'
dst_file = 'test_allFea.csv'

plusFea(file1,file2,file3,dst_file)