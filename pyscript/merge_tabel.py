import pandas
import openpyxl
import os
from Bio import SeqIO
def read_table(path):
    table=pandas.read_csv(path,comment="#",sep="\t",dtype=str)
    return table

def get_seq(path):
    with open(path,'r')as handle:
        record=SeqIO.parse(handle,format="fasta")
        seq_list=list(record)
        return seq_list

def create_new_df(seq_list):
    df1=None
    for i in seq_list:
        df=pandas.DataFrame([[i.id,i.seq.__str__()]],columns=["id","seq"])
        df1=pandas.concat([df1,df])
    return df1
def merge(table,df1):
    if len(table)==len(df1):
        t1 = pandas.merge(table, df1, left_on="contig", right_on="id", how='left')
        t1=t1.drop('id',axis=1)
        # df=pandas.concat([table,df1],axis=1)
    return t1

if __name__=="__main__":
    df0 = None
    for i in os.listdir("/Users/dasiweida/proj/T_extract"):
        try:
            table=read_table('/Users/dasiweida/proj/T_py_Pangenome/ICE_result/ser_lsr/{}.MGE.txt'.format(i))
            seq_list=get_seq('/Users/dasiweida/proj/T_extract/{}/merge.fasta'.format(i))
            new_df=create_new_df(seq_list)
            writer=merge(table,new_df)
            # excel=pandas.ExcelWriter('/Users/dasiweida/proj/T_extract/AD1/AD1')
            df0 = pandas.concat([df0, writer], axis=0)
            # excel.save()
        except FileNotFoundError:
            print("not find {} file or it is empty".format(i))
    df0.to_excel("/Users/dasiweida/proj/MERGE.xlsx")