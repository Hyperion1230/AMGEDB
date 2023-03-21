import pandas as pd
import os,re,glob

def parse_txt(path):
    name=re.search(".*/(.*?).MGE",path).group(1)
    handle=pd.read_csv(path,comment="#",sep="\t")
    ME_num=len(handle[handle["class"].str.contains("ME")])
    ICE_num=len(handle)-ME_num
    record=name+"\t"+str(ICE_num)+"\t"+str(ME_num)+"\n"
    return record

record="sample_name"+"\t"+"ICE"+"\t"+"ME"+'\n'
for path in glob.glob(input('路径')):
    record=record+parse_txt(path)
with open("./statistic.tsv","w") as handle:
    handle.write(record)
print("done")
