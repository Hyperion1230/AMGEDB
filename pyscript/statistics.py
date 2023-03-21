import pandas as pd
import os,re,glob

def parse_txt(path):
    name=re.search(".*/(.*?).MGE",path).group(1)
    with pd.read_csv(path,commit="#",sep="\t") as handle:
        ME_num=len(handle[handle["class"].str.contains("ME")])
        ICE_num=len(handle)-ME_num
        record=name+"\t"+ICE_num+"\t"+ME_num+"\n"
    return record

record="sample_name"+"\t"+"ICE"+"\t"+"ME"
for path in glob.glob(input('路径')):
    record=parse_txt(path)
with open("./statistic.tsv","w") as handle:
    handle.write(record)
print("done")
