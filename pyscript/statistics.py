import pandas as pd
import os, re, glob, time
from functools import wraps


def timer(func):
    @wraps(func)
    def inner(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        end = time.time()
        print("消耗时间为{}秒".format(end - start))
        return result

    return inner


@timer
def parse_txt(path):
    name = re.search(".*/(.*?).MGE", path).group(1)
    handle = pd.read_csv(path, comment="#", sep="\t")
    ME_num = len(handle[handle["class"].str.contains("ME")])
    ICE_num = len(handle) - ME_num
    record = name + "\t" + str(ICE_num) + "\t" + str(ME_num) + "\n"
    return record


def getPATHlist(Path1):
    pathlist = os.listdir(Path1)  # MAC适用
    pathlist.sort()
    pathlist.pop(0)
    #############################去除MAC的.DS_store
    PathList = list(map(lambda x: '/Users/dasiweida/proj/ICE/' + str(x), pathlist))
    return PathList


def merge(path2):
    Name = os.path.basename(path2).split('.')[0]
    tb = pd.read_csv(path2, comment="#", sep="\t")
    count = len(tb)
    sampleNameList = [[Name] for _ in range(count)]
    pd1 = pd.DataFrame(data=sampleNameList, columns=['sampleName'])
    ResutlTable = pd.concat([pd1, tb], axis=1)
    return ResutlTable

@timer
def summary():
    tbr=None
    l1=getPATHlist('/Users/dasiweida/proj/ICE')
    for i in l1:
        for i2 in os.listdir(i):
            i2=str(i)+'/'+str(i2)
            tb1=merge(i2)
            tbr=pd.concat([tbr,tb1])
    tbr.reset_index()
    return tbr
if __name__ == '__main__':
    while 0:
        record = "sample_name" + "\t" + "ICE" + "\t" + "ME" + '\n'
        for path in glob.glob(input('路径')):
            record = record + parse_txt(path)
        with open("./statistic.tsv", "w") as handle:
            handle.write(record)
        print("done")

    tb=summary()
    tb.to_excel("/Users/dasiweida/proj/total.xlsx")