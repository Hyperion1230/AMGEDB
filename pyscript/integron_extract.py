import pandas
import os
import functools
import time
import logging
import traceback
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(filename)s [%(levelname)s] %(message)s',
                    datefmt='%a %d %b %Y %H:%M:%S',
                    # filename='my.log',
                    # filemode='w'
                    )
def timer(func):
    functools.wraps(func)
    def inner(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        end = time.time()
        logging.info("消耗时间为{}秒".format(end - start))
        return result
    return inner
@timer
def parse_table(path,name):
    sampleName=name
    try:
        table_integron=pandas.read_csv(path,comment="#",sep="\t")
    except:
        logging.warning(traceback.print_exc())
        return "cont"
    table_pure=table_integron[table_integron['ID_integron']!="ID_integron"]
    finalTable = table_pure[table_pure["type"] == 'complete']#筛选出complete的数据
    nameList=[sampleName for _ in range(len(finalTable))]
    nameTable=pandas.DataFrame(data=nameList,columns=["sample_name"])#构建一个以样本名字为数据的列
    tableOutput=pandas.concat([finalTable,nameTable],axis=1)#合并两个表格
    return tableOutput

if __name__ == '__main__':
    paths=os.listdir("/gss2/home_new/liujx02/integron/07_integron")
    table=pandas.DataFrame()
    for sampleName in paths:
        logging.info("start sample {} parse ...".format(sampleName))
        absPath="/gss2/home_new/liujx02/integron/07_integron/"+sampleName+"/final.contigs.integrons"#简易的输入integreon数据表
        table_son=parse_table(absPath,sampleName)
    if table_son == "cont":
            logging.warning("empty table !!")
            continue
        logging.info("The complete INTEGRON is available in {}".format(table_son.groupby("ID_replicon").ngroups))
        table=pandas.concat([table,table_son])
    table.to_csv("final_integron.txv",sep='\t')
