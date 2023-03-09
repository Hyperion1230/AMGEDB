import os,time,traceback
import json,re
def download_data(url, outpath):
    try:
        t1 = time.time()
        # global filename  # ！！！！！！！！！！！需要测试！！！！！！！！！！！！
        # filename = wget.download(url, out=outpath)  # 这里还需要测试一次文件名的格式！！！！！！！！
        os.system("wget -P {} {}".format(outpath,url))
        t2 = time.time()
        print("Done!\nThe download time is {} seconds".format(t2 - t1))
    except:
        print("Download Fail!,The eorro message is")
        traceback.print_exc()

def get_url(taxonomy, path):
    handle = open(path, 'r', )
    json_obj = json.load(handle)
    for record in json_obj:
        t=0
        if re.search(taxonomy, record['species']) != None:
            t=1
            print(record['species'])
            break
        # else:
        #     # print(record['species'])
        #     continue
    if t==1:
        url_anno = "https://ngdc.cncb.ac.cn/propan/downloadData/annotationsMatrix/" + record['panMc']
        url_pro = "https://ngdc.cncb.ac.cn/propan/downloadData/proteinSequence/" + record['seqpro']
        url_list = [url_anno, url_pro]
        return url_list
    else:
        return