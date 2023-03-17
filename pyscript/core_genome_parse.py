import  os, json, re, time,io, numpy ,pdb,glob,warnings,tqdm
import traceback
import argparse
import pandas
import tempfile
import Download_panGenome as DP

##############初始化##################
parser = argparse.ArgumentParser(description='A program used to analyze the core genes of MGE', epilog="version   0.1.0")
parser.add_argument('--url', '-u', dest='url', type=str, nargs="+", help='url table')
# parser.add_argument('--temporary files', '-tmp', dest='tmp', type=str, help='diamond temporart file local', default=".")不再需要临时文件夹了
parser.add_argument('--hmmout', '-ho', dest='hmmout', type=str, help='hmmout path',required=1)
parser.add_argument('--out', '-o', dest='out', type=str, help='out dir', default=".")
parser.add_argument('--faa', '-f', dest='faa', type=str, help='prodigal output faa file dir(no faa file)',required=1)
parser.add_argument('--threads', '-t', dest='threads', type=int,help='Maximum number of threads available to a program(Default to maximum)', default=0)
parser.add_argument('--diamond_index_path','-dip', dest='Diamondindexpath', type=str, help='Diamond indexpath path',required=1)
# parser.add_argument('--faapath', dest='faapath', type=str, help='faapath path',required=1)
parser.add_argument('--taxa',  dest='taxa', type=str, help='taxa path',required=1)
parser.add_argument('--annopath','-ap',dest='annoPath',type=str, help="pangenome anno data dir",required=1)
parser.add_argument('--macsydir','-mp',dest='macsy',type=str,help="macsy result dir")
parser.add_argument('--name','-n', dest='name', type=str, help='out name',required=1)





####################
class MGE:
    left=0
    right=0
    Class=""
    def __init__(self,l,r,C):
        self.left=l
        self.right=r
        self.Class=C


def read_taxonomy(path):
    def a (x):
        if str(x)=="nan":
            return
        else:
            out=x.split(';')
            out=x[1:]
            return str(out)
    handle = pandas.read_csv(path, sep='\t', header=None)
    table = handle.drop([0, 2, 3, 4, 5, 6], axis=1)
    table.rename(columns={1: "contig", 7: 'taxanomy'}, inplace=True)  # 行标题改名
    table=table.set_index("contig")
    table["taxanomy"] = table["taxanomy"].str.split(";", expand=True)[2]
    # table=table.drop(table[table["taxanomy"].str.contains("NA",na=True)])
    # table=table.dropna(axis=0, subset=["taxanomy"])
    # table = table.drop([table["taxanomy"].str.contains("NA", na=True)])
    t2=table.to_dict()#限速步骤
    dict_value=list(t2["taxanomy"].values())
    dict_key=list(t2["taxanomy"].keys())
    dict_value2=list(map(a,dict_value))
    t2={"contig":dict_key,
        "taxa":dict_value2}
    t2=pandas.DataFrame(t2)
    # t2=dict(zip(dict_key,dict_value2))
    return t2


def parser_pangenome(outpath):

    handle = pandas.read_csv(outpath)
    df = handle.iloc[:, 0:2]
    df2 = df[df["Pangenome Class"] == "core"]
    ls = df2["Clusters"].to_list()
    return set(ls)


def diamond(taxa,contig):#返回的是一个pandas的df对象
    #pdb.set_trace()
    with tempfile.TemporaryDirectory() as contig_dir:#临时文件不再保留
        # indexpath="/gss1/home/liujx02/tangyj/proj_MGEdatabase/Gao_data/Pangenome/diamond_index/"
        indexpath=diamond_index
        # faapath="/gss1/home/liujx02/tangyj/proj_MGEdatabase/Gao_data/prodigal/{}.megahit.faa".format(args.name)
        faapath=faaPath+"/{}.megahit.faa".format(args.name)
        # extract_contig=os.system("seqkit grep -r -p '{}' {} -o {}.{}.tmp.faa".format(re.search("k.*_.*_",contig).group(1),faapath,outpath,contig))
        os.system("seqkit grep -r -p '{}_' {} -o {}/{}.tmp.faa".format(contig,faapath,contig_dir,contig))
        index=os.listdir(indexpath)
        # index=list(map(lambda x:os.path.basename(x)),getlist)

        faaname_gene = taxa
        for i in index:
            if faaname_gene[0] in i:
                if faaname_gene[1] in i:
                    global taxa_targe
                    taxa_targe=i[:-5]
                    os.system('diamond blastp -q {}/{}.tmp.faa -d {} -o {}/tmp.dimout -e 0.00673 --very-sensitive --quiet'.format(contig_dir,contig, indexpath+"/"+taxa_targe, contig_dir))
        # os.popen('diamond makedb --in {} -d ./{} -p {}'.format(faaname,args.out,args.threads))
        # if args.tmp == '.':
            # msg1 = os.popen('diamond makedb --in {} -d {} -p {}'.format(faaname, filename, args.threads))
        # os.system('diamond blastp -q {} -d {} -o {} -e 0.00673 --very-sensitive'.format(faa_path, taxa, outpath))
        # else:
        #     msg1 = os.popen('diamond makedb --in {} -d {} -p {}'.format(faaname, filename, args.threads))
        #     msg2 = os.popen('diamond blastp -q {} -d {} -o {} -e 0.00673 --very-sensitive'.format(faa_path, filename, outpath))
                    try:
                        dim_table=pandas.read_csv(contig_dir+"/tmp.dimout",sep='\t',header=None).iloc[:,[0,1,6,7]]
                    except (pandas.errors.EmptyDataError):
                        print("No corresponding Pangenome is matched, skip !!")
                        return "NO"
                    except:
                        traceback.print_exc()
                        return "NO"

                    dim_table=dim_table.rename(columns={0: "contig", 1: 'Cluster'})
                    dim_table["Cluster"] = dim_table.Cluster.map(lambda x: re.search(".*_\d+", x).group(0))
                    return dim_table,taxa_targe
        return "NO"
def extracth(hmmout):#提取hmm文件中比对上的数据
    #这个部分比较不稳定，会受到hmmout文件格式的影响
    #但是目前情况下运行还是比较稳定的
    hmmoutExtract=os.popen('sed -n {},99999999999p {}'.format(15, hmmout)).read()#如果使用的是PFAM的hmm模型就要改成17，其他的用15
    finalLine = hmmoutExtract.find('Domain', 1)
    fileIO=io.StringIO(hmmoutExtract[:finalLine])
    try:
        warnings.filterwarnings("ignore")
        table=pandas.read_csv(fileIO,'\s+',error_bad_lines=False,index_col=False).dropna().T.T#errorbadlines之后的版本将被废弃
        table=table[table['-------'].astype('float')<0.00673]
        list1=numpy.array(table.iloc[:,8:9])
        list2=list1.reshape(len(list1),)
        return list2
    except:
        traceback.print_exc()
        return 'no'

def group_taxa(list,kaijudict):# abandon def
    for item in list:
        taxamony_1=kaijudict[item]
        if taxamony_1=="NA":
            continue
        else:
            df=pandas.DataFrame(kaijudict)
def sort_core(annopath,diamond_out):
    # annopath="/gss1/home/liujx02/tangyj/proj_MGEdatabase/Gao_data/Pangenome/annotation/csv/*."
    annopath=annoPath+"/*."
    annopath=glob.glob(annopath+diamond_out[1]+"*")
    core_set=parser_pangenome(annopath[0])
    df1=diamond_out[0][diamond_out[0]['Cluster'].isin(core_set)]
    getlocal=set(df1['contig'])
    contigset = list(map(lambda x: eval(x.split("_")[2]), getlocal))
    return set(contigset)#maybe use set() can sort the list?

def get_local(local_list , re_local_num , T_local):
    tag='complete'
    for local in local_list:# Judge position
        if local < eval(re_local_num):
            left=local
        elif local == eval(re_local_num):
            tag='incomplete'
            right=local
            break
        else:
            right=local
            break
    if len(T_local)==0:MGE_class='ME';T1="NONE"#if haven't T4SS system contig
    for T in T_local:
        if eval(T)>=left and eval(T)<=right:#T4SS sys inside nocore gene island
            MGE_class="ICE"
            T1=T
            break
        else:
            MGE_class="ME"
            T1=T
    return [left,right,tag,MGE_class,T1]

def macsy_parse(path):
    df = pandas.read_csv(path, header=0, comment='#', sep='\t')
    df = df[df['gene_name'].str.contains("T4SS")]
    ls1=df["hit_id"].to_list()
    return set(ls1)

if __name__=="__main__":
    ######参数处理#######
    args = parser.parse_args()
    diamond_index = args.Diamondindexpath
    faaPath = args.faa
    annoPath = args.annoPath
    macsyPath=args.macsy
    hmmoutPath = args.hmmout
    kaijuPath = args.taxa
    #收集所有contig的物种信息
    while 0:
        taxalist=set()
        for path in os.listdir():
            taxa=read_taxonomy('/Users/dasiweida/Downloads/MFM16_3_kaiju_names.txt')
            ls=list(filter(None,list(taxa.values())))
            taxalist=taxalist|set(ls)
        for taxamony in taxalist:
            url=DP.get_url(taxamony,'/Users/dasiweida/Downloads/downloadData.json')
            if url==None:
                continue
            DP.download_data(url[1],"/Users/dasiweida/Downloads/downloadtest")
            DP.download_data(url[0], "/Users/dasiweida/Downloads/downloadtest")

    ######比对########
    while 1:
        taxaPath=kaijuPath+"/{}_kaiju_names.txt".format(args.name)
        hmmout=hmmoutPath+"/{}.megahit.faa_{}.hmm.hmmout".format(args.name,os.path.basename(args.hmmout))
        # taxa = read_taxonomy("/gss1/home/liujx02/tangyj/proj_MGEdatabase/Gao_data/Species_annotation/kaijuname/{}_kaiju_names.txt".format(args.name))
        taxa = read_taxonomy(taxaPath)
        # contig_list_raw = extracth("/gss1/home/liujx02/tangyj/proj_MGEdatabase/Gao_data/hmmout/ser_lsr/{}.megahit.faa_ser_lsr.hmm.hmmout".format(args.name))
        contig_list_raw = extracth(hmmout)
        contig_list=list(map(lambda x:(re.search(".*?_\d+",x).group(0)),contig_list_raw))
        # contig_list=list(map(lambda x:(re.search(".*?__\d+",x).group(0)),contig_list_raw))#for some special data，e.g Yin_data
        taxa2=taxa[taxa['contig'].isin(contig_list)]
        # T4SS_list=macsy_parse("/gss1/home/liujx02/tangyj/proj_MGEdatabase/Gao_data/macsy/{}.megahit.faa/all_systems.tsv".format(args.name))
        T4SS_list = macsy_parse(macsyPath+"/{}.megahit.faa/all_systems.tsv".format(args.name))
        if os.path.exists(args.out):
            pass
        else:
            os.makedirs(args.out)
        if os.path.exists("{}/{}.MGE.txt".format(args.out,args.name)):
            break
        M1 = open("{}/{}.MGE.txt".format(args.out,args.name), "w")
        record='#ice_identification——result：\ncontig\tspecies\tcore1\trecombinase\tconj\tcore2\tintegrity\tclass\n'
        dict1=zip(taxa2['contig'],taxa2["taxa"])
        for tup in tqdm.tqdm(dict1,total=len(taxa2)):
            # print(tup[0])
            # re_conitg=re.search(".*?_\d+",tup[0]).group(0)
            if tup[1]=="NA" or tup[1]==None:
                continue
            else:
                try:
                    # pdb.set_trace()
                    T4SS_local = []
                    t0=0
                    for T4SS_contig in T4SS_list:
                        # pdb.set_trace()
                        if tup[0]+"_" in T4SS_contig:
                            # print(tup[0],T4SS_contig)
                            T4SS_local.append(T4SS_contig.split("_")[2])
                    vel=tup[1].split(" ")
                    df=diamond(vel,tup[0])
                    if df[0].__str__()=="NO":continue
                    for i in contig_list:
                        if tup[0] in i:
                            local=sort_core(vel,df)
                            break
                    for contig in contig_list_raw:
                        if tup[0] in contig:
                            re_local=contig.split("_")[2]
                    if len(local)==0:continue

                    # pdb.set_trace()

                    if eval(re_local)>local[0] and eval(re_local) < local[-1]:
                        # pdb.set_trace()
                        local_pair=get_local(local,re_local,T4SS_local.sort())
                        record=record+tup[0]+"\t"+tup[1]+"\t"+str(local_pair[0])+"\t"+re_local+"\t"+local_pair[4]+"\t"+str(local_pair[1])+"\t"+local_pair[2]+"\t"+local_pair[3]+'\n'
                        # print(record)

                except:
                    print(tup)
                    traceback.print_exc()
                    break
        M1.write(record)
        print(record)
        M1.close()
        print("statistics：\nhmm_contig:{}\nmacsy:{}\nkaiju_record:{}\n".format(len(contig_list),len(T4SS_list),len(taxa)))
        print('done')
        exit()

