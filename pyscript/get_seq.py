import pandas,os,re,io,pdb,time
import tempfile,traceback
import argparse
import logging

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(filename)s [%(levelname)s] %(message)s',
                    datefmt='%a %d %b %Y %H:%M:%S',
                    # filename='my.log',
                    # filemode='a+'
                    )
from functools import wraps
parser = argparse.ArgumentParser(description='It is used to extract contig fragments corresponding to MGE', epilog="version   0.1.0")
parser.add_argument('--table', '-t', dest='table', type=str, help='mge table',required=True)
parser.add_argument('--contig_files', '-c', dest='contig', type=str, help='contig file local',required=True)
parser.add_argument('--faa_file', '-f', dest='faa', type=str, help='faa path',required=True)
parser.add_argument('--out', '-o', dest='out', type=str, help='out dir', default=".")

def timer(func):
    @wraps(func)
    def inner(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        end = time.time()
        logging.info("消耗时间为{}秒".format(end - start))
        return result

    return inner
def parse(table,faa,contig,outpath):
    try:
        df = pandas.read_csv(table, comment="#", sep="\t", dtype=str)###path
    except pandas.errors.EmptyDataError:
        logging.warning("the MGE table {} is empty!!".format(os.path.basename(table)))
        exit()
    for item in df.itertuples():
        list1=get_seq_local(item[3],item[6],faa,item[1])
        get_orf_file(item[3],item[6],faa,contig,outpath)
        logging.info("left:{}<------->reigt:{},contig_name:{}".format(list1[1],list1[2],item[1]))
        get_seq(list1[1],list1[2],item[1],contig,outpath)
def get_orf_file(left,right,faa_path,contig,outpath):
    with tempfile.TemporaryDirectory() as orf_dir:
        name_list=""
        for i in range(left-1,right):
            name_list=name_list+str(contig)+"_"+str(i)+"\t"
        a = open(orf_dir + "/orf_list.txt", "w")
        a.write(name_list)
        a.close()
        faaOutPath=outpath+"/"+contig+".faa"
        os.popen("seqkit grep -f {}/orf_list.txt {}/protein.faa -o {}".format(orf_dir,faa_path,faaOutPath)).read()

def get_seq_local(left,right,faa_path,contig):
    try:
        table=os.popen("cat {} |grep '{}'".format(faa_path,contig)).read()
        tableIO = io.StringIO(table)
        df= pandas.read_csv(tableIO,sep="#",header=None,dtype=str)
        seq = [df.iloc[0, 0], df.iloc[eval(left)-1, 1], df.iloc[eval(right)-1, 2]]#此步骤保证了核心基因被剔除在外
        seq[0] = re.search(">.*_.*_", seq[0]).group()[:-1]
        return seq
    except pandas.errors.EmptyDataError:
        logging.warning("the file is {} and {}".format(faa_path, contig))
        logging.warning("the {} MGE table is empty!!".format(re.search("(PRJ.*)/pro.*",faa_path).group(0)))
        exit()
    except:
        traceback.print_exc()
        logging.error("model {} have something error".format("get_seq_local"))
        logging.error("the file is {} and {}".format(faa_path,contig))
        pdb.set_trace()
        exit()

def get_seq(start,over,contig,contig_path,outpath):
    if eval(start)>eval(over):
        start,over=over,start
    with tempfile.TemporaryDirectory() as seq_dir:
        a=open(seq_dir+"/list.txt","w")
        a.write(contig)
        a.close()
        logging.info("seqkit grep -f {}/list.txt {}/final.contig.fa |seqkit subseq -r {}:{} -o {}/{}.fasta".format(seq_dir,contig_path,eval(start),eval(over),outpath,contig))
        t=os.popen("seqkit grep -f {}/list.txt {}/final.contigs.fa |seqkit subseq -r {}:{} -o {}/{}.fasta".format(seq_dir,contig_path,eval(start),eval(over),outpath,contig)).read()
        # time.sleep(0.5)
    return

if __name__ =="__main__":
    args = parser.parse_args()
    if os.path.exists(args.out):
        pass
    else:
        os.makedirs(args.out)
    parse(args.table,args.faa,args.contig,args.out)
