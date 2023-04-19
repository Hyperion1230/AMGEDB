import pandas,os,re,io,pdb,time
import tempfile,traceback
import argparse
import logging

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(filename)s [%(levelname)s] %(message)s',
                    datefmt='%a %d %b %Y %H:%M:%S',
                    filename='my.log',
                    filemode='w')

from functools import wraps
parser = argparse.ArgumentParser(description='It is used to extract contig fragments corresponding to MGE', epilog="version   0.1.0")
parser.add_argument('--table', '-t', dest='table', type=str, help='mge table',required=1)
parser.add_argument('--contig_files', '-c', dest='contig', type=str, help='contig file local',required=1)
parser.add_argument('--faa_file', '-f', dest='faa', type=str, help='faa path',required=1)
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
    df = pandas.read_csv(table, comment="#", sep="\t", dtype=str)###path
    for item in df.itertuples():
        list1=get_seq_local(item[4],item[7],faa,item[2])
        logging.info("left:{}<------->reigt:{},contig_name:{}".format(list1[1],list1[2],item[1]))
        get_seq(list1[1],list1[2],item[1],contig,outpath)

def get_seq_local(left,right,faa_path,contig):
    try:
        table=os.popen("cat {} |grep '{}'".format(faa_path,contig)).read()
        tableIO = io.StringIO(table)
        df= pandas.read_csv(tableIO,sep="#",header=None,dtype=str)
        seq = [df.iloc[0, 0], df.iloc[eval(left)-1, 1], df.iloc[eval(right)-1, 2]]
        seq[0] = re.search(">.*_.*_", seq[0]).group()[:-1]
        return seq
    except:

        logging.error("model {} have something error".format("get_seq_local"))
        logging.error("the file is {} and {}".format(faa_path,contig))
        pdb.set_trace()

def get_seq(start,over,contig,contig_path,outpath):
    if eval(start)>eval(over):
        start,over=over,start
    with tempfile.TemporaryDirectory() as seq_dir:
        a=open(seq_dir+"/list.txt","w")
        a.write(contig)
        a.close()
        logging.info("seqkit grep -f {}/list.txt {}/final.contig_1000.fa |seqkit subseq -r {}:{} -o {}/{}.fasta".format(seq_dir,contig_path,eval(start),eval(over),outpath,contig))
        t=os.popen("seqkit grep -f {}/list.txt {}/final.contigs_1000.fa |seqkit subseq -r {}:{} -o {}/{}.fasta".format(seq_dir,contig_path,eval(start),eval(over),outpath,contig)).read()
        # time.sleep(0.5)
    return

if __name__ =="__main__":
    args = parser.parse_args()
    if os.path.exists(args.out):
        pass
    else:
        os.makedirs(args.out)
    parse(args.table,args.faa,args.contig,args.out)
