import pandas,os,re,io,pdb,time
import tempfile,traceback
def parse():
    df = pandas.read_csv('/Users/dasiweida/PycharmProjects/pangenome/test_data/MCM16_2.MGE.txt', comment="#", sep="\t", dtype=str)###path
    for item in df.itertuples():
        list1=get_seq_local(item[3],item[6],"/Users/dasiweida/PycharmProjects/pangenome/test_data/MCM16_2.megahit.faa",item[1])
        get_seq(list1[1],list1[2],item[1],'/Users/dasiweida/PycharmProjects/pangenome/test_data/final.contigs_1000.fa',"/Users/dasiweida/PycharmProjects/pangenome/test_data")

def get_seq_local(left,right,faa_path,contig):
    # faa_path = faaPath + "/{}.megahit.faa".format(args.name)
    table=os.popen("cat {} |grep '{}'".format(faa_path,contig)).read()
    tableIO = io.StringIO(table)
    df= pandas.read_csv(tableIO,sep="#",header=None,dtype=str)
    seq = [df.iloc[0, 0], df.iloc[eval(left), 1], df.iloc[eval(right), 2]]
    seq[0] = re.search(">.*_.*_", seq[0]).group()[:-1]
    return seq
def get_seq(start,over,contig,contig_path,outpath):
    if eval(start)>eval(over):
        start,over=over,start
    with tempfile.TemporaryDirectory() as seq_dir:
        a=open(seq_dir+"/list.txt","w")
        a.write(contig)
        a.close()
        print("seqkit grep -f {}/list.txt {} |seqkit subseq -r {}:{} -o {}/{}.fasta".format(seq_dir,contig_path,eval(start),eval(over),outpath,contig))
        t=os.popen("seqkit grep -f {}/list.txt {} |seqkit subseq -r {}:{} -o {}/{}.fasta".format(seq_dir,contig_path,eval(start),eval(over),outpath,contig))
        time.sleep(1)
    return

if __name__ =="__main__":
    parse()