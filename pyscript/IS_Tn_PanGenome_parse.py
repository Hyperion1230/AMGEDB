import core_genome_parse as CGP
import argparse
import os,pdb
import re
from tqdm import tqdm
parser = argparse.ArgumentParser(description='A program used to analyze the core genes of MGE', epilog="version   0.1.0")
parser.add_argument('--url', '-u', dest='url', type=str, nargs="+", help='url table')
# parser.add_argument('--temporary files', '-tmp', dest='tmp', type=str, help='diamond temporart file local', default=".")不再需要临时文件夹了
parser.add_argument('--hmmout', '-ho', dest='hmmout', type=str, help='hmmout path',required=True)
parser.add_argument('--out', '-o', dest='out', type=str, help='out dir', default=".")
parser.add_argument('--faa', '-f', dest='faa', type=str, help='prodigal output faa file dir(no faa file)',required=1)
parser.add_argument('--threads', '-t', dest='threads', type=int,help='Maximum number of threads available to a program(Default to maximum)', default=0)
parser.add_argument('--diamond_index_path','-dip', dest='Diamondindexpath', type=str, help='Diamond indexpath path',required=1)
# parser.add_argument('--faapath', dest='faapath', type=str, help='faapath path',required=1)
parser.add_argument('--taxa',  dest='taxa', type=str, help='taxa path',required=1)
parser.add_argument('--annopath','-ap',dest='annoPath',type=str, help="pangenome anno data dir",required=1)
parser.add_argument('--name','-n', dest='name', type=str, help='out name',required=1)

args = parser.parse_args()

######参数处理#######
diamond_index=args.Diamondindexpath
faaPath=args.faa
annoPath=args.annoPath
# macsyPath=args.macsy
hmmoutPath=args.hmmout
kaijuPath=args.taxa

####################
taxaPath = kaijuPath + "/{}_kaiju_names.txt".format(args.name)
hmmout = hmmoutPath + "/{}.megahit.faa_{}.hmm.hmmout".format(args.name, os.path.basename(args.hmmout))
taxa = CGP.read_taxonomy(taxaPath)
contig_list_raw = CGP.extracth(hmmout)
contig_list=list(map(lambda x:(re.search(".*?_\d+",x).group(0)),contig_list_raw))
pdb.set_trace()
taxa2=taxa[taxa['contig'].isin(contig_list)]
if os.path.exists(args.out):
    pass
else:
    os.makedirs(args.out)
if os.path.exists("{}/{}.MGE.txt".format(args.out, args.name)):
    exit()
with open("{}/{}.MGE.txt".format(args.out,args.name), "w") as handle:
    record = '#ice_identification——result：\ncontig\tspecies\tcore1\trecombinase\tcore2\tintegrity\tclass\n'
    dict1 = zip(taxa2['contig'], taxa2["taxa"])
    for tup in tqdm(dict1, total=len(taxa2)):
        if tup[1] == "NA" or tup[1] == None:
            continue
        else:
            vel = tup[1].split(" ")
            df = CGP.diamond(vel, tup[0])
            if df.__str__() == "NO": continue
            for i in contig_list:
                if tup[0] in i:
                    locals = CGP.sort_core(vel, df)
                    break
                for contig in contig_list_raw:
                    if tup[0] in contig:
                        re_local = contig.split("_")[2]
                if len(locals) == 0: continue
                if eval(re_local) > locals[0] and eval(re_local) < locals[-1]:
                    tag = 'complete'
                    for local in locals:
                        if local < eval(re_local):
                            left = local
                        elif local == eval(re_local):
                            tag = 'incomplete'
                            right = local
                            break
                        else:
                            right = local
                            break
                record = record + tup[0] + "\t" + tup[1] + "\t" + str(left) + "\t" + re_local + "\t" + str(right) + "\t" + tag + "\t" + "RI" + '\n'
print("done")
