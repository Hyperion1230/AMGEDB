########################################
###THIS SCRIPT will be used to assign taxonomy to assembled contigs 
# https://github.com/bioinformatics-centre/kaiju

#call gcc
export GCC_HOME=/gss1/biosoft/gcc-8.2.0
export PATH=${GCC_HOME}/bin:${PATH}
export LD_LIBRARY_PATH=${GCC_HOME}/lib:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${GCC_HOME}/lib64:${LD_LIBRARY_PATH}
export MANPATH=${GCC_HOME}/share/man:${MANPATH}


base_dir=/gss1/home/liujx02/tangyj/proj_MGEdatabase/Gao_data
in_dir=$base_dir/contig
out_dir=$base_dir/Species_annotation
s_dir=/gss1/home/liujx01/program/kaiju-master/bin
kaiju_db=/gss1/home/liujx01/database/kaijudb_2022
mkdir $out_dir

echo "START AT: "; date
#identify the sequences and file names 
# the input file would be the filtered contigs from step 10; the ARG-containing contigs 
for file in $in_dir/*megahit
do

        file1=$file
        filename=$(basename "$file1")
        samplename=${filename%.megahit}
		#samplename=$(echo $filename | cut -f 1-2 -d "_")


# nr: NCBI BLAST non-redundant protein database "nr", only Archaea, bacteria, and viruses

$s_dir/kaiju -t $kaiju_db/nodes.dmp \
-f $kaiju_db/kaiju_db_nr.fmi \
-z 48 -v -s 65 -a greedy -e 5 \
-i $file1/final.contigs_1000.fa -o $out_dir/${samplename}_kaiju.out \

# add taxa names to output names 
$s_dir/kaiju-addTaxonNames -t $kaiju_db/nodes.dmp -n $kaiju_db/names.dmp \
-i $out_dir/${samplename}_kaiju.out -o $out_dir/${samplename}_kaiju_names.txt -r family,genus,species


done

# classfication summary 
# species
$s_dir/kaiju2table -t $kaiju_db/nodes.dmp -n $kaiju_db/names.dmp \
-r species -o $out_dir/kaiju_summary_species.tsv $out_dir/*_kaiju.out

# genes
$s_dir/kaiju2table -t $kaiju_db/nodes.dmp -n $kaiju_db/names.dmp \
-r genus -o $out_dir/kaiju_summary_genus.tsv $out_dir/*_kaiju.out

# family
$s_dir/kaiju2table -t $kaiju_db/nodes.dmp -n $kaiju_db/names.dmp \
-r family -o $out_dir/kaiju_summary_family.tsv $out_dir/*_kaiju.out

echo "DONE AT: "; date





