## A custom megaphage identified pipline

```
###01 MGV virus identification process
cd ~ && mkdir -p 01_MGV_viral_analysis && cd 01_MGV_viral_analysis
cp -r ~/MGV/viral_detection_pipeline ./

## Locate input metagenomic assembly
#cp ../combined.all.fasta ./viral_detection_pipeline/input/combined.all.fna
## Call viral genes
cd viral_detection_pipeline
prodigal -i $seq -a input/combined.all.faa -d input/combined.all.ffn -p meta -f gff > input/combined.all.gff

# Run HMMER on imgvr and pfam databases
hmmsearch -Z 1 --cpu 32 --noali --tblout output/imgvr.out input/imgvr.hmm input/combined.all.faa
hmmsearch -Z 1 --cut_tc --cpu 32 --noali  --tblout output/pfam.out input/pfam.hmm input/combined.all.faa

#Count genes hitting viral and microbial marker genes
python count_hmm_hits.py input/combined.all.fna input/combined.all.faa output/imgvr.out output/pfam.out > output/hmm_hits.tsv

# Run VirFinder
Rscript virfinder.R combined.all.fna output/combined.all.virfinder.tsv

#Quantify strand switch rate of genes
python strand_switch.py input/combined.all.fna input/combined.all.faa > output/strand_switch.tsv

#Create master table of sequence features
python master_table.py output/hmm_hits.tsv output/virfinder.tsv output/strand_switch.tsv > output/master_table.tsv
#Predict viral contigs
python viral_classify.py \
--features output/master_table.tsv \
--in_base input/combined.all \
--out_base output/combined.all
mv output/combined.all.fna ../../combined.all.mgv-check.fasta 


###02 Sequence statistic length, select 200 kb sequences
less combined.all.mgv-check.fasta | seqkit fx2tab --length --name > combined.all.mgv-check.length
less combined.all.mgv-check.length | awk -F '\t' '{if($2>=200000) print $0}' > combined.all.200k.length
less combined.all.200k.length | awk -F '\t' '{print $1}' > combined.all.200k.IDs
seqkit grep -j 100 -f combined.all.200k.IDs ../combined.all.mgv-check.fasta -o combined.all.200k.IDs.fasta


###03 checkv Assessment of virus integrity and removal of some low-quality areas
source activate checkv
Database=~/Database/checkv_softdb/checkv-db-v1.0
checkv end_to_end -t 90 -d ${Database} combined.all.200k.IDs.fasta combined.all.200k_checkv

cat combined.all.200k_checkv/proviruses.fna |sed 's/^>/>Proviruses_/g' >> combined.all.200k.checkv.fna
cat combined.all.200k_checkv/viruses.fna | sed 's/^>/>Viruses_/g' >> combined.all.200k.checkv.fna
cat combined.all.200k_checkv/quality_summary.tsv >> combined.all.200k.checkv.quality.summary.tsv

###04 Access to genetic codes and proteins
## >100kb >5%
## GC content evaluation
python ~/scripts/custom_scripts/gc_content.py combined.all.200k.checkv.fna 50 GC_up_50.IDs.fasta GC_low_50.IDs.fasta
python ~/scripts/custom_scripts/genetic_codes_evalution.py combined.all.200k.checkv.fna > combined.all.200k.checkv.gff
less combined.all.200k.checkv.gff | grep "11$" > combined.all.200k.11.IDs 
less combined.all.200k.checkv.gff | grep "11$" > combined.all.200k.15.IDs 
less combined.all.200k.checkv.gff | grep "11$" > combined.all.200k.90.IDs 
less combined.all.200k.checkv.gff | grep "11$" > combined.all.200k.91.IDs # 0
seqkit -f combined.all.200k.11.IDs ./combined.all.200k.checkv.fna -o combined.all.200k.11.fasta
seqkit -f combined.all.200k.15.IDs ./combined.all.200k.checkv.fna -o combined.all.200k.15.fasta
seqkit -f combined.all.200k.90.IDs ./combined.all.200k.checkv.fna -o combined.all.200k.90.fasta
seqkit -f combined.all.200k.91.IDs ./combined.all.200k.checkv.fna -o combined.all.200k.91.fasta

~/soft/prodigal -i combined.all.200k.11.fasta -a combined.all.200k.11.faa -d combined.all.200k.11.fna -g 11 -f gff -o combined.all.200k.11.gff
~/soft/prodigal -i combined.all.200k.15.fasta -a combined.all.200k.15.faa -d combined.all.200k.15.fna -g 15 -f gff -o combined.all.200k.15.gff
~/soft/prodigal -i combined.all.200k.90.fasta -a combined.all.200k.90.faa -d combined.all.200k.90.fna -g 90 -f gff -o combined.all.200k.90.gff

cat combined.all.200k.11.fasta combined.all.200k.15.fasta combined.all.200k.90.fasta > final.combined.all.200k.fasta
cat combined.all.200k.11.faa combined.all.200k.15.faa combined.all.200k.90.faa > final.combined.all.200k.protein.fasta
cat combined.all.200k.11.fna combined.all.200k.15.fna combined.all.200k.90.fna > final.combined.all.200k.nucl.fasta
cat combined.all.200k.11.gff combined.all.200k.15.gff combined.all.200k.90.gff > final.combined.all.200k.gff

### 05 Identification of the megaophage genome
#05_01:Function
#1 custom terminase database, Structural proteins of megaphages are difficult to annotate
cd ~ && mkdir -p /Function/diamond_db
cd Function
## diamond blastp 50 scores 
source activate diamond
diamond makedb --in ~/Database/terminase-and-major-capsid_custom_database/new.terminase.all.fasta -d diamond_db/terminase.all -p 50
diamond blastp -d diamond_db/terminase.all -q ../final.combined.all.200k.protein.fasta -o custom.terminase.matches.m8 -f 6 -p 90 --min-score 50 -e 1e-5
python get_diamond_best.py custom.terminase.matches.m8 > 01_diamond.function.info 

#2 pfama tigrfam vogdb
less 01_diamond.function.info | awk -F '\t' '{print $1}' > 1.IDs
less ../final.combined.all.200k.protein.fasta | grep ">" | awk -F ' ' '{print $1}' | sed 's/>//g' > 2.IDs
less 1.IDs 2.IDs | sort | uniq -u > eles.IDs # 3524580
seqkit grep -f eles.IDs ../final.combined.all.200k.protein.fasta -o eles.IDs.fasta
rm 1.IDs 2.IDs eles.IDs

# pfama
hmmsearch -Z 1 --cpu 100 --tblout eles.pfama.tb.csv -E 0.00001 ~/Database/pfam/Pfam-A.hmm eles.IDs.fasta 
## tigrfam
hmmsearch -Z 1 --cpu 100 --tblout eles.tigrfam.tb.csv -E 0.00001 ~/Database/TIGRFAM/all.trigrfam.HMM eles.IDs.fasta
## vogdb
hmmsearch -Z 1 --cpu 100 --tblout eles.vogdb.tb.csv -E 0.00001 ~/Database/VOGDB/all.vogdb.hmm eles.IDs.fasta

cat eles.pfama.tb.csv eles.tigrfam.tb.csv eles.vogdb.tb.csv > else.all.tb.csv
python get_hmmsearch_best.py else.all.tb.csv > 02_diamond.function.info
cat 01_diamond.function.info 02_diamond.function.info > final.function.info

#3 filter megaphages based on special proteins
## (capsid\|tail\|terminase\|spike\|holin\|portal\|baseplate\|virus\|viral\|phage\|prohead\|tape measure\|head\|DNA packaging\|base plate\|virion\|tapemeasure neck\|bacteriophage\|prophage\|T4\|p22) # select special proeteins
## (tRNA synthetase, tRNA synthase, ribosomal protein, preprotein translocase, DNA gyrase subunit A) # remove proteins associated with prokaryotes
less final.function.info | grep -i "capsid\|tail\|terminase\|spike\|holin\|portal\|baseplate\|terminase\|prohead\|tape measure\|head\|DNA packaging\|base plate\|virion\|tapemeasure neck\|bacteriophage\|p22" > tmp.prepare 
less tmp.prepare | grep -v "GET4\|ABC_choline_bnd\|AMT4_dom_C\|antiphage_ZorA_1\|antiphage_ZorA_2\|antiphage_ZorA_3\|antiphage_ZorA_4\|Cadherin_tail\|CDRT4\|Choline_bind_1\|Choline_bind_2\|Choline_bind_3\|choline_CutC\|Choline_kinase\|Choline_kin_N\|Choline_sulf_C\|Choline_transpo\|CT47\|Filament_head\|Forkhead\|Forkhead_N\|Headcase\|Keratin_2_head\|Keratin_2_tail\|Myosin_head\|Myosin_tail_1\|myxo_SS_tail\|omphalotin_tail\|Por_Secre_tail\|IFT43\|IFT46_B_C\|Lac_bphage_repr\|Mpv17_PMP22\|onco_T4SS_CagA\|Ost4\|RAC_head\|Stork_head\|REFSEQ@Gp22.23\|REFSEQ@gp22,@hypothetical@protein\|sat4\|sp|P22409|TRA_STRLI@Protein@Tra\|sp|P22877|INTR_SACER@Integrase\|sp|P71276|RT48_ECOLX@Retron@Ec48@reverse@transcriptase\|sp|Q80DT4|PROF_CWPXG@Profilin\|Spt4\|Spt46\|T4BSS_DotH_IcmK\|T4BSS_DotI_IcmL\|T4bSS_IcmS\|T4_deiodinase\|T4SS\|T4SS_assoc_EirA\|T4SS_CagC\|T4SS_ComB10\|T4SS-DNA_transf\|T4SS_DotA\|T4SS_effec_BID\|T4SS_IcmE_DotG\|T4SS_lipo_DotD\|T4SS_pilin\|T4SS_pilin_TrwL\|T4SS_TraI\|T4SS_VirB10\|tail_200_repeat\|tail-700\|Uso1_p115_head\|WIAG-tail\|YadA_head\|nucleocapsid\|REFSEQ@antiholin\|sp|P03705|HOLIN_LAMBD@Antiholin" | awk -F '\t' '{print $1}' | grep -v "Proviruses_" | sort -u | awk -F '_' '{$NF=""; print $0}' | sed 's/ /_/g' | sed 's/_$//g' | awk -F '_' '{$NF=""; print $0}' | sed 's/ /_/g' | sed 's/_$//g' | sort | uniq -c | sed 's/^ *//g' | sort -n | awk -F ' ' '{if($1>=2) print $2}' > tmp.combined.all.200k.by-function.phage.IDs
less final.function.info | grep -i "tRNA synthetase\|tRNA synthase\|ribosomal protein\|preprotein translocase\|DNA gyrase subunit A" >> test.info # 0
cat tmp.combined.all.200k.by-function.phage.IDs | sort -u > combined.all.200k.by-function.phage.IDs

#05_02:Host predict
mkdir -p blastn_index
makeblastdb -in final.combined.all.200k.fasta -dbtype nucl -out blastn_index/combined.all.200k
blastn -query Database/db_uhgg_gtdb_meta_spacer.fa -out ./spacer.blastn.result.out -num_threads 90 \
-db blastn_index/combined.all.200k \
-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' \
-max_target_seqs 10000000 -dust no -word_size 8 -evalue 10
python get_blastn_best.py spacer.blastn.result.out > spacer.blastn.result # cov=100, 0/1 mismatch
less spacer.blastn.result | awk -F '\t' '{print $1}' | sort -u > combined.all.200k.by-spacer.phage.IDs

#05_03:Phamer predict
python ~/soft/PhaMer/PhaMer_preprocessing.py --contigs final.combined.all.200k.fasta --midfolder ./phamer_predict --threads 100 --dbdir Database/PhaMer &
less phamer_predict/phamer.out.csv | grep ",phage," | awk -F ',' '{print $1}' > combined.all.200k.by-phamer.phage.IDs

# Phages predict
## 1. combine the results from function and host prediction
less combined.all.200k.by-spacer.phage.IDs combined.all.200k.by-function.phage.IDs | sort | uniq -d > combined.all.200k.by-host-And-function.phage.IDs
## 2. Add Phamer result
less combined.all.200k.by-phamer.phage.IDs combined.all.200k.by-host-And-function.phage.IDs | sort -u > combined.all.200k.by-host-And-function-And-phamer.phage.IDs 
```














 