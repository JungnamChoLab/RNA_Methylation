#Here use HS Col-0 as an example, DENA was used to find m6A sites along transcripts
python3 DENA/step4_predict/LSTM_extract.py get_pos --fasta TAIR10_cdna_onsen1  --motif 'RRACH' --output ./candidate_predict_pos_onsen.txt
cd ../nanopore_DRS_m6A/1.data/HS-Col-0/20220803_0857_X3_FAT94668_ec441eab/fastq_pass/
cat FAT94668_pass_c80d3cc8_* > pass.fq
cd ..
conda activate DENA
multi_to_single_fast5 -t 20 -i fast5_pass -s single_reads
mv single_reads fast5_pass
f5="~/nanopore_DRS_m6A/1.data/HS-Col-0/20220803_0857_X3_FAT94668_ec441eab/fast5_pass/"
fq="~/nanopore_DRS_m6A/1.data/HS-Col-0/20220803_0857_X3_FAT94668_ec441eab/fastq_pass/pass.fq"
file="hs_col0"
tombo resquiggle --rna --processes 20 --corrected-group RawGenomeCorrected_001 --basecall-group Basecall_1D_000 --include-event-stdev --overwrite --ignore-read-locks ${f5}/single_reads/ TAIR10_cdna_onsen1
minimap2/minimap2 -ax map-ont -L -p 0 -N 10 TAIR10_cdna_onsen1 $fq | samtools view -bh -F 2324 | samtools sort -O bam > 2.mapping/${file}_basecalls.bam
samtools index 2.mapping/${file}_basecalls.bam
mkdir 3.dena/${file}.tmp
cd 3.dena/${file}.tmp
python ~/nanopore_DRS_m6A/DENA/step4_predict/LSTM_extract.py predict --processes 20 --fast5 ${f5}/single_reads/ --corr_grp "RawGenomeCorrected_001" --bam ~/nanopore_DRS_m6A/2.mapping/${file}_basecalls.bam --sites ~/nanopore_DRS_m6A/candidate_predict_pos_onsen.txt --label 'unknown' --windows 2 2
cd ../
python ~/nanopore_DRS_m6A/DENA/step4_predict/LSTM_predict.py -i ${file}.tmp -m ~/nanopore_DRS_m6A/model -o ${file}_m6a_onsen -p ${file} -d

#Use Nanocount get TPM value as expression level
NanoCount --extra_tx_info -i ../2.mapping/hs_col0_basecalls.bam -o col0_counts.tsv -b hs_col0.bam  --max_dist_3_prime -1

#use R package "GenomicFeatures" to convert transcriptome coordinate to genome coordinate
library("GenomicFeatures")
Txdb_gtf <- makeTxDbFromGFF("~/genome/TAIR10_GFF3_genes_exons.gtf")
GRList <- exonsBy(Txdb_gtf, by = "tx")
names(GRList ) <- id2name(Txdb_gtf, "tx")
tx_coor <- read.delim("hs_col0.tsv",header=F)
x  <- GRanges(tx_coor$V1, IRanges(tx_coor$V2,tx_coor$V2+1))
xx <- mapFromTranscripts(x,GRList,ignore.strand = FALSE)
hs_col0 <- cbind(as.data.frame(xx),tx_coor)
write.table(hs_col0,"hs_col0_genome.tsv",sep="\t",quote=F,row.names = F)

#Bedgraph file was generated from m6a sites whose m6a ratio above 0.1 and reads number above 10, which can be shown in IGV
awk '{if ($5=="+")print $1"\t"$2"\t"$5"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13;if ($5=="-")print $1"\t"$3"\t"$5"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13}' hs_col0_genome.tsv | sort -k1,1 -k2,2n > hs_col0_genomecor.tsv
awk '$9>=0.1 && $8 >=10' hs_col0_genomecor.tsv |awk '{print "Chr"$1"\t"$2-1"\t"$2"\t"$4"\t"$5"\t"$3"\t"$6"\t"$7"\t"$8"\t"$9}'|sort -k1,1 -k2,2n -k9,9nr |awk '!a[$2]++{print}'| sed 's#Chr##' > hs_col0_m6a_unique.bed
awk -v OFS="\t" '{print $1,$2,$3,$10}' hs_col0_m6a_unique.bed > hs_col0_m6a_unique.bedgraph

#m6A level was calculated by average m6A raito from modified sites
hs_col0 <- read.delim("hs_col0_genomecor.tsv",header=F)
hs_col0_dena_f <- aggregate(hs_col0[hs_col0$V9>=0.1&hs_col0$V8>=10,]$V9, by=list(trans=hs_col0[hs_col0$V9>=0.1&hs_col0$V8>=10,]$V4),mean)
