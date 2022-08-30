###All RNA-seq and SG-RNA-seq follows the same process
mkdir 0_fastqc 3_trim 4_hisat2 5_stringtie 6_featurecounts
cd 2.cleandata;ls | sed /name/d > name
cp name ../3_trim;cp name ../4_hisat2;cp name ../5_stringtie;

###check data completeness and move data;
for i in `cat name`;do cd $i; md5sum -c MD5.txt >> ../../alllog 2>&1;mv *.fq.gz ..;cd ..;done

###use FastQC to check data quality
cd 2.cleandata
~/softwares/FastQC/fastqc -t 10 -o ../0_fastqc *.fq.gz;cd ../3_trim;

###If there is some adapters, use Trimmomatic to trim again. If the reads are very clean, this step was skipped
cp ~/softwares/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa .
cd 3_trim
for i in `cat name`;do
echo $i;
java -jar ~/softwares/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 10 -phred33 \
../2.cleandata/${i}_clean_R1.fq.gz ../2.cleandata/${i}_clean_R2.fq.gz -baseout ${i}.fq.gz \
ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:1 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50 >> trim_record 2>&1;
done

###use hisat2 to map the transcriptome to TAIR10 genome 
cd ../4_hisat2;
for i in `cat name`;do 
echo $i >> mapping_record;
hisat2 -p 28 --dta --rg-id $i --rg SM:$i --rna-strandness RF --fr -x ~/genome/tair10seq \
-1 ../3_trim/${i}_1P.fq.gz -2 ../3_trim/${i}_2P.fq.gz -U ../3_trim/${i}_1U.fq.gz,../3_trim/${i}_2U.fq.gz \
-S ${i}.sam >> mapping_record 2>&1;
samtools view -bSF 4 ${i}.sam > ${i}.bam;samtools sort -@ 15 -o ${i}.sorted.bam ${i}.bam;
samtools index ${i}.sorted.bam;
rm ${i}.sam;
done
for i in `cat name`;do rm ${i}.bam;done
for i in `cat name`;do bamCoverage -p 10 -bs 1 -b ${i}.sorted.bam -o ${i}.bw;done

###use stringtie to get the genes FPKM
cd ../5_stringtie;
for i in `cat name`;do 
echo $i;
stringtie ../4_hisat2/${i}.sorted.bam --rf -e -B -A ${i}.gene.abundance -G ~/genome/TAIR10_GFF3_genes_exons.gtf \
-p 10 -o ${i}/${i}.gtf;sort ${i}.gene.abundance -o ${i}.gene.abundance;
done
ls */*.gtf > sample_lst.txt ### edit sample_lst.txt ###%s#\(\S\+\)\/#\1 \1/#
python prepDE.py -i sample_lst.txt

###gene_count_matrix.csv was used for DESeq2 DEG analysis
g <- read.csv("gene_count_matrix.csv",row.names = "gene_id")
colnames(g)
g <- g[,c(1:6,10:12,7:9)]
colnames(g)
colnames(g)<- c("Col-0 rep1","Col-0 rep2","Col-0 rep3","SG Col-0 rep1","SG Col-0 rep2","SG Col-0 rep3","alkbh9b-1 rep1",
"alkbh9b-1 rep2","alkbh9b-1 rep3","SG alkbh9b-1 rep1","SG alkbh9b-1 rep2","SG alkbh9b-1 rep3")
countData <- as.matrix(g)
library("DESeq2")
library("ggplot2")
#rep(c(1,2),3,each=3)
colData <- data.frame("genotype"=rep(c("WT","Mutant"),1,each=6),
                      "condition"=rep(c("Total","SG"),2,each=3))
colData
rownames(colData) <- colnames(countData)
colData
all(rownames(colData) %in% colnames(countData))
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~genotype+condition)
dds$genotype
dds$condition
dds$genotype <- factor(dds$genotype, levels = c("WT","Mutant"))
dds$condition <- factor(dds$condition, levels = c("Total","SG"))
dds$genotype
dds$condition
dds <- DESeq(dds)
#Drawing sample matrix and PCA to check replicates
rld <- rlog(dds)
vsd <- vst(dds)
sampleDists <- dist( t( assay(rld) ) )
#sampleDists <- dist( t( assay(vsd) ) )
sampleDistMatrix <- as.matrix( sampleDists )
library( pheatmap )
library( "RColorBrewer" )
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap( sampleDistMatrix,
          clustering_distance_rows = sampleDists,
          clustering_distance_cols = sampleDists,
          col = colours, border_color = NA )
plotPCA(rld,intgroup=c("condition","genotype"))
