####lingwang@cemps.ac.cn

###mapping m6A transcriptome data
for i in `cat sample`;do
echo $i;
hisat2 -p 12 --dta --rg-id $i --rg SM:$i --rna-strandness RF --fr -x ~/genome/tair10seq -1 ../clean_data/${i}_1P.gz -2 ../clean_data/${i}_2P.gz -S ${i}.sam;rm ${i}.sam;
samtools view -bS ${i}.sam > ${i}.bam;samtools sort -@ 10 ${i}.bam ${i}.sorted;samtools index ${i}.sorted.bam
done

###use MACS2 call m6A peak
macs2 callpeak -t ../hisat2/w0f2_ip.bam -c ../hisat2/w0f2_input.bam -f BAM -g 65084214 -p 5e-2 --nomodel --extsize 50 -n w0f2_pe5_ex50
macs2 callpeak -t ../hisat2/w0f1_ip.bam -c ../hisat2/w0f1_input.bam -f BAM -g 65084214 -p 5e-2 --nomodel --extsize 50 -n w0f1_pe5_ex50
macs2 callpeak -t ../hisat2/w3f2_ip.bam -c ../hisat2/w3f2_input.bam -f BAM -g 65084214 -p 5e-2 --nomodel --extsize 50 -n w3f2_pe5_ex50
macs2 callpeak -t ../hisat2/w3f1_ip.bam -c ../hisat2/w3f1_input.bam -f BAM -g 65084214 -p 5e-2 --nomodel --extsize 50 -n w3f1_pe5_ex50
macs2 callpeak -t ../hisat2/w0l2_ip.bam -c ../hisat2/w0l2_input.bam -f BAM -g 65084214 -p 5e-2 --nomodel --extsize 50 -n w0l2_pe5_ex50
macs2 callpeak -t ../hisat2/w0l1_ip.bam -c ../hisat2/w0l1_input.bam -f BAM -g 65084214 -p 5e-2 --nomodel --extsize 50 -n w0l1_pe5_ex50
macs2 callpeak -t ../hisat2/w3l2_ip.bam -c ../hisat2/w3l2_input.bam -f BAM -g 65084214 -p 5e-2 --nomodel --extsize 50 -n w3l2_pe5_ex50
macs2 callpeak -t ../hisat2/w3l1_ip.bam -c ../hisat2/w3l1_input.bam -f BAM -g 65084214 -p 5e-2 --nomodel --extsize 50 -n w3l1_pe5_ex50
macs2 callpeak -t ../hisat2/m0f2_ip.bam -c ../hisat2/m0f2_input.bam -f BAM -g 65084214 -p 5e-2 --nomodel --extsize 50 -n m0f2_pe5_ex50
macs2 callpeak -t ../hisat2/m0f1_ip.bam -c ../hisat2/m0f1_input.bam -f BAM -g 65084214 -p 5e-2 --nomodel --extsize 50 -n m0f1_pe5_ex50
macs2 callpeak -t ../hisat2/m0l2_ip.bam -c ../hisat2/m0l2_input.bam -f BAM -g 65084214 -p 5e-2 --nomodel --extsize 50 -n m0l2_pe5_ex50
macs2 callpeak -t ../hisat2/m0l1_ip.bam -c ../hisat2/m0l1_input.bam -f BAM -g 65084214 -p 5e-2 --nomodel --extsize 50 -n m0l1_pe5_ex50
macs2 callpeak -t ../hisat2/m3f2_ip.bam -c ../hisat2/m3f2_input.bam -f BAM -g 65084214 -p 5e-2 --nomodel --extsize 50 -n m3f2_pe5_ex50
macs2 callpeak -t ../hisat2/m3f1_ip.bam -c ../hisat2/m3f1_input.bam -f BAM -g 65084214 -p 5e-2 --nomodel --extsize 50 -n m3f1_pe5_ex50
macs2 callpeak -t ../hisat2/m3l2_ip.bam -c ../hisat2/m3l2_input.bam -f BAM -g 65084214 -p 5e-2 --nomodel --extsize 50 -n m3l2_pe5_ex50
macs2 callpeak -t ../hisat2/m3l1_ip.bam -c ../hisat2/m3l1_input.bam -f BAM -g 65084214 -p 5e-2 --nomodel --extsize 50 -n m3l1_pe5_ex50

###m6A peaks from two replicates were considered valid
intersectBed -a w0f1_pe5_ex50_peaks.narrowPeak -b w0f2_pe5_ex50_peaks.narrowPeak > w0f.bed
intersectBed -a w3f1_pe5_ex50_peaks.narrowPeak -b w3f2_pe5_ex50_peaks.narrowPeak > w3f.bed
intersectBed -a w3l1_pe5_ex50_peaks.narrowPeak -b w3l2_pe5_ex50_peaks.narrowPeak > w3l.bed
intersectBed -a w0l1_pe5_ex50_peaks.narrowPeak -b w0l2_pe5_ex50_peaks.narrowPeak > w0l.bed
intersectBed -a m0f1_pe5_ex50_peaks.narrowPeak -b m0f2_pe5_ex50_peaks.narrowPeak > m0f.bed
intersectBed -a m3f1_pe5_ex50_peaks.narrowPeak -b m3f2_pe5_ex50_peaks.narrowPeak > m3f.bed
intersectBed -a m3l1_pe5_ex50_peaks.narrowPeak -b m3l2_pe5_ex50_peaks.narrowPeak > m3l.bed
intersectBed -a m0l1_pe5_ex50_peaks.narrowPeak -b m0l2_pe5_ex50_peaks.narrowPeak > m0l.bed

###use R package "ChIPseeker" to annoate peaks to genes
BiocManager::install("ChIPseeker")
library("ChIPseeker")
library("GenomicFeatures")
Txdb_gtf <- makeTxDbFromGFF("TAIR10_GFF3_genes_exons.gtf")
seqlevels(Txdb_gtf)
files <- list.files(".", pattern= "bed", full.names=T)
peakAnnoList <- lapply(files, annotatePeak,
                       TxDb=Txdb_gtf,tssRegion = c(-500,500),overlap="all",
                       genomicAnnotationPriority=c("5UTR","3UTR","Exon","Intron","Promoter","Downstream","Intergenic"))
names(peakAnnoList)<- c("Mutant0hflower","Mutant0hleaf","Mutant3hflower","Mutant3hleaf","WT0hflower","WT0hleaf","WT3hflower","WT3hleaf")
plotAnnoBar(peakAnnoList)
write.table(as.data.frame(peakAnnoList$Mutant0hflower@anno),"Mutant_0h_flower_peak_annotation.txt",sep="\t",quote=F,row.names = F)
write.table(as.data.frame(peakAnnoList$Mutant3hflower@anno),"Mutant_3h_flower_peak_annotation.txt",sep="\t",quote=F,row.names = F)
write.table(as.data.frame(peakAnnoList$Mutant0hleaf@anno),"Mutant_0h_leaf_peak_annotation.txt",sep="\t",quote=F,row.names = F)
write.table(as.data.frame(peakAnnoList$Mutant3hleaf@anno),"Mutant_3h_leaf_peak_annotation.txt",sep="\t",quote=F,row.names = F)
write.table(as.data.frame(peakAnnoList$WT0hflower@anno),"WT_0h_flower_peak_annotation.txt",sep="\t",quote=F,row.names = F)
write.table(as.data.frame(peakAnnoList$WT3hflower@anno),"WT_3h_flower_peak_annotation.txt",sep="\t",quote=F,row.names = F)
write.table(as.data.frame(peakAnnoList$WT0hleaf@anno),"WT_0h_leaf_peak_annotation.txt",sep="\t",quote=F,row.names = F)
write.table(as.data.frame(peakAnnoList$WT3hleaf@anno),"WT_3h_leaf_peak_annotation.txt",sep="\t",quote=F,row.names = F)

###use R package "Guitar" to check m6A peaks along mRNA
BiocManager::install("Guitar")
library("Guitar")
library("GenomicFeatures")
Txdb_gtf <- makeTxDbFromGFF("TAIR10_GFF3_genes_exons.gtf")
files <- list.files(".", pattern= "bed", full.names=T)
extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")
import.narrowPeak <- function(...) {
  import(..., format="BED", extraCols=extraCols_narrowPeak)
}
#WT flower leaf
fileswt<-files[5:8]
stGRangelist <- list()
for (i in 1:length(fileswt)){
  stGRangelist[[i]] <- blocks(import.narrowPeak(fileswt[[i]]))
}
p <- GuitarPlot(txTxdb = Txdb_gtf,
           stGRangeLists = stGRangelist,
           pltTxType="mrna",
           enableCI=FALSE,
           headOrtail = FALSE,
           stGroupName=c("Flower CS","Leaf CS","Flower HS","Leaf HS")
           )
colors=c("red","pink","black","grey")
p+ theme_classic()+scale_color_manual(values=colors)+
  theme(axis.text=element_text(size=10,color="black"))+
  xlab( label = "")
#leaf
stGRangelist <- list()
files <- list.files(".", pattern= "bed", full.names=T)
filesleaf<-files[c(2,6)]
filesleaf
for (i in 1:length(filesleaf)){
  stGRangelist[[i]] <- blocks(import.narrowPeak(filesleaf[[i]]))
}
pl <- GuitarPlot(txTxdb = Txdb_gtf,
                stGRangeLists = stGRangelist,
                pltTxType="mrna",
                enableCI=FALSE,
                headOrtail = FALSE,
                stGroupName=c("10b-1 CS","Col-0 CS")
)
colors=c("blue","black")
pl+ theme_classic()+scale_color_manual(values=colors)+
  theme(axis.text=element_text(size=10,color="black"))+
  xlab( label = "")
  









