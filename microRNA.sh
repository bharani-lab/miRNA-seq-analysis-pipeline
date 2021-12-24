#Sequence_filter_below16-above30nt :
for prefix in $(ls *.fastq | uniq)
do
/home/kathirvel/Documents/Glaucoma/RNA/Tools/bbmap/reformat.sh ${prefix} out=${prefix%.fastq}.fastq  minlen=16 maxlen=30
done
#Mapping:
for prefix in $(ls *.fastq | uniq)
do
STAR --runThreadN 16 --runMode alignReads --genomeDir /home/kathirvel/Documents/Glaucoma/miRNA/STAR/Indexed_Reference/ --readFilesIn ${prefix} --outFileNamePrefix ${prefix%.fastq} --outSAMstrandField All --outSAMtype BAM SortedByCoordinate
done
#Quantification
for prefix in $(ls *.bam | uniq)
do
featureCounts -t miRNA -g Name -O -s 1 -M -a /home/kathirvel/Desktop/Cat_TB/hsa.gff3 -o ${prefix%.bam}.csv ${prefix}
done
#Move files
mkdir Mapping Quantification Normalisation DiffExp
mv *.bam *.out *.tab Mapping
mv *.csv Quantification
cd Quantification
#post Quantification filter >5 
	#remove unwanted coloumns 
for prefix in $(ls *sortedByCoord.out.csv | uniq); do awk 'NR>1 {print $1,$7}' ${prefix} > ${prefix}.csv; done
rm *sortedByCoord.out.csv
	#Combine all csv files into one
paste *Aligned.sortedByCoord.out.csv.csv > Counts.csv 
awk '{print $1,$2,$4,$6,$8,$10,$12,$14,$16,$18,$20,$22,$24,$26,$28,$30,$32,$34,$36,$38,$40,$42,$44,$46,$48,$50,$52,$54,$56,$58,$60}' Counts.csv > Counts_1.csv
rm Counts.csv
	#Remove miRNAs <5 counts
awk '$2>5 || $3>5 || $4>5 || $5>5 || $6>5 || $7>5 || $8>5 || $9>5 || $10>5 || $11>5 || $12>5 || $13>5 || $14>5 || $15>5 || $16>5 || $17>5 || $18>5 || $19>5 || $20>5 || $21>5 {print}' Counts_1.csv > Counts.csv
cd Normalisation
#DE Analysis
#Differential expression analysis
#edgeR + Quantile + FeatureCounts
R
setwd("/home/kathirvel/Normalisation/")
library(readr)
Counts <- read_csv("~Normalisation/Counts_>5.csv")
GeneNames <- Counts$Geneid
Counts <- Counts[, c(-1)]
row.names(Counts) <- Names
#Condition file creation
file <- c('Etoh_Counts1','Etoh_Counts2','Etoh_Counts3','Etoh_Counts4','Dex_Counts1','Dex_Counts2','Dex_Counts3','Dex_Counts4')
condition <- c('1EtoH','1EtoH','1EtoH','1EtoH','2Dex','2Dex','2Dex','2Dex')
Info <- data.frame(file, condition)
library(edgeR)
library(HTSFilter)
d <-Counts
d <- DGEList(counts = d, group = Info$condition)
TMM <- calcNormFactors(d, method="TMM")
d <- estimateCommonDisp(TMM)
d <- estimateTagwiseDisp(d)
et <- exactTest(d)
p <- et$table$PValue
adj.p <- p.adjust(p, method = "BH")
res <- cbind(id=rownames(et$table), et$table, adj.p, threshold=-p)
Down <- res[res$PValue<0.05,]
DEG <- Down[order(Down$logFC), ]
write.csv(DEG, "/home/kathirvel/DiffExp/DEG_NR.csv")
#Make a basic volcano plot
with(res, plot(logFC, -log10(PValue), pch=20, xlim=c(-6,10)))
with(subset(res,PValue<.05 & abs(logFC)>2), points(logFC, -log10(PValue), pch=20, col="red"))
library(calibrate)
with(subset(res,adj.p<.05 & abs(logFC)>2), textxy(logFC, -log10(PValue), labs=id, cex=.6))
# Add colored points: red if padj<0.05, orange of log2FC>2, green if both)
library(EnhancedVolcano)
res <- res[, c(-1)]
EnhancedVolcano(res, lab = rownames(res), x = 'logFC', y = 'PValue', xlim = c(-7, 7), ylim = c(0, 25), FCcutoff = 2, pointSize = 2, labSize = 3)




	

















