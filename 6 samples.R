1- Downloading Data using SRA-toolkit
# create the project directory
mkdir -p bulk_rna_raw
cd bulk_rna_raw

# download SRA files (Single_end)
prefetch SRR3617193 SRR3617194 SRR3617201 SRR3617202 SRR3617206 SRR3617207

# Convert to FASTQ
fasterq-dump SRR3617193 
fasterq-dump SRR3617194 
fasterq-dump  SRR3617201
fasterq-dump  SRR3617202
fasterq-dump SRR3617206 
fasterq-dump SRR3617207

2- Quality Control using Fastqc
# create the output directory for qc files
mkdir -p qc
# run fastqc to all fastq files
fastqc -t 4 -o qc *.fastq
# -t 4 use 4 cpu threads, -o to get output reports in qc directory
# make sure you are in qc directory
cd qc


3- Trimming using trimmomatic
# 1. Create a folder to keep things organized
mkdir -p trimmed_fastq

# 2. Run the loop for your specific SRR IDs
for id in SRR3617193 SRR3617194 SRR3617201 SRR3617202 SRR3617206 SRR3617207; do
  echo "Trimming sample: $id"
  trimmomatic SE -threads 4 \
  "${id}.sra.fastq" \
  "trimmed_fastq/${id}_trimmed.fastq.gz" \
  ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done

4- Quality control for the trimming files
fastqc -t 4 -o qct *fastq.gz*
  
5- Indexing & Alignment
#create durectory for genome index
mkdir -p genome_index
# Download the Primary Assembly
wget https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Unzip it (HISAT2-build needs it unzipped)
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Rename it to something simple
mv Homo_sapiens.GRCh38.dna.primary_assembly.fa genome.fasta

#Build the Genome Index (Do this once)
hisat2-build -p 8 genome.fasta genome

# Create directory for results
mkdir -p mapped_bam

# Align, Sort, and Index in a single "pipe" to save disk space
hisat2 -p 8 -x genome -U "$f" | samtools sort -@ 4 -o "mapped_bam/${id}.bam"
    
    samtools index "mapped_bam/${id}.bam"
    
    echo "Sample $id is finished!"
done

6-Generate the gene count matrix for all 6 samples 

# Preparation: The Annotation File (GTF)
wget https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz

gunzip Homo_sapiens.GRCh38.110.gtf.gz

7- Quantification
featureCounts -T 8 \
    -a Homo_sapiens.GRCh38.110.gtf \
    -o counts_matrix.txt \
    -g gene_id \
    -t exon \
    *.bam
    
# Cleaning the Matrix for R/DESeq2

cut -f 1,7,8,9,10,11,12 counts_matrix.txt > final_counts_for_R.txt

sed -i 's/mapped_bam\///g; s/.bam//g' final_counts.txt

tail -n +2 final_counts.txt | tr '\t' ',' > final_counts.csv

8- PCA, DESeq analysis, Gene Ontology (GO) and Kegg

library(dplyr)
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(EnhancedVolcano)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GEOquery)
metadata <- read.csv("SraRunTable.csv", row.names = 1)
data <- read.csv("final_counts.csv", row.names = 1)
metadata_filtered <- metadata[, c(16, 29, 34)]
all(colnames(data)==row.names(metadata_filtered))
data_round <- round(data)
metadata_filtered <- metadata_filtered[!is.na(metadata_filtered$health_state), ]
metadata_filtered$health_state <- factor(metadata_filtered$health_state)

dds <- DESeqDataSetFromMatrix(countData = data_round,
                              colData = metadata_filtered,
                              design = ~health_state)

#filtering dds
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep, ]
dds <- DESeq(dds)
vsd <- vst(dds, blind = TRUE)
pcaData <- plotPCA(vsd, intgroup = "health_state", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(x = PC1, y = PC2, color = health_state)) + 
  geom_point(size = 3) +                                          
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal() +
  ggtitle("PCA of RNA-Seq Samples (VST Transformed)")

res <- lfcShrink(dds, coef = 2, type = "apeglm")
head(res[order(res$pvalue), ])
summary(res)

plotMA(res, ylim=c(-5,5))

res_sig <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1.0), ] 
res_sig <- na.omit(res_sig)

# 1. Convert the DESeqResults object to a standard data frame

library(EnhancedVolcano)
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1,
                title = 'Disease_Compare')
top_genes <- head(rownames(res_sig[order(res_sig$padj), ]), 50)

vsd <- vst(dds, blind=TRUE)

df <- as.data.frame(colData(dds)[, c("health_state")])
rownames(df) <- colnames(vsd)
colnames(df) <- "health_state"


pheatmap(assay(vsd)[top_genes, ], 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = TRUE, 
         annotation_col = df,  # Using the 'df' we just created
         scale = "row",
         main = "DEG Expression across Samples")

library(org.Hs.eg.db)
library(clusterProfiler) # You'll need this shortly

entrez_ids <- mapIds(org.Hs.eg.db,
                     keys = top_genes,
                     column = "ENTREZID",
                     keytype = "ENSEMBL",
                     multiVals = "first")
clean_entrez_ids <- na.omit(entrez_ids) 
clean_entrez_ids <- unique(clean_entrez_ids) 
ego <- enrichGO(gene = clean_entrez_ids,
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                readable = TRUE)
dotplot(ego)

## To color genes by log2 fold changes, we need to extract the log2 fold changes from our results table creating a named vector
OE_foldchanges <- res_sig$log2FoldChange
names(OE_foldchanges) <- res_sig$gene

#The cnetplot.png provides a map of which specific genes are driving the enrichment of these pathways.
cnetplot(ego, 
         categorySize="pvalue", 
         showCategory = 5, 
         foldChange=OE_foldchanges, 
         vertex.label.font=6)

#To perform GSEA analysis of KEGG gene sets,

# Rerun KEGG analysis with pvalueCutoff = 1.0 using the CLEAN list
kegg_results_all <- enrichKEGG(gene = clean_entrez_ids,
                               organism = 'hsa', 
                               pvalueCutoff = 1.0) # Set to 1.0
# Print the top pathways to see if they are relevant
print(head(kegg_results_all))

# Plot the top 10 enriched pathways
dotplot(kegg_results_all, 
        showCategory = 10, 
        title = "Top 10 Enriched KEGG Pathways (p-value < 1.0)")