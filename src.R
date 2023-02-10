# Configurations --------------------------------------------------------------
graphics.off() # clear plots
rm(list = ls()) # clear global environment
cat("\014") # clear console
set.seed(1) # seed random number generator

# Set working directory -------------------------------------------------------
setwd("/Users/swafa/Desktop/DGP 486-0/Assignment 1")

# Load libraries --------------------------------------------------------------
library(DESeq2)
library(pheatmap)
library(ggplot2)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(KEGGREST)
library(pcaExplorer)
library(limma)

# Load and process sample information -----------------------------------------
SampleInfo = read.csv("data/GSE120024_series_matrix.csv",
                      stringsAsFactors=FALSE,
                      header=FALSE)
SampleInfo = t(SampleInfo[, 1:15]) # get mRNA but not m6A enriched RNA

# Transform sample information into a data frame
SampleInfo = as.data.frame(data.matrix(SampleInfo))

# Define row and column names
rownames(SampleInfo) = NULL
colnames(SampleInfo) = c("ID", "GeoNum", "InputType", "Group", "Sex", "Batch")

# Process strings
SampleInfo$ID = gsub(" ", "-", fixed=TRUE, SampleInfo$ID)
SampleInfo$Batch = gsub("batch: ", "", fixed=TRUE, SampleInfo$Batch)
SampleInfo$InputType = "mRNA"
SampleInfo$Group = gsub("group: ", "", fixed=TRUE, SampleInfo$Group)
SampleInfo$Group = gsub("type II diabetes", "T2D", fixed=TRUE, SampleInfo$Group)
SampleInfo$Sex = gsub("gender: ", "", fixed=TRUE, SampleInfo$Sex)

# Add normalized age
age = c(45, 47, 60, 65, 68, 26, 57, 58, 57, 52, 47, 52, 52, 49, 65)
SampleInfo$NormAge = (age - mean(age)) / mean(age)

# Set reference levels for factors
SampleInfo$Group = factor(SampleInfo$Group, levels=c("control", "T2D"))
SampleInfo$Sex = factor(SampleInfo$Sex, levels=c("M", "F"))
SampleInfo$Batch = factor(SampleInfo$Batch, levels=c("A", "B", "C"))

# Load and process gene counts ------------------------------------------------
countData = read.table("data/refGene_counts.txt", 
                       sep="\t",
                       stringsAsFactors=FALSE)

# Remove Geneid
geneID = countData$V1
geneID = geneID[-1]

# Remove the first 6 columns and column names
colName = countData[1, ]
colName = colName[-c(1, 2, 3, 4, 5, 6)]
countData = countData[-1, -c(1, 2, 3, 4, 5, 6)]

# Transform gene counts into a data frame
countData = as.data.frame(data.matrix(countData))

# Define row names
rownames(countData) = geneID

# Process column name strings
colName = gsub("1Bam/", "", fixed=TRUE, colName)
colName = gsub(".bam", "-input", fixed=TRUE, colName)
colName = gsub("l0", "l", fixed=TRUE, colName)
colName = gsub("D0", "D", fixed=TRUE, colName)
colName = gsub("Ctrl", "Ctl", fixed=TRUE, colName)

# Define row and column names
colnames(countData) = colName

# Match SampleInfo and countData ----------------------------------------------

length(unique(intersect(colName, SampleInfo$ID)))
countData = countData[, match(SampleInfo$ID, colnames(countData))]
identical(SampleInfo$ID, colnames(countData))

# Data filtering --------------------------------------------------------------

countData = countData[rowSums(countData) > 0, ] # remove genes with 0 count

# Construct DESeqDataSet ------------------------------------------------------

# Initialize variables
cts = countData
coldata = SampleInfo
myDesign = ~ Batch + Sex + NormAge + Group

# Check cts and coldata are in the same order
identical(coldata$ID, colnames(cts))

# Construct a DESeqDataSet
dds = DESeqDataSetFromMatrix(countData=cts, colData=coldata, design=myDesign)

# Filter genes with 10+ count
keep = rowSums(counts(dds)) >= 10
dds = dds[keep, ]

# Perform differential expression analysis ------------------------------------
dds = DESeq(dds)
res = results(dds)
summary(res)
res = as.data.frame(res)
res = na.omit(res)
resOrdered = res[order(res$padj), ]
head(resOrdered)

# Data transformation ---------------------------------------------------------
vsd = vst(dds, blind=FALSE)
mat = assay(vsd)
mm = model.matrix(~Group, colData(vsd))
mat = limma::removeBatchEffect(mat, batch=vsd$Batch, design=mm)
assay(vsd) = mat

#vsd = vst(dds, blind=FALSE)
#rld = rlog(dds, blind=FALSE)
#ntd = normTransform(dds)
#normCount = counts(dds, normalized=TRUE)

# Visualization and outlier detection -----------------------------------------
p0 = plotPCA(vsd, intgroup="Batch") + 
  labs(col="Batch") + 
  theme(legend.position = "bottom") + 
  geom_point(size=5) +
  theme(text=element_text(size=20))
png("results/pca_batch.png", height = 8, width = 8, units = "in", res = 1200)
p0
dev.off()

p1 = plotPCA(vsd, intgroup="Sex") +
  labs(col="Sex") +
  theme(legend.position = "bottom") +
  geom_point(size=5) +
  theme(text=element_text(size=20))
png("results/pca_sex.png", height = 8, width = 8, units = "in", res = 1200)
p1
dev.off()

p2 = plotPCA(vsd, intgroup="NormAge") + 
  labs(col="NormAge") + 
  geom_point(size=5) +
  theme(text=element_text(size=20))
png("results/pca_age.png", height = 8, width = 8, units = "in", res = 1200)
p2
dev.off()

p3 = plotPCA(vsd, intgroup="Group") +
  labs(col="Group") + 
  theme(legend.position = "bottom") +
  geom_point(size=5) +
  theme(text=element_text(size=20))
png("results/pca_group.png", height = 8, width = 8, units = "in", res = 1200)
p3
dev.off()

df = as.data.frame(colData(dds)[, c("Group", "Batch", "Sex", "NormAge")])
p4 = pheatmap(assay(vsd), 
              cluster_rows=TRUE, 
              show_rownames=FALSE,
              cluster_cols=TRUE,
              annotation_col=df,
              fontsize=11.5)
png("results/full_heatmap.png", height = 12, width = 8, units = "in", res = 1200)
p4
dev.off()

select = rownames(res[order(res$padj, decreasing=FALSE), ])[1:10]
df = as.data.frame(colData(dds)[, c("Group", "Batch", "Sex", "NormAge")])
p5 = pheatmap(assay(vsd[select, ]), 
         cluster_rows=TRUE, 
         show_rownames=TRUE,
         cluster_cols=TRUE,
         annotation_col=df,
         fontsize=12.5)
png("results/de_heatmap.png", height = 12, width = 8, units = "in", res = 1200)
p5
dev.off()

# KEGG pathway analysis -------------------------------------------------------

search_kegg_organism('hsa', by='kegg_code')
homo = search_kegg_organism('Homo sapiens', by='scientific_name')

candi_gene = geneID[res$padj<0.05]
gene.id = bitr(candi_gene, fromType="SYMBOL",
                toType=c("ENTREZID"),
                OrgDb=org.Hs.eg.db)
gene_input = gene.id$ENTREZID

write.csv(gene_input, "genes.csv", row.names=T)

kk = enrichKEGG(gene=gene_input,
                 organism='hsa',
                 pvalueCutoff=0.05,
                 minGSSize=15,
                 maxGSSize=500,
                 use_internal_data=F)

df_kegg = kk@result[1:sum(kk@result$pvalue < 0.05), ]

df_kegg$GeneRatio = sapply(df_kegg$GeneRatio, function(x) eval(parse(text = x)))

p6 = ggplot(df_kegg, aes(x=GeneRatio, y=Description)) + 
  geom_point(aes(size=GeneRatio, color=pvalue)) +
  theme_bw(base_size=14) +
  scale_y_discrete(limits=rev) +
  scale_colour_gradient(limits=c(0, 0.05), low="blue", high="red") +
  ylab(NULL) +
  ggtitle("KEGG pathway analysis") +
  theme(text=element_text(size=20))
png("results/kegg.png", height = 12, width = 16, units = "in", res = 1200)
p6
dev.off()

candi_gene = geneID[res$pvalue<0.05]
gene.id = bitr(candi_gene, fromType="SYMBOL",
               toType=c("ENTREZID"),
               OrgDb=org.Hs.eg.db)
gene_input = gene.id$ENTREZID

ngenes = dim(gene.id)[1]
gene_idx = matrix(NA)
gene_eids = matrix(NA)
for (i in 1:ngenes){
  temp = gene.id$SYMBOL[i]
  if (any(temp==rownames(res))){
    gene_idx = rbind(gene_idx, which(temp == rownames(res)))
    gene_eids = rbind(gene_eids, gene.id$ENTREZID[i])
  }
}
gene_idx = na.omit(gene_idx)
gene_eids = na.omit(gene_eids)

original_gene_list = res$log2FoldChange[gene_idx]
names(original_gene_list) = gene_eids
gene_list = na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing=TRUE)
gse <- gseGO(geneList=gene_list, 
             ont="ALL",
             keyType="ENTREZID",
             minGSSize=15,
             maxGSSize=500,
             pvalueCutoff=1,
             OrgDb=org.Hs.eg.db)

p7 = dotplot(gse, split=".sign", showCategory=10, title="GO enrichment analysis")
#p7 = dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
png("results/go.png", height = 12, width = 16, units = "in", res = 1200)
p7
dev.off()
