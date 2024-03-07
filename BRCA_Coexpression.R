setwd("E:\\biyesheji_data_set\\brca")
library(biomaRt)
library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads() 

WGCNA_matrix <- read.table("Merge_RNA_seq_FPKM.txt",header=T,
                           comment.char = "",
                           check.names=F)
TCGA_BRCA_RNA_seq= as.data.frame(t(WGCNA_matrix[,-1]))
names(TCGA_BRCA_RNA_seq) = WGCNA_matrix$Tag
rownames(TCGA_BRCA_RNA_seq) = names(WGCNA_matrix[,-1])
TCGA_BRCA_Clinical <- read.table("Clinical.BCR Biotab")
col1_clinical <- TCGA_BRCA_Clinical[,1]
rowname_clinical <- col1_clinical[-1]
row1_clinical <- TCGA_BRCA_Clinical[1,]
colname_clinical <- row1_clinical[-1]
datclinical <- TCGA_BRCA_Clinical[-1,-1]
names(datclinical) <- colname_clinical
rownames(datclinical) <- rowname_clinical
library(org.Hs.eg.db)

# 获取基因矩阵的行名（假设为 gene_ids）
gene_ids <- names(TCGA_BRCA_RNA_seq)
# 执行转换
gene_symbols <- mapIds(org.Hs.eg.db, keys = gene_ids, column = "SYMBOL", keytype = "ENSEMBL")

# 将转换后的基因符号替换为基因矩阵的行名
names(TCGA_BRCA_RNA_seq) <- gene_symbols

for (i in nrow(datclinical):1)
{
  if(datclinical[i,][1]!="Positive" && datclinical[i,][1]!="Negative" || datclinical[i,][2]!="Positive" && datclinical[i,][2]!="Negative" || datclinical[i,][3]!="Positive" && datclinical[i,][3]!="Negative")
  {
    datclinical <- datclinical[-i, ]
  }
}

set.seed(123)

# 从矩阵中随机抽取100个样本
sample_indices <- sample(1:nrow(datclinical), 100, replace = FALSE)

# 保留抽取的100个样本，删除其他多余的样本
datclinical <- datclinical[sample_indices, ]
TCGA_BRCA_RNA_seq = TCGA_BRCA_RNA_seq[match(row.names(datclinical), gsub("-01$", "", rownames(TCGA_BRCA_RNA_seq))),]
TCGA_BRCA_RNA_seq <- na.omit(TCGA_BRCA_RNA_seq)
datclinical = datclinical[match(row.names(TCGA_BRCA_RNA_seq), paste0(rownames(datclinical), "-01")),]
datclinical <- na.omit(datclinical)
TCGA_BRCA_RNA_seq = t(TCGA_BRCA_RNA_seq)
m.vars=apply(TCGA_BRCA_RNA_seq,1,var)
TCGA_BRCA_RNA_seq.upper = TCGA_BRCA_RNA_seq[which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.25))[4]),]
TCGA_BRCA_RNA_seq2 = as.data.frame(t(TCGA_BRCA_RNA_seq.upper))

datclinical=datclinical[match(row.names(TCGA_BRCA_RNA_seq2),paste0(row.names(datclinical),'-01')),]
#trainDt=as.matrix(cbind(ifelse(datclinical[,1]=='Positive',0,1),#将阴性的样本标记为1
#                        ifelse(datclinical[,2]=='Positive',0,1),#将阴性的样本标记为1
#                        ifelse(datclinical[,3]=='Positive',0,1),#将阴性的样本标记为1
#                        ifelse(datclinical[,1]=='Positive'&datclinical[,2]=='Positive'&datclinical[,3]=='Negative',0,ifelse(datclinical[,1]=='Positive'&datclinical[,2]=='Positive'&datclinical[,3]=='Positive',1,ifelse(datclinical[,1]=='Negative'&datclinical[,2]=='Negative'&datclinical[,3]=='Positive',2,ifelse(datclinical[,1]=='Negative'&datclinical[,2]=='Negative'&datclinical[,3]=='Negative',3,NA))))))

n <- nrow(datclinical)
trainDt <- matrix(NA, nrow = n, ncol = 4)

for (i in 1:n) {
  trainDt[i, 1] <- ifelse(datclinical[i, 1] == 'Positive', 0, 1)
  trainDt[i, 2] <- ifelse(datclinical[i, 2] == 'Positive', 0, 1)
  trainDt[i, 3] <- ifelse(datclinical[i, 3] == 'Positive', 0, 1)
  
  if (datclinical[i, 1] == 'Positive' && datclinical[i, 2] == 'Positive') {
    if (datclinical[i, 3] == 'Negative') {
      trainDt[i, 4] <- 0
    } else if (datclinical[i, 3] == 'Positive') {
      trainDt[i, 4] <- 1
    }
  } else if (datclinical[i, 1] == 'Negative' && datclinical[i, 2] == 'Negative') {
    if (datclinical[i, 3] == 'Positive') {
      trainDt[i, 4] <- 2
    } else if (datclinical[i, 3] == 'Negative') {
      trainDt[i, 4] <- 3
    }
  }
}
column_names = c("ER_NEGA", "PR_NEGA", "HER2_NEGA","Category")
colnames(trainDt) = column_names
row.names(trainDt) = row.names(datclinical)
trainDt <- na.omit(trainDt)
TCGA_BRCA_RNA_seq2 = TCGA_BRCA_RNA_seq2[match(row.names(trainDt), gsub("-01$", "", rownames(TCGA_BRCA_RNA_seq2))),]
gsg = goodSamplesGenes(TCGA_BRCA_RNA_seq2, verbose = 3)
gsg$allOK
sampleTree = hclust(dist(TCGA_BRCA_RNA_seq2), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers"
     , sub="", xlab="")
clust = cutreeStatic(sampleTree, cutHeight = 80000, minSize = 10)
table(clust)
keepSamples = (clust==1)
TCGA_BRCA_RNA_seq2 = TCGA_BRCA_RNA_seq2[keepSamples, ]
nGenes = ncol(TCGA_BRCA_RNA_seq2)
nSamples = nrow(TCGA_BRCA_RNA_seq2)
save(datExpr2, file = "FPKM-01-dataInput.RData")
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(TCGA_BRCA_RNA_seq2, powerVector = powers, verbose = 5)
##画图##
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

pow=4

cor <- WGCNA::cor
net = blockwiseModules(TCGA_BRCA_RNA_seq2, power = pow, maxBlockSize = 7000,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "FPKM-TOM",
                       verbose = 3)
cor<-stats::cor
table(net$colors)
# open a graphics window
#sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    groupLabels = c("Module colors", 
                                    "GS.weight"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


moduleLabelsAutomatic = net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
moduleColorsFemale = moduleColorsAutomatic
MEs0 = moduleEigengenes(TCGA_BRCA_RNA_seq2, moduleColorsFemale)$eigengenes
MEsFemale = orderMEs(MEs0)

MEsFemale = MEsFemale[match(row.names(trainDt), gsub("-01$", "", rownames(MEsFemale))),]
# 检查trainDt是否还包含NA
sum(is.na(trainDt))
modTraitCor = cor(MEsFemale, trainDt, use = "p")
modTraitP = corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

colnames(modTraitCor) = column_names
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(trainDt), yLabels = names(MEsFemale), 
               ySymbols = colnames(moduleColorsFemale), colorLabels = FALSE, colors = greenWhiteRed(50), 
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1)
               , main = paste("Module-trait relationships"))

TCGA_BRCA_RNA_seq2 = TCGA_BRCA_RNA_seq2[match(row.names(trainDt), gsub("-01$", "", rownames(TCGA_BRCA_RNA_seq2))),]
modTraitCor = cor(MEsFemale, TCGA_BRCA_RNA_seq2, use = "p")
modTraitP = corPvalueStudent(modTraitCor, nSamples)
MEyellow=modTraitCor[which(row.names(modTraitCor)=='MEyellow'),]
head(MEyellow[order(-MEyellow)])

TOM = TOMsimilarityFromExpr(TCGA_BRCA_RNA_seq2, power = pow);
probes = names(TCGA_BRCA_RNA_seq2)
mc='yellow'
mcInds=which(match(moduleColorsAutomatic, gsub('^ME','',mc))==1)
modProbes=probes[mcInds]
modTOM = TOM[mcInds, mcInds];

dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("edges-", mc, ".txt", sep=""),
                               nodeFile = paste("nodes-", mc, ".txt", sep=""),
                               weighted = TRUE,
                               threshold = median(modTOM),
                               nodeNames = modProbes,
                               #altNodeNames = modGenes,
                               nodeAttr = moduleColorsAutomatic[mcInds])
