#!/usr/bin/env  Rscript
library(docopt)
"Usage: deseq2.r  -i <file> -o <dir> -s <file> [--Rlib <dir>]
Options:
    -i, --input <file>              输入文件，样本reads count原始表达矩阵。 第一列是基因、转录本的ID，后面的所有列可以是样本，也可以是注释。 样本的表达量不能有缺失，缺失应该设为0。矩阵分隔符是'\\t'
    -o, --output <dir>              结果输出目录
    -s, --samplefile <file>	    样本分组数据，第一列为样本名，第二列为分组，后面不限 。
    --Rlib <dir>                    R包路径  [default: /home/genesky/software/r/3.5.1/lib64/R/library] " -> doc

opts                     <- docopt(doc, version = '对样本进行聚类 \n          李澎鹏\n')
input                    <- opts$input
output_dir               <- opts$output
samplefile               <- opts$samplefile
Rlib                     <- opts$Rlib
.libPaths(Rlib)

# 部分信息提前准备

if (!file.exists(output_dir)) dir.create(output_dir,recursive=TRUE)

library("DESeq2")
library("BiocParallel")
require(ggplot2)
require(ggthemes)
library("RColorBrewer")
library('pheatmap')
library("corrplot")
register(MulticoreParam(32))
condFile = samplefile
coldata <- read.delim( condFile, header=TRUE,stringsAsFactors=TRUE )


exampleFile = input
cts <- read.delim( exampleFile, header=TRUE, stringsAsFactors=TRUE, row.names=1  ) 

rownames(coldata)<- coldata$Sample
anno_col<- subset(coldata, select = c(Condition))
cts <- read.delim( exampleFile, header=TRUE, stringsAsFactors=TRUE, row.names=1  ) 
anno_row<- subset(coldata, select = c(Condition))
need<- colnames(cts)%in%rownames(coldata)

cts<- cts[,need]
need2 <-rownames(coldata)%in%colnames(cts)
coldata <- coldata[need2,]
cts<- cts[,rownames(coldata)]
summary(cts)
#输出所有表达量

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ Condition)
dds$Condition <- factor(dds$Condition)

dds1 <- DESeq(dds)
baseMeanPerLvl <- sapply( levels(dds1$Condition), function(lvl) rowMeans( counts(dds1,normalized=TRUE)[,dds1$Condition == lvl, drop=F] ) )
h<- as.data.frame(baseMeanPerLvl)
data=cbind(GeneID=row.names(h), h)
write.table( data, row.names=FALSE,file=paste0(output_dir,"/EachSampleFPKM.tsv"),quote=FALSE,sep="\t")
Data = read.table(paste0(output_dir,"/EachSampleFPKM.tsv"),header=TRUE, stringsAsFactors=TRUE  )

rownames(Data)<-Data[,1]
Data = Data[,-1]
corr <- cor(Data)
pdf(paste0(output_dir,"/Correlate.pdf"),width=20,height=20)
corrplot.mixed( corr,lower = "number", upper = "pie" )
dev.off()


#画聚类图

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ Condition)
dds <- dds[ rowSums(counts(dds)) > 1, ]


############聚类分析，首先对表达count进行标准化处理，之后进行聚类**这一步时间较长
rld <- vst(dds, blind=FALSE)
###############样本聚类图
cor<-cor(assay(rld), method='pearson', use='pairwise.complete.obs')

sampleDists <- as.dist(1-cor)
sampleDistMatrix <- as.matrix(sampleDists)
colnames(sampleDistMatrix) <- colnames(rld)
colors <- colorRampPalette(colors = c("green","yellow","red")) (255)
pheatmap(sampleDistMatrix,
		clustering_distance_rows=sampleDists,
		clustering_distance_cols=sampleDists,
		annotation_col  = anno_col,
		annotation_row  = anno_row,
		col=colors,cellwidth = 15, cellheight = 12, fontsize = 8,filename = paste0(output_dir,"/clustering_sample.pdf") 
		)
		
pdf(file=paste0(output_dir, "/PCA.pdf"))
plotPCA(rld, intgroup=c("Condition"))
dev.off()


