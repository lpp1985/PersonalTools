#!/usr/bin/env Rscript 
library(docopt)
library(plyr)
"Usage: clusterProfiler.R  [options]  INPUT OUTPUTDIR
Options:
   -s --srgan=srgan    the organism, the value in <hsa, mmu, rno> [default: hsa]

Arguments:
   INPUT       the name of input files
   OUTPUTDIR   the output dir name" -> doc

opts   <- docopt(doc)
input  <- opts$INPUT
prefix <- opts$OUTPUTDIR 
organ  <- opts$s

dir.create(prefix, showWarnings = FALSE,recursive = TRUE)


library(clusterProfiler)
library(topGO)
library(rlist)
library(tidyr)
db       <- list(
'hsa' = "org.Hs.eg.db",
'mmu' = 'org.Mm.eg.db',
'rno' = 'org.Rn.eg.db',
'eco' = 'org.EcK12.eg.db'

)

gene     <- read.delim(input, stringsAsFactors=F,header=T)$Name

all_data <- read.delim(input, stringsAsFactors=F,header=T)


eg <- bitr(gene, fromType="ENSEMBL", toType="ENTREZID", OrgDb=db[[organ]])
symbol <- bitr(gene, fromType="ENSEMBL", toType="SYMBOL", OrgDb=db[[organ]])
all_gene <- read.delim(input, stringsAsFactors=F,header=T)
colnames(symbol)[1]<-"Name"
colnames(eg)[1] <- "Name"

if (!("SYMBOL" %in% colnames(all_gene))) {
	merge(all_gene,symbol,by="Name")->all_gene


	
}

if (!("ENTREZID" %in% colnames(all_gene))) {
	merge(all_gene,eg,by="Name")->all_gene


	
}

write.table(all_gene, input, sep="\t", quote=F, row.names=F)

gene     <- read.delim(input, stringsAsFactors=F,header=T)$SYMBOL

if (  ( "Status"   %in% colnames(all_gene) )  ){
	all_data <- all_gene[,c("SYMBOL","Status")]
	write.table( all_data, paste0(input,".status" ),  sep="\t", quote=F, row.names=F)
}
eg <- bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb=db[[organ]])

gene     <- eg[[2]]

ego_cc <- enrichGO(gene          = gene,
                   OrgDb=db[[organ]],
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)



ego_bp <- enrichGO(gene          = gene,
                   OrgDb=db[[organ]],
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)


ego_mf <- enrichGO(gene          = gene,
                   OrgDb=db[[organ]],
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

cc <- as.data.frame(ego_cc)
if(nrow(cc) >= 1) cc$Type <- "CC"

bp <- as.data.frame(ego_bp)
if(nrow(bp) >= 1) bp$Type <- "BP"

mf <- as.data.frame(ego_mf)
if(nrow(mf) >= 1) mf$Type <- "MF"

dir.create(prefix, showWarnings = FALSE)

go_enrichment   <- paste(prefix, "go_enrichment.xls", sep="/")
go_enrichmentanno   <- paste(prefix, "go_enrichmentWithAnnotation.xls", sep="/")



kegg_enrichment <- paste(prefix, "kegg_enrichment.xls", sep="/")



go_bp_pdf          <- paste(prefix, "go.bp.pdf", sep="/")
go_cc_pdf          <- paste(prefix, "go.cc.pdf", sep="/")
go_mf_pdf          <- paste(prefix, "go.mf.pdf", sep="/")


pdf(go_bp_pdf)
plotGOgraph(ego_bp)
dev.off()

pdf(go_cc_pdf)
plotGOgraph(ego_cc)
dev.off()

pdf(go_mf_pdf)
plotGOgraph(ego_mf)
dev.off()



res <- rbind(bp, cc, mf)

res_raw <- res
res_raw$geneID <- unlist(lapply(1:length(res_raw$geneID), function(t){  paste(eg[[2]][eg[[1]] %in% unlist(strsplit(res_raw$geneID[t], "/"))], collapse = "/")  }))
res_raw<-as.data.frame(res_raw)
rename( res_raw,c("geneID"="ENTREZID")    ) -> res_raw
rename( res_raw,c("pvalue"="EnrichedPvalue")    ) -> res_raw

res_raw %>% separate_rows(ENTREZID, sep = "/") ->longtext
all_data <- read.delim(input, stringsAsFactors=F,header=T)
merge(  longtext,all_data,by="ENTREZID") ->output_long

if(nrow(res) >= 1){
	
	write.table(res, go_enrichment, sep="\t", quote=F, row.names=F)

	write.table(output_long, go_enrichmentanno, sep="\t", quote=F, row.names=F)


}
ego <- res
types       <- unique(ego$Type)
DataForPlot <- list()
for(j in 1:length(types)){
    DataForPlot[[j]] <-  ego[grep(types[j],ego$Type),]
    # 展示最显著的topN个
    TopN <- 10
    if(nrow(DataForPlot[[j]]) < TopN){

        DataForPlot[[j]] <- DataForPlot[[j]][c(1:nrow(DataForPlot[[j]])), c(2,5,9,10)]
        DataForPlot[[j]] <- DataForPlot[[j]][order(DataForPlot[[j]]$Count, decreasing = T), ]
    }else{

        DataForPlot[[j]]=DataForPlot[[j]][c(1:TopN),c(2,5,9,10)]
        DataForPlot[[j]]=DataForPlot[[j]][order(DataForPlot[[j]]$Count,decreasing = T),]
    }
                
}

GO <- list.rbind(DataForPlot)     # rbind all elements in a list by row
str(GO)

color  <- colorRampPalette(c("light blue","pink","light green"))(length(unique(GO$Type)))
color1 <- rep(color,as.numeric(table(GO$Type)))

#########plot go_barplot#####################
go_bar <- paste(prefix, "GO_barplot.pdf", sep="/")
color  <- colorRampPalette(c("Goldenrod1","Tomato1","MediumPurple4"))(length(unique(GO$Type)))
color1 <- rep(color,as.numeric(table(GO$Type)))

pdf(go_bar, width=15, height=8)
if(max(nchar(GO$Description)) >= 100){
	layout(matrix(c(3,2,1), nrow = 1, byrow = T),widths = c(1.2,0.1,1))
}else{
	layout(matrix(c(3,2,1), nrow = 1, byrow = T),widths = c(0.8,0.1,1))
}

par(mar = c(5,0.1, 1, 3))
gobar <- barplot(GO$Count, plot=T, 
        cex.lab=1, las=1, ps=0.5, border = F,
        xlim=c(0,max(GO$Count)*1.25),
        #xaxt= "n", 
        cex.names=0.8, axis.lty=1, 
        axes=TRUE,col=color1,horiz=T,
        mgp=c(0,-0.5,-1))
text(cex = 1, y = gobar, x = GO$Count+max(GO$Count)/35, lab=c(GO$Count))
#axis(side = 1, at =c(0,max(GO$Count)*1.25),labels = c("0","Gene Number"))
title(xlab = "Gene Number",line = 1, cex.lab = 1.2)
legend("topright", legend=unique(GO$Type),bty="n",fill=color,cex = 1.2)
y1 = par("usr")[3]
y2 = par("usr")[4]

par(mar = c(5,0, 1, 0))
p = GO$pvalue * -1
p_bar = barplot(p,horiz=T,col="PaleGreen3",border = NA, xaxt= "n")
max = format(max(GO$pvalue),scientific=TRUE,digit=2)
axis(side = 1, line = -1, at = c(0,min(p)), labels = F)
title(xlab = "P value",line = 0.91, cex.lab = 1.2)

par(mar = c(5,3, 1, 1))
a = plot(1:5,ylim = c(y1,y2),type = "n", xaxs = "i", yaxs = "i", axes=F, xlab="",ylab="")
for(i in 1:length(GO$Description)){
    text(x = 5 ,y = gobar[i], labels = GO$Description[i], cex = 1.4, adj = 1)
}
dev.off()

# # pdf(go_bar, width=15, height=8)
# # par(mar = c(3,max(nchar(GO$Description)) * 0.37, 1, 2)) 
# # gobar <- barplot(GO$Count, plot=T, names.arg=GO$Description, 
# # 		cex.lab=1, las=1, ps=0.5, 
# # 		xlim=c(0,max(GO$Count)*1.25), 
# # 		cex.names=0.8, axis.lty=1, 
# # 		xlab="Gene Number",axes=TRUE,col=color1,horiz=T,
# # 		mgp=c(1.95,0.55,0))
# # text(cex = 0.8, y = gobar, x = GO$Count+max(GO$Count)/35, lab=c(GO$Count))
# # legend("topright", legend=unique(GO$Type),bty="n",fill=color)
# # dev.off()

#########plot kegg#####################
all_data <- read.delim(input, stringsAsFactors=F,header=T)
kk <- enrichKEGG(gene         = gene,
                 organism = organ,
                 pvalueCutoff = 0.05,
		 qvalueCutoff = 0.05)
kk_enrich <- as.data.frame(kk)
kk_enrich$geneID <- unlist(lapply(1:length(kk_enrich$geneID), function(t){  paste(eg[[1]][eg[[2]] %in% unlist(strsplit(kk_enrich$geneID[t], "/"))], collapse = "/")  }))

kk_enrich_raw <- as.data.frame(kk)
rename(kk_enrich_raw,c("geneID"="ENTREZID") )->kk_enrich_raw


rename (kk_enrich,c("ID"="PathwayID")   )-> kk_enrich
rename (kk_enrich_raw,c("ID"="PathwayID")   )-> kk_enrich_raw




kk_enrich_raw <- rename(   kk_enrich_raw,c("pvalue"="EnrichedPvalue") )
kk_enrich_raw %>% separate_rows(ENTREZID, sep = "/") ->longtext

kegg_enrichmenglong <- paste(prefix, "kegg_enrichmentWithAnnotation.xls", sep="/")

kegg_pdf <- paste(prefix, "kegg_dotplot.pdf", sep="/")
if(nrow(kk_enrich) >= 1){
	pdf(kegg_pdf,height=7,width=12)
	p <-dotplot(kk)
	print(p)
	dev.off()
	merge(  longtext,all_data,by="ENTREZID") ->output_long
	output_long<- output_long[order(output_long$PathwayID),]
	write.table(output_long, kegg_enrichmenglong, sep="\t", quote=F, row.names=F)
	write.table(kk_enrich, kegg_enrichment, sep="\t", quote=F, row.names=F)
}


###############gsea analysis############
geneList <- unique(all_data[,c("ENTREZID","log2FoldChange")] )
glist <- geneList[,2]
names(glist) <- as.character(geneList[,1])

glist <- sort(glist,decreasing = T)

 
kk <- gseKEGG(gene         =glist  ,
				organism = organ,	
				pvalueCutoff =0.05,
				pAdjustMethod     = "BH",
				)

kk_enrich <- as.data.frame(kk)
kk_enrich_raw <- as.data.frame(kk)
rename(kk_enrich,c("core_enrichment"="geneID") ) -> kk_enrich
rename (kk_enrich,c("ID"="PathwayID")   )-> kk_enrich

rename(kk_enrich_raw,c("core_enrichment"="ENTREZID") )->kk_enrich_raw
rename(kk_enrich_raw,c("ID"="PathwayID")   )-> kk_enrich_raw
kk_enrich_raw <- rename(   kk_enrich_raw,c("pvalue"="EnrichedPvalue") )


kk_enrich$geneID <- unlist(lapply(1:length(kk_enrich$geneID), function(t){  paste(eg[[1]][eg[[2]] %in% unlist(strsplit(kk_enrich$geneID[t], "/"))], collapse = "/")  }))
kegg_pdf <- paste(prefix, "kegg_gsea_dotplot.pdf", sep="/")
kegg_enrichment <- paste(prefix, "kegg_gsea_enrichment.xls", sep="/")
kegg_enrichmenglong <- paste(prefix, "kegg_gsea_enrichmentWithAnnotation.xls", sep="/")



if(nrow(kk_enrich) >= 1){
	pdf(kegg_pdf,height=7,width=12)
	p <-dotplot(kk)
	print(p)
	dev.off()
	
	write.table(kk_enrich, kegg_enrichment, sep="\t", quote=F, row.names=F)
	kk_enrich_raw %>% separate_rows(ENTREZID, sep = "/") ->longtext
	merge(  longtext,all_data,by="ENTREZID") ->output_long
	output_long<- output_long[order(output_long$PathwayID),]
	write.table(output_long, kegg_enrichmenglong, sep="\t", quote=F, row.names=F)
	
}


