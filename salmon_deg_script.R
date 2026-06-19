################################################### SALMON
library(GenomicFeatures)
library(tximport)
library(txdbmaker)
suppressPackageStartupMessages(library(sva))

working_dir<-'~/storage/AnTrD_personal_task/salmon/'
setwd(working_dir)
getwd()

# Data loading
metadata_salmon<-read.table('salmon.tsv',sep='\t',header = T)

metadata_salmon <- metadata_salmon %>%
  filter(!Sample %in% c("SRR13388751", "SRR13388736"))

metadata_salmon$Group <- factor(metadata_salmon$Group, levels = c("young", "old")) 

tx2gene <- read.table("~/storage/AnTrD_personal_task/salmon/tx2genes.tsv",sep='\t',header=T)
tx2gene$TXNAME<-gsub('\\..*','',tx2gene$TXNAME)
tx2gene$GENEID<-gsub('\\..*','',tx2gene$GENEID)
tx2gene %>% head()

txi <- tximport(metadata_salmon$SalmonPath, type="salmon",
                tx2gene=tx2gene, countsFromAbundance = "lengthScaledTPM",
                ignoreTxVersion = TRUE)

colnames(txi$counts) <- metadata_salmon$Sample
colnames(txi$abundance) <- metadata_salmon$Sample

#Normalization in DESeq2 MRN (~RLE)
#Create DESeq2 object
ddsalmon <- DESeqDataSetFromTximport(txi = txi,
                                colData = metadata_salmon,
                                design = ~ Sex + Group)

#Estimate size factors for normalization
ddsalmon <- estimateSizeFactors(ddsalmon)
range(ddsalmon$sizeFactor) #большой разброс это плохо
hist(ddsalmon$sizeFactor,breaks = 30)
boxplot(ddsalmon$sizeFactor)

#filter out genes via edgeR package (Add min.count to filter noisy genes)
idx <-edgeR::filterByExpr(counts(ddsalmon,normalized=T), 
                          model.matrix(~Sex + Group, metadata_salmon),
                          min.count = 10)
ddsalmon <- ddsalmon[idx,]
print(ddsalmon)

############################################### SVA

norm_counts <- counts(ddsalmon, normalized=TRUE)
mod <- model.matrix(~ Sex + Group, data = metadata_salmon)
mod0 <- model.matrix(~ Sex, data = metadata_salmon)
norm_counts <- norm_counts[rowSums(norm_counts) > 0,]
fit <- svaseq(norm_counts, mod=mod, mod0=mod0) # n.sv=2)
print(paste("Автоматически найдено компонент:", fit$n.sv))

colData(ddsalmon)$SV1 <- fit$sv[,1]
colData(ddsalmon)$SV2 <- fit$sv[,2]

design(ddsalmon) <- ~Sex + SV1 + SV2 + Group

################################################## EDA

counts(ddsalmon,normalized=FALSE)%>%head
counts(ddsalmon,normalized=TRUE)%>%head

#Transform normalized counts with VST
vsd <- varianceStabilizingTransformation(ddsalmon)
vsd <- assay(vsd)

#Calculate correlation between samples with PCC and SCC
PCC_vsd<-cor(vsd,method = 'pearson')
SCC_vsd<-cor(vsd,method = 'spearman')

#Correlations are identical for SCC, but not for PCC (VST vs raw)!
pheatmap(PCC_vsd, main = "PCC_vsd")
pheatmap(SCC_vsd, main = "SCC_vsd")


#Run PCA
#Choose top N variable genes
ntop<-(vsd%>%apply(1,sd))%>%sort%>%rev%>%head(500)%>%names
#ntop<-rownames(vsd)

pca<-prcomp(t(vsd[ntop,]))
pca_df<-data.frame(PC1=pca$x[,1],
                   PC2=pca$x[,2],
                   group=metadata_salmon$Group,
                   sample_id=metadata_salmon$Sample,
                   shape=metadata_salmon$Sex)

pca_SN<-ggplot(pca_df,aes(x = PC1, y = PC2, color = group,shape = shape,label=sample_id)) +
  geom_text(hjust=0, vjust=0,show.legend = F,size=3)+
  xlab(paste0('PC1 ',summary(pca)$importance[2,1]))+ylab(paste0('PC2 ',summary(pca)$importance[2,2]))+
  geom_point(size = 2) + #scale_color_manual(values = brewer.pal(4,'Set1')[1:2]) +#stat_ellipse()+
  theme_bw() + coord_cartesian(clip = "off") + xlim(-25, 70) 

pca_SN
dev.off()

################################################ remove batch effect

#Remove batch-effect from vsd for visualization
mm <- model.matrix(~Group, colData(ddsalmon))
#vsd_nobatch<-removeBatchEffect(vsd,design = mm,batch = sampleBatch)

##Remove batch with continuous variable (e.g. Svaseq's surogate variable(sv))
#vsd_nobatch<-removeBatchEffect(vsd,design = mm,covariates = colData(dds)$sv1)
##Remove batch with continuous and categorical variables (e.g. Svaseq's sv + known batch)
vsd_nobatch<-removeBatchEffect(vsd,design = mm,
                               batch = colData(ddsalmon)$Sex,
                               covariates = colData(ddsalmon)$SV1)

#Choose top N variable genes
#ntop<-rownames(vsd_nobatch)
ntop<-(vsd_nobatch%>%apply(1,sd))%>%sort%>%rev%>%head(500)%>%names

pca<-prcomp(t(vsd_nobatch[ntop,]))
pca_df<-data.frame(PC1=pca$x[,1],
                   PC2=pca$x[,2],
                   condition=metadata_salmon$Group,
                   batch=metadata_salmon$Sex,
                   sample_id=metadata_salmon$Sample)

ggplot(pca_df,aes(x = PC1, y = PC2, color = condition,label=sample_id)) +
  geom_text(hjust=0, vjust=0,show.legend = F,size=3)+
  xlab(paste0('PC1 ',summary(pca)$importance[2,1]))+ylab(paste0('PC2 ',summary(pca)$importance[2,2]))+
  geom_point(size = 2) + xlim(-25, 70) +
  theme_bw()+
  ggplot(pca_df,aes(x = PC1, y = PC2, color = batch,label=sample_id)) +
  geom_text(hjust=0, vjust=0,show.legend = F,size=3)+
  xlab(paste0('PC1 ',summary(pca)$importance[2,1]))+ylab(paste0('PC2 ',summary(pca)$importance[2,2]))+
  geom_point(size = 2) +xlim(-25, 70) +
  theme_bw()


################################################# DESEQ

dds_analysis <- DESeq(ddsalmon)
#Plot dispersion estimates from DESeq2
plotDispEsts(dds_analysis, main="Dispersion plot")

resultsNames(dds_analysis)
res <- results(dds_analysis, name="Group_old_vs_young")
summary(res)
as.data.frame(res)%>%subset(!is.na(padj))

y_o<-as.data.frame(res)%>%subset(!is.na(padj))

#Obtain gene symbols from ENS
gene_names <- data.frame(symbol=mapIds(org.Hs.eg.db, sub("\\..*$", "", rownames(y_o)), "SYMBOL", "ENSEMBL"),
                         entrez=mapIds(org.Hs.eg.db, sub("\\..*$", "", rownames(y_o)), "ENTREZID", "ENSEMBL"),
                         row.names = y_o%>%rownames)%>%mutate(across(.fns=as.character))

y_o<-cbind(gene_names,y_o)

top_genes_salmon <- y_o %>% 
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  mutate(Direction = ifelse(log2FoldChange > 0, "Up", "Down")) %>%
  arrange(desc(log2FoldChange), padj)

top_genes_salmon %>% 
 filter(Direction == "Down") %>%
  nrow()

top_genes_salmon %>% 
  dplyr::select(symbol, log2FoldChange, stat, padj) %>% View

View(top_genes_salmon)

y_o %>%
  as.data.frame() %>%
  rownames_to_column(var = "ens_id") %>%
  write_csv("../results/all_genes_salmon.csv")

getwd()
write_csv(top_genes_salmon, "../results/top_genes_salmon.csv")

################################################################ VOLCANO

EnhancedVolcano<-EnhancedVolcano(y_o, lab = y_o$symbol, title = 'Young vs Old', subtitle = '',
                                 legendLabels = c('NS','Log2FC', 'padj', 
                                                  'Log2FC padj'), 
                                 legendPosition = 'top', x = 'log2FoldChange', y = 'padj', 
                                 xlab = bquote(~Log[2]~ 'fold change'),
                                 col = c('grey30', brewer.pal(9,'Set1')[2],'gold',brewer.pal(9,'Set1')[1]),
                                 pCutoff = 0.05, FCcutoff = 1,
                                 labSize = 4, #legendIconSize = 12, 
                                 #pointSize = 4, legendLabSize = 30, axisLabSize = 30, captionLabSize = 25,
                                 gridlines.major = F, gridlines.minor = F)

#png(filename=paste0(DE_dir,DE_name,"_Volkanoplot.png"),1100,1000,res = 175)
ggsave("../results/EnhancedVolcano_salmon.png", plot = EnhancedVolcano, width = 10, height = 5)
plot(EnhancedVolcano)

#######################################################

#Create gct object for Phantasus
counts<-counts(ddsalmon, normalized=T)
counts <- counts[which(rownames(counts) %in% rownames(y_o)),]
counts <- counts[match(rownames(y_o),rownames(counts)),]

es <- ExpressionSet(counts)
es$condition <- colData(ddsalmon)$Group
es$batch <- colData(ddsalmon)$Sex

head(es)
fData(es)<-y_o%>%as.data.frame
fData(es)$Row.names <- NULL

head(fData(es))
pData(es)
#Top 15 up- and downregulated genes
DEGs<-arrange(fData(es),-stat)%>%{rbind(head(.,20),tail(.,20))}
DEGs<-merge(DEGs, counts, by=0, all=F)%>%arrange(-stat)
DEGs<-dplyr::select(DEGs,c('symbol',tail(colnames(DEGs),ncol(counts))))
DEGs <- DEGs %>% filter(!is.na(symbol) & symbol != "")
DEGs$symbol <- make.unique(DEGs$symbol)
DEGs<-DEGs%>%column_to_rownames('symbol')

#png(filename=paste0(DE_dir,DE_name,"_top20_pheatmap.png"), width = 1600 , height = 1600,res=135)
top20_salmon <-pheatmap(DEGs,
                 annotation_col = 
                   colData(ddsalmon)[, c("Group"),drop=F]%>%as.data.frame(),
                 legend = T,
                 color = rev(RColorBrewer::brewer.pal(9,"RdYlBu")),
                 main = '', cluster_rows = F, cluster_cols = F,
                 scale = 'row')

ggsave("../results/top20_salmon.png", plot = top20_salmon, width = 10, height = 5)
#dev.off()

####################################################### GSEA
#######################################################

p_spec <- "Homo_sapiens"
p_cat<-'H'
res <- y_o%>%arrange(stat)
head(res)

expression <- res[,c("symbol", "stat")]
expression <- expression[!duplicated(expression$symbol),]
expression <- subset(expression, symbol!= "<NA>")

head(expression)

ranks <- setNames(expression$stat, expression$symbol)

head(ranks)
pathways = msigdbr(species = "Homo sapiens", category = "H") %>% split(x = .$gene_symbol, f = .$gs_name)
fgseaRes <- fgseaMultilevel(pathways = pathways,
                            stats = ranks,
                            eps=0.0,
                            minSize=15,
                            maxSize=500)%>%arrange(-NES)

topPathwaysUp <- fgseaRes[(ES > 0)&(pval<0.1)&(abs(NES)>1.00)][head(rev(order(NES)), n=20), pathway]
topPathwaysDown <- fgseaRes[(ES < 0)&(pval<0.1)&(abs(NES)>1.00)][head(order(NES), n=20), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

#fwrite(fgseaRes, file=paste0(gsub(' ','_',p_spec),'_',p_cat,'.tsv'), sep="\t", sep2=c("", " ", ""))

plotGseaTable(pathways[topPathways], ranks, fgseaRes, gseaParam = 0.5)
png(filename=paste0(gsub(' ','_',p_spec),'_',p_cat,"_GseaTable.png"),
    width = 2000, height = 2000,res = 135)
plot(plot)
dev.off()


##########################################################################


fgseaRes$leadingEdge<-fgseaRes%>%pull(leadingEdge)%>%lapply(function(x)paste(x,collapse = ' '))%>%unlist
list_gsea<-list(fgseaRes%>%as.data.frame)
names(list_gsea)[1]<-p_cat

gsea_plots<-lapply(names(list_gsea),function(name){
  
  pathway<-name
  data<-list_gsea[[pathway]]
  data<-data%>%subset(padj<0.05)
  
  pos=sum(data$NES>0)
  neg=sum(data$NES<0)
  
  n_top<-5
  
  if(pos<n_top&neg<n_top){data<-data%>%{rbind(head(.,pos),tail(.,neg))}
  }else if(pos<n_top){data<-data%>%{rbind(head(.,pos),tail(.,n_top))}
  }else if(neg<n_top){data<-data%>%{rbind(head(.,n_top),tail(.,neg))}
  }else if(pos>=n_top&neg>=n_top){data<-data%>%{rbind(head(.,n_top),tail(.,n_top))}
  }else{data<-NULL}
  
  if (pos>0&neg>0){pallete_use='RdYlBu'
  }else if (pos==0){pallete_use='Blues'
  }else{pallete_use='Reds'}
  
  
  
  data$pathway<-data$pathway%>%gsub('_',' ',.)
  long<-data$pathway%>%str_split(' ')%>%lapply(length)>2
  
  data$pathway[long]<-data$pathway[long]%>%lapply(function(x){
    str_spl_x<-x%>%str_split(' ')%>%unlist
    len_test<-str_spl_x%>%length
    if (len_test%%2==0){
      x<-paste0(
        paste(str_spl_x[1:(len_test/2)],collapse = ' '),'\n',
        paste(str_spl_x[-(1:(len_test/2))],collapse = ' '))
      return(x)
    }else{
      x<-paste0(
        paste(str_spl_x[1:(len_test%/%2)],collapse = ' '),'\n',
        paste(str_spl_x[-(1:(len_test%/%2))],collapse = ' '))
      return(x)
    }})%>%unlist
  
  data$pathway<-factor(data$pathway,levels=data%>%arrange(NES)%>%pull(1))
  plt <- ggplot(data)+
    geom_col(aes(NES, pathway, fill = NES),width = 0.8) + 
    scale_fill_distiller(palette = pallete_use)+
    geom_text(
      data = data%>%subset(NES>0),
      aes(0+0.05, y = pathway, label = pathway),
      hjust = 0,
      size = 3,
      #nudge_x = 0.3,
      lineheight = .75,
      colour = "black",
      family = "Econ Sans Cnd"
    )+
    geom_text(
      data = data %>% subset(NES < 0),
      # Привязываем к нулю и отступаем чуть-чуть влево (-0.05)
      aes(0 - 0.05, y = pathway, label = pathway),
      hjust = 1, # Выравнивание по правому краю текста
      size = 3,
      lineheight = .75,
      colour = "black",
      family = "Econ Sans Cnd")+xlim(min(data$NES) - 1, max(data$NES) + 1)+
    theme_bw()+theme(axis.text.y=element_blank(),
                     text = element_text(size=15),
                     axis.ticks.y=element_blank(),
                     panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
  
  #png(file=paste0(p_spec,'_',pathway,'_GSEA.png'),1100,800,res=200)
  plot(plt)
  dev.off()
  return(plt)
})
names(gsea_plots)<-names(list_gsea)

library(Cairo)
lapply(names(gsea_plots)[1],function(pathway){
  ggsave(
    filename = paste0(p_spec, "_", pathway, "_GSEA.pdf"),
    plot = gsea_plots[[pathway]],
    device = CairoPDF,
    width = 5.5,
    height = 4,
    units = "in"
  )
})
