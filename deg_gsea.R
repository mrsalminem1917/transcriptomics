######################################################### Uploading data
chooseCRANmirror()
BiocManager::install('EnhancedVolcano')
suppressPackageStartupMessages(library(EnhancedVolcano))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(patchwork))

########################

setwd("/home/student/storage/AnTrD_personal_task/")
metadata <- read.table('metadata.tsv', header = TRUE)
rownames(metadata)<-metadata$Sample
raw_counts <- lapply(metadata$FeaturePath, function(path){
  count_mtx <- read.table(path, sep='\t', header=TRUE, row.names = 'Geneid') %>%
    dplyr::select(6)
  return(count_mtx)
}) %>% bind_cols()
colnames(raw_counts) <- metadata$Sample
identical(colnames(raw_counts),rownames(metadata))


###### REMOVE OUTLIERS

raw_counts <- raw_counts[, !(colnames(raw_counts) %in% c("SRR13388751"))] #"SRR13388736"
metadata <- metadata %>%
  filter(!Sample %in% c("SRR13388751"))

######################

metadata$Group <- factor(metadata$Group, levels = c("young", "old")) 
metadata$Sex <- factor(metadata$Sex, levels = unique(metadata$Sex))
#metadata$Splice <- c(7184292, 15102220, 6288151, 9141049, 13762893, 7410666, 11548654, 10753811, 7602493, 12801877)
#metadata$Depth <- c(44186692, 49659198, 35204411, 53877644, 113363405, 43964127, 56014356, 116732027, 36377529, 53696405)
#metadata$Batch <- metadata$Splice/metadata$Depth
#metadata$Dev <- metadata$Splice/metadata$Depth
#metadata$Batch <- c("low", "medium", "low", "medium", "high", "low", "medium", "high", "low", "medium")
#metadata$Batch <- factor(metadata$Batch, levels = c("low", "medium", "high"))
#metadata$Batch <- as.numeric(scale(metadata$Batch))
#metadata$Batch <- cut(metadata$Batch, breaks = 2, labels = c("low", "high"))

############################################################# DESEQ2

design_matrix <- ~Sex + Group

dds <- DESeqDataSetFromMatrix(countData = raw_counts, 
                              colData = metadata, 
                              design = design_matrix)
dds <- estimateSizeFactors(dds) #RLE normalization

# Проверка распределения (разброс в пределах 0.5-2.0 — это норма)
print(range(dds$sizeFactor))
hist(dds$sizeFactor, breaks = 10, main="Size Factors Distribution")

idx <- filterByExpr(counts(dds, normalized=TRUE),
                   design = design_matrix,
                   group = metadata$Group,
                   min.count = 10)
#idx <- rowSums(raw_counts >= 5 ) >= 5
dds <- dds[idx, ]


################################################################### SVA
#https://biodatascience.github.io/compbio/dist/sva.html

norm_counts <- counts(dds, normalized=TRUE)
mod <- model.matrix(~ Sex + Group, data = metadata)
mod0 <- model.matrix(~ Sex, data = metadata)
norm_counts <- norm_counts[rowSums(norm_counts) > 0,]
fit <- svaseq(norm_counts, mod=mod, mod0=mod0)
print(paste("Автоматически найдено компонент:", fit$n.sv))

colData(dds)$SV1 <- fit$sv[,1]
colData(dds)$SV2 <- fit$sv[,2]

design(dds) <- ~ Sex + SV1 + SV2 + Group


##########################
vsd <- varianceStabilizingTransformation(dds)
vsd <- assay(vsd)

ntop<-(vsd%>%apply(1,sd))%>%sort%>%rev%>%head(500)%>%names

pca<-prcomp(t(vsd[ntop,]))
pca_df <- data.frame(PC1 = pca$x[, 1],
                     PC2 = pca$x[, 2],
                     group = metadata$Group,
                     sample_id = metadata$Sample,
                     #batch = metadata$Batch,
                     sex = metadata$Sex)

#pca_FC<-ggplot(pca_df,aes(x = PC1, y = PC2, color = group ,label=sample_id)) +
  #geom_text(hjust=0, vjust=0,show.legend = F,size=3)+
  #xlab(paste0('PC1 ',summary(pca)$importance[2,1]))+ylab(paste0('PC2 ',summary(pca)$importance[2,2]))+
  #geom_point(size = 3) + #scale_color_manual(values = brewer.pal(4,'Set1')[1:2]) +#stat_ellipse()+
  #theme_bw() + coord_cartesian(clip = "off") #+ xlim(-90, 120) 
#pca_FC

p1 <- ggplot(pca_df,aes(x = PC1, y = PC2, color = sex,label=sample_id)) +
  geom_text(hjust=0, vjust=0,show.legend = F,size=3)+
  xlab(paste0('PC1 ',summary(pca)$importance[2,1]))+ylab(paste0('PC2 ',summary(pca)$importance[2,2]))+
  geom_point(size = 2) + xlim(-60, 60) 
  theme_bw()
p2 <- ggplot(pca_df,aes(x = PC1, y = PC2, color = group,label=sample_id)) +
  geom_text(hjust=0, vjust=0,show.legend = F,size=3)+
  xlab(paste0('PC1 ',summary(pca)$importance[2,1]))+ylab(paste0('PC2 ',summary(pca)$importance[2,2]))+
  geom_point(size = 2) + xlim(-60, 60) 
  theme_bw()

pca_graph <- p1+p2
ggsave("./results/pca_graph_without_remove_batch_effect.png", plot = pca_graph, width = 10, height = 5)

#########################################################
#design_matrix <- model.matrix(~ Group, data = metadata)

# С помощью limma временно "вычитаем" эффекты пола и плавающего качества библиотек ТОЛЬКО для красивого графика

#mat_corrected <- removeBatchEffect(
 # vsd, 
 # batch = colData(dds)$Sex, 
  #covariates = colData(dds)$SV1,
  #design = design_matrix)

# Строим новый MDS на очищенной матрице
#plotMDS(mat_corrected, col=ifelse(colData(dds)$Group == "young", "cyan4", "violetred"))

############################################################### PCA REMOVE BATCH EFFECT
#Remove batch-effect from vsd for visualization
mm <- model.matrix(~Group, colData(dds))
#vsd_nobatch<-removeBatchEffect(vsd,design = mm, batch = colData(dds)$Batch)

##Remove batch with continuous variable (e.g. Svaseq's surogate variable(sv))
#vsd_nobatch<-removeBatchEffect(vsd,design = mm,covariates = colData(dds)$sv1)
##Remove batch with continuous and categorical variables (e.g. Svaseq's sv + known batch)
vsd_nobatch<-removeBatchEffect(vsd,design = mm,batch = colData(dds)$Sex,covariates = colData(dds)$SV1) 

#Choose top N variable genes
#ntop<-rownames(vsd_nobatch)
ntop<-(vsd_nobatch%>%apply(1,sd))%>%sort%>%rev%>%head(500)%>%names

pca<-prcomp(t(vsd_nobatch[ntop,]))
pca_df<-data.frame(PC1=pca$x[,1],
                   PC2=pca$x[,2],
                   condition=metadata$Group,
                   batch=metadata$Sex,
                   sample_id=metadata$Sample)

wobatch <- ggplot(pca_df,aes(x = PC1, y = PC2, color = condition,label=sample_id)) +
  geom_text(hjust=0, vjust=0,show.legend = F,size=3)+
  xlab(paste0('PC1 ',summary(pca)$importance[2,1]))+ylab(paste0('PC2 ',summary(pca)$importance[2,2]))+
  geom_point(size = 2) + xlim(-60, 60) +
  theme_bw()+
  ggplot(pca_df,aes(x = PC1, y = PC2, color = batch,label=sample_id)) +
  geom_text(hjust=0, vjust=0,show.legend = F,size=3)+
  xlab(paste0('PC1 ',summary(pca)$importance[2,1]))+ylab(paste0('PC2 ',summary(pca)$importance[2,2]))+
  geom_point(size = 2) + xlim(-60, 60) +
  theme_bw()

ggsave("./results/pca_graph_wo_batch_effect.png", plot = wobatch, width = 10, height = 5)

################################################

# Classical MDS
# N rows (objects) x p columns (variables)
# each row identified by a unique row name
d <- cor(vsd_nobatch[ntop,],method = 'sp')-1 # Spearman cor distances between the samples
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
fit # view results

# plot solution
mds_df<-data.frame(MDS1=fit$points[,1],
                   MDS2=fit$points[,2],
                   condition=metadata$Group,
                   batch=metadata$Sex,
                   sample_id=metadata$Sample)

mds <- ggplot(mds_df,aes(x = MDS1, y = MDS2, color = condition,label=sample_id)) +
  geom_text(hjust=0, vjust=0,show.legend = F,size=3)+
  xlab('MDS1')+ylab('MDS2')+
  geom_point(size = 2) + xlim(-0.5, 0.5) +
  theme_bw()+
  ggplot(mds_df,aes(x = MDS1, y = MDS2, color = batch,label=sample_id)) +
  geom_text(hjust=0, vjust=0,show.legend = F,size=3)+
  xlab('MDS1')+ylab('MDS2')+
  geom_point(size = 2)+ xlim(-0.5, 0.5) +
  theme_bw()

ggsave("./results/mds_wo_batch_effect.png", plot = mds, width = 10, height = 5)


######################### CORR
PCC_vsd<-cor(vsd,method = 'pearson')
SCC_vsd<-cor(vsd,method = 'spearman')
pheatmap(PCC_vsd, main = "PCC_vsd")
spearman <- pheatmap(SCC_vsd, main = "SCC_vsd")

ggsave("./results/spearman.png", plot = spearman, width = 10, height = 5)

############################

calc_CV <- function(x) {sd(x) / mean(x)}
CV_dat <- apply(vsd,1,calc_CV)

# show the highest CD genes and their CV values
sort(CV_dat, decreasing = TRUE)[1:5]
#>   CAMK2N1    AKR1C1      AQP2     GDF15       REN 
#> 0.5886006 0.5114973 0.4607206 0.4196469 0.4193216

# Identify genes in the top 3rd of the CV values
GOI <- names(CV_dat)[CV_dat > quantile(CV_dat, 0.75)]

#Plot heatmap with top 25% variable genes
largeheatmap <- pheatmap(vsd[GOI,],
         scale = "row", 
         show_rownames = F, show_colnames = T,
         border_color = NA,
         clustering_method = "average",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         breaks = seq(-3, 3, 0.05),
         #color = colorRampPalette(brewer.pal(8,'RdYlBu')%>%rev)(120),
         annotation_col = 
           colData(dds)[, c("Group", "Sex"),drop=F]%>%as.data.frame())

ggsave("./results/largeheatmap.png", plot = largeheatmap, width = 10, height = 5)

####################################################
#Create gct object for Phantasus
counts<-counts(dds, normalized=T)
counts <- counts[which(rownames(counts) %in% rownames(Group_young_vs_old)),]
counts <- counts[match(rownames(Group_young_vs_old),rownames(counts)),]

es <- ExpressionSet(counts)
es$condition <- colData(dds)$Group
es$batch <- colData(dds)$Sex

head(es)
fData(es)<-Group_young_vs_old%>%as.data.frame
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
top20 <-pheatmap(DEGs,
         annotation_col = 
           colData(dds)[, c("Group"),drop=F]%>%as.data.frame(),
         legend = T,
         color = rev(RColorBrewer::brewer.pal(9,"RdYlBu")),
         main = '', cluster_rows = F, cluster_cols = F,
         scale = 'row')

ggsave("./results/top20.png", plot = top20, width = 10, height = 5)
#dev.off()

########################### MDS

d <- dist(t(vsd))
fit <- MASS::isoMDS(d, k=2) # k is the number of dim
mds_df<-data.frame(MDS1=fit$points[,1],
                   MDS2=fit$points[,2],
                   group=metadata$Group,
                   sample_id=metadata$Sample)
ggplot(mds_df,aes(x = MDS1, y = MDS2, color = group,label=sample_id)) +
  geom_text(hjust=0, vjust=0,show.legend = F,size=3)+
  xlab('MDS1')+ylab('MDS2')+
  geom_point(size = 2) + #scale_color_manual(values = brewer.pal(4,'Set1')[1:2]) +#stat_ellipse()+
  theme_bw()

######################################################### DESEQ WALD

design(dds)
dds <- DESeq(dds, test = 'Wald')
plotDispEsts(dds, main="Dispersion plot")
resultsNames(dds)
res_wald <- results(dds, name="Group_old_vs_young")
Group_young_vs_old <- as.data.frame(res_wald) %>% subset(!is.na(padj))
gene_names <- data.frame(
  symbol = mapIds(org.Hs.eg.db, sub("\\..*$", "", rownames(Group_young_vs_old)), "SYMBOL", "ENSEMBL"),
  entrez = mapIds(org.Hs.eg.db, sub("\\..*$", "", rownames(Group_young_vs_old)), "ENTREZID", "ENSEMBL"),
  row.names = rownames(Group_young_vs_old)) %>% mutate(across(.fns=as.character))
Group_young_vs_old <- cbind(gene_names, Group_young_vs_old)
top_genes <- Group_young_vs_old %>% 
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  arrange(padj)

View(top_genes)

as.vector(top_genes$symbol[!is.na(top_genes$symbol)])

top_genes<- Group_young_vs_old %>% 
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  mutate(Direction = ifelse(log2FoldChange > 0, "Up", "Down")) %>%
  arrange(desc(log2FoldChange), padj)

top_genes %>% 
  filter(Direction == "Up") %>%
  nrow()

top_genes %>% 
  dplyr::select(symbol, log2FoldChange, stat, padj) %>% View

View(top_genes)




getwd()
write_csv(top_genes, "./results/top_exon_level.csv")
write_csv(Group_young_vs_old, "./results/all_deg_gene_level.csv")

######################################################### VOLCANO
EnhancedVolcano<-EnhancedVolcano(Group_young_vs_old, lab = Group_young_vs_old$symbol, title = 'Young vs Old', subtitle = '',
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
ggsave("./results/EnhancedVolcano_exon_level.png", plot = EnhancedVolcano, width = 10, height = 5)
plot(EnhancedVolcano)
#dev.off()


###########################################################





# Создаем списки генов (замените ddeg_list_2 на ваш второй список генов)
genes_salmon <- top_genes_salmon[!is.na(top_genes_salmon$symbol), "symbol"]
genes <- top_genes[!is.na(top_genes$symbol), "symbol"]

c("CDKN2B", "FAM83B", "LG1", "CFAP61", "SKAP2", "C12orf75", "CRIM1", "OXCT1",
  "NR2F2", "LRP1B", "ENAM", "LPP", "ESRRG", "EDA2R", "KCNQ5", "FAM117B",
  "ZNF844") %>% cat(., sep = "\n")

down <- c(
  "LDHA", "TPM1", "MYLPF", "SLN", "TPM1", "FAM49A", "MYLK4", "TNNT3",
  "TBC1D1", "MYL1", "HCFC1R1", "MSS51", "TUBA8", "CASQ1") %>% cat(., sep = "\n")

article <- c(up, down)

###############################################


getwd()
setwd("/home/student/storage/AnTrD_personal_task/")

top_genes <- read.csv("./results/top_exon_level.csv")
top_genes_salmon <- read.csv("./results/top_genes_salmon.csv")
View(top_genes)

top_genes %>% 
  mutate(Direction = ifelse(log2FoldChange > 0, "Up", "Down")) %>%
  arrange(desc(log2FoldChange), padj) -> top_genes

top_genes %>% 
  filter(Direction == "Down") %>%
  nrow()

top_genes %>% 
  filter(!is.na(symbol)) %>%
  dplyr::select(symbol) %>%
  filter(Direction == "Down") %>%
  pull %>%cat(., sep = "\n")


# Объединяем в один именованный список
venn_list <- list(
  "Salmon + DESeq2" = genes_salmon,
  "FeatureCounts + DESeq2"  = genes,
  "article" = article)


install.packages("ggvenn")
library(ggvenn)

venna <- ggvenn(
  venn_list, 
  fill_color = c("red", "green", "blue"), # Цвета кругов
  stroke_size = 0.5,                         # Толщина границ
  set_name_size = 4,                         # Размер шрифта названий списков
  text_size = 4                              # Размер цифр внутри
)

shared_genes <- intersect(genes_salmon, article)
shared_genes <- intersect(intersect(genes, article), genes_salmon)
getwd()
ggsave("../results/venna.png", plot = venna, width = 10, height = 5)


top_genes_salmon %>%
  filter(symbol %in% shared_genes) %>%
  arrange(-stat) %>% 
  dplyr::select(symbol, log2FoldChange, stat, padj) %>% View






############################## ORA
BiocManager::install("clusterProfiler")
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.Hs.eg.db))

msig_hallmarks <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)

top_genes %>% 
  filter(!is.na(symbol)) %>%
  filter(Direction == "Down") %>%
  dplyr::select(symbol) %>%
  pull -> query_down


hallmark_down_res <- enricher(
  gene = query_down, 
  TERM2GENE = msig_hallmarks, 
  universe = as.vector(na.omit(Group_young_vs_old$symbol)))

dotplot(go_up_res, showCategory = 15, title = "Upregulated GO Biological Processes")
packageVersion("ggplot2") 

library(devtools)
install_version("ggplot2", version = "3.5.0", repos = "http://cran.r-project.org")

#############################################


library(decoupleR)
library(msigdbr)
library(tidyverse)

# Очищаем данные от NA в символах генов
df_clean <- Group_young_vs_old %>% 
  filter(!is.na(symbol)) %>% 
  distinct(symbol, .keep_all = TRUE)

msig_h <- msigdbr(species = "Homo sapiens", category = "H") %>%
  #dplyr::select(gs_name, gene_symbol) %>%
  #rename(source = gs_name, target = gene_symbol) %>%
  distinct(gs_name, gene_symbol)

colnames(msigdbr(species = "Homo sapiens", category = "H"))

# Создаем матрицу: 1 столбец со stat, строки — символы генов
mat <- as.matrix(df_clean$stat)
rownames(mat) <- df_clean$symbol
colnames(mat) <- "Group_young_vs_old"  # Название вашего контраста

# 2. Запускаем ORA
# n_up = 100 берет топ-100 генов с наибольшим stat (активированные)
# n_bottom = 100 берет топ-100 генов с наименьшим stat (подавленные)
ora_h_res <- run_ora(
  mat = mat,
  network = msig_h,
  .source = "gs_name",
  .target = "gene_symbol",
  n_up = 100,
  n_bottom = 100,
  n_background = nrow(mat), # наш очищенный universe
  minsize = 5)

# Смотрим результаты
ora_h_res

