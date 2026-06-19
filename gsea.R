####################################################GSEA
chooseCRANmirror()
BiocManager::install("fgsea")
install.packages('babelgene')
install.packages('https://cran.r-project.org/src/contrib/Archive/msigdbr/msigdbr_7.5.1.tar.gz',
                 repos=NULL, type = 'source')

suppressPackageStartupMessages(library(msigdbr))
suppressPackageStartupMessages(library(fgsea))
library(data.table)

getwd()
setwd("/home/student/storage/AnTrD_personal_task/")
Group_young_vs_old <- read.csv("./results/all_deg_gene_level.csv")

#Generating Ranks
p_spec <- "Homo_sapiens"
p_cat<-'H'
res <- Group_young_vs_old%>%arrange(-stat)
head(res)


Group_young_vs_old %>%
  filter(symbol %in% pathways$HALLMARK_OXIDATIVE_PHOSPHORYLATION)

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

fwrite(fgseaRes, file=paste0(gsub(' ','_',p_spec),'_',p_cat,'.tsv'), sep="\t", sep2=c("", " ", ""))

plot<-plotGseaTable(pathways[topPathways], ranks, fgseaRes, gseaParam = 0.5)
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
  
  
  png(file=paste0(p_spec,'_',pathway,'_GSEA.png'),1100,800,res=200)
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




