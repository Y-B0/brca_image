##enviroment
{
  library(do)
  library(gtools)
  library(stringr)
  library(magrittr)
  library(ComplexHeatmap)
  library(FactoMineR)
  library(factoextra)
  library(magrittr)
  library(glmnet)
  library(pROC)
  library(survminer)
  library(survival)
  library(corpcor)
  library(Hmisc)
  library(corrplot)
  library(ggsci)
  library(ROC)
  library(plotROC)
  library(lmQCM)
  library(forestplot)
  library(TCGAbiolinks)
  library(tidyverse)
  library(clusterProfiler) 
  library(RIdeogram)
  library(reshape2)
  library(ggpubr)
  library(gtools)
  library(XML)
  library(RCurl)
  library(httr)
  library(xml2)
  library(estimate)
  library(rms)
  library(ggbeeswarm)
  library(ConsensusClusterPlus)
  library(ggcorrplot)
  library(psych)
  source("D:/project/my_function/go_demo.R")
}

##basic data process
{
  # TCGA.BRCA.htseq_counts.tsv.gz <- read.delim("rawdata/TCGA-BRCA.htseq_counts.tsv.gz.tsv")
  # gencode.v22.annotation <- read.delim("rawdata/gencode.v22.annotation.gene")
  # TCGA.BRCA.htseq_counts.tsv.gz<-TCGA.BRCA.htseq_counts.tsv.gz[match(gencode.v22.annotation$id,TCGA.BRCA.htseq_counts.tsv.gz$Ensembl_ID),]
  # TCGA.BRCA.htseq_counts.tsv.gz$Ensembl_ID<-gencode.v22.annotation$gene
  # BRCA.im.select <- read.delim("BRCA.im.select.txt", header=T,row.names = 1)
  # #
  # colnames(TCGA.BRCA.htseq_counts.tsv.gz)[1]<-"Gene"
  # TCGA.BRCA.htseq_counts.tsv.gz<-TCGA.BRCA.htseq_counts.tsv.gz[!duplicated(TCGA.BRCA.htseq_counts.tsv.gz$Gene),]
  # colnames(TCGA.BRCA.htseq_counts.tsv.gz)<-str_sub(colnames(TCGA.BRCA.htseq_counts.tsv.gz),1,15)
  # TCGA.BRCA.htseq_counts.tsv.gz<-TCGA.BRCA.htseq_counts.tsv.gz[,c("Gene",intersect(colnames(BRCA.im.select),colnames(TCGA.BRCA.htseq_counts.tsv.gz)))]
  # rownames(TCGA.BRCA.htseq_counts.tsv.gz)<-TCGA.BRCA.htseq_counts.tsv.gz$Gene
  # TCGA.BRCA.htseq_counts.tsv.gz<-TCGA.BRCA.htseq_counts.tsv.gz[,-1]
  # TCGA.BRCA.htseq_counts.tsv.gz<-t(TCGA.BRCA.htseq_counts.tsv.gz)
  # TCGA.BRCA.htseq_counts.tsv.gz<-scale(TCGA.BRCA.htseq_counts.tsv.gz,center = T,scale = T)
  # rownames(TCGA.BRCA.htseq_counts.tsv.gz)<-str_sub(rownames(TCGA.BRCA.htseq_counts.tsv.gz),1,12)
  # save(TCGA.BRCA.htseq_counts.tsv.gz,file = "raw_genecount.RData")
  # protein_coding.hg38 <- read.delim("protein_coding.hg38.position", header=FALSE)
  # TCGA.BRCA.htseq_counts.tsv.gz<-TCGA.BRCA.htseq_counts.tsv.gz[TCGA.BRCA.htseq_counts.tsv.gz$Gene%in%protein_coding.hg38$V4,]
  # TCGA.BRCA.htseq_counts.tsv.gz<-TCGA.BRCA.htseq_counts.tsv.gz[rowSums(TCGA.BRCA.htseq_counts.tsv.gz[,-1]>=1)>=ncol(TCGA.BRCA.htseq_counts.tsv.gz[,-1])/2,]
  # write.table(TCGA.BRCA.htseq_counts.tsv.gz,file = "genecount.txt",sep = "\t",col.names = T,row.names = F,quote = F)
  # 
  
  gene <- read.csv("coexp_modules.csv", header=F,row.names = 1)%>%t()
  colnames(gene)<-paste("Module",1:ncol(gene),sep = "")
  
  ##pca read in
  pca <- read.csv("coexp_eigengene_matrix.csv",row.names = 1) %>%t(.)
  colnames(pca)<-paste("Module",1:ncol(pca),sep = "")
  rownames(pca)<-str_sub(rownames(pca),1,12)
  ###correcte process
  pca[,-c(2,4,5,6,7,9,10,11,12,18)]<--pca[,-c(2,4,5,6,7,9,10,11,12,18)]
  
  ##mRNA read in
  BRCA.mrna.select <- read.csv("datamatrix_2.csv",row.names = 1)
  BRCA.mrna.select<-scale(BRCA.mrna.select,center = T,scale = T)
  BRCA.mrna.select<-as.data.frame(BRCA.mrna.select)
  colnames(BRCA.mrna.select)<-str_sub(colnames(BRCA.mrna.select),1,12)
  BRCA.mrna.select<-t(BRCA.mrna.select)
  
  ##image feature read in
  BRCA.im.select <- read.delim("BRCA.im.select.txt", header=T,row.names = 1)
  colnames(BRCA.im.select)<-str_sub(colnames(BRCA.im.select),1,12)
  BRCA.im.select<-t(BRCA.im.select)%>%as.data.frame(.)
  BRCA.im.select<-BRCA.im.select[,grep("^Rm|^Gm|^Bm",grep("1|5",colnames(BRCA.im.select),value = T),value = T,invert = T)]
  BRCA.im.select<-scale(BRCA.im.select,center = T,scale = T)
  
  ##clinic  
  clinic <- read.delim("rawdata/TCGA-BRCA.survival.tsv")
  clinic$X_PATIENT<-gsub("-",".",clinic$X_PATIENT)
  colnames(clinic)<-c("sample","vital_status","sample1","days_to_last_follow_up")
  clinic<-clinic[clinic$days_to_last_follow_up<3000,]
  clinic<-clinic[!duplicated(clinic$sample1),]
  
  ##subtype
  subtype<-as.data.frame(TCGAquery_subtype(tumor = "brca"))
  subtype<-subtype[,c("patient","BRCA_Subtype_PAM50")]
  subtype$patient<-gsub("-",".",subtype$patient)
  
  ##stage
  load("basic_data.RData")
  clinical <- read.delim("rawdata/clinical.tsv")
  clinical<-clinical[,c("case_submitter_id","ajcc_pathologic_stage","age_at_index","ethnicity","gender","race")]
  clinical<-clinical[!duplicated(clinical),]
  clinical$ajcc_pathologic_stage<-Replace(clinical$ajcc_pathologic_stage,pattern = c("Stage IA:Stage I","Stage IB:Stage I","Stage IIA:Stage II","Stage IIB:Stage II","Stage IIIA:Stage III","Stage IIIB:Stage III","Stage IIIC:Stage III"))
  clinical$case_submitter_id<-gsub("-","\\.",clinical$case_submitter_id)
  clinical<-clinical[match(clinic$sample1,clinical$case_submitter_id),]
  other_clinic<-clinical;rm(clinical)
  clinic$sample<-clinical$ajcc_pathologic_stage
  colnames(clinic)[1]<-"Stage"
}
save.image("basic_data.RData")

## cor map (image,model,tumor purity)
{
  #load("raw_genecount.RData")
  #TCGA.BRCA.htseq_counts.tsv.gz<-t(TCGA.BRCA.htseq_counts.tsv.gz)
  #write.table(TCGA.BRCA.htseq_counts.tsv.gz,file = "expr.txt",row.names = T,col.names = T,sep = "\t",quote = F)
  #filterCommonGenes(input.f = "expr.txt",output.f = "Estimate_gene.gct",id = "GeneSymbol")
  #estimateScore(input.ds = "Estimate_gene.gct",output.ds = "Estimate_score.gct",platform = "illumina")
  scores=read.table("Estimate_score.gct",skip = 2,header = T)
  rownames(scores)=scores[,1]
  scores=t(scores[,3:ncol(scores)])%>%as.data.frame(.)
  scores$TumorPurity=cos(0.6049872018+0.0001467884*scores$ESTIMATEScore)
  scores<-scores[,-3]
  rownames(scores)<-str_sub(rownames(scores),1,12)
  load("basic_data.RData")
  BRCA.im.select<-BRCA.im.select[,paste(c("Area","Major","Minor","Ratio","DistMin","DistMax","DistMean"),rep(c(1,5),each=7),sep = "")]
  colnames(BRCA.im.select)<-c("Small_Nucleus_Area","Small_Major_Axis","Small_Minor_Axis","Small_Aspect_Ratio",
                              "Small_Min_Distance","Small_Max_Distance","Small_Mean_Distance",
                              "Large_Nucleus_Area","Large_Major_Axis","Large_Minor_Axis",
                              "Large_Aspect_Ratio","Large_Min_Distance","Large_Max_Distance","Large_Mean_Distance")
  BRCA.im.select<-BRCA.im.select[match(rownames(BRCA.mrna.select),rownames(BRCA.im.select)),]
  ##total18:18
  cor<-cbind(BRCA.im.select,pca,scores)
  cor<-rcorr(as.matrix(cor),type = "spearman")
  
  ggcorrplot(cor$r[15:63,rev(c(1:14,61:63))],p.mat = cor$p[15:63,rev(c(1:14,61:63))],hc.order = F,outline.color = "white",
             ggtheme = theme_bw(),sig.level = 0.05,insig = "blank",lab_size = 4,digits = 1,
             colors = c("#6D9EC1","white","#E46726"),lab = T)+
    theme(axis.text.y = element_text(color = "red"),axis.text.x = element_text(color = "red"),axis.text.x.top = element_text(hjust = 0,vjust = 0))+
    scale_x_discrete(position = "top")
}

## enrich
{
  load("basic_data.RData")
  enrich_list<-list()
  cytoband_list<-list()
  gene<-as.data.frame(gene)
  cytoband_gmt <- read.gmt("D:/project/msigdb_v7.4_files_to_download_locally/msigdb_v7.4_GMTs/c1.all.v7.4.entrez.gmt")
  for (i in colnames(gene)) {
    key<-i
    genelist<-gene[,key]
    genelist<-bitr(genelist,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
    try({
      #dir.create(path=paste("enrich/",key,sep = ""),recursive = T)
      setwd(paste("enrich/",key,sep = ""))
      #tmp<-go_demo(genesymbol = genelist$SYMBOL,plot = F)
      #enrich_list[[key]]<-tmp
      cytoband<-enricher(genelist$ENTREZID,TERM2GENE = cytoband_gmt)
      cytoband_list[[key]]<-cytoband@result
      write.table(as.data.frame(cytoband@result),file = "cytoband.txt",row.names = F,col.names = T,quote = F,sep = "\t")
    })
    setwd("../../")
  }
  save(enrich_list,file="enrich/enrich_model.RData")
  save(cytoband_list,file = "enrich/cytoband_module.RData")
  ##process enrich result
  load("enrich/enrich_model.RData")
  enrich_list<-enrich_list[c("Module3","Module4","Module7","Module11","Module31")]
  for (i in 1:length(enrich_list)) {
    enrich_list[[i]]$method<-names(enrich_list[i])
    enrich_list[[i]]$model<-names(enrich_list[i])
    enrich_list[[i]]<-enrich_list[[i]][enrich_list[[i]]$ONTOLOGY=="BP",]
    enrich_list[[i]]<-enrich_list[[i]][enrich_list[[i]]$pvalue<0.05&enrich_list[[i]]$qvalue<0.05,]
    enrich_list[[i]]<-head(enrich_list[[i]],n=20)
  }
  enrich_list<-Reduce(rbind,enrich_list)
  ##go 10:10
  enrich_list_go<-enrich_list
  enrich_list_go$method<-Replace(enrich_list_go$method,pattern = c("Module31:Immune system","Module3:Immune system","Module4:Extracellular matrix","Module7:Immune system","Module11:Vascular development"))
  enrich_list_go$method<-factor(enrich_list_go$method,levels = c("Immune system","Extracellular matrix","Vascular development","Antigen presentation"))
  enrich_list_go$model<-factor(enrich_list_go$model,levels = unique(enrich_list_go$model))
  enrich_list_go<-enrich_list_go[order(enrich_list_go$method,enrich_list_go$model),]
  
  enrich_list_go$Description<-factor(enrich_list_go$Description,level=unique(enrich_list_go$Description))
  enrich_list_go$ratio<-unlist(lapply(enrich_list_go$GeneRatio, function(x){eval(parse(text = x))}))
  
  ggplot()+geom_point(aes(colour=enrich_list_go$p.adjust,x=enrich_list_go$model,y = enrich_list_go$Description,size=enrich_list_go$ratio))+theme_bw()+
    xlab("")+ylab("")+scale_y_discrete(labels=function(x)str_sub(x,1,70))+scale_colour_gsea(name="P-value")+scale_size(name="GeneRatio")+
    theme(text = element_text(size = 12),axis.text.y = element_text(color = rep(c("#1F77B4B2","#FF7F0EB2","#2CA02CB2"),time=c(49,20,20))))
  
  ##cytoband 
  cytoband<-read.delim("cytoBand.txt", header=FALSE)%>%as.data.frame(.)
  cytoband$V4<-str_split(cytoband$V4,"\\.",simplify = T)[,1]
  cytoband$key<-paste(cytoband$V1,cytoband$V4,sep = "")
  cytoband<-data.frame(chr=tapply(cytoband[,1], cytoband$key,unique),start=tapply(cytoband[,2], cytoband$key, min),end=tapply(cytoband[,3], cytoband$key, max))
  
  load("enrich/cytoband_module.RData")
  enrich_list_cytoband<-cytoband_list
  enrich_list_cytoband<-enrich_list_cytoband[c("Module3","Module4","Module7","Module11","Module31")]
  for (i in 1:length(enrich_list_cytoband)) {
    enrich_list_cytoband[[i]]$method<-names(enrich_list_cytoband[i])
    enrich_list_cytoband[[i]]$model<-names(enrich_list_cytoband[i])
    enrich_list_cytoband[[i]]<-enrich_list_cytoband[[i]][enrich_list_cytoband[[i]]$pvalue<0.05&enrich_list_cytoband[[i]]$qvalue<0.05,]
  }
  enrich_list_cytoband<-Reduce(rbind,enrich_list_cytoband)
  enrich_list_cytoband<-enrich_list_cytoband[mixedorder(enrich_list_cytoband$model,enrich_list_cytoband$ID,decreasing = F),]
  cytoband<-cytoband[match(enrich_list_cytoband$ID,rownames(cytoband)),]
  enrich_list_cytoband<-cbind(enrich_list_cytoband,cytoband)
  cp<-enrich_list_cytoband
  enrich_list_cytoband<-split(enrich_list_cytoband,enrich_list_cytoband$ID)%>%lapply(.,function(x){x$Count<-max(x$Count);x})%>%Reduce(rbind,.)
  ##plot
  enrich_list_cytoband$chr<-gsub("chr","",enrich_list_cytoband$chr)
  enrich_list_cytoband$Shape<-enrich_list_cytoband$model%>%Replace(.,pattern = c("Module31:circle","Module3:circle","Module4:circle","Module7:circle"))
  enrich_list_cytoband$color<-enrich_list_cytoband$model%>%Replace(.,pattern = c("Module31:ff7f00","Module3:0093df","Module4:33a02c","Module7:6a3d9a"))
  label<-enrich_list_cytoband[,c("model","Shape","chr","start","end","color")]
  label$model<-factor(label$model,levels = c("Module3","Module4","Module7","Module31"))
  colnames(label)<-c("Type","Shape","Chr","Start","End","color")
  overlaid<-enrich_list_cytoband[,c("chr","start","end","GeneRatio")]
  colnames(overlaid)<-c("Chr","Start","End","Value")
  overlaid$Value<-unlist(lapply(overlaid$Value, function(x){eval(parse(text = x))}))
  label$Type<-factor(label$Type,levels = c("Module3","Module4","Module7","Module31"))
  label<-label[order(label$Type),]
  ideogram(karyotype = human_karyotype,label = label,label_type = "marker",overlaid = overlaid)
  convertSVG("chromosome.svg",file = "figure/cytoband.pdf",device = "pdf",width=9,height = 12)
  #"1:8ECFC9","2:FFBE7A","3:FA7F6F","4:82B0D2","5:BEB8DC"
  
  ggplot()+geom_point(aes(x=enrich_list_cytoband$model,y = enrich_list_cytoband$Description,color=enrich_list_cytoband$model,size=enrich_list_cytoband$Count))+theme_bw()+
    xlab("")+ylab("")+scale_color_d3(name="")+
    theme(text = element_text(size = 12),axis.text.y = element_text(colour = rep(c("#1F77B4B2","#FF7F0EB2","#2CA02CB2"),each=20)))
  
}

## TF plot 16:5
{
  tf_list<-list()
  for (model in c("module3","module4","module7","module11","module31")) {
    TFs <- read.delim(paste("tf/",model,".tsv",sep = ""))%>%head(.,10)
    TFs$Library<-model
    tf_list[[model]]<-separate_rows(TFs, Overlapping_Genes,sep = ",")%>%as.data.frame(.)
  }
  tf_list<-Reduce(rbind,tf_list)
  tf_list$Library<-gsub("m","M",tf_list$Library)
  tmp<-as.data.frame(table(tf_list$TF,tf_list$Library))%>%.[.$Freq!=0,]
  tmp$Var2<-factor(tmp$Var2,level=c("Module3","Module4","Module7","Module11","Module31"))
  tmp<-tmp[mixedorder(tmp$Var2),]
  tmp$proportion<-tmp$Freq/rep(c(166,139,85,54,17),each=10);tmp$proportion<-round(tmp$proportion*100,0)
  tmp$Var1<-as.character(tmp$Var1);tmp$Var2<-as.character(tmp$Var2)
  tmp$proportion<-tmp$proportion*0.01
  tmp[tmp$Var1%in%tmp$Var1[duplicated(tmp$Var1)],"Var1"]<-paste(tmp[tmp$Var1%in%tmp$Var1[duplicated(tmp$Var1)],"Var2"],tmp[tmp$Var1%in%tmp$Var1[duplicated(tmp$Var1)],"Var1"],sep = "_")
  tmp$Var2<-factor(tmp$Var2,level=c("Module3","Module4","Module7","Module11","Module31"))
  tmp<-tmp[order(tmp$Var2,tmp$proportion,decreasing = F),]
  tmp$Var1<-factor(tmp$Var1,levels = tmp$Var1)
  ggplot(tmp)+geom_bar(aes(x=Var1,fill=Var2,y=proportion),position = "dodge",width = 0.6,stat = "identity")+theme_bw()+ylab("Rate")+xlab("")+
    theme(text=element_text(size=20),axis.text.x = element_text(colour = rep(pal_futurama()(5),each=10),angle = 45,vjust = 1,hjust = 1))+scale_fill_futurama(name="")+
    geom_text(aes(x = Var1, y= 1,label = Freq),size=3.5,check_overlap = T)
}

## hubgene and TF expression cor with image feature
{
  load("basic_data.RData")
  load("raw_genecount.RData")
  BRCA.mrna.select<-BRCA.mrna.select[intersect(clinic$sample1,rownames(BRCA.mrna.select)),]
  BRCA.im.select<-BRCA.im.select[match(rownames(BRCA.mrna.select),rownames(BRCA.im.select)),paste(c("Area","Major","Minor","Ratio","DistMin","DistMax","DistMean"),rep(c(1,5),each=7),sep = "")]
  BRCA.im.select<-BRCA.im.select[intersect(clinic$sample1,rownames(BRCA.mrna.select)),paste(c("Area","Major","Minor","Ratio","DistMin","DistMax","DistMean"),rep(c(1,5),each=7),sep = "")]
  colnames(BRCA.im.select)<-c("Small_Nucleus_Area","Small_Major_Axis","Small_Minor_Axis","Small_Aspect_Ratio",
                              "Small_Min_Distance","Small_Max_Distance","Small_Mean_Distance",
                              "Large_Nucleus_Area","Large_Major_Axis","Large_Minor_Axis",
                              "Large_Aspect_Ratio","Large_Min_Distance","Large_Max_Distance","Large_Mean_Distance")
  
  tmp<-list()
  for (model in c("module3","module4","module7","module11","module31")) {
    TFs <- read.delim(paste("tf/",model,".tsv",sep = ""))%>%head(.,10)
    TFs$Library<-"TF"
    hubgene<-read.table(paste("ppi/",model,"_cytohubba.csv",sep = ""),skip = 1,sep = ",",header = T)[,2]%>%as.data.frame(.)
    hubgene$Library<-"Gene"
    
    tmp[[model]]<-data.frame(name=c(hubgene$.,TFs$TF),type=c(hubgene$Library,TFs$Library),model=model)
    g<-tmp[[model]]
    g<-g[order(g$type),]
    data<-TCGA.BRCA.htseq_counts.tsv.gz[rownames(BRCA.mrna.select),g$name]
    
    ###hubgenes & TFs    10:7.8
    {
      test1<-cbind(TCGA.BRCA.htseq_counts.tsv.gz[rownames(BRCA.mrna.select),g$name[g$type=="TF"]],BRCA.im.select)
      test2<-cbind(TCGA.BRCA.htseq_counts.tsv.gz[rownames(BRCA.mrna.select),g$name[g$type=="Gene"]],BRCA.im.select)
      cor<-corr.test(as.matrix(test1),as.matrix(test2[,rev(colnames(test2))]),method = "spearman")
      p<-ggcorrplot(cor$r,p.mat = cor$p,hc.order = F,outline.color = "white",
                    ggtheme = theme_bw(),sig.level = 0.05,insig = "blank",lab_size = 3,digits = 1,
                    colors = c("#6D9EC1","white","#E46726"),lab = T)+
        theme(axis.text.x = element_text(color = rep(c("#C24294","red"),c(10,14))),axis.text.y = element_text(color = rep(c("red","#4EADEA"),c(14,10))),axis.text.x.top = element_text(hjust = 0,vjust = 0))+
        scale_x_discrete(position = "top")
      ggsave(paste("figure/",model,"_image_gene_heatmap.pdf",sep = ""),p,width = 10,height = 7.8)
    }
  }
  ##all gene cor
  {
    load("basic_data.RData")
    load("raw_genecount.RData")
    BRCA.mrna.select<-BRCA.mrna.select[intersect(clinic$sample1,rownames(BRCA.mrna.select)),]
    BRCA.im.select<-BRCA.im.select[match(rownames(BRCA.mrna.select),rownames(BRCA.im.select)),paste(c("Area","Major","Minor","Ratio","DistMin","DistMax","DistMean"),rep(c(1,5),each=7),sep = "")]
    BRCA.im.select<-BRCA.im.select[intersect(clinic$sample1,rownames(BRCA.mrna.select)),paste(c("Area","Major","Minor","Ratio","DistMin","DistMax","DistMean"),rep(c(1,5),each=7),sep = "")]
    colnames(BRCA.im.select)<-c("Small_Nucleus_Area","Small_Major_Axis","Small_Minor_Axis","Small_Aspect_Ratio",
                                "Small_Min_Distance","Small_Max_Distance","Small_Mean_Distance",
                                "Large_Nucleus_Area","Large_Major_Axis","Large_Minor_Axis",
                                "Large_Aspect_Ratio","Large_Min_Distance","Large_Max_Distance","Large_Mean_Distance")
    
    coexp_modules <- read.csv("coexp_modules.csv", header=FALSE,row.names = 1)%>%t(.)%>%as.data.frame(.)
    colnames(coexp_modules)<-paste("module",1:46,sep = "")
    for (model in colnames(coexp_modules)) {
      
      g<-coexp_modules[,model]
      g<-g[!(is.na(g) | g=="")]
      data<-TCGA.BRCA.htseq_counts.tsv.gz[rownames(BRCA.mrna.select),g]
      
      test1<-BRCA.im.select
      test2<-TCGA.BRCA.htseq_counts.tsv.gz[rownames(BRCA.mrna.select),head(g,50)]
      cor<-corr.test(as.matrix(test2),as.matrix(test1[,rev(colnames(test1))]),method = "spearman")
      p<-ggcorrplot(cor$r,hc.order = F,outline.color = "white",
                    ggtheme = theme_bw(),sig.level = 0.05,insig = "blank",
                    colors = c("#6D9EC1","white","#E46726"))+
        theme(axis.text.y = element_text(color = "red"),axis.text.x = element_text(color = "#4EADEA"),axis.text.x.top = element_text(hjust = 0,vjust = 0))+
        scale_x_discrete(position = "top")
      ggsave(paste("figure/all_gene_cor/",model,"_image_gene_heatmap.pdf",sep = ""),p,width = 12,height = 10)
      
    }
  }
}

## regulator network
{
  for (model in c("module3","module4","module7","module11","module31")) {
    library(tidyr)
    TFs <- read.delim(paste("tf/",model,".tsv",sep = ""))%>%.[,c("TF","Overlapping_Genes")]%>%head(.,10)
    TFs<-separate_rows(TFs, Overlapping_Genes,sep = ",")%>%as.data.frame(.)
    TFs$geneType<-"TF"
    TFs<-as.data.frame(TFs)
    write.table(TFs,file = paste("tf/",model,"_regulate_network_allgene.txt",sep = ""),quote = F,sep = "\t",row.names = F,col.names = T)
    type<-data.frame(node=c(unique(TFs$TF[!(TFs$TF%in%TFs$Overlapping_Genes)]),
                            unique(TFs$Overlapping_Genes[!(TFs$Overlapping_Genes%in%TFs$TF)]),
                            intersect(TFs$TF,TFs$Overlapping_Genes)),
                     type=c(rep("TF",length(unique(TFs$TF[!(TFs$TF%in%TFs$Overlapping_Genes)]))),
                            rep("gene",length(unique(TFs$Overlapping_Genes[!(TFs$Overlapping_Genes%in%TFs$TF)]))),
                            rep("Both",length(intersect(TFs$TF,TFs$Overlapping_Genes)))))
    write.table(type,file = paste("tf/",model,"_regulate_network_type_allgene.txt",sep = ""),quote = F,sep = "\t",row.names = F,col.names = T)
    
  }

}

##predict model
{
  setwd("D:/project/brca5/")
  
  load("basic_data.RData")
  load("raw_genecount.RData")
  tmp<-list()
  for (model in c("module3","module4","module7","module11","module31")) {
    TFs <- read.delim(paste("tf/",model,".tsv",sep = ""))%>%head(.,10)
    TFs$Library<-"TF"
    hubgene<-read.table(paste("ppi/",model,"_cytohubba.csv",sep = ""),skip = 1,sep = ",",header = T)[,2]%>%as.data.frame(.)
    hubgene$Library<-"Gene"
    
    tmp[[model]]<-data.frame(name=c(hubgene$.,TFs$TF),type=c(hubgene$Library,TFs$Library),model=model)
  }
  tmp<-Reduce(rbind,tmp)
  
  clinic<-clinic[clinic$sample1%in%intersect(clinic$sample1,rownames(BRCA.mrna.select)),]
  BRCA.mrna.select<-TCGA.BRCA.htseq_counts.tsv.gz;rm(TCGA.BRCA.htseq_counts.tsv.gz)
  BRCA.im.select<-BRCA.im.select[match(clinic$sample1,rownames(BRCA.im.select)),c("DistMax1","DistMin1","DistMean1","Area5","Major5","Minor5","Ratio5","DistMax5","DistMin5","DistMean5")]
  BRCA.mrna.select<-BRCA.mrna.select[match(clinic$sample1,rownames(BRCA.mrna.select)),unique(tmp$name)]
  
  
  comb<-lapply(1:5, function(x){a<-combn(c("module3","module4","module7","module11","module31"),m = x,simplify = T)})
  
  method<-lapply(1:3, function(x){a<-combn(c("Gene","TF","image"),m = x,simplify = F)})
  method<-Reduce(c,method)
  
  
  years<-list()
  roc<-list()
  roc_total<-data.frame()
  coef_total<-list()
  total<-list()
  plot<-list()
  raw_data<-list()
  riskscore_plot<-list()
  surv_plot<-list()
  coef<-list()
  riskscore_data<-list()
  model_list<-list()
  for (y in 365*c(10)) {
    for (k in 1:(length(method))) {
      for (i in 1:length(comb)) {
        for (j in 1:ncol(comb[[i]])) {
          m<-comb[[i]][,j]
          ###select TF|Gene
          
          if ("image"==method[[k]][1]) {
            data<-cbind(BRCA.im.select,clinic[,c("vital_status","days_to_last_follow_up")])
          }else if (("image"%in%method[[k]])) {
            matrix<-BRCA.mrna.select[,intersect(c(tmp$name[tmp$model%in%m&tmp$type%in%method[[k]]]),colnames(BRCA.mrna.select))]
            data<-cbind(matrix,BRCA.im.select,clinic[,c("vital_status","days_to_last_follow_up")])
          }else{
            matrix<-BRCA.mrna.select[,intersect(c(tmp$name[tmp$model%in%m&tmp$type%in%method[[k]]]),colnames(BRCA.mrna.select))]
            data<-cbind(matrix,clinic[,c("vital_status","days_to_last_follow_up")])
          }
          ### survival days control
          data<-data[data$days_to_last_follow_up<y,]
          ###module train 
          alpha_fit<-cv.glmnet(data.matrix(data[,1:(ncol(data)-2)]),data$vital_status,alpha = 1,family = "binomial", type.measure="auc")
          model_list[[as.character(y)]][[paste(m,collapse = "_")]][[paste(method[[k]],collapse = "_")]]<-alpha_fit
          coef[[paste(m,collapse = "&")]]<-as.matrix(coef(alpha_fit,s=alpha_fit$lambda.min))
          
          riskscore<-data.frame(score=(as.matrix(data[,1:(ncol(data)-2)])%*%coef[[paste(m,collapse = "&")]][2:nrow(coef[[paste(m,collapse = "&")]]),]),sample=rownames(data))
          riskscore<-riskscore[order(riskscore$score),]
          riskscore$score<-riskscore$score-median(riskscore$score)
          riskscore$sample<-factor(riskscore$sample,levels = riskscore$sample)
          riskscore$plot_group<-ifelse(riskscore$score>0,"High","Low")
          
          
          clinic1<-clinic[match(riskscore$sample,clinic$sample1),]
          riskscore<-cbind(riskscore,clinic1[,c("vital_status","days_to_last_follow_up")])
          cutoff<-median(riskscore$score)
          
          
          riskscore$group<-ifelse(riskscore[,1]>=cutoff,"High","Low")
          riskscore$plot_group<-factor(riskscore$plot_group,levels = c("Low","High"))
          riskscore_plot[[as.character(y)]][[paste(m,collapse = "_")]][[paste(method[[k]],collapse = "_")]]<-ggplot()+geom_point(aes(y=riskscore$score,x=riskscore$sample,color=riskscore$plot_group))+theme_bw()+
            theme(text = element_text(size=20),axis.text.x=element_blank(),axis.ticks.x=element_blank())+xlab("Samples")+ylab("Riskscore")+scale_color_aaas(name="Group")
          
          
          riskscore<-riskscore[order(riskscore$score),]
          riskscore_data[[as.character(y)]][[paste(m,collapse = "_")]][[paste(method[[k]],collapse = "_")]]<-riskscore
          sfit1<-survfit(Surv(days_to_last_follow_up,
                              vital_status)~group,data = riskscore)
          surv_plot[[as.character(y)]][[paste(m,collapse = "_")]][[paste(method[[k]],collapse = "_")]]<-ggsurvplot(sfit1,
                                                                                                                   data = riskscore,
                                                                                                                   legend.title = NULL,
                                                                                                                   conf.int = T,
                                                                                                                   pval = T,
                                                                                                                   title = "Riskscore",
                                                                                                                   risk.table = T,
                                                                                                                   palette = c("#EE0000FF","#3B4992FF"),
                                                                                                                   risk.table.col = 'strata',
                                                                                                                   tables.height = 0.3,
                                                                                                                   xlab = "Time In Days",
                                                                                                                   pval.method = TRUE,
                                                                                                                   surv.median.line = "hv")
          
          
          data<-data[match(riskscore$sample,rownames(data)),]
          data$riskscore.score<-riskscore$score
          
          #predict
          basicplot<-ggplot(data,aes(d = vital_status, m = riskscore.score)) + geom_roc()+style_roc()
          auc<-calc_auc(basicplot)$AUC
          roc[[paste(m,collapse = "&")]]<-c(paste(m,collapse = "&"),auc)
          roc_total[(nrow(roc_total)+1),1:3]<-c(paste(m,collapse = "_"),paste(method[[k]],collapse = "_"),auc)
          plot[[as.character(y)]][[paste(m,collapse = "_")]][[paste(method[[k]],collapse = "_")]]<-basicplot+annotate("text",x = .75, y = .25,
                                                                                                                      label = paste("AUC =", round(auc, 2)))+theme(text = element_text(size=20))
          raw_data[[as.character(y)]][[paste(m,collapse = "_")]][[paste(method[[k]],collapse = "_")]]<-data
          coef_total[[as.character(y)]][[paste(m,collapse = "_")]][[paste(method[[k]],collapse = "_")]]<-as.matrix(coef(alpha_fit,s=alpha_fit$lambda.min))
        }
      }
      total[[k]]<-as.data.frame(t(as.data.frame(roc)))
    }
    test<-Reduce(cbind,total)%>%.[,colnames(.)=="V2"]
    colnames(test)<-unlist(lapply(method,function(x){x<-paste(x,collapse = "_")}))
    years[[as.character(y)]]<-as.data.frame(test)
  }
  save(riskscore,file = "riskscore.RData")
  load("alter_predict_final.RData")
  cof<-data.frame(coef_total[["3650"]][["module3_module7_module11_module31"]][["Gene_TF_image"]])
  cof$gene<-rownames(cof)
  cof<-cof[cof$s1!=0,]
  cof<-cof[order(cof$s1,decreasing = T),]
  #write.table(years[["3650"]],file = "Table2_auc.xls",quote = F,sep = "\t",row.names = T,col.names = T)
  #write.table(cof,file = "Table3_coef.xls",quote = F,sep = "\t",row.names = F,col.names = T)
  
  riskscore_plot[["3650"]][["module3_module7_module11_module31"]][["Gene_TF_image"]]
  
  surv_plot[["3650"]][["module3_module7_module11_module31"]][["Gene_TF_image"]]
  
  plot[["3650"]][["module3_module7_module11_module31"]][["Gene_TF_image"]]
  
  riskscore<-riskscore_data[["3650"]][["module3_module7_module11_module31"]][["Gene_TF_image"]]
  save(riskscore,file = "riskscore.RData")
  rm();gc()
}

##after riskscore
{
  
  load("riskscore.RData")
  load("basic_data.RData")
  BRCA.im.select<-BRCA.im.select[,paste(c("Area","Major","Minor","Ratio","DistMin","DistMax","DistMean"),rep(c(1,5),each=7),sep = "")]
  colnames(BRCA.im.select)<-c("Small_Nucleus_Area","Small_Major_Axis","Small_Minor_Axis","Small_Aspect_Ratio",
                              "Small_Min_Distance","Small_Max_Distance","Small_Mean_Distance",
                              "Large_Nucleus_Area","Large_Major_Axis","Large_Minor_Axis",
                              "Large_Aspect_Ratio","Large_Min_Distance","Large_Max_Distance","Large_Mean_Distance")
  riskscore<-riskscore[order(riskscore$group),]
  riskscore_cp<-riskscore[order(riskscore$group,decreasing = T),]
  riskscore_cp<-riskscore_cp[c(1:300,643:942),]
  BRCA.im.select<-BRCA.im.select[match(riskscore$sample,rownames(BRCA.im.select)),]
  BRCA.im.select_cp<-BRCA.im.select[match(riskscore_cp$sample,rownames(BRCA.im.select)),]
  ## feature wilcox.test
  test<-cbind(riskscore,BRCA.im.select)
  test<-test[,6:20]
  tmp<-melt(test,id.vars = "group",variable.name = "Feature",value.name = "Value")
  tmp$group<-factor(tmp$group,levels = c("Low","High"))
  #12.55:6
  ggviolin(tmp,x="Feature",y="Value",
           color = "group", 
           fill="NA",
           palette ="jco", 
           alpha=0.7, 
           linetype=1,
           trim=F,
           size=1,
           width=1, 
           draw_quantiles = 0.5, 
           add=c("jitter","mean"), 
           add.params = list(color="group",alpha=0.6),
           error.plot = "errorbar",
           position=position_dodge(1),
           xlab="Morphological features",
           ylab="Value",
           title="The distribution of morphological features",
           ggtheme=theme_bw())+ 
    stat_compare_means(aes(group = group),method = "wilcox.test",label="p.signif")+
    theme(plot.title = element_text(hjust = 0.5),text = element_text(size=20),axis.text.x = element_text(angle = 90,hjust = 1,vjust = 1))+
    scale_color_aaas(name="Riskscore group")
  
  other_clinic_cp<-other_clinic
  subtype<-subtype[match(other_clinic$case_submitter_id,subtype$patient),]
  subtype$BRCA_Subtype_PAM50<-gsub("Normal","Normal like",subtype$BRCA_Subtype_PAM50)
  other_clinic_cp$ethnicity<-subtype$BRCA_Subtype_PAM50
  colnames(other_clinic_cp)[4]<-"Subtype"
  other_clinic<-other_clinic[match(riskscore_cp$sample,other_clinic$case_submitter_id),]
  subtype_cp<-subtype
  subtype<-subtype[match(riskscore_cp$sample,subtype$patient),]
  
  group_row <- factor(rep(c("Low","High"),times=c(300,300)))
  group_col <- factor(rep(c("Small","Large"),times=c(7,7)))
  top_anno<-HeatmapAnnotation(Group=riskscore_cp$group,Subtype=subtype$BRCA_Subtype_PAM50,Stage=other_clinic$ajcc_pathologic_stage,
                              col = list(Subtype=c("Basal"="#2F7FC1","Her2"="#96C37D",
                                                   "LumA"="#F3D266","LumB"="#D8383A","Normal like"="#A9B8C6"),
                                         Group=c("High"="#EE0000FF","Low"="#3B4992FF"),
                                         Stage=c("Stage I"="#2F7FC1","Stage II"="#96C37D","Stage III"="#F3D266","Stage IV"="#D8383A")
                              ),which = "row")
  #heatmap 8:5
  Heatmap(as.matrix(BRCA.im.select_cp),show_row_names = F,show_column_names = T,cluster_columns = T,cluster_rows = F,
          left_annotation = top_anno,row_split = group_row,column_split = group_col, name = "Value",
          col = circlize::colorRamp2(c(-4,-1,0,1,4),c('#008837','#a6dba0','#f7f7f7','#c2a5cf','#7b3294')))
  
  
  other_clinic_cp<-other_clinic_cp[match(riskscore$sample,other_clinic_cp$case_submitter_id),]
  data<-merge(other_clinic_cp,riskscore,by.x = "case_submitter_id",by.y = "sample")
  colnames(data)<-c("Sample","Stage","Age","Subtype","Gender","Race","Riskscore","Plot_group","Status","Time","Group")
  data$Stage<-as.factor(data$Stage)%>%as.numeric()
  data$Race<-as.factor(data$Race)%>%as.numeric()
  data$Gender<-as.factor(data$Gender)%>%as.numeric()
  data$Group<-as.factor(data$Group)%>%as.numeric()
  data$Subtype<-as.factor(data$Subtype)%>%as.numeric()
  
  ##single cox
  covariates <- c("Stage", "Subtype", "Age",  "Gender", "Race", "Riskscore")
  univ_formulas <- sapply(covariates,function(x) as.formula(paste('Surv(Time, Status)~', x)))
  univ_models <- lapply( univ_formulas, function(x){coxph(x, data = data)})
  univ_results <- lapply(univ_models,
                         function(x){
                           x <- summary(x)
                           p.value <- signif(x$waldtest["pvalue"], digits = 2)
                           wald.test <- signif(x$waldtest["test"], digits = 2)
                           beta <- signif(x$coef[1], digits = 2);
                           HR <- signif(x$coef[2], digits = 2);
                           HR.confint.lower <- signif(x$conf.int[ ,"lower .95"], 2)
                           HR.confint.upper <- signif(x$conf.int[ ,"upper .95"], 2)
                           HR <- paste0(HR, " (",
                                        HR.confint.lower, "-", HR.confint.upper, ") ")
                           res <- c(beta, HR, wald.test, p.value,HR.confint.lower,HR.confint.upper)
                           names(res) <- c("beta", "HR(95% CI for HR)", "wald.test",
                                           "p.value","HR.confint.lower","HR.confint.upper")
                           return(res)
                         })
  
  res <- as.data.frame(t(as.data.frame(univ_results, check.name = F)))
  
  out_multi<-t(res)
  out_multi<-t(out_multi)%>%as.data.frame(.)
  out_multi$HR<-str_split(out_multi$`HR(95% CI for HR)`," ",simplify = T)[,1]
  tabletext <- cbind(c(NA,"Risk Factors",rownames(out_multi)),
                     c(NA,"P value",out_multi$p.value),
                     c(NA,"Hazard Ratio(95% CI)",out_multi$`HR(95% CI for HR)`))
  forestplot(labeltext=tabletext, 
             graph.pos=3,  
             col=fpColors(box="#D55E00", lines="#CC79A7", zero = "gray50"),
             mean=c(NA,NA,out_multi$HR),
             lower=c(NA,NA,out_multi$HR.confint.lower), #95%
             upper=c(NA,NA,out_multi$HR.confint.upper), #95%
             boxsize=0.2,lwd.ci=2,  
             ci.vertices.height = 0.2,ci.vertices=TRUE,
             zero=1,lwd.zero=1,     
             colgap=unit(7,"mm"),    
             xticks = c(0,0.5, 1,1.5,2), 
             lwd.xaxis=1,            
             lineheight = unit(2.6,"cm"), 
             graphwidth = unit(.4,"npc"), 
             cex=0.9, fn.ci_norm = fpDrawCircleCI, 
             hrzl_lines=list("2" = gpar(lwd=4, col="black"),
                             "3" = gpar(lwd=2, col="black")),
             mar=unit(rep(0.5, times = 4), "cm"),
             
             txt_gp=fpTxtGp(label=gpar(cex=1.2),
                            ticks=gpar(cex=1.2),
                            xlab=gpar(cex = 1.25),
                            title=gpar(cex = 1.2)),
             xlab="Hazard Ratio",
             title = "Univariate cox regression analysis")
  res$p.value<-as.numeric(res$p.value)
  write.table(res,file = "single_cox_result.txt",sep = "\t",row.names = T,col.names = T)
  res_filter<-res[res$p.value<0.05,]
  
  ##multi cox
  res.cox <- coxph(as.formula(paste("Surv(Time, Status) ~ ",paste(rownames(res_filter),collapse = "+"),collapse = "")),
                   data = data)
  x <- summary(res.cox)
  pvalue=signif(as.matrix(x$coefficients)[,5],2)
  HR=signif(as.matrix(x$coefficients)[,2],2)
  low=signif(x$conf.int[,3],2)
  high=signif(x$conf.int[,4],2)
  multi_res=data.frame(p.value=pvalue,
                       HR=paste(HR," (",low,"-",high,")",sep=""),
                       stringsAsFactors = F
  )
  colnames(multi_res)<-c("p.value","HR(95% CI for HR)")
  res<-res[,c("p.value","HR(95% CI for HR)")]
  multi_res<-multi_res[match(rownames(res),rownames(multi_res)),]
  tmp<-cbind(res,multi_res)
  write.table(tmp,file = "multi_cox_result.txt",sep = "\t",row.names = T,col.names = T)
  
}

## hubgene and TF expression 14:5
{
  load("riskscore.RData")
  load("basic_data.RData")
  load("raw_genecount.RData")
  clinic<-clinic[match(riskscore$sample,clinic$sample1),]
  for (i in c("Module3","Module4","Module7","Module11","Module31")) {
    model=i
    hubgene<-read.table(paste("ppi/",model,"_cytohubba.csv",sep = ""),skip = 1,sep = ",",header = T)[,2]
    TFs <- read.delim(paste("tf/",model,".tsv",sep = ""))%>%.[,c("TF","Overlapping_Genes")]%>%head(.,10)%>%.$TF
    if (i=="Module3") {
      TFs<-TFs[c(1:6,8:10,7)]
    }
    BRCA.mrna.select<-TCGA.BRCA.htseq_counts.tsv.gz
    BRCA.mrna.select<-BRCA.mrna.select[match(riskscore$sample,rownames(BRCA.mrna.select)),c(hubgene,TFs)]
    
    BRCA.mrna.select<-BRCA.mrna.select[,rev(colnames(BRCA.mrna.select))]
    data<-cbind(BRCA.mrna.select,riskscore)
    data<-melt(data,id.vars = colnames(riskscore),variable.name = "Gene",value.name = "Value")
    data$group<-factor(data$group,levels = c("Low","High"))
    
    if (i=="Module3") {
      p<-ggplot(data=data)+geom_boxplot(aes(x = Gene,y = Value,fill=group))+theme_bw()+xlab("")+ylab("Expression")+
        stat_compare_means(aes(x=Gene,y=Value,group=group),label="p.signif",label.y = 4)+ggtitle(model)+coord_cartesian(ylim = c(-4,4))+
        theme(text = element_text(size = 20),axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,colour = c("#41AB5D",rep(c("#DD3497","#3690C0"),c(9,10)))))+
        scale_fill_aaas(name="Riskscore group")
      ggsave(filename = paste("./figure/",model,"_expression.pdf",collapse = "",sep = ""),p,width = 14,height = 5)
    }else{
      p<-ggplot(data=data)+geom_boxplot(aes(x = Gene,y = Value,fill=group))+theme_bw()+xlab("")+ylab("Expression")+
        stat_compare_means(aes(x=Gene,y=Value,group=group),label="p.signif",label.y = 4)+ggtitle(model)+coord_cartesian(ylim = c(-4,4))+
        theme(text = element_text(size = 20),axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,colour = rep(c("#DD3497","#3690C0"),each=10)))+
        scale_fill_aaas(name="Riskscore group")
      ggsave(filename = paste("./figure/",model,"_expression.pdf",collapse = "",sep = ""),p,width = 14,height = 5)
    }
  }
}


##verification
###cptac
load("basic_data.RData")
rm(list=ls()[!ls()=="gene"])
colnames(gene)<-str_to_lower(colnames(gene))
cptac<- read.delim("rawdata/CPTAC2_Breast_Prospective_Collection_BI_Proteome.tmt10.tsv", row.names=1)
cptac[is.na(cptac)]<-0
cptac<-cptac[4:nrow(cptac),1:(ncol(cptac)-6)]
cptac<-cptac[,grep("*.Unshared.Log.Ratio",colnames(cptac),invert = T)]
cptac<-t(cptac)

###protein
protein<-read.delim("rawdata/new.protein.txt", row.names=1)
protein<-t(protein)

load("basic_data.RData")
load("raw_genecount.RData")
BRCA.im.select<-BRCA.im.select[match(rownames(BRCA.mrna.select),rownames(BRCA.im.select)),paste(c("Area","Major","Minor","Ratio","DistMin","DistMax","DistMean"),rep(c(1,5),each=7),sep = "")]
BRCA.im.select<-BRCA.im.select[intersect(clinic$sample1,rownames(BRCA.mrna.select)),paste(c("Area","Major","Minor","Ratio","DistMin","DistMax","DistMean"),rep(c(1,5),each=7),sep = "")]
colnames(BRCA.im.select)<-c("Small_Nucleus_Area","Small_Major_Axis","Small_Minor_Axis","Small_Aspect_Ratio",
                            "Small_Min_Distance","Small_Max_Distance","Small_Mean_Distance",
                            "Large_Nucleus_Area","Large_Major_Axis","Large_Minor_Axis",
                            "Large_Aspect_Ratio","Large_Min_Distance","Large_Max_Distance","Large_Mean_Distance")
protein<-protein[intersect(rownames(protein),rownames(BRCA.im.select)),]
BRCA.im.select<-BRCA.im.select[intersect(rownames(protein),rownames(BRCA.im.select)),]
pca<-pca[intersect(rownames(protein),rownames(BRCA.im.select)),]
colnames(pca)<-str_to_lower(colnames(pca))
tmp<-list()
for (model in c("module3","module4","module7","module11","module31")) {
  TFs <- read.delim(paste("tf/",model,".tsv",sep = ""))%>%head(.,10)
  TFs$Library<-"TF"
  hubgene<-read.table(paste("ppi/",model,"_cytohubba.csv",sep = ""),skip = 1,sep = ",",header = T)[,2]%>%as.data.frame(.)
  hubgene$Library<-"Gene"
  
  tmp[[model]]<-data.frame(name=c(hubgene$.,TFs$TF),type=c(hubgene$Library,TFs$Library),model=model)
  g<-tmp[[model]]
  g<-g[order(g$type),]
  data<-protein[rownames(protein),intersect(g$name,colnames(protein))]
  
  ###hubgenes & TFs    10:7.8
  {
    test1<-cbind(protein[rownames(protein),intersect(g$name[g$type=="TF"],colnames(protein))],BRCA.im.select,pca[,model])
    colnames(test1)[ncol(test1)]<-"eigengene"
    test2<-cbind(protein[rownames(protein),intersect(g$name[g$type=="Gene"],colnames(protein))],BRCA.im.select,pca[,model])
    colnames(test2)[ncol(test2)]<-"eigengene"
    cor<-corr.test(as.matrix(test1),as.matrix(test2[,rev(colnames(test2))]),method = "spearman")
    p<-ggcorrplot(cor$r,p.mat = cor$p,hc.order = F,outline.color = "white",
                  ggtheme = theme_bw(),sig.level = 0.05,insig = "blank",lab_size = 3,digits = 1,
                  colors = c("#6D9EC1","white","#E46726"),lab = T)+
      theme(axis.text.x = element_text(color = rep(c("#C24294","red"),c((ncol(test1)-15),15))),axis.text.y = element_text(color = rep(c("red","#4EADEA"),c(15,(ncol(test2)-15)))),axis.text.x.top = element_text(hjust = 0,vjust = 0))+
      scale_x_discrete(position = "top")
    ggsave(paste("verification/",model,"_image_protein_heatmap.pdf",sep = ""),p,width = 10,height = 7.8)
  }
}

### gse dataset
GSE196666 <- read.delim("rawdata/GSE196666_STAR.counts.WDR5.txt", row.names=1)
GSE196666<-GSE196666[!duplicated(GSE196666$hgnc),]
rownames(GSE196666)<-GSE196666$hgnc;GSE196666<-GSE196666[,-ncol(GSE196666)]
GSE196666<-log2(GSE196666+1)
GSE196666<-t(GSE196666)

tmp<-list()
for (model in c("module3","module4","module7","module11","module31")) {
  TFs <- read.delim(paste("tf/",model,".tsv",sep = ""))%>%head(.,10)
  TFs$Library<-"TF"
  hubgene<-read.table(paste("ppi/",model,"_cytohubba.csv",sep = ""),skip = 1,sep = ",",header = T)[,2]%>%as.data.frame(.)
  hubgene$Library<-"Gene"
  
  tmp[[model]]<-data.frame(name=c(hubgene$.,TFs$TF),type=c(hubgene$Library,TFs$Library),model=model)
  g<-tmp[[model]]
  g<-g[order(g$type),]
  
  ###hubgenes & TFs    10:7.8
  {
    test1<-GSE196666[,intersect(g$name[g$type=="TF"],colnames(GSE196666))]
    
    test2<-GSE196666[,intersect(g$name[g$type=="Gene"],colnames(GSE196666))]
    cor<-corr.test(as.matrix(test1),as.matrix(test2[,rev(colnames(test2))]),method = "spearman")
    p<-ggcorrplot(cor$r,p.mat = cor$p,hc.order = F,outline.color = "white",
                  ggtheme = theme_bw(),sig.level = 0.05,insig = "blank",lab_size = 3,digits = 1,
                  colors = c("#6D9EC1","white","#E46726"),lab = T)+
      theme(axis.text.x = element_text(color = "#C24294"),axis.text.y = element_text(color = "#4EADEA"),axis.text.x.top = element_text(hjust = 0,vjust = 0))+
      scale_x_discrete(position = "top")
    ggsave(paste("verification/",model,"_gse196666_gene_heatmap.pdf",sep = ""),p,width = 10,height = 7.8)
  }
}

### gse dataset
GSE212143 <- read_excel("rawdata/GSE212143_15_breast_cancer_cell_lines_RNAseq.xlsx")
GSE212143<-GSE212143[!duplicated(GSE212143$gene_short_name),]
GSE212143<-as.data.frame(GSE212143)
rownames(GSE212143)<-GSE212143$gene_short_name
GSE212143<-GSE212143[,-c(1:2)]
GSE212143<-t(GSE212143)
GSE212143<-apply(GSE212143, 2,as.numeric)
tmp<-list()
for (model in c("module3","module4","module7","module11","module31")) {
  TFs <- read.delim(paste("tf/",model,".tsv",sep = ""))%>%head(.,10)
  TFs$Library<-"TF"
  hubgene<-read.table(paste("ppi/",model,"_cytohubba.csv",sep = ""),skip = 1,sep = ",",header = T)[,2]%>%as.data.frame(.)
  hubgene$Library<-"Gene"
  
  tmp[[model]]<-data.frame(name=c(hubgene$.,TFs$TF),type=c(hubgene$Library,TFs$Library),model=model)
  g<-tmp[[model]]
  g<-g[order(g$type),]
  
  ###hubgenes & TFs    10:7.8
  {
    test1<-GSE212143[,intersect(g$name[g$type=="TF"],colnames(GSE212143))]
    
    test2<-GSE212143[,intersect(g$name[g$type=="Gene"],colnames(GSE212143))]
    cor<-corr.test(as.matrix(test1),as.matrix(test2[,rev(colnames(test2))]),method = "spearman")
    p<-ggcorrplot(cor$r,p.mat = cor$p,hc.order = F,outline.color = "white",
                  ggtheme = theme_bw(),sig.level = 0.05,insig = "blank",lab_size = 3,digits = 1,
                  colors = c("#6D9EC1","white","#E46726"),lab = T)+
      theme(axis.text.x = element_text(color = "#C24294"),axis.text.y = element_text(color = "#4EADEA"),axis.text.x.top = element_text(hjust = 0,vjust = 0))+
      scale_x_discrete(position = "top")
    ggsave(paste("verification/",model,"_gse212143_gene_heatmap.pdf",sep = ""),p,width = 10,height = 7.8)
  }
}

### gse dataset
GSE211729 <- read.delim("rawdata/GSE211729_series_matrix.txt")
rownames(GSE211729)<-GSE211729$Sample_source_name_ch1
GSE211729<-GSE211729[,grep("Breast.Tumor.*",colnames(GSE211729))]
colnames(GSE211729)<-GSE211729[1,]
GSE211729<-GSE211729[-1,]
rownames(GSE211729)<-gsub("_at","",rownames(GSE211729))

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
symbol<-getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
              filters = 'ensembl_gene_id',
              values = rownames(GSE211729), 
              mart = ensembl)
symbol<-symbol[!duplicated(symbol$hgnc_symbol),]
symbol<-symbol[match(intersect(rownames(GSE211729),symbol$ensembl_gene_id),symbol$ensembl_gene_id),]
GSE211729<-GSE211729[match(intersect(rownames(GSE211729),symbol$ensembl_gene_id),rownames(GSE211729)),]
rownames(GSE211729)<-symbol$hgnc_symbol
GSE211729<-t(GSE211729)
GSE211729<-apply(GSE211729, 2, as.numeric)
tmp<-list()
for (model in c("module3","module4","module7","module11","module31")) {
  TFs <- read.delim(paste("tf/",model,".tsv",sep = ""))%>%head(.,10)
  TFs$Library<-"TF"
  hubgene<-read.table(paste("ppi/",model,"_cytohubba.csv",sep = ""),skip = 1,sep = ",",header = T)[,2]%>%as.data.frame(.)
  hubgene$Library<-"Gene"
  
  tmp[[model]]<-data.frame(name=c(hubgene$.,TFs$TF),type=c(hubgene$Library,TFs$Library),model=model)
  g<-tmp[[model]]
  g<-g[order(g$type),]
  
  ###hubgenes & TFs    10:7.8
  {
    test1<-GSE211729[,intersect(g$name[g$type=="TF"],colnames(GSE211729))]
    
    test2<-GSE211729[,intersect(g$name[g$type=="Gene"],colnames(GSE211729))]
    cor<-corr.test(as.matrix(test1),as.matrix(test2[,rev(colnames(test2))]),method = "spearman")
    p<-ggcorrplot(cor$r,p.mat = cor$p,hc.order = F,outline.color = "white",
                  ggtheme = theme_bw(),sig.level = 0.05,insig = "blank",lab_size = 3,digits = 1,
                  colors = c("#6D9EC1","white","#E46726"),lab = T)+
      theme(axis.text.x = element_text(color = "#C24294"),axis.text.y = element_text(color = "#4EADEA"),axis.text.x.top = element_text(hjust = 0,vjust = 0))+
      scale_x_discrete(position = "top")
    ggsave(paste("verification/",model,"_gse211729_gene_heatmap.pdf",sep = ""),p,width = 10,height = 7.8)
  }
}
