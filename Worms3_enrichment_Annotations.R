library('ggplot2')
library('gplots')
library('plyr')
library('reshape2')
library('tidyr')
library('quantreg')
library('splitstackshape')
library(gtable)
library(gridExtra)
library(xlsx)
library('gtools')


theme_set(theme_light())

theme_update(panel.background = element_rect(colour = "black"),
             axis.text = element_text(colour = "black"))

setwd("~/Projects/2015-Metformin/Worms")

capFirst <- function(s) {
  paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "")
}

GOsplit<-function(x,rtrn="GO"){
  if (grepl('~',x)) {
    sep='~'
  } else if (grepl(':',x)) {
    sep=':'
  } else {
    if (rtrn=="GO"){
      return(NA)
    } else if (rtrn=="Term"){
      return(x)
    }
  }
  splt<-strsplit(x,split=sep)[[1]]
  res<-list()
  if (length(splt)==2){
    res$GO<-as.character(splt[1])
    res$Term<-as.character(splt[2])
  } else {
    res$GO<-NA
    res$Term<-as.character(splt[1])
  }
  if (rtrn=="GO"){
    return(res$GO)
  } else if (rtrn=="Term"){
    return(res$Term)
  }
  
}

read.GO<-function(flnm,filter=TRUE,pval='Benjamini'){
  GO<-read.table(flnm,sep='\t',quote = '"',header = TRUE,stringsAsFactors=FALSE)
  as.data.frame(GO)
  #rownames(GO)<-str(rownames(GO))
  GO$ID<-as.character(lapply(GO$Term,GOsplit,rtrn="GO"))
  GO$Term<-as.character(lapply(GO$Term,GOsplit,rtrn="Term"))
  if (pval=='Benjamini'){
    GO$adj_pval<-GO$Benjamini
  } else if (pval=='Bonferoni'){
    GO$adj_pval<-GO$Bonferoni
  } else if (pval=='FDR'){
    GO$adj_pval<-GO$FDR
  } else if (pval=='PValue'){
    GO$adj_pval<-GO$PValue
  } else {
    GO$adj_pval<-GO$PValue
  }
  
  GO[GO$Category=="KEGG_PATHWAY",]$Category<-"KP"
  GO[GO$Category=="GOTERM_BP_FAT",]$Category<-"BP"
  GO[GO$Category=="GOTERM_CC_FAT",]$Category<-"CC"
  GO[GO$Category=="GOTERM_MF_FAT",]$Category<-"MF"
  if (filter){
    GO<-subset(GO,Category %in% c("BP","MF","CC"))
  }
  
  return(GO)
}


read.Annot<-function(flnm,filter=TRUE){
  dar<-read.table(flnm,sep='\t',quote = '"',header = TRUE,stringsAsFactors=FALSE)
  dar<-dar[,c('ID',"Gene.Name","Species","OFFICIAL_GENE_SYMBOL","GOTERM_BP_FAT","GOTERM_CC_FAT","GOTERM_MF_FAT","KEGG_PATHWAY")]
  dar<-dar[,c('ID',"Gene.Name","Species","GOTERM_BP_FAT","GOTERM_CC_FAT","GOTERM_MF_FAT","KEGG_PATHWAY")]
  darm<-melt(dar,id=c('ID',"Gene.Name","Species"),variable.name = 'Category',value.name='Term')
  darmm<-cSplit(darm, "Term", sep = ",", direction = "long")
  darmm<-data.frame(darmm)
  darmm$Term<-as.character(darmm$Term)
  darmm$Category<-as.character(darmm$Category)
  darmm$Term_ID<-as.character(lapply(darmm$Term,GOsplit,rtrn="GO"))
  darmm$Term<-as.character(lapply(darmm$Term,GOsplit,rtrn="Term"))
  darmm[darmm$Category=="KEGG_PATHWAY",]$Category<-"KP"
  darmm[darmm$Category=="GOTERM_BP_FAT",]$Category<-"BP"
  darmm[darmm$Category=="GOTERM_CC_FAT",]$Category<-"CC"
  darmm[darmm$Category=="GOTERM_MF_FAT",]$Category<-"MF"
  darmm$Strain<-ifelse(darmm$Category=='KP',substr(darmm$Term_ID,1,3),'eco')
  return(darmm)
}




#Annotations

ecor<-read.table('../EColi_annotation/EcoCyc_Patwhays.tsv',sep='\t',header=TRUE,stringsAsFactors = FALSE)
ecom<- cSplit(ecor, "Genes", sep = ";", direction = "long")
ecocyc<-rename(ecom,c('Genes'='Gene','Pathway'='Term'))
ecocyc<-as.data.frame(ecocyc)
ecocyc$Link<-as.character(ecocyc$Link)
ecocyc$Term_ID<-separate(data = ecocyc, col = Link, into = c("Base", "Mid",'Term_ID'), sep = "\\=")$Term_ID
ecocyc$Category<-'EcoCyc'
ecocyc$Link<-NULL

#Get KEGG pathway annotations
keggr<-read.table('../EColi_annotation/Gene_KeggEco.csv',sep=',',
                  header=TRUE,stringsAsFactors = TRUE,colClasses=c(rep("factor",6)))
keggc<-subset(keggr,! is.na(KeggPathID))
patc<-data.frame(table(keggc$KeggPathID))
keggp<-merge(keggc,patc, by.x='KeggPathID',by.y='Var1',all.x=TRUE)
keggp$X<-NULL
keggp$Category<-'KEGG'
keggm<-rename(keggp,c('PathDescription'='Term','Freq'='Size','KeggPathID'='Term_ID'))

#Get GO annotaations
goECr<-read.table('../EColi_annotation/Gene_GO.csv',sep=',',
                  header=TRUE,stringsAsFactors = TRUE)#,colClasses=c(rep("factor",6))
goEC<-subset(goECr,!is.na(GO_ID))
goEC$X<-NULL
goc<-data.frame(table(goEC$GO_ID))
goECm<-merge(goEC,goc,by.x='GO_ID',by.y='Var1',all.x=TRUE)
goECm$Category<-paste('GO',as.character(goECm$Ontology),sep='_')
goECm<-rename(goECm,c('GO_ID'='Term_ID','Freq'='Size'))


#Gene set sizes need to be estimated for whole E coli, screened Keio and MIC>5

allannot<-merge(keggm[,c('Term','Size','Gene','Term_ID','Category')],ecocyc,all.x=TRUE,all.y=TRUE)
allannot<-merge(goECm[,c('Term','Size','Gene','Term_ID','Category')],allannot,all.x=TRUE,all.y=TRUE)
allannot$ID<-NULL
allannot<-rename(allannot,c('Size'='Size_EC'))
allannot<-allannot[,c('Gene','Category','Term_ID','Term','Size_EC')]


allannot$Term<-capFirst(allannot$Term)
allannot[allannot$Term=='NANA','Term']<-NA

#Final Allannot
write.csv(allannot,paste(ddir,'/Enrichement_all_terms.csv',sep=''),row.names = FALSE)
