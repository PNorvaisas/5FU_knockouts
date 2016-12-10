library('ggplot2')
library('gplots')
library('plyr')
library('reshape2')
library('tidyr')
#library('Vennerable')
#library('GOplot')
library('quantreg')
#library('KEGGREST')
#library('devtools')
library('splitstackshape')

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


odir<-'Figures_v2'
ddir<-'Data_v2'

qbacmicq<-read.table(paste(ddir,'/MICs_and_bacterial_growth-Clean.csv',sep=''),sep=',',header=TRUE)
qbacmicq<-qbacmicq[,!colnames(qbacmicq) %in% c('X','JW_id','ECK','X0','X1','X2.5','X5','bno','EG','GI','MOPS_24hr','MOPS_48hr')]


ecor<-read.table('../EColi_annotation/EcoCyc_Patwhays.tsv',sep='\t',header=TRUE,stringsAsFactors = FALSE)
ecom<- cSplit(ecor, "Genes", sep = ";", direction = "long")
ecocyc<-rename(ecom,c('Genes'='Gene','Pathway'='Term'))
ecocyc<-as.data.frame(ecocyc)
ecocyc$Link<-as.character(ecocyc$Link)
ecocyc$Term_ID<-separate(data = ecocyc, col = Link, into = c("Base", "Mid",'Term_ID'), sep = "\\=")$Term_ID
ecocyc$Category<-'EcoCyc'
ecocyc$Link<-NULL


danr<-read.Annot('../EColi_annotation/EColi_Keio_DAVID_GO.txt')
dan<-subset(danr,Strain=='eco' & Term_ID!='NA')
dan$Strain<-NULL
dan$Species<-NULL
dan$Gene.Name<-NULL
danm<-merge.data.frame(dan,qbacmicq[,c('Gene','UniACC')],by.x='ID',by.y='UniACC',all.x=TRUE)

#ecocycm<-merge.data.frame(ecocyc,danm[,c('Gene','Gene.Name')],by='Gene')


allannot<-merge(ecocyc,danm,all.x=TRUE,all.y=TRUE)
allannot$ID<-NULL
allannot<-rename(allannot,c('Size'='Set_size'))
allannot<-allannot[,c('Gene','Category','Term_ID','Term','Set_size')]



circm<-merge(qbacmicq,allannot, by.x='Gene',by.y='Gene',all.x = TRUE)

circm$Term<-as.factor(circm$Term)
circm$Category<-as.factor(circm$Category)


circms<-ddply(subset(circm,MIC>=5),.(Category,Term_ID),summarise,
               MIC_avg_sum=mean(MIC,na.rm = TRUE),MIC_sd_sum=sd(MIC,na.rm = TRUE),
               MIC_med_sum=median(MIC,na.rm = TRUE),
               MIC_05_sum=quantile(MIC,probs=c(0.05),na.rm = TRUE),MIC_95_sum=quantile(MIC,probs=c(0.95),na.rm = TRUE),
               MIC_25_sum=quantile(MIC,probs=c(0.25),na.rm = TRUE),MIC_75_sum=quantile(MIC,probs=c(0.75),na.rm = TRUE),
               NGM_D_avg_sum=mean(NGM_D,na.rm = TRUE),NGM_D_sd_sum=sd(NGM_D,na.rm = TRUE),
               NGM_D_med_sum=median(NGM_D,na.rm = TRUE),
               NGM_D_05_sum=quantile(NGM_D,probs=c(0.05),na.rm = TRUE),NGM_D_95_sum=quantile(NGM_D,probs=c(0.95),na.rm = TRUE),
               NGM_D_25_sum=quantile(NGM_D,probs=c(0.25),na.rm = TRUE),NGM_D_75_sum=quantile(NGM_D,probs=c(0.75),na.rm = TRUE),
               NGM_C_avg_sum=mean(NGM_C,na.rm = TRUE),NGM_C_sd_sum=sd(NGM_C,na.rm = TRUE),
               NGM_C_med_sum=median(NGM_C,na.rm = TRUE),
               NGM_C_05_sum=quantile(NGM_C,probs=c(0.05),na.rm = TRUE),NGM_C_95_sum=quantile(NGM_C,probs=c(0.95),na.rm = TRUE),
               NGM_C_25_sum=quantile(NGM_C,probs=c(0.25),na.rm = TRUE),NGM_C_75_sum=quantile(NGM_C,probs=c(0.75),na.rm = TRUE),
               LB_sum=mean(LB_22hr,na.rm = TRUE),LB_sd_sum=sd(LB_22hr,na.rm = TRUE),
               LB_med_sum=median(LB_22hr,na.rm = TRUE),
               LB_05_sum=quantile(LB_22hr,probs=c(0.05),na.rm = TRUE),LB_95_sum=quantile(LB_22hr,probs=c(0.95),na.rm = TRUE),
               LB_25_sum=quantile(LB_22hr,probs=c(0.25),na.rm = TRUE),LB_75_sum=quantile(LB_22hr,probs=c(0.75),na.rm = TRUE),
               gcount=as.numeric(length(unique(Gene))) )



enr<-merge(circm,circms,by=c('Category','Term_ID'),all.x=TRUE)
enrnames<-read.table('Data_v2/Enrichment_dataset_ref.csv',sep=',')
enrnames<-rename(enrnames,c('V2'='Name'))
enrnames$Name<-as.character(enrnames$Name)
enrnames<-enrnames[enrnames$Name!='x','Name']

try(if(length(intersect(colnames(enr),enrnames))!=44)
  stop("Possible change in enrichment dataset!"))
enr<-enr[,colnames(enr)[c(1:2,14:15,44,3:13,16:43)]]
enr[enr=='NaN']<-NA


write.csv(allannot,paste(ddir,'/All_enrichment_terms.csv',sep=''))
write.csv(circms,paste(ddir,'/Distributions_for_terms.csv',sep=''))
write.csv(enr,paste(ddir,'/Complete_data_term_distribution.csv',sep=''))


theme_set(theme_light())

mtres<-1
gtres<-4
mttres<-10

Slevel=c(0)
SAbr=c('all')
Stat<-data.frame(slevel=Slevel,sabr=SAbr)


TName=c('BP','KP','EcoCyc')
TTitle=c('Biological process', 'KEGG pathway','EcoCyc pathway')
TAbr=c('BP','KEGG-PWY','EcoCyc-PWY')
Types<-data.frame(tname=TName,tabr=TAbr,ttitle=TTitle)

Media=c('LB_22hr','NGM_D','NGM_C')
MTitle=c('LB 22hr growth','NGM + 100uM 5FU growth','NGM growth')
MAbr=c('LB','NGM-5FU','NGM')
MOrder=c('LB_med_sum','NGM_D_med_sum','NGM_C_med_sum')
MQvars=c('LB','NGM_D','NGM_C')
MScale=c(1,0.3,0.3)
Medias<-data.frame(media=Media,mtitle=MTitle,mabr=MAbr,morder=MOrder,mscale=MScale,mqvars=MQvars)

##Growth distributions
for(m in 1:nrow(Medias)) {
  mr<-Medias[m,]
  for(s in 1:nrow(Stat)) {
    sr <- Stat[s,]
    for(t in 1:nrow(Types)) {
      tr <- Types[t,]
      sel<-subset(enr,Category==as.character(tr$tname) & !is.na(enr[,as.character(mr$media)])  & MIC>mtres & gcount>gtres & MIC_avg_sum>mttres)
      #& logp_sum>2 & logp_sum>2  & logp_sum>2  & logp_sum>2
      gtitle<-paste('Enriched E. Coli MG1655 ',tr$ttitle,' - ',mr$mtitle,' distribution, ',sr$sabr,sep = '')
      fname<-paste(odir,'/Enrichment/Bacterial_growth/',tr$tabr,'_BG_',mr$mabr,'_',sr$sabr,'.pdf',sep = '')
      print(gtitle)
      print(fname)
      termBG<-ggplot(sel,
                     aes(y=sel[,as.character(mr$media)],x=reorder(sel$Term,sel[,as.character(mr$morder)]),
                         fill=gcount,color=MIC_avg_sum))+
        geom_boxplot(position='identity')+
        coord_flip()+labs(fill='Number of genes',color='MIC average')+
        scale_fill_gradient( high="red",low='white')+
        scale_color_gradient(limits=c(0,100), high="blue",low='gray')+
        ylab('OD')+ylim(0,mr$mscale)+
        xlab(tr$ttitle)+ggtitle(gtitle)
      #termBG
      cairo_pdf(fname,width=9,height=9)
      print(termBG)
      dev.off()
      #scale_color_gradient(limits=c(0,8), high="blue",low='gray')+
      #,
      #ggsave(plot=termBG,file=fname,width=9,height=9)
      #dev.copy2pdf(device=cairo_pdf,file=fname,width=12,height=12)
    }
  }
}

##MIC distributions
for(s in 1:nrow(Stat)) {
  sr <- Stat[s,]
  for(t in 1:nrow(Types)) {
    tr <- Types[t,]
    sel<-subset(enr,Category==as.character(tr$tname) & MIC>mtres & gcount>gtres& MIC_avg_sum>mttres)
    #& logp_sum>2 & logp_sum>2  & logp_sum>2  & logp_sum>2
    gtitle<-paste('Enriched E. Coli MG1655 ',tr$ttitle,' - MIC distribution, ',sr$sabr,sep = '')
    fname<-paste(odir,'/Enrichment/MIC/',tr$tabr,'_MIC_',sr$sabr,'.pdf',sep = '')
    print(gtitle)
    print(fname)
    termMIC<-ggplot(sel,
                    aes(y=MIC,x=reorder(Term,MIC_med_sum),fill=gcount,color=MIC_avg_sum))+
      geom_boxplot(position='identity')+coord_flip()+
      labs(fill='Number of genes',color='MIC average')+
      scale_fill_gradient( high="red",low='white')+
      scale_color_gradient(limits=c(0,100), high="blue",low='gray')+
      ylab('MIC [5FU], uM')+ylim(0,100)+
      xlab(tr$ttitle)+ggtitle(gtitle)
    cairo_pdf(fname,width=9,height=9)
    print(termMIC)
    dev.off()
    #      
  }
}



Media2=c('LB_22hr','NGM_D','NGM_C')
MTitle2=c('LB 22hr growth','NGM + 100uM 5FU growth','NGM growth')
MAbr2=c('LB','NGM-5FU','NGM')
MOrder2=c('LB_med_sum','NGM_D_med_sum','NGM_C_med_sum')
MQvars2=c('LB','NGM_D','NGM_C')
MScale2=c(1,0.2,0.2)
Medias2<-data.frame(media=Media2,mtitle=MTitle2,mabr=MAbr2,morder=MOrder2,mscale=MScale2,mqvars=MQvars2)

##MIC vs Media
for(m in 1:nrow(Medias2)) {
  mr<-Medias2[m,]
  for(s in 1:nrow(Stat)) {
    sr <- Stat[s,]
    for(t in 1:nrow(Types)) {
      tr <- Types[t,]
#       if (mr$mqvars=='LB') {
#         sel<-subset(circms,category==as.character(tr$tname) & strain=='eco' & logp_sum>sr$slevel)
#       } else {
#         sel<-subset(circms,category==as.character(tr$tname) & strain=='eco' & logp_sum>sr$slevel)
#       }
      sel<-subset(enr,Category==as.character(tr$tname) & !is.na(enr[,paste(mr$mqvars,'_med_sum',sep='')] ) &MIC>mtres & gcount>gtres& MIC_avg_sum>mttres)
      #& logp_sum>2 & logp_sum>2  & logp_sum>2  & logp_sum>2
      gtitle<-paste('Enriched E. Coli MG1655 ',tr$ttitle,': MIC',' - ',mr$mtitle,' distribution, ',sr$sabr,sep = '')
      fname<-paste(odir,'/Enrichment/MICvsOD/',tr$tabr,'_MIC_vs_BG_',mr$mabr,'_',sr$sabr,'.pdf',sep = '')
      print(gtitle)
      print(fname)
      kpcor<-ggplot(sel,
                    aes(x=MIC_med_sum,y=sel[,paste(as.character(mr$mqvars),'_med_sum',sep='')],color=gcount))+
        geom_point()+
        geom_errorbarh(aes(xmax=MIC_75_sum,
                                      xmin=MIC_25_sum),height=mr$mscale*0.0125,alpha=0.5)+
        geom_errorbar(aes(ymax=sel[,paste(mr$mqvars,'_75_sum',sep='')],
                                     ymin=sel[,paste(mr$mqvars,'_25_sum',sep='')]),width=1,alpha=0.5)+
        geom_text(aes(label=Term),color='black',hjust=-0.1, vjust=-0.5,size=3)+
        scale_color_gradient(limits=c(0,8), high="blue",low='gray')+
        labs(color='Number of genes')+
        xlab(expression(paste('MIC [5FU], ',mu,'M')))+ylab('OD')+ylim(0,mr$mscale)+xlim(0,100)+
        ggtitle(gtitle)
      cairo_pdf(fname,width=9,height=9)
      print(kpcor)
      dev.off()
      #aes(size=gcount)
      #ggsave(plot=termBG,file=fname,width=9,height=9)
      #dev.copy2pdf(device=cairo_pdf,file=fname,width=12,height=12)
    }
  }
}

fitbac<-lm(NGM_D ~ NGM_C,data=qbacmicq)
confint(fitbac,'NGM_C',level=0.95)

for(s in 1:nrow(Stat)) {
  sr <- Stat[s,]
  for(t in 1:nrow(Types)) {
    tr <- Types[t,]
        
    sel<-subset(enr,Category==as.character(tr$tname) & MIC>mtres & gcount>gtres& MIC_avg_sum>mttres)
    #& logp_sum>2 & logp_sum>2  & logp_sum>2  & logp_sum>2
    gtitle<-paste('Enriched E. Coli MG1655 ',tr$ttitle,': NGM - NGM + 100um 5FU correlation, ',sr$sabr,sep = '')
    fname<-paste(odir,'/Enrichment/NGMvsNGM-5FU/',tr$tabr,'_NGM_vs_NGM-5FU_',sr$sabr,'.pdf',sep = '')
    print(gtitle)
    print(fname)
    kpcor<-ggplot(sel,aes(x=NGM_C_med_sum,y=NGM_D_med_sum,color=MIC_avg_sum))+
      geom_point(aes(size=gcount))+
      geom_abline(aes(fill='1:1'),intercept=0,slope=1,alpha=0.5,color='grey',linetype='longdash')+
      geom_abline(aes(fill='Trend for all knockouts'),intercept=coef(fitbac)[[1]],
                  slope=coef(fitbac)[[2]],alpha=0.5,color='red')+
      geom_errorbarh(aes(xmax=NGM_C_75_sum,xmin=NGM_C_25_sum),height=0.15*0.0125,alpha=0.5)+
      geom_errorbar(aes(ymax=NGM_D_75_sum,ymin=NGM_D_25_sum),width=0.2*0.0125,alpha=0.5)+
      geom_text(aes(label=Term),color='black',hjust=-0.1, vjust=-0.5,size=3)+
      labs(size='Number of genes',color='MIC average')+
      xlab('OD, NGM - Control')+
      ylab(expression(paste('OD, NGM + 100',mu,'M 5FU',sep='')))+
      ylim(0,0.2)+xlim(0,0.2)+
      ggtitle(gtitle)
    cairo_pdf(fname,width=9,height=9)
    print(kpcor)
    dev.off()
    #      scale_color_gradient(limits=c(0,8), high="blue",low='gray')+ ,fill='Guides'
  }
}



#ggsave(plot=termBG,file=fname,width=9,height=9)
#dev.copy2pdf(device=cairo_pdf,file=fname,width=12,height=12)


# 
# circgo<-subset(circm,category %in% c('KP','BP'))
# 
# GOBar(subset(circgo, category == 'KP'))
# GOBar(circgo, display = 'multiple')
# GOBar(subset(circgo,adj_pval>0.01), display = 'multiple')
# GOBubble(circgo, title = 'Bubble plot', color = c('deeppink', 'blue', 'chartreuse1'), display = 'multiple', labels = 4)
# GOBubble(circgo)
# circgo<-subset(circgo,adj_pval<0.01 & count>4)
# GOCircle(subset(circgo, category == 'BP'))
# #dev.copy2pdf(device=cairo_pdf,file="GO/GO_AUCnorm_All.pdf",width=12,height=8)
# 
# BP<-unique(subset(circgo, category == 'BP' & count>4)[,c('ID','term','adj_pval','zscore')])$ID
# CC<-unique(subset(circgo, category == 'CC' & count>4)[,c('ID','term','adj_pval','zscore')])$ID
# GOCircle(circgo,nsub=c(CC[1],BP[1:9]))
# 
# GOCircle(subset(circgo,zscore<0 & category == 'BP'))
# #dev.copy2pdf(device=cairo_pdf,file="GO/GO_AUCnorm_All.pdf",width=12,height=8)
# GOCircle(subset(circgo,zscore>0 & category == 'BP'))
# #dev.copy2pdf(device=cairo_pdf,file="GO/GO_Peak_UP.pdf",width=12,height=8)
# 
# #circ<-subset(circ,category=='BP')
# 
# 
# MIC_top<-unique(circ[with(circ,order(logFC,decreasing=TRUE)),]$genes)[1:10]
# genes_bottom<-unique(circ[with(circ,order(logFC,decreasing=FALSE)),]$genes)[1:10]
# 
# 
# BP_top<-as.character(unique(subset(circgo, category == 'BP' & count>4)[,c('ID','term','adj_pval','zscore')]$term))
# KP_top<-as.character(unique(subset(circgo, category == 'KP' & count>4)[,c('ID','term','adj_pval','zscore')]$term))
# 
# 
# chord <- chord_dat(circgo, process=KP_top[1:10])
# chord <- chord_dat(circgo, process=BP_top[1:10])
# GOChord(chord, space = 0.02, gene.space = 0.25, gene.size = 5)
# #dev.copy2pdf(device=cairo_pdf,file="GO/GO_AUCnorm_All_circle.pdf",width=16,height=16)
# 
# GOCluster(circ, KP_top[1:4], clust.by = 'term', term.width = 2)
# 
# 


