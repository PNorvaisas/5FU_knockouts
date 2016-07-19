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

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
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




odir<-'Figures_final'
ddir<-'Data_final'

bacmic<-read.table(paste(ddir,'/MICs_and_bacterial_growth-Complete.csv',sep=''),sep=',',header=TRUE)
bacmic<-bacmic[,!colnames(bacmic) %in% c('X','JW_id','ECK','X0','X1','X2.5','X5','bno','EG','GI','MOPS_24hr','MOPS_48hr')]


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


circmr<-merge(bacmic,allannot, by.x='Gene',by.y='Gene',all.x = TRUE)
circmr$Term<-as.factor(circmr$Term)
circmr$Category<-as.factor(circmr$Category)

#Get gene counts for screened Keio library
circsumfull<-ddply(circmr,.(Category,Term_ID),summarise, Size_KeioS=as.numeric(length(unique(Gene))) )

circm<-merge(circmr,circsumfull,by=c('Category','Term_ID'),all.x=TRUE)

unique(allannot$Category)


theme_set(theme_light())

Slevel=c(0)
SAbr=c('all')
Stat<-data.frame(slevel=Slevel,sabr=SAbr)

# Old lists
# TName=c('BP','KP','EcoCyc')
# TTitle=c('Biological process', 'KEGG pathway','EcoCyc pathway')
# TAbr=c('BP','KEGG-PWY','EcoCyc-PWY')
# Types<-data.frame(tname=TName,tabr=TAbr,ttitle=TTitle)


TName=c('KEGG','EcoCyc','GO_BP')
TTitle=c('KEGG pathway','EcoCyc pathway','GO Biological process')
TAbr=c('KEGG-PWY','EcoCyc-PWY','GO_BP')
Types<-data.frame(tname=TName,tabr=TAbr,ttitle=TTitle)

Thrnames<-c('Coverage in E. coli','Coverage in valid Keio','Score in whole E. coli','Score in screened Keio','Score in whole E. coli (median)','Score in screened Keio (median)')
Thrcols<-c('Coverage_EC','Coverage_KeioS','Score_EC','Score_KeioS','Score_EC_Med','Score_KeioS_Med')
Thraxis<-c('Coverage','Coverage','Score','Score','Score','Score')
Thrabr<-c('CovEC','CovKS','Score_EC','Score_KS','Score_ECmed','Score_KSmed')
Thresholds<-data.frame(names=Thrnames,cols=Thrcols,axis=Thraxis,abr=Thrabr)


Media=c('LB_22hr','OD_T_Mean','OD_C_Mean')
MTitle=c('LB 22hr growth','NGM + 100uM 5FU growth','NGM growth')
MAbr=c('LB','NGM-5FU','NGM')
MOrder=c('LB_S_Med','OD_T_Mean_S_Med','OD_C_Mean_S_Med')
MQvars=c('LB','OD_T_Mean','OD_C_Mean')
MScale=c(1,0.3,0.3)
Medias<-data.frame(media=Media,mtitle=MTitle,mabr=MAbr,morder=MOrder,mscale=MScale,mqvars=MQvars)


#Enrichment for genes with MIC>smth
#micthres<-5
#MIC thres in figures
mtres<-1
#Number of genes thres
gtres<-3
#MIC med thres
mttres<-5


micthres<-5

circmc<-circm[, -grep("_SD", colnames(circm))]
circmc<-circmc[, -grep("_pval", colnames(circmc))]
circmc<-circmc[,!colnames(circmc) %in% c('N')]

circmelt<-melt(subset(circmc,MIC>micthres),id=c('Category','Term_ID','Term','Gene','Plate','Well','Size_EC','Size_KeioS'),
               variable.name = 'Measure',value.name='Value')

circmstat<-ddply(circmelt,.(Category,Term_ID,Term,Size_EC,Size_KeioS,Measure),summarise,
              Size_MICo5=as.numeric(length(unique(Gene))),
              S_Mean=mean(Value,na.rm = TRUE),S_SD=sd(Value,na.rm = TRUE),
              S_Med=median(Value,na.rm = TRUE),
              S_Q05=quantile(Value,probs=c(0.05),na.rm = TRUE),S_Q95=quantile(Value,probs=c(0.95),na.rm = TRUE),
              S_Q25=quantile(Value,probs=c(0.25),na.rm = TRUE),S_Q75=quantile(Value,probs=c(0.75),na.rm = TRUE))


circmstatm<-melt(circmstat,measure.vars = c('S_Mean','S_SD','S_Med','S_Q05','S_Q25','S_Q75','S_Q95'),
              variable.name = 'Stat',value.name='Value')

circms<-dcast(circmstatm,Category+Term_ID+Term+Size_EC+Size_KeioS+Size_MICo5~ Measure+Stat,value.var = 'Value')



circms$Coverage_EC<-circms$Size_MICo5/circms$Size_EC
circms$Coverage_KeioS<-circms$Size_MICo5/circms$Size_KeioS
circms$Score_EC<-circms$Coverage_EC*circms$MIC_S_Mean
circms$Score_KeioS<-circms$Coverage_KeioS*circms$MIC_S_Mean
circms$Score_EC_Med<-circms$Coverage_EC*circms$MIC_S_Med
circms$Score_KeioS_Med<-circms$Coverage_KeioS*circms$MIC_S_Med






enr<-merge(circm,circms[,! colnames(circms) %in% c('Size_EC','Size_KeioS')],by=c('Category','Term_ID','Term'),all.x=TRUE)
enr[enr=='NaN']<-NA
dim(enr)



circms_re<-subset(circms,!is.na(Category))
circms_re<-circms_re[, -grep("logOD", colnames(circms_re))]
circms_re<-circms_re[, -grep("WTDiff", colnames(circms_re))]
#circms_re<-circms_re[, -grep("CTWTDiff", colnames(circms_re))]
circms_re<-circms_re[, -grep("CTODDiff", colnames(circms_re))]

circms_re<-circms_re[,c(1:6,
                        match(c('Coverage_EC'),colnames(circms_re)):match(c('Score_KeioS_Med'),colnames(circms_re)),
                        7: (match(c('Coverage_EC'),colnames(circms_re))-1) )]

annotexpl<-read.table(paste(ddir,'/Enrichment_distribution_column_explanation.csv',sep=''),
                    sep=',',quote = '"',header = TRUE,stringsAsFactors=FALSE)




explrm<-subset(annotexpl,Column %in% colnames(circms_re))
explrm<-explrm[match(colnames(circms_re),explrm$Column),]
write.csv(circms_re,
          paste(ddir,'/Enrichment_Distributions_for_terms_MICo',micthres,'.csv',sep=''),
          row.names = FALSE)
write.xlsx(explrm, file=paste(ddir,'/Enrichment_Distributions_for_terms_MICo',micthres,'.xlsx',sep=''),
           sheetName="Readme",row.names = FALSE,showNA=FALSE)
write.xlsx(circms_re, file=paste(ddir,'/Enrichment_Distributions_for_terms_MICo',micthres,'.xlsx',sep=''),
           sheetName="Data", append=TRUE,row.names = FALSE,showNA=FALSE)




#write.csv(enr,paste(ddir,'/Enrichment_Complete_data_term_distribution_MICo',micthres,'.csv',sep=''))


#Create necessary directories
envar<-paste('Enrichment_MICo',micthres,sep='')
dir.create(paste(odir,'/',envar,sep=''), showWarnings = TRUE, recursive = FALSE, mode = "0777")

micname<-expression(paste('Worm MIC\n median [5FU], ',mu,'M',sep=''))
mic<-expression(paste('Worm MIC [5FU], ',mu,'M',sep=''))

##MIC distributions and coverage
print('MIC distributions and coverage')
dir.create(paste(odir,'/',envar,'/MIC_coverage_and_distribution',sep=''), showWarnings = TRUE, recursive = FALSE, mode = "0777")
for(s in 1:nrow(Stat)) {
  sr <- Stat[s,]
  for(t in 1:nrow(Types)) {
    tr <- Types[t,]
    for(m in 1:nrow(Thresholds)) {
      mr<-Thresholds[m,]
      sel<-subset(enr,Category==as.character(tr$tname) & MIC>mtres & Size_MICo5>gtres&
                    MIC_S_Med>mttres &!is.na(Term) &
                    ! Term %in% c('Metabolic pathways',
                                  'Biosynthesis of secondary metabolites',
                                  'Microbial metabolism in diverse environments'))
      selcov<-subset(circms,Category==as.character(tr$tname) & Size_MICo5>gtres&
                    MIC_S_Med>mttres &!is.na(Term) &
                    ! Term %in% c('Metabolic pathways',
                                  'Biosynthesis of secondary metabolites',
                                  'Microbial metabolism in diverse environments'))
      
      gtitle<-paste('Enriched E. Coli MG1655 ',tr$ttitle,' - MIC distribution, ',sr$sabr,' MIC>',micthres,sep = '')
      fname<-paste(odir,'/',envar,'/MIC_coverage_and_distribution/',tr$tabr,'_MIC_',sr$sabr,'_MICo',micthres,'_',mr$abr,'.pdf',sep = '')
      #print(gtitle)
      
      termMICdist<-ggplot(sel,
                          aes(y=MIC,x=reorder(sel$Term,sel[,as.character(mr$cols)]),fill=Size_MICo5,color=MIC_S_Med))+
        geom_boxplot(position='identity')+
        coord_flip()+
        labs(fill='Number of genes',color='MIC median')+
        scale_fill_gradient( high="red",low='white')+
        #scale_color_gradient(limits=c(0,100), high="blue",low='gray')+
        scale_colour_gradientn(colours = c('black','orange','red'),
                               breaks=c(5,25,50),limits=c(5,50),guide='legend',name=micname)+
        ylab(mic)+ylim(0,100)+
        xlab(tr$ttitle)+ggtitle(gtitle)+theme(legend.position="left")
      gtitle2<-paste(mr$names,' MIC>',micthres,sep = '')
      termMICcov<-ggplot(selcov,
                         aes(y=selcov[,as.character(mr$cols)],x=reorder(selcov$Term,selcov[,as.character(mr$cols)])))+
        geom_bar(stat = "identity",fill='red',alpha=0.8)+
        coord_flip()+
        labs(fill='Number of genes')+
        scale_fill_gradient( high="red",low='gray')+
        ylab(mr$axis)+
        xlab(tr$ttitle)+ggtitle(gtitle2)+theme(axis.text.y=element_blank(),
                                               axis.title.y=element_blank(),
                                               axis.ticks.y=element_blank(),
                                               legend.position="none")
      #termMICcov
      #fill=Size_MICo5
      gt1 <- ggplot_gtable(ggplot_build(termMICdist))
      gt2 <- ggplot_gtable(ggplot_build(termMICcov))
      #newWidth = unit.pmax(gt1$widths[2:3], gt2$widths[2:3])
      #gt1$widths[2:3] = as.list(newWidth)
      #gt2$widths[2:3] = as.list(newWidth)
      
      
      print(fname)
      cairo_pdf(fname,width=12,height=9)
      grid.arrange(gt1, gt2, ncol=2,widths=c(3,1))
      #multiplot(termMICdist, termMICcov,  cols=2)
      dev.off()
    }
  }
}

# 
# p + theme(axis.line=element_blank(),axis.text.x=element_blank(),
#           axis.text.y=element_blank(),axis.ticks=element_blank(),
#           axis.title.x=element_blank(),
#           axis.title.y=element_blank(),legend.position="none",
#           panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
#           panel.grid.minor=element_blank(),plot.background=element_blank())

##Coverage
print('Coverage')
dir.create(paste(odir,'/',envar,'/Coverage',sep=''), showWarnings = TRUE, recursive = FALSE, mode = "0777")
for(m in 1:nrow(Thresholds)) {
  mr<-Thresholds[m,]
  for(s in 1:nrow(Stat)) {
    sr <- Stat[s,]
    for(t in 1:nrow(Types)) {
      tr <- Types[t,]
      sel<-subset(circms,Category==as.character(tr$tname)  &
                    Size_MICo5>gtres & ! is.na(Term) &
                    ! Term %in% c('Metabolic pathways',
                                  'Biosynthesis of secondary metabolites',
                                  'Microbial metabolism in diverse environments'))

      gtitle<-paste('Enriched E. Coli MG1655 ',tr$ttitle,' - ',mr$names,' MIC>',micthres,sep = '')
      fname<-paste(odir,'/',envar,'/Coverage/',tr$tabr,'_',mr$cols,'_MICo',micthres,'.pdf',sep = '')
      print(gtitle)
      print(fname)
      termBG<-ggplot(sel,
                     aes(y=sel[,as.character(mr$cols)],x=reorder(sel$Term,sel[,as.character(mr$cols)]),
                         fill=Size_MICo5))+
        geom_bar(stat = "identity")+
        coord_flip()+
        labs(fill='Number of genes')+
        scale_fill_gradient( high="red",low='grey')+
        ylab(mr$axis)+#ylim(0,mr$mscale)+
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
print('MIC distributions')
dir.create(paste(odir,'/',envar,'/MIC',sep=''), showWarnings = TRUE, recursive = FALSE, mode = "0777")
for(s in 1:nrow(Stat)) {
  sr <- Stat[s,]
  for(t in 1:nrow(Types)) {
    tr <- Types[t,]
    sel<-subset(enr,Category==as.character(tr$tname) & MIC>mtres & Size_MICo5>gtres&
                  MIC_S_Mean>mttres &!is.na(Term) &
                  ! Term %in% c('Metabolic pathways',
                                'Biosynthesis of secondary metabolites',
                                'Microbial metabolism in diverse environments'))
    #& logp_sum>2 & logp_sum>2  & logp_sum>2  & logp_sum>2
    gtitle<-paste('Enriched E. Coli MG1655 ',tr$ttitle,' - MIC distribution, ',sr$sabr,' MIC>',micthres,sep = '')
    fname<-paste(odir,'/',envar,'/MIC/',tr$tabr,'_MIC_',sr$sabr,'_MICo',micthres,'.pdf',sep = '')
    print(gtitle)
    print(fname)
    termMIC<-ggplot(sel,
                    aes(y=MIC,x=reorder(Term,MIC_S_Med),fill=Size_MICo5,color=MIC_S_Mean))+
      geom_boxplot(position='identity')+coord_flip()+
      scale_colour_gradientn(colours = c('black','orange','red'),
                             breaks=c(5,25,50),limits=c(5,50),guide='legend',name=micname)+
      labs(fill='Number of genes',color='MIC average')+
      scale_fill_gradient( high="red",low='white')+
      #scale_color_gradient(limits=c(0,100), high="blue",low='gray')+
      ylab(mic)+ylim(0,100)+
      xlab(tr$ttitle)+ggtitle(gtitle)
    cairo_pdf(fname,width=9,height=9)
    print(termMIC)
    dev.off()
    #      
  }
}





##Growth distributions
print('Growth distributions')
dir.create(paste(odir,'/',envar,'/Bacterial_growth',sep=''), showWarnings = TRUE, recursive = FALSE, mode = "0777")
for(m in 1:nrow(Medias)) {
  mr<-Medias[m,]
  for(s in 1:nrow(Stat)) {
    sr <- Stat[s,]
    for(t in 1:nrow(Types)) {
      tr <- Types[t,]
      #print(mr$media)
      #print(tr$tname)
      sel<-subset(enr,Category==as.character(tr$tname) & !is.na(enr[,as.character(mr$media)])  &
                    MIC>mtres & Size_MICo5>gtres & MIC_S_Mean>mttres &
                    !is.na(Term) &! Term %in% c('Metabolic pathways',
                                                 'Biosynthesis of secondary metabolites',
                                                 'Microbial metabolism in diverse environments'))
      gtitle<-paste('Enriched E. Coli MG1655 ',tr$ttitle,' - ',mr$mtitle,' distribution, ',sr$sabr,' MIC>',micthres,sep = '')
      fname<-paste(odir,'/',envar,'/Bacterial_growth/',tr$tabr,'_BG_',mr$mabr,'_',sr$sabr,'_MICo',micthres,'.pdf',sep = '')
      print(gtitle)
      print(fname)
      termBG<-ggplot(sel,aes(y=sel[,as.character(mr$media)],
                             x=reorder(sel$Term,sel[,as.character(mr$morder)]),
                             fill=Size_MICo5,color=MIC_S_Mean))+
        geom_boxplot(position='identity')+
        coord_flip()+labs(fill='Number of genes',color='MIC average')+
        scale_fill_gradient( high="red",low='white')+
        scale_color_gradient(limits=c(0,100), high="blue",low='gray')+
        ylab('OD')+ylim(0,mr$mscale)+
        xlab(tr$ttitle)+ggtitle(gtitle)
      termBG
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


Media2=c('LB_22hr','OD_T_Mean','OD_C_Mean')
MTitle2=c('LB 22hr growth','NGM + 100uM 5FU growth','NGM growth')
MAbr2=c('LB','NGM-5FU','NGM')
MOrder2=c('LB_S_Med','OD_T_Mean_S_Med','OD_C_Mean_S_Med')
MQvars2=c('LB','OD_T_Mean','OD_C_Mean')
MScale2=c(1,0.2,0.4)
Medias2<-data.frame(media=Media2,mtitle=MTitle2,mabr=MAbr2,morder=MOrder2,mscale=MScale2,mqvars=MQvars2)

##MIC vs Media
dir.create(paste(odir,'/',envar,'/MICvsOD',sep=''), showWarnings = TRUE, recursive = FALSE, mode = "0777")
print('MIC vs Media')
for(m in 1:nrow(Medias2)) {
  mr<-Medias2[m,]
  for(s in 1:nrow(Stat)) {
    sr <- Stat[s,]
    for(t in 1:nrow(Types)) {
      tr <- Types[t,]
      sel<-subset(enr,Category==as.character(tr$tname) &
                    !is.na(enr[,paste(mr$mqvars,'_S_Med',sep='')] ) &
                    MIC>mtres & Size_MICo5>gtres& MIC_S_Mean>mttres)

      gtitle<-paste('Enriched E. Coli MG1655 ',tr$ttitle,': MIC',' - ',mr$mtitle,' distribution, ',sr$sabr,' MIC>',micthres,sep = '')
      fname<-paste(odir,'/',envar,'/MICvsOD/',tr$tabr,'_MIC_vs_BG_',mr$mabr,'_',sr$sabr,'_MICo',micthres,'.pdf',sep = '')
      print(gtitle)
      print(fname)
      kpcor<-ggplot(sel,
                    aes(x=MIC_S_Med,y=sel[,paste(as.character(mr$mqvars),'_S_Med',sep='')],color=Size_MICo5))+
        geom_point()+
        geom_errorbarh(aes(xmax=MIC_S_Q75,
                           xmin=MIC_S_Q25),height=mr$mscale*0.0125,alpha=0.5)+
        geom_errorbar(aes(ymax=sel[,paste(mr$mqvars,'_S_Q75',sep='')],
                          ymin=sel[,paste(mr$mqvars,'_S_Q25',sep='')]),width=1,alpha=0.5)+
        geom_text(aes(label=Term),color='black',hjust=-0.1, vjust=-0.5,size=3)+
        scale_color_gradient(limits=c(0,8), high="blue",low='gray')+
        labs(color='Number of genes')+
        xlab(mic)+ylab('OD')+ylim(0,mr$mscale)+xlim(0,100)+
        ggtitle(gtitle)
      cairo_pdf(fname,width=9,height=9)
      print(kpcor)
      dev.off()
    }
  }
}

fitbac<-lm(OD_T_Mean ~ OD_C_Mean,data=bacmic)
confint(fitbac,'OD_C_Mean',level=0.95)

dir.create(paste(odir,'/',envar,'/NGMvsNGM-5FU',sep=''), showWarnings = TRUE, recursive = FALSE, mode = "0777")
print('NGMvsNGM-5FU')
for(s in 1:nrow(Stat)) {
  sr <- Stat[s,]
  for(t in 1:nrow(Types)) {
    tr <- Types[t,]
    
    sel<-subset(enr,Category==as.character(tr$tname) & MIC>mtres & Size_MICo5>gtres & MIC_S_Mean>mttres)
    #& logp_sum>2 & logp_sum>2  & logp_sum>2  & logp_sum>2
    gtitle<-paste('Enriched E. Coli MG1655 ',tr$ttitle,': NGM - NGM + 100um 5FU correlation, ',sr$sabr,' MIC>',micthres,sep = '')
    fname<-paste(odir,'/',envar,'/NGMvsNGM-5FU/',tr$tabr,'_NGM_vs_NGM-5FU_',sr$sabr,'_MICo',micthres,'.pdf',sep = '')
    print(gtitle)
    print(fname)
    kpcor<-ggplot(sel,aes(x=OD_C_Mean_S_Med,y=OD_T_Mean_S_Med,color=MIC_S_Mean))+
      geom_point(aes(size=Size_MICo5))+
      geom_abline(aes(fill='1:1'),intercept=0,slope=1,alpha=0.5,color='grey',linetype='longdash')+
      geom_abline(aes(fill='Trend for all knockouts'),intercept=coef(fitbac)[[1]],
                  slope=coef(fitbac)[[2]],alpha=0.5,color='red')+
      geom_errorbarh(aes(xmax=OD_C_Mean_S_Q75,xmin=OD_C_Mean_S_Q25),height=0.15*0.0125,alpha=0.5)+
      geom_errorbar(aes(ymax=OD_T_Mean_S_Q75,ymin=OD_T_Mean_S_Q25),width=0.2*0.0125,alpha=0.5)+
      geom_text(aes(label=Term),color='black',hjust=-0.1, vjust=-0.5,size=3)+
      labs(size='Number of genes',color='MIC average')+
      xlab('OD, NGM - Control')+
      ylab(expression(paste('OD, NGM + 100',mu,'M 5FU',sep='')))+
      ylim(0,0.4)+xlim(0,0.4)+
      ggtitle(gtitle)
    cairo_pdf(fname,width=9,height=9)
    print(kpcor)
    dev.off()
    #      scale_color_gradient(limits=c(0,8), high="blue",low='gray')+ ,fill='Guides'
  }
}
  
#}




