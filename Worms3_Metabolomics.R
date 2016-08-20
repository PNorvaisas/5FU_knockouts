library('ggplot2')
library('gplots')
library('plyr')
library('reshape2')
library(tidyr)
library(BSDA)
library(xlsx)
library(gridExtra)

library("cluster")
library(ggbiplot)
library(contrast)

theme_set(theme_light())

setwd("~/Projects/2015-Metformin/Worms")

ddir<-'Data_final'
odir<-'Figures_final/Metabolomics'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")


getinfo<-function(cof) {
  df<-data.frame(cof)
  df$Comparisons<-rownames(df)
  return(df)
}

.simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1, 1)), substring(s, 2),
        sep = "", collapse = " ")
}

read.and.clean<-function(fname,sheet='') {
  print(fname)
  #print(grepl('.xlsx',fname))
  #print(grepl('.csv',fname))
  if (grepl('xlsx',fname)) {
    dft<-read.xlsx2(fname,sheetName = sheet,
                    stringsAsFactors = FALSE,
                    header=FALSE)
  } else if (grepl('csv',fname)) {
    dft<-read.csv(fname,sep=',',quote = "\"",stringsAsFactors = FALSE,
                  header=FALSE)
    dft[1,which(dft[1,] == "")] <- "Metabolite Set"
  }
  #print(length(colnames(dft)))
  
  colnames(dft)<-dft[1,]
  df<-dft[-1,]
  
  return(df)
}

#Old
fls<-c('Metabolomics/078_bacteria_report_no_outlier.xlsx','Metabolomics/078_celegans_report.xlsx')
orgs<-c('E.coli','C.elegans')
orgsel<-data.frame('Organism'=orgs,'File'=fls)
orgsel<-apply(orgsel,2,as.character)

#QEA, ' (QEA)',
datn<-c('Diff')
dat<-c(' (Diff)')
datatp<-data.frame('Type'=datn,'Code'=dat)
datatp<-apply(datatp,2,as.character)

tp<-c('5FU/Control','Arabinose/Control','Arabinose+5FU/ArabinoseControl')
shtec<-c('C vs A','B vs A','D vs B')
shtce<-c('W vs Z','X vs Z','Y vs X')
colsel<-data.frame('Type'=tp,'E.coli'=shtec,'C.elegans'=shtce)
colsel<-apply(colsel,2,as.character)


alldat<-list()
#for (dt in 1:nrow(datatp)) {
dt=1
alldif<-data.frame()
datatype='Diff'#datatp[dt,'Type']
datacode=' (Diff)'#datatp[dt,'Code']
for (or in 1:nrow(orgsel)) {
  organism<-orgsel[or,'Organism']
  fl<-orgsel[or,'File']
  for (cs in 1:nrow(colsel)) {
    colcode<-colsel[cs,organism]
    coltype<-colsel[cs,'Type']
    print(fl)
    print(colcode)
    df<-read.and.clean(as.character(fl),paste(colcode,datacode,sep=''))
    df$Type<-coltype
    df$Organism<-organism
    
    if (dim(all)[1]==0) {
      alldif<-df
    } else {
      alldif<-merge(alldif,df,all=TRUE)
    }
    
  }
}
alldat[['Diff']]<-alldif





hclustfunc <- function(x) hclust(x, method="complete")
# #distfunc <- function(x) dist(x,method="euclidean")
# # Try using daisy GOWER function 
distfunc <- function(x) daisy(x,metric="gower")
reorderfun_mean = function(d,w) { reorder(d, w, agglo.FUN = mean) }
reorderfun_median = function(d,w) { reorder(d, w, agglo.FUN = median) }

gyr <- colorRampPalette(c("gray","yellow", "red"))(n = 32)
gyrs<-colorRampPalette(c("gray90","steelblue1","blue",'cyan'))(n = 5)



#QEA data

enr<-c('QEA','QEAKegg')
enfil<-c('_msea_qea_result.csv','_pathway_results.csv')
enscale<-c("FoldEnrichment","Impact")
enscnm<-c("Fold Enrichment","Impact")
enrichments<-data.frame('Type'=enr,'File'=enfil,'Scale'=enscale,'Name'=enscnm)
enrichments<-apply(enrichments,2,as.character)



fls<-c('Metabolomics/Bacteria_enrichment/','Metabolomics/Worms_enrichment/')
orgs<-c('E.coli','C.elegans')
orgsel<-data.frame('Organism'=orgs,'File'=fls)
orgsel<-apply(orgsel,2,as.character)

tp<-c('5FU/Control','Arabinose/Control','Arabinose+5FU/Control','Arabinose+5FU/Arabinose','Arabinose+5FU/5FU')
shtec<-c('5FUvsControl','LAravsControl','LAra5FUvsControl','LAra5FUvsLAra','LAra5FUvs5FU')
colsel<-data.frame('Type'=tp,'Abr'=shtec)
colsel<-apply(colsel,2,as.character)

for (eni in 1:nrow(enrichments)) {
  entype<-enrichments[eni,'Type']
  enfl<-enrichments[eni,'File']
  ensc<-enrichments[eni,'Scale']
  ensn<-enrichments[eni,'Name']
  #alldat<-list()
  allqea<-data.frame()
  for (or in 1:nrow(orgsel)) {
    organism<-orgsel[or,'Organism']
    fl<-orgsel[or,'File']
    for (cs in 1:nrow(colsel)) {
      colcode<-colsel[cs,'Abr']
      coltype<-colsel[cs,'Type']
      fln<-paste(fl,colcode,enfl,sep='')
      #print(fln)
      #print(colcode)
      df<-read.and.clean(as.character(fln))
      df$Type<-coltype
      df$Organism<-organism
      #print(colnames(df))
      if (dim(allqea)[1]==0) {
        allqea<-df
      } else {
        allqea<-merge(allqea,df,all=TRUE)
      }
      
    }
  }
  alldat[[entype]]<-allqea
  
  stats<-alldat[[entype]]
  
  stats[,2:8]<-apply(stats[,2:8],2,as.numeric)
  
  stats$Organism<-as.factor(stats$Organism)
  stats$Significant<-ifelse(stats$`Raw p`<0.05,TRUE,FALSE)
  stats$Significant<-factor(stats$Significant,levels = c(TRUE,FALSE),labels=c('True','False'))
  stats$logp<--log(stats$`Raw p`,10)
  
  pathexclude<-c("Intracellular Signalling Through Adenosine Receptor A2A And Adenosine | Intracellular Signalling Through Adenosine Receptor A2B And Adenosine",
                 "Intracellular Signalling Through Fsh Receptor And Follicle Stimulating Hormone | Intracellular Signalling Through Lhcgr Receptor And Luteinizing Hormone/Choriogonadotropin",
                 "INTRACELLULAR SIGNALLING THROUGH HISTAMINE H2 RECEPTOR AND HISTAMINE",
                 "INTRACELLULAR SIGNALLING THROUGH PGD2 RECEPTOR AND PROSTAGLANDIN D2",
                 "INTRACELLULAR SIGNALLING THROUGH PROSTACYCLIN RECEPTOR AND PROSTACYCLIN",
                 "BILE ACID BIOSYNTHESIS",
                 "VASOPRESSIN REGULATION OF WATER HOMEOSTASIS",
                 "CORTICOTROPIN ACTIVATION OF CORTISOL PRODUCTION",
                 "CATECHOLAMINE BIOSYNTHESIS",
                 "INSULIN SIGNALLING",
                 "Novobiocin biosynthesis",
                 "Streptomycin biosynthesis",
                 "EXCITATORY NEURAL SIGNALLING THROUGH 5-HTR 4 AND  SEROTONIN | EXCITATORY NEURAL SIGNALLING THROUGH 5-HTR 7 AND  SEROTONIN  | EXCITATORY NEURAL SIGNALLING THROUGH 5-HTR 6 AND SEROTONIN ")
  
  pathexclude2<-toupper(pathexclude)
  stats<-subset(stats,!`Metabolite Set` %in% pathexclude)
  stats<-subset(stats,!`Metabolite Set` %in% pathexclude2)
  
  if (entype=='QEA') {
    stats$`-log(p)`<-NULL
    stats$FoldEnrichment<-stats$`Statistic Q`/stats$`Expected Q`
    
    stats$`Metabolite Set`<-tolower(stats$`Metabolite Set`)
    stats$`Metabolite Set`<-paste(toupper(substr(stats$`Metabolite Set`, 1, 1)), substr(stats$`Metabolite Set`, 2, nchar(stats$`Metabolite Set`)), sep="")
    
    
    metsubs<-c(" b6 "=" B6 ","Rna"="RNA"," coa "=" CoA ")
    
    for (i in names(metsubs)) {
      stats$`Metabolite Set`<-gsub(i,metsubs[[i]],stats$`Metabolite Set`)
    }
  }
  #
  
  write.csv(stats,paste(ddir,'/Metabolomics_All_',entype,'.csv',sep=''),row.names = FALSE)
  
  print("Data written!")
  cstats<-stats
  metabr<-c(" metabolism"=" met."," biosynthesis"=" biosyn."," degradation"=" deg.")
  
  for (i in names(metabr)) {
    cstats$`Metabolite Set`<-gsub(i,metabr[[i]],cstats$`Metabolite Set`)
  }
  
  print("Generating overview plots!")
  cstatss<-subset(cstats,`Raw p`<0.05)
  overs<-ggplot(cstatss,
         aes_string(y=ensc,x="`Metabolite Set`"))+#,fill=-log10(P.value)
    geom_bar(stat = "identity",fill='red')+
    coord_flip()+
    xlab('Metabolite set')+
    ylab(ensn)+
    facet_grid(Type~Organism)
#   dev.copy2pdf(device=cairo_pdf,
#                file=paste(odir,"/Metabolomics_",entype,"_Significant.pdf",sep=''),
#                width=12,height=20)
  cairo_pdf(paste(odir,"/Metabolomics_",entype,"_Significant.pdf",sep=''),width=12,height=20)
  print(overs)
  dev.off()
  
  overa<-ggplot(cstats,
         aes_string(y=ensc,x="`Metabolite Set`",alpha="Significant"))+#,fill=-log10(P.value)
    geom_bar(stat = "identity",fill='red')+
    scale_alpha_discrete(range=c(1,0.2))+
    coord_flip()+
    xlab('Metabolite set')+
    ylab(ensn)+
    facet_grid(Type~Organism)
#   dev.copy2pdf(device=cairo_pdf,
#                file=paste(odir,"/Metabolomics_",entype,"_All.pdf",sep=''),
#                width=12,height=25)
  cairo_pdf(paste(odir,"/Metabolomics_",entype,"_All.pdf",sep=''),width=12,height=25)
  print(overa)
  dev.off()
  
  
  tpcol<-c('5FU/Control','Arabinose/Control','Arabinose+5FU/Control','Arabinose+5FU/Arabinose','Arabinose+5FU/5FU')
  tpabr<-c('5FUvsC','LAravsC','LAra5FUvsC','LAra5FUvsLAra','LAra5FUvs5FU')
  types<-data.frame(column=tpcol,abreviation=tpabr)
  
  
  print("Generating individual plots!")
  for(t in 1:nrow(types)){
    tp<-types[t,]
    sel<-subset(cstats,Type==tp$column )
    fname=paste(odir,"/Metabolomics_",entype,"_",tp$abreviation,"_All.pdf",sep='')
    fac<-ggplot(sel,
                aes_string(y=ensc,x="`Metabolite Set`",alpha="Significant"))+
      geom_bar(stat = "identity",fill='red')+
      scale_alpha_discrete(range=c(1,0.2))+
      #     scale_fill_gradientn(colours = c('orange','red'),
      #                          breaks=c(1.3,2),limits=c(1.3,7))+
      coord_flip()+
      ggtitle(tp$column)+
      xlab('Metabolite set')+
      ylab(ensn)+
      facet_grid(.~Organism)
    fac
    
    print(fname)
    cairo_pdf(fname,width=12,height=9)
    print(fac)
    #grid.arrange(gt1, gt2, ncol=2,widths=c(3,1))
    #multiplot(termMICdist, termMICcov,  cols=2)
    dev.off()
    
  }
  
  
  for(t in 1:nrow(types)){
    tp<-types[t,]
    sel<-subset(cstats,Type==tp$column & `Raw p`<0.05)
    fname=paste(odir,"/Metabolomics_",entype,"_",tp$abreviation,"_Significant.pdf",sep='')
    fac<-ggplot(sel,
                aes(y=ensc,x="`Metabolite Set`"))+
      geom_bar(stat = "identity",fill='red')+
      #     scale_fill_gradientn(colours = c('orange','red'),
      #                          breaks=c(1.3,2),limits=c(1.3,7))+
      coord_flip()+
      ggtitle(tp$column)+
      xlab('Metabolite set')+
      ylab(ensn)+
      facet_grid(.~Organism)
    fac
    
    print(fname)
    cairo_pdf(fname,width=12,height=9)
    print(fac)
    #grid.arrange(gt1, gt2, ncol=2,widths=c(3,1))
    #multiplot(termMICdist, termMICcov,  cols=2)
    dev.off()
    
  }
  
  print("Generating heatmaps!")
  #Heatmap for Enrichment
  #has to be removed do to problems with data type
  statsc<-cstats[,!colnames(cstats) %in% c('Significant')]
  
  
  #statscf<-statsc
  
  statsm<-melt(statsc,id.vars = c('Metabolite Set','Total Cmpd','Hits','Organism','Type'),
               value.name = 'Value',variable.name = 'Measure', na.rm = FALSE)
  statsm$Value<-as.numeric(statsm$Value)
  
  
  statsmf<-subset(statsm,Measure=='logp')
  statsall<-dcast(statsmf,`Metabolite Set`~Organism+Type+Measure,value.var = 'Value')
  

#   hexclude<-c('Metabolite Set')
#   cnm<-setdiff(colnames(statsall),hexclude)
  
  if (entype=='QEA'){
    cnm<-c("C.elegans_5FU/Control_logp",
           "C.elegans_Arabinose/Control_logp",
           "C.elegans_Arabinose+5FU/Control_logp",
           "E.coli_Arabinose/Control_logp")
#     ,
#     "C.elegans_Arabinose+5FU/5FU_logp",
#     "C.elegans_Arabinose+5FU/Arabinose_logp"
#     "E.coli_Arabinose+5FU/5FU_logp",
#     "E.coli_Arabinose+5FU/Control_logp"
  } else if (entype=='QEAKegg') {
    cnm<-c("C.elegans_5FU/Control_logp",
           "C.elegans_Arabinose/Control_logp",
           "C.elegans_Arabinose+5FU/Control_logp",
           "E.coli_5FU/Control_logp",
           "E.coli_Arabinose/Control_logp",
           "E.coli_Arabinose+5FU/Control_logp")
#     "C.elegans_Arabinose+5FU/5FU_logp",
#     "C.elegans_Arabinose+5FU/Arabinose_logp",
#     "E.coli_Arabinose+5FU/5FU_logp",
#     "E.coli_Arabinose+5FU/Arabinose_logp"
  }
  
  
  hdata<-statsall[!apply(statsall[,cnm], 1, function(x)(all(x< -log(0.05,10),na.rm=TRUE))),]
  rnm<-hdata$`Metabolite Set`
  hdata<-hdata[,cnm]
  
  
  hdatafill<-hdata
  
  
  hdatafill[is.na(hdatafill)]<-0
  
  
  hmap<-heatmap.2(as.matrix(hdatafill),key=TRUE,Colv=FALSE,trace='none',labRow=rnm,
                  labCol=cnm,col=gyrs,
                  xlab='Comparison',Rowv=TRUE,
                  dendrogram="row",scale="none",na.color="white",
                  cexRow=0.8,cexCol=0.5,margin=c(8,16),
                  lwid=c(0.2,0.8),symkey=FALSE,
                  reorderfun=reorderfun_mean)
  
  
  #cairo_pdf(paste(odir,'/Metabolomics_',entype,'_heatmap.pdf',sep = ''),width=8,height=15)
  hmapr<-heatmap.2(as.matrix(hdata),key=TRUE,Colv=FALSE,trace='none',labRow=rnm,
            labCol=cnm,col=gyrs,
            xlab='Comparison',Rowv=hmap$rowDendrogram,
            dendrogram='row',scale="none",na.color="white",
            cexRow=0.7,cexCol=0.7,margin=c(16,16),
            lwid=c(0.2,0.8),symkey=FALSE,
            reorderfun=reorderfun_mean,
            breaks=c(0,-log(0.05,10),2,3,4,5))
  #dev.off()
  dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Metabolomics_',entype,'_heatmap.pdf',sep = ''),
               width=8,height=17)
  
  
}



#Data for comparison

diffs<-alldat$Diff
diffm<-melt(diffs,id.vars = c('ID','Name','Type','Organism'),variable.name = 'Stats',value.name = 'Value')

dif<-dcast(diffm,ID+Name~Organism+Type+Stats,value.var = 'Value')
dif[,3:37]<-apply(dif[,3:37],2,as.numeric)


#Raw data

#Read raw data
bac<-read.and.clean('Metabolomics/078_bacteria_report_no_outlier.xlsx','Data (normalized)')
worm<-read.and.clean('Metabolomics/078_celegans_report.xlsx','Data (normalized)')



gchange<-c('A'='Control','B'='LAra','C'='5FU','D'='LAra5FU','Z'='Control','X'='LAra','W'='5FU','Y'='LAra5FU')
bac$Group<-revalue(bac$Group, gchange)
worm$Group<-revalue(worm$Group, gchange)
bac$ID<-paste(bac$Group,bac$Replicate,sep='-')
worm$ID<-paste(worm$Group,worm$Replicate,sep='-')


exclude<-c('Taurochenodesoxycholic acid','Chenodeoxycholic acid','Taurocholic acid','Glycocholic acid')
bacc<-bac[,!colnames(bac) %in% exclude]
wormc<-worm[,!colnames(worm) %in% exclude]

bacw<-bacc[,!colnames(bacc) %in% c('Replicate','5-FU','L-Arabinose')]
wormw<-wormc[,!colnames(wormc) %in% c('Replicate','5-FU','L-Arabinose')]

write.csv(bacw,paste('Metabolomics','/Bacteria_raw_clean.csv',sep=''),row.names = FALSE)
write.csv(wormw,paste('Metabolomics','/Worms_raw_clean.csv',sep=''),row.names = FALSE)

#Change to bac and worm if all data needs to be analysed
bacm<-melt(bacc,id.vars = c('ID','Group','Replicate','5-FU','L-Arabinose'),variable.name = 'Metabolite',value.name = 'Conc')
wormm<-melt(wormc,id.vars = c('ID','Group','Replicate','5-FU','L-Arabinose'),variable.name = 'Metabolite',value.name = 'Conc')
bacm$Organism<-'E.coli'
wormm$Organism<-'C.elegans'

mets<-merge(bacm,wormm,all.x = TRUE,all.y=TRUE)
mets$Organism<-factor(mets$Organism,levels=c('E.coli','C.elegans'))
mets$Metabolite<-as.factor(mets$Metabolite)
mets$`5-FU`<-as.factor(mets$`5-FU`)
mets$`L-Arabinose`<-as.factor(mets$`L-Arabinose`)
mets$Conc<-as.numeric(mets$Conc)
mets$logConc<-log(mets$Conc,2)

write.csv(mets,paste(ddir,'/Metabolomics_Clean_Data.csv',sep=''),row.names = FALSE)





#Linear modeling
################
#Data for comparisons

allmets<-unique(as.character(mets$Metabolite))

mets$GroupG<-paste(mets$`5-FU`,mets$`L-Arabinose`,sep='')
metsg<-dcast(mets,Metabolite+Replicate+`5-FU`+`L-Arabinose` ~Organism,mean,value.var = c('logConc'),fill = as.numeric(NA),drop=FALSE)


metslm<-dcast(mets,Replicate+`5-FU`+`L-Arabinose`+Organism ~Metabolite,mean,value.var = c('logConc'),fill = as.numeric(NA),drop=FALSE)
metslm$`L-Arabinose`<-relevel(metslm$`L-Arabinose`, ref = "0mM")
metslm$`5-FU`<-relevel(metslm$`5-FU`, ref = "0uM")

ecolimiss<-as.character(unique(subset(metsg,is.na(E.coli) & !is.na(C.elegans) & Replicate!=2)$Metabolite))
celegansmiss<-as.character(unique(subset(metsg,is.na(C.elegans) & !is.na(E.coli))$Metabolite))

#PCA analysis

# 
# library(devtools)
# library(ggbiplot)
# install_github("ggbiplot", "vqv")
#General
metslm$Group<-paste(metslm$Organism,'_5-FU:',metslm$`5-FU`,'_L-Arabinose:',metslm$`L-Arabinose`,sep='')

metslms<-subset(metslm,!(Replicate==2 & Organism=='E.coli' & `5-FU`=='1uM' & `L-Arabinose`=='10mM'))
metsfilt<-metslms[,!colnames(metslms) %in% c(ecolimiss,celegansmiss)]

pca.group<-metsfilt$Group
pca.dat<-metsfilt[,!colnames(metsfilt) %in% c('Replicate','Organism','5-FU','L-Arabinose','Group')]

ir.pca <- prcomp(pca.dat,
                 center = TRUE,
                 scale. = TRUE,
                 na.action=na.omit) 
summary(ir.pca)
plot(ir.pca,type='l')
generalpca <- ggbiplot(ir.pca, obs.scale = 1,
                       var.scale = 1,
                       groups = pca.group,
                       ellipse = TRUE,
                       circle = TRUE,
                       var.axes = 0)+
  scale_color_discrete(name = '')+
  ylim(-10,10)+
  theme(legend.direction = 'vertical',legend.position = 'right')
generalpca

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Metabolomics_PCA.pdf",sep=''),
             width=12,height=9)

#C.elegans
metslms<-subset(metslm,Organism=='C.elegans')
metsfilt<-metslms[,!colnames(metslms) %in% c(celegansmiss)]

pca.group<-metsfilt$Group
pca.dat<-metsfilt[,!colnames(metsfilt) %in% c('Replicate','Organism','5-FU','L-Arabinose','Group')]

#PCa by default centers and scales
ir.pca <- prcomp(pca.dat,
                 center = TRUE,
                 scale. = TRUE,
                 na.action=na.omit) 
summary(ir.pca)
plot(ir.pca,type='l')
celeganspca <- ggbiplot(ir.pca, obs.scale = 1,
                       var.scale = 1,
                       groups = pca.group,
                       ellipse = TRUE,
                       circle =TRUE,
                       var.axes = 0)+
  ylim(-10,10)+
  scale_color_discrete(name = '')+
  theme(legend.direction = 'vertical',legend.position = 'right')
celeganspca

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Metabolomics_PCA_Celegans.pdf",sep=''),
             width=12,height=9)

#E.coli
metslms<-subset(metslm,Organism=='E.coli' & !(Replicate==2 & Organism=='E.coli' & `5-FU`=='1uM' & `L-Arabinose`=='10mM'))
metsfilt<-metslms[,!colnames(metslms) %in% c(ecolimiss)]

pca.group<-metsfilt$Group
pca.dat<-metsfilt[,!colnames(metsfilt) %in% c('Replicate','Organism','5-FU','L-Arabinose','Group')]
pca.dat.c<-scale(pca.dat,scale=T,center=T)

#Centering and scaling by default
ir.pca <- prcomp(pca.dat,
                 center = TRUE,
                 scale. = TRUE,
                 na.action=na.omit) 

summary(ir.pca)

# 
# ir.pca.noc <- prcomp(pca.dat.c,
#                  center = FALSE,
#                  scale. = FALSE,
#                  na.action=na.omit) 
# summary(ir.pca.noc)
# plot(ir.pca.noc,type='l')
# 
# 
# ecolipcac <- ggbiplot(ir.pca.noc, obs.scale = 1,
#                      var.scale = 1,
#                      groups = pca.group,
#                      ellipse = TRUE,
#                      circle = TRUE,
#                      var.axes = 0)+
#   ylim(-10,10)+
#   scale_color_discrete(name = '')+
#   theme(legend.direction = 'vertical',legend.position = 'right')
# ecolipcac


ecolipca <- ggbiplot(ir.pca, obs.scale = 1,
                        var.scale = 1,
                        groups = pca.group,
                        ellipse = TRUE,
                        circle = TRUE,
                        var.axes = 0)+
  ylim(-10,10)+
  scale_color_discrete(name = '')+
  theme(legend.direction = 'vertical',legend.position = 'right')
ecolipca

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Metabolomics_PCA_Ecoli.pdf",sep=''),
             width=12,height=9)


#Linear modeling
###############
#General comparison

model<-lm(logConc~`5-FU`*`L-Arabinose`*Organism,data=mets,na.action =na.exclude)#
result<-summary(model)
result



#Full comparison

cleanmets<-setdiff(allmets,union(ecolimiss,celegansmiss))

metstring<-paste(cleanmets,collapse='`,`')#[1:10]

formula<-as.formula(paste("cbind(`",metstring,"`) ~ `5-FU`*`L-Arabinose`*Organism",sep=''))
model<-lm(formula,data=metslm,na.action =na.exclude)#
result<-summary(model)




lresultsgen<-ldply(result, function(x) getinfo(x$coefficients))
lresultsgen$Type<-'Interorganism'

#E.coli specific
metslmeco=subset(metslm,Organism=='E.coli')


ecolimets<-setdiff(allmets,ecolimiss)
metstring<-paste(ecolimets,collapse='`,`')#[1:10]

formula<-as.formula(paste("cbind(`",metstring,"`) ~ `5-FU`*`L-Arabinose`",sep=''))
model<-lm(formula,data=metslmeco,na.action =na.exclude)#
result<-summary(model)

#Get main results
lresultseco<-ldply(result, function(x) getinfo(x$coefficients))
lresultseco$Type<-'E.coli'


formula<-as.formula(paste("cbind(`",metstring,"`) ~as.factor(`5-FU`:`L-Arabinose`)",sep=''))
model<-lm(formula,data=metslmeco,na.action =na.exclude)#
result<-summary(model)

#Get additional results
lresultsecob<-ldply(result, function(x) getinfo(x$coefficients))
lresultsecob$Type<-'E.coli'


#C.elegans specific
metslmcel=subset(metslm,Organism=='C.elegans')

celegansmets<-setdiff(allmets,celegansmiss)
metstring<-paste(celegansmets,collapse='`,`')#[1:10]

formula<-as.formula(paste("cbind(`",metstring,"`) ~ `5-FU`*`L-Arabinose`",sep=''))
model<-lm(formula,data=metslmcel,na.action =na.exclude)#
result<-summary(model)

#Get main results
lresultscel<-ldply(result, function(x) getinfo(x$coefficients))
lresultscel$Type<-'C.elegans'

formula<-as.formula(paste("cbind(`",metstring,"`) ~ as.factor(`5-FU`:`L-Arabinose`)",sep=''))
model<-lm(formula,data=metslmcel,na.action =na.exclude)#
result<-summary(model)

#Get additional results
lresultscelb<-ldply(result, function(x) getinfo(x$coefficients))
lresultscelb$Type<-'C.elegans'

#Get bonus results together

comptb<-merge(lresultsecob,lresultscelb,all = TRUE)

#Get main results together
compt<-merge(lresultseco,lresultscel,all = TRUE)

#Add general results
allcomp<-merge(compt,lresultsgen,all=TRUE)

#Add additional results
allcompfull<-merge(allcomp,comptb,all=TRUE)


allcompfull$Metabolite<-gsub('Response ','',allcompfull$.id)
allcompfull$.id<-NULL
allcompfull<-rename(allcompfull,c("Pr...t.."='p.value',"Std..Error"='SD','Estimate'='Difference'))


compchange<-c("`5-FU`1uM"="5FU vs Control",
              "`L-Arabinose`10mM"="L-Arabinose vs Control",
              "as.factor(`5-FU`:`L-Arabinose`)1uM:10mM"="5FU & L-Arabinose vs Control",
              "`5-FU`1uM:`L-Arabinose`10mM"="Interaction: 5FU & L-Arabinose",
              "OrganismC.elegans"="C.elegans vs E.coli")
allcompfull$Comparisons<-revalue(allcompfull$Comparisons,compchange)






#Melt everything
# allcompm<-melt(allcompclean,id.vars = c('Metabolite','Type','Comparisons'),
#                variable.name = 'Stat',value.name = 'Value')
# 
# fullstatt<-dcast(allcompm,.id~Type+Comparisons+Stat,value.var = 'Value',drop=TRUE,fill=as.numeric(NA))
# fullstat<-merge(fullstatt,dif,by.x=c('Metabolite'),by.y=c('Name'))

#fullstat<-merge(fullstata,bstatt,by='Metabolite')

# 
# #C.elegans 5-FU
# ggplot(fullstat,aes(x=fullstat$'C.elegans_5FU/Control_Difference',
#                     y=fullstat$'C.elegans_`5-FU`1uM_Estimate'))+
#   geom_point()
# #C. elegans L-Arabinose
# ggplot(fullstat,aes(x=fullstat$`C.elegans_Arabinose/Control_Difference`,
#                     y=fullstat$"C.elegans_`L-Arabinose`10mM_Estimate"))+
#   geom_point()
# 
# #E. coli 5-FU
# ggplot(fullstat,aes(x=fullstat$'E.coli_5FU/Control_Difference',
#                     y=fullstat$'E.coli_`5-FU`1uM_Estimate'))+
#   geom_point()
# #E.coli L-Arabinose
# ggplot(fullstat,aes(x=fullstat$`E.coli_Arabinose/Control_Difference`,
#                     y=fullstat$"E.coli_`L-Arabinose`10mM_Estimate"))+
#   geom_point()
# 

selcol<-c("C.elegans_5FU vs Control",
          "C.elegans_L-Arabinose vs Control",
          "C.elegans_5FU & L-Arabinose vs Control",
          "E.coli_5FU vs Control",
          "E.coli_L-Arabinose vs Control",
          "E.coli_5FU & L-Arabinose vs Control")



cleannames<-gsub('_',' ',selcol)

allcompclean<-subset(allcompfull,Type!='Interorganism' & Comparisons!='(Intercept)')

allcompclean$FDR<-p.adjust(allcompclean$p.value,method = 'fdr')


allstat<-dcast(allcompclean,Metabolite~Type+Comparisons,value.var = 'Difference',drop=FALSE,fill=as.numeric(NA))
allstat.pval<-dcast(allcompclean,Metabolite~Type+Comparisons,value.var = 'p.value',drop=FALSE,fill=as.numeric(NA))
allstat.fdr<-dcast(allcompclean,Metabolite~Type+Comparisons,value.var = 'FDR',drop=FALSE,fill=as.numeric(NA))


#Draw heatmap


#allstatsel<-allstat#subset(allstat,!is.na(SD))
Names<-allstat$Metabolite


statsh<-apply(allstat[,selcol], 2, as.numeric )
statsh.p<-apply(allstat.pval[,selcol], 2, as.numeric )
statsh.fdr<-apply(allstat.fdr[,selcol], 2, as.numeric )

# 
# hclustfunc <- function(x) hclust(x, method="complete")
# # #distfunc <- function(x) dist(x,method="euclidean")
# # # Try using daisy GOWER function 
#distfunc <- function(x) daisy(x,metric="gower")
#distfunc<-function(x) as.dist((1-cor(t(x)))/2)

#Prepare data for heatmap
nstep<-14

bgg <- colorRampPalette(c("blue", "gray", "green"))(n = nstep)
cs<-colorRampPalette(c("cornsilk"))(n = 1)

colscale<-c(cs,bgg)

brkst<-seq(-3.5,3.5,by=7/nstep)
brks<-c(-6,brkst)


reorderfun_mean = function(d,w) { reorder(d, w, agglo.FUN = mean) }
reorderfun_median = function(d,w) { reorder(d, w, agglo.FUN = median) }


statsfill<-statsh

statsfill[is.na(statsfill)]<-0


statsh.pc<-statsh

statsh.fc<-statsh

statsh.pc[statsh.p>0.05]<- -5
statsh.fc[statsh.fdr>0.05]<- -5




hmap<-heatmap.2(as.matrix(statsfill),key=TRUE,Colv=FALSE,trace='none',labRow=Names,
          labCol=cleannames,col=bgg,
          xlab='Comparison',Rowv=TRUE,
          dendrogram="row",scale="none",na.color="plum",
          cexRow=0.4,cexCol=1,margin=c(8,16),
          lwid=c(0.2,0.8),symkey=TRUE,
          reorderfun=reorderfun_mean)



#Rowv=dendrogram

heatmap.2(as.matrix(statsh),key=TRUE,Colv=FALSE,trace='none',labRow=Names,
          labCol=cleannames,col=bgg,
          xlab='Comparison',Rowv=hmap$rowDendrogram,
          dendrogram='row',scale="none",na.color="white",
          cexRow=0.4,cexCol=0.4,margin=c(16,16),
          lwid=c(0.2,0.8),symkey=FALSE,
          reorderfun=reorderfun_mean,breaks=brkst)

dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Metabolomics_Stats_All.pdf',sep = ''),
             width=6,height=15)

heatmap.2(as.matrix(statsh.pc),key=TRUE,Colv=FALSE,trace='none',labRow=Names,
          labCol=cleannames,col=colscale,
          xlab='Comparison',Rowv=hmap$rowDendrogram,
          dendrogram='row',scale="none",na.color="white",
          cexRow=0.4,cexCol=0.4,margin=c(16,16),
          lwid=c(0.2,0.8),symkey=FALSE,
          reorderfun=reorderfun_mean,breaks=brks)

dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Metabolomics_Stats_Significant_p.pdf',sep = ''),
             width=6,height=15)

heatmap.2(as.matrix(statsh.fc),key=TRUE,Colv=FALSE,trace='none',labRow=Names,
          labCol=cleannames,col=colscale,
          xlab='Comparison',Rowv=hmap$rowDendrogram,
          dendrogram='row',scale="none",na.color="white",
          cexRow=0.4,cexCol=0.4,margin=c(16,16),
          lwid=c(0.2,0.8),symkey=FALSE,
          reorderfun=reorderfun_mean,breaks=brks)

dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Metabolomics_Stats_Significant_FDR.pdf',sep = ''),
             width=6,height=15)

