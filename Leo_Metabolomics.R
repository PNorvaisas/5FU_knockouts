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
theme_update(panel.background = element_rect(colour = "black"),
             axis.text = element_text(colour = "black"))

setwd("~/Projects/2015-Metformin/Worms")

ddir<-'Data_final'
odir<-'Figures_final/Nucleotide_metabolomics'
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

#setwd("~/Projects/2015-Metformin/Worms/Leo_metabolomics")

ecoliraw1<-read.and.clean('Leo_metabolomics/Bacteria MS nucleotide data_BW, ndk, pyrE_13_09_16 v2.xlsx','Clean')
yjjg<-read.and.clean('Leo_metabolomics/Bacteria MS nucleotide data_BW, yjjg (batch 5)_22_09_16.xlsx','Clean')


ecoliraw<-merge(ecoliraw1,yjjg,by='Compound',all=TRUE)
#ecoliraw<-ecoliraw1


ecoliM<-melt(ecoliraw,id.vars = c('Class','Type','Group','Compound'),
             value.name = 'Conc',variable.name = 'Experiment')
ecoliM$Experiment<-as.character(ecoliM$Experiment)
ecoliM$Conc<-as.numeric(ecoliM$Conc) 
ecoliM$LogConc<-log(ecoliM$Conc,2)
ecoliM$Compound<-trimws(ecoliM$Compound)

ecoliM$Replicate<-substr(ecoliM$Experiment,nchar(ecoliM$Experiment),nchar(ecoliM$Experiment))
ecoliM$Strain<-substr(ecoliM$Experiment,1,nchar(ecoliM$Experiment)-1)
head(ecoliM)


ecoliM$Class<-as.factor(ecoliM$Class)
ecoliM$Type<-as.factor(ecoliM$Type)
ecoliM$Strain<-as.factor(ecoliM$Strain)
ecoliM$Compound<-as.factor(ecoliM$Compound)


ecoliM$Batch<-1
ecoliM$Batch<-ifelse(ecoliM$Strain=='yjjG',2,ecoliM$Batch)
ecoliM$Batch<-ifelse(ecoliM$Replicate %in% c(5,6,7,8),2,ecoliM$Batch)


#Filter data
#ecoliM<-subset(ecoliM,Batch==1)

#ecoliM<-subset(ecoliM,Batch==2)


#General comparison
#Strain
fitStrain<-lm(LogConc~Strain,ecoliM)
summary(fitStrain)

#No significant differences in general compounds levels

#Group
fitGroup<-lm(LogConc~0+Group,ecoliM)
summary(fitGroup)

ecoliM$Group<-relevel(ecoliM$Group, ref = "A")
fitGroup<-lm(LogConc~Group,ecoliM)
summary(fitGroup)


#Class
fitClass<-lm(LogConc~0+Class,ecoliM)
summary(fitClass)


ecoliM$Class<-relevel(ecoliM$Class, ref = "Nucleobases")
fitClass<-lm(LogConc~Class,ecoliM)
summary(fitClass)


#Type
fitType<-lm(LogConc~0+Type,ecoliM)
summary(fitType)

ecoliM$Type<-relevel(ecoliM$Type, ref = "Nucleobases")
fitType<-lm(LogConc~Type,ecoliM)
summary(fitType)


fitSClass<-lm(LogConc~Strain*Class,ecoliM)
summary(fitSClass)

fitSType<-lm(LogConc~Strain*Type,ecoliM)
summary(fitSType)

fitSGroup<-lm(LogConc~Strain*Group,ecoliM)
summary(fitSGroup)





metslm<-dcast(ecoliM,Batch+Replicate+Strain+Experiment~Compound,mean,value.var = c('LogConc'),fill = as.numeric(NA),drop=TRUE)
cols<-setdiff(colnames(metslm),c('Batch','Replicate','Strain','Experiment'))
metslmw<-metslm[,c('Experiment','Strain','Batch',cols)]


write.csv(metslmw,paste(ddir,'/Nucleotide_metabolomics_Ecoli_Raw_LogConc.csv',sep=''),row.names = FALSE)

#Find compounds with missing values
miss<-apply(metslm, 2, function(x) any(is.na(x)))
missing<-names(miss[miss==TRUE])

pca.group<-metslm$Strain

pca.dat<-metslm[,!colnames(metslm) %in% c('Replicate','Strain','Experiment','Batch',missing)]

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
  theme(legend.direction = 'vertical',legend.position = 'right')
generalpca

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Nucleotide_metabolomics_PCA.pdf",sep=''),
             width=12,height=9)




allmets<-unique(as.character(ecoliM$Compound))
cleanmets<-setdiff(allmets,missing)
metstring<-paste(cleanmets,collapse='`,`')#[1:10]

formula<-as.formula(paste("cbind(`",metstring,"`) ~ Strain",sep=''))
batches<-unique(metslm$Batch)
genresults<-data.frame()
references<-unique(metslm$Strain)
#Results collected in each batch separately




for (bat in batches) {
  for (rf in references) {
    metslmsel<-subset(metslm,Batch==bat)
    metslmsel$Strain<-relevel(metslmsel$Strain,ref = rf)
    model<-lm(formula,metslmsel)#,na.action =na.exclude
    result<-summary(model)
    #print(result)
    resultseco<-ldply(result, function(x) getinfo(x$coefficients))
    resultseco$Reference<-rf
    if (dim(genresults)[[1]]==0) {
      genresults<-resultseco
    } else {
      genresults<-merge(genresults,resultseco,all=TRUE)
    }
  }
}

misresults<-data.frame()
for (met in missing) {
  for (bat in batches) {
    for (rf in references) {
      print(met)
      metslmsel<-subset(metslm,Batch==bat)
      metslmsel$Strain<-relevel(metslmsel$Strain,ref = rf)
      formulamet<-as.formula(paste("`",met,"` ~ Strain",sep=''))
      misres<-summary(lm(formulamet,metslmsel))
      missum<-getinfo(misres$coefficients)
      missum$.id<-paste('Response ',met,sep='')
      missum$Reference<-rf
      if (dim(misresults)[[1]]==0) {
        misresults<-missum
      } else {
        misresults<-merge(misresults,missum,all=TRUE)
      }
    }
  }
}

allresults<-genresults#merge(genresults,misresults,all=TRUE)


allresults$Metabolite<-gsub('Response ','',allresults$.id)
allresults$.id<-NULL
allresults<-rename(allresults,c("Pr...t.."='p.value',
                                "Std..Error"='SE',
                                'Estimate'='Difference'))


allresults$Comparisons<-gsub('Strain','',allresults$Comparisons)

fullinfo<-merge(ecoliraw[,1:4],allresults,by.x = 'Compound',by.y='Metabolite')
fullinfo$FoldChange<-ifelse(fullinfo$Difference>0,2^fullinfo$Difference,-(2^(-fullinfo$Difference)))

fullinfo$FDR<-p.adjust(fullinfo$p.value,method = 'fdr')
fullinfo<-rename(fullinfo,c("Comparisons"="Target"))
fullinfo<-rename(fullinfo,c('Difference'='logFC'))

fullinfo$Comparison<-paste(fullinfo$Target,'/',fullinfo$Reference,sep='')

write.csv(fullinfo,paste(ddir,'/Nucleotide_metabolomics_Ecoli_Statistics.csv',sep=''),row.names = FALSE)

multcomp<-subset(fullinfo,Comparison %in% c('pyrE/BW','ndk/BW','ndk/pyrE'))
multcomp$FDR<-p.adjust(multcomp$p.value,method = 'fdr')


multcompm<-melt(multcomp,id.vars = c('Class','Type','Group','Compound',
                                     'Comparison','Target','Reference'),
                variable.name = 'Stat',value.name = 'Value')

multcompcast<-dcast(multcompm,Type+Class+Group+Compound~Comparison+Stat,mean,value.var = c('Value'),
                    fill = as.numeric(NA),drop=TRUE)

maincomp<-'ndk/pyrE logFC'
brks<-c(-4,-2,0,2,4)

gradcols<-c('blue','blue','black','red','red')

erralpha<-1
errcolor<-'grey80'

ggplot(multcompcast,aes(y=`ndk/BW_logFC`,x=`pyrE/BW_logFC`))+
  geom_vline(xintercept = 0,color='red',alpha=0.5)+
  geom_hline(yintercept = 0,color='red',alpha=0.5)+
  geom_abline(aes(slope=1,intercept=0),color='grey',linetype='longdash',size=0.5)+
  geom_errorbar(aes(ymin=`ndk/BW_logFC`-`ndk/BW_SE`,ymax=`ndk/BW_logFC`+`ndk/BW_SE`),alpha=erralpha,color=errcolor,width=0)+
  geom_errorbarh(aes(xmin=`pyrE/BW_logFC`-`pyrE/BW_SE`,xmax=`pyrE/BW_logFC`+`pyrE/BW_SE`),alpha=erralpha,color=errcolor,height=0)+
  geom_point(aes(colour=`ndk/pyrE_logFC`),alpha=0.9,size=3)+
  scale_x_continuous(breaks=seq(-3,6,by=1),limits=c(-3,6))+
  scale_y_continuous(breaks=seq(-3,6,by=1),limits=c(-3,6))+
  ylab('ndk/BW logFC')+
  xlab('pyrE/BW logFC')+
  geom_text(aes(label=ifelse(`ndk/pyrE_p.value`<0.05 & `ndk/BW_p.value`<0.05,as.character(Compound),'')),
            hjust=1, vjust=-0.5,size=3,colour = "red")+
  scale_colour_gradientn(colours = gradcols,
                         breaks=brks,limits=c(-5,5),name=maincomp)
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/ndk_pyrE_comparison.pdf",sep = ''),
             width=6,height=4)



multcompcastw<-multcompcast[, -grep("_t.value", colnames(multcompcast))]
multcompcastw<-multcompcastw[, -grep("_FoldChange", colnames(multcompcastw))]

write.csv(multcompcastw,paste(ddir,'/Nucleoside_metabolomics_ndk_Comparisons.csv',sep=''),row.names = FALSE)




selresults<-subset(fullinfo,Comparison %in% c('ndk/BW','pyrE/BW','yjjG/BW','ndk/pyrE'))
selresults$FDR<-p.adjust(selresults$p.value,method = 'fdr')




selresultsw<-selresults




colnames(selresultsw)
selresultsw<-selresultsw[,c('Type','Group','Class','Comparison','Compound','FoldChange','logFC','SE','t.value','p.value','FDR')]
colnames(selresultsw)



Metsexpl<-read.table(paste(ddir,'/Metabolomics_nucleotides_statistics_explanation.csv',sep=''),
                     sep=',',quote = '"',header = TRUE,stringsAsFactors=FALSE)

explst<-subset(Metsexpl,Column %in% colnames(selresultsw))

explst<-explst[match(colnames(selresultsw),explst$Column),]


write.xlsx2(explst, file=paste(ddir,'/Metabolomics_nucleotides_Statistics.xlsx',sep=''), sheetName="Readme",row.names = FALSE,showNA=FALSE)
write.xlsx2(selresultsw, file=paste(ddir,'/Metabolomics_nucleotides_Statistics.xlsx',sep=''), sheetName="Data", append=TRUE,row.names = FALSE)#showNA=FALSE


# #Update manually
# write.xlsx2(explst, file='/Users/Povilas/Projects/B-D-H paper/Without Chemical screen/tables/Scott et al_Table S6.xlsx',
#             sheetName="Nucleoside_metabolomics_Statistics_Readme",row.names = FALSE,showNA=FALSE,append=TRUE)
# write.xlsx2(selresultsw, file='/Users/Povilas/Projects/B-D-H paper/Without Chemical screen/tables/Scott et al_Table S6.xlsx',
#             sheetName="Nucleodie_metabolomics_Statistics", append=TRUE,row.names = FALSE)#showNA=FALSE
# # 



#Volcano plot

txtsize<-5
nudgex<-0.1
nudgey<-0.3
erralpha<-1
errcolor<-'grey80'
baralpha<-0.2
barcolor<-'red'


straincols<-c('red','darkgreen','#89CFF0')



Nucp<-ggplot(subset(selresults,Comparison!='ndk/pyrE'),aes(x=logFC,y=-log(p.value,10),color=Comparison))+
  geom_hline(yintercept = -log(0.05,10),color=barcolor,alpha=baralpha)+
  geom_vline(xintercept = -1,color=barcolor,alpha=baralpha)+
  geom_vline(xintercept = 1,color=barcolor,alpha=baralpha)+
  labs(color='Strain')+
  xlab('logFC')+
  xlim(-3,6)+
  scale_color_manual(values=straincols)+
  scale_y_continuous(breaks=seq(0,10,by=2),limits=c(0,10.5))+
  theme(panel.grid.minor = element_blank())
Nucp+geom_errorbarh(aes(xmin=logFC-SE,xmax=logFC+SE),color=errcolor)+
  geom_point(size=3)
dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Volcano_Significant_p.value_clean.pdf',sep = ''),
             width=9,height=5,useDingbats=FALSE)
Nucp+geom_text(aes(label=ifelse(p.value<0.05 & abs(logFC)>1,as.character(Compound),'')),
                 nudge_x=nudgex,nudge_y=nudgey,size=txtsize)
dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Volcano_Significant_p.value_clean_text.pdf',sep = ''),
             width=9,height=5,useDingbats=FALSE)



Nucp+geom_errorbarh(aes(xmin=logFC-SE,xmax=logFC+SE),color=errcolor)+
  geom_point(size=3)+
  geom_text(aes(label=ifelse(p.value<0.05 & abs(logFC)>1,as.character(Compound),'')),
               nudge_x=nudgex,nudge_y=nudgey,size=txtsize)
dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Volcano_Significant_p.value_Full.pdf',sep = ''),
             width=9,height=5,useDingbats=FALSE)




#Heatmaps raw data prep
metsstat<-ddply(ecoliM,Compound+Strain~Compound,summarise,Mean=mean(LogConc),SD=sd(LogConc))

metsfullstat<-merge(ecoliraw[,1:4],metsstat,by.x = 'Compound',by.y='Compound')

#allstat<-dcast(metsfullstat,Class+Type+Compound~Strain,mean,value.var = c('Mean'),fill = as.numeric(NA),drop=TRUE)
allstat<-dcast(ecoliM,Class+Type+Compound~Strain+Replicate,mean,value.var = c('LogConc'),fill = as.numeric(NA),drop=TRUE)



allstats<-allstat[order(allstat$Class,allstat$Type),]
Names<-allstats[,'Compound']

#Draw heatmap


nstep<-14

bgg <- colorRampPalette(c("blue", "gray90", "red"))(n = nstep)

colscale<-bgg#c(cs,bgg)

brkst<-seq(-3.5,3.5,by=7/nstep)
brks<-c(-6,brkst)


selcol<-c('BW_1','BW_2','BW_3','BW_4',
          'BW_5','BW_6','BW_7','BW_8',
          'ndk_1','ndk_2','ndk_3','ndk_4',
          'pyrE_1','pyrE_2','pyrE_4',
          'yjjG_1','yjjG_2','yjjG_3','yjjG_4')


cleannames<-selcol

statsh<-allstats[,selcol]

reorderfun_mean = function(d,w) { reorder(d, w, agglo.FUN = mean) }
reorderfun_median = function(d,w) { reorder(d, w, agglo.FUN = median) }

statsfill<-statsh
statsfill[is.na(statsfill)]<-0


statshf<-statsh[,selcol]
statsfill<-statshf
statsfill[is.na(statsfill)]<-0



hmap<-heatmap.2(as.matrix(statshf),key=TRUE,Colv=TRUE,trace='none',labRow=Names,
                labCol=cleannames,col=bgg,
                xlab='Comparison',Rowv=TRUE,
                dendrogram="both",scale="none",na.color="white",
                cexRow=1,cexCol=1,margin=c(8,16),
                lwid=c(0.2,0.8),symkey=TRUE)

dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Metabolomics_Heatmap.pdf',sep = ''),
             width=12,height=15)






