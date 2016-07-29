library('ggplot2')
library('gplots')
library('plyr')
library('reshape2')
library(tidyr)
library(BSDA)
library(xlsx)
library(gridExtra)
theme_set(theme_light())

odir<-'Figures_final/Metabolomics'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")
ddir<-'Data_final'

bac<-read.table('Metabolomics/Bacteria_log_normalised.csv',sep=',',stringsAsFactors = FALSE)
worm<-read.table('Metabolomics/Worms_log_normalised.csv',sep=',',stringsAsFactors = FALSE)

#Enrichment
bacsCT<-read.table('Metabolomics/Bacteria_C-A_QEA.csv',sep=',',stringsAsFactors = FALSE,quote = '"',header = TRUE)
wormsCT<-read.table('Metabolomics/Worms_W-Z_QEA.csv',sep=',',stringsAsFactors = FALSE,quote = '"',header = TRUE)

bacsCT$Type<-'5FU/Control'
wormsCT$Type<-'5FU/Control'
bacsCT$Organism<-'E.coli'
wormsCT$Organism<-'C.elegans'

bacsA<-read.table('Metabolomics/Bacteria_B-A_QEA.csv',sep=',',stringsAsFactors = FALSE,quote = '"',header = TRUE)
wormsA<-read.table('Metabolomics/Worms_X-Z_QEA.csv',sep=',',stringsAsFactors = FALSE,quote = '"',header = TRUE)

bacsA$Type<-'Arabinose/Control'
wormsA$Type<-'Arabinose/Control'
bacsA$Organism<-'E.coli'
wormsA$Organism<-'C.elegans'


bacsACT<-read.table('Metabolomics/Bacteria_D-B_QEA.csv',sep=',',stringsAsFactors = FALSE,quote = '"',header = TRUE)
wormsACT<-read.table('Metabolomics/Worms_Y-X_QEA.csv',sep=',',stringsAsFactors = FALSE,quote = '"',header = TRUE)

bacsACT$Type<-'Arabinose+5FU/Control'
wormsACT$Type<-'Arabinose+5FU/Control'
bacsACT$Organism<-'E.coli'
wormsACT$Organism<-'C.elegans'



statsCT<-merge(bacsCT,wormsCT,all.x = TRUE,all.y=TRUE)
statsA<-merge(bacsA,wormsA,all.x = TRUE,all.y=TRUE)
statsACT<-merge(bacsACT,wormsACT,all.x = TRUE,all.y=TRUE)

stats1<-merge(statsCT,statsACT,all.x = TRUE,all.y=TRUE)
stats<-merge(stats1,statsA,all.x = TRUE,all.y=TRUE)

stats$FoldEnrichment<-stats$Statistic/stats$Expected
stats$Organism<-as.factor(stats$Organism)
stats$Significant<-ifelse(stats$P.value<0.05,TRUE,FALSE)
stats$Significant<-factor(stats$Significant,levels = c(TRUE,FALSE),labels=c('True','False'))

write.csv(stats,paste(ddir,'/Metabolomics_All_QEA.csv',sep=''),row.names = FALSE)


statsm<-melt(stats,id.vars = c('Metabolite.Set','Total','Organism','Type','Metabolites..Total.'),
            value.name = 'Value',variable.name = 'Measure')
statsm$Value<-as.numeric(statsm$Value)

statsall<-dcast(statsm,Metabolite.Set+Total+Type~Measure+Organism,value.var = 'Value')



ggplot(subset(stats,P.value<0.05),
       aes(y=FoldEnrichment,x=reorder(Metabolite.Set,FoldEnrichment)))+#,fill=-log10(P.value)
  geom_bar(stat = "identity",fill='red')+
#   scale_fill_gradientn(colours = c('black','orange','red','red'),
#                          breaks=c(1,2,3,4,5,6),limits=c(1,6),guide='legend')+
  coord_flip()+
  xlab('Metabolite set')+
  ylab('Fold enrichment')+
  facet_grid(Type~Organism)
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Metabolomics_QEA_significant.pdf",sep=''),
             width=12,height=18)

cstats<-subset(stats,!Metabolite.Set %in% c("Intracellular Signalling Through Adenosine Receptor A2A And Adenosine | Intracellular Signalling Through Adenosine Receptor A2B And Adenosine",
                                             "Intracellular Signalling Through Fsh Receptor And Follicle Stimulating Hormone | Intracellular Signalling Through Lhcgr Receptor And Luteinizing Hormone/Choriogonadotropin",
                                            "Excitatory Neural Signalling Through 5-Htr 4 And Serotonin | Excitatory Neural Signalling Through 5-Htr 7 And Serotonin | Excitatory Neural Signalling Through 5-Htr 6 And Serotonin"))
ggplot(cstats,
       aes(y=FoldEnrichment,x=reorder(Metabolite.Set,FoldEnrichment),alpha=Significant))+#,fill=-log10(P.value)
  geom_bar(stat = "identity",fill='red')+
  scale_alpha_discrete(range=c(1,0.2))+
  coord_flip()+
  xlab('Metabolite set')+
  ylab('Fold enrichment')+
  facet_grid(Type~Organism)
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Metabolomics_QEA.pdf",sep=''),
             width=12,height=20)



ggplot(subset(stats,Organism=='C.elegans' & Type=='5FU/Control' & P.value<0.05),
       aes(y=FoldEnrichment,x=reorder(Metabolite.Set,FoldEnrichment),fill=-log10(P.value)))+#,fill=-log10(P.value)
  geom_bar(stat = "identity")+
  scale_fill_gradientn(colours = c('red','orange','yellow'),
                       breaks=c(1.3,2,3,4,5,6,7),
                       limits=c(1.3,7))+
  coord_flip()+
  labs(fill='-log10(p-value)')+
  ggtitle('C. elegans 5FU/Control')+
  xlab('Metabolite set')+
  ylab('Fold enrichment')
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Metabolomics_QEA_5FUvsControl_Celegans_significant.pdf",sep=''),
             width=6,height=7)




tpcol<-c('5FU/Control','Arabinose/Control','Arabinose+5FU/Control')
tpabr<-c('5FUvsC','L-ArabvsC','L-Arab5FUvsC')
types<-data.frame(column=tpcol,abreviation=tpabr)

for(t in 1:nrow(types)){
  tp<-types[t,]
  sel<-subset(cstats,Type==tp$column )
  fname=paste(odir,"/Metabolomics_QEA_",tp$abreviation,"_All.pdf",sep='')
  fac<-ggplot(sel,
         aes(y=FoldEnrichment,x=reorder(Metabolite.Set,FoldEnrichment),alpha=Significant))+
    geom_bar(stat = "identity",fill='red')+
    scale_alpha_discrete(range=c(1,0.2))+
#     scale_fill_gradientn(colours = c('orange','red'),
#                          breaks=c(1.3,2),limits=c(1.3,7))+
    coord_flip()+
    ggtitle(tp$column)+
    xlab('Metabolite set')+
    ylab('Fold enrichment')+
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
  sel<-subset(cstats,Type==tp$column & P.value<0.05)
  fname=paste(odir,"/Metabolomics_QEA_",tp$abreviation,"_significant.pdf",sep='')
  fac<-ggplot(sel,
              aes(y=FoldEnrichment,x=reorder(Metabolite.Set,FoldEnrichment)))+
    geom_bar(stat = "identity",fill='red')+
    #     scale_fill_gradientn(colours = c('orange','red'),
    #                          breaks=c(1.3,2),limits=c(1.3,7))+
    coord_flip()+
    ggtitle(tp$column)+
    xlab('Metabolite set')+
    ylab('Fold enrichment')+
    facet_grid(.~Organism)
  fac
  
  print(fname)
  cairo_pdf(fname,width=12,height=9)
  print(fac)
  #grid.arrange(gt1, gt2, ncol=2,widths=c(3,1))
  #multiplot(termMICdist, termMICcov,  cols=2)
  dev.off()
  
}



colnames(bac)<-bac[1,]
colnames(worm)<-worm[1,]

bac<-bac[-c(1),]
worm<-worm[-c(1),]

bacm<-melt(bac,id.vars = c('ID','Group','Replicate','5-FU','L-Arabinose'),variable.name = 'Metabolite',value.name = 'logConc')
wormm<-melt(worm,id.vars = c('ID','Group','Replicate','5-FU','L-Arabinose'),variable.name = 'Metabolite',value.name = 'logConc')
bacm$Organism<-'E.coli'
wormm$Organism<-'C.elegans'

mets<-merge(bacm,wormm,all.x = TRUE,all.y=TRUE)
mets$Organism<-factor(mets$Organism,levels=c('E.coli','C.elegans'))
mets$Metabolite<-as.factor(mets$Metabolite)
mets$`5-FU`<-as.factor(mets$`5-FU`)
mets$`L-Arabinose`<-as.factor(mets$`L-Arabinose`)
mets$logConc<-as.numeric(mets$logConc)


write.csv(mets,paste(ddir,'/Metabolomics_All_Data_logNormalised.csv',sep=''),row.names = FALSE)


metssum<-ddply(mets, .(Group,`5-FU`,`L-Arabinose`,Organism,Metabolite), summarise,
               Mean=mean(logConc,na.rm = TRUE),
               SD=sd(logConc,na.rm = TRUE))

metsumm<-melt(metssum,id.vars = c('Group','5-FU','L-Arabinose','Organism','Metabolite'),
              variable.name = 'Stats',value.name = 'Value')

metO<-dcast(metsumm,`5-FU`+`L-Arabinose`+Metabolite~Organism+Stats,value.var = 'Value')


metssumA<-ddply(mets, .(Group,`5-FU`,Organism,Metabolite), summarise,
               Mean=mean(logConc,na.rm = TRUE),
               SD=sd(logConc,na.rm = TRUE))


metsdm<-melt(mets,measure.vars = c('5-FU','L-Arabinose','Organism'))


ggplot(metO,aes(x=E.coli_Mean,y=C.elegans_Mean,color=`5-FU`,shape=`L-Arabinose`))+geom_point()

#Linear modeling


#mets <- within(mets, Organism <- relevel(Organism, ref = 1))
#contrasts(mets$Metabolite)<-contr.sum(length(unique(mets$Metabolite)))

#:Metabolite

selmets<-subset(mets,Organism=='E.coli' & `L-Arabinose`=='0mM')
model<-lm(logConc ~`5-FU`*Metabolite,data=selmets)#+`L-Arabinose`+`L-Arabinose`+Organism+
summary(model)

modelc<-data.frame(smodel$coefficients)

plot(model)

ggplot(selmets,aes(x=logConc))+geom_histogram()

ggplot(selmets,aes(x=logConc,fill=`5-FU`))+geom_histogram(position='identity',alpha=0.5)

