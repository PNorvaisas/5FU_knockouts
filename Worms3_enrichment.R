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






odir<-'Figures_final'
ddir<-'Data_final'

bacmic<-read.table(paste(ddir,'/MICs_and_bacterial_growth-Complete.csv',sep=''),sep=',',header=TRUE)
bacmic<-bacmic[,!colnames(bacmic) %in% c('X','JW_id','ECK','X0','X1','X2.5','X5','bno','EG','GI','MOPS_24hr','MOPS_48hr')]



#Read Annotations
allannot<-read.table(paste(ddir,'/Enrichement_all_terms.csv',sep=''),sep=',',header=TRUE)

circmr<-merge(bacmic,allannot, by.x='Gene',by.y='Gene',all.x = TRUE)
circmr$Term<-as.factor(circmr$Term)
circmr$Category<-as.factor(circmr$Category)

#Get gene counts for screened Keio library
circsummic<-ddply(circmr,.(Category,Term_ID),summarise, Size_KeioS=as.numeric(length(unique(Gene))) )


#circsumfull<-merge(circsummic,circsumbac,by=c('Category','Term_ID'),all.x=TRUE,all.y=TRUE)



circm<-merge(circmr,circsummic,by=c('Category','Term_ID'),all.x=TRUE)

circm<-subset(circm,Category %in% c('KEGG','EcoCyc') & Term!='' & ! is.na(Term))

unique(allannot$Category)
unique(circm$Category)

unique(circm$Term)


Slevel=c(0)
SAbr=c('all')
Stat<-data.frame(slevel=Slevel,sabr=SAbr)

# Old lists
# TName=c('BP','KP','EcoCyc')
# TTitle=c('Biological process', 'KEGG pathway','EcoCyc pathway')
# TAbr=c('BP','KEGG-PWY','EcoCyc-PWY')
# Types<-data.frame(tname=TName,tabr=TAbr,ttitle=TTitle)


TName=c('KEGG','EcoCyc')#,'EcoCyc','GO_BP'
TTitle=c('KEGG pathway','EcoCyc pathway')#,'EcoCyc pathway','GO Biological process'
TAbr=c('KEGG-PWY','EcoCyc-PWY')#,'EcoCyc-PWY','GO_BP'
Types<-data.frame(tname=TName,tabr=TAbr,ttitle=TTitle)

Thrnames<-c('Coverage in E. coli by MIC','Coverage in valid Keio by MIC','Coverage by KO Interacting','Coverage by KO Synergistic','Coverage by KO Antagonistic',
            'Coverage by KO Interacting (FDR)','Coverage by KO Synergistic (FDR)','Coverage by KO Antagonistic (FDR)')
Thrcols<-c('Coverage_EC','Coverage_KeioS','Coverage_BacInteracting','Coverage_BacSynergistic','Coverage_BacAntagonistic',
           'Coverage_BacInteracting_FDR','Coverage_BacSynergistic_FDR','Coverage_BacAntagonistic_FDR')
Thraxis<-c('Coverage','Coverage','Coverage','Coverage','Coverage','Coverage','Coverage','Coverage')
Thrabr<-c('CovEC','CovKS','CovBacInt','CovBacSyn','CovBacAnt',
          'CovBacInt_FDR','CovBacSyn_FDR','CovBacAnt_FDR')
Thresholds<-data.frame(names=Thrnames,cols=Thrcols,axis=Thraxis,abr=Thrabr)


Media=c('LB_22hr','T_OD_Mean','C_OD_Mean')
MTitle=c('LB 22hr growth','NGM + 100uM 5FU growth','NGM growth')
MAbr=c('LB','NGM-5FU','NGM')
MOrder=c('LB_S_Med','T_OD_Mean_S_Med','C_OD_Mean_S_Med')
MQvars=c('LB','T_OD_Mean','C_OD_Mean')
MScale=c(1,0.3,0.3)
Medias<-data.frame(media=Media,mtitle=MTitle,mabr=MAbr,morder=MOrder,mscale=MScale,mqvars=MQvars)



#Enricchment threshold should be set to 5
micthres<-5


#Clean for MIC enrichment
circmc<-circm[, -grep("_SD", colnames(circm))]
circmc<-circmc[, -grep("_p.value", colnames(circmc))]
circmc<-circmc[, -grep("_t.value", colnames(circmc))]
circmc<-circmc[, -grep("_SE", colnames(circmc))]
circmc<-circmc[, -grep("_FDR", colnames(circmc))]
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

circmst<-dcast(circmstatm,Category+Term_ID+Term+Size_EC+Size_KeioS+Size_MICo5~ Measure+Stat,value.var = 'Value')


#
# circAllm<-melt(subset(circm,MIC>micthres & GT_p.value<0.05),id=c('Category','Term_ID','Term','Gene','Plate','Well','Size_EC','Size_KeioS'),
#                variable.name = 'Measure',value.name='Value')



# circmBAll<-subset(circm,MIC>micthres & GT_p.value<0.05)
# circmBAll$Interaction<-'BacInteracting'
# circmBSyn<-subset(circm,MIC>micthres & GT_Interaction < 0 & GT_p.value<0.05)
# circmBSyn$Interaction<-'BacSynergistic'
# circmBAnt<-subset(circm,MIC>micthres & GT_Interaction > 0 & GT_p.value<0.05)
# circmBAnt$Interaction<-'BacAntagonistic'

circmBAll<-subset(circm,MIC>micthres & GT_FDR < 0.05)
circmBAll$Interaction<-'BacInteracting'
circmBSyn<-subset(circm,MIC>micthres & GT_Interaction < 0 & GT_FDR<0.05)
circmBSyn$Interaction<-'BacSynergistic'
circmBAnt<-subset(circm,MIC>micthres & GT_Interaction > 0 & GT_FDR<0.05)
circmBAnt$Interaction<-'BacAntagonistic'

circmBSA_p<-merge(circmBSyn,circmBAnt,all=TRUE)
circmB<-merge(circmBSA_p,circmBAll,all=TRUE)



circmBStat<-ddply(circmB,.(Category,Term_ID,Interaction),summarise,
                    Size=as.numeric(length(unique(Gene))))


circmBStatM<-melt(circmBStat,measure.vars = c('Size'),
                 variable.name = 'Stat',value.name='Value')
circmBs<-dcast(circmBStatM,Category+Term_ID ~ Stat+Interaction,value.var = 'Value')

circms<-merge(circmst,circmBs,by=c('Category','Term_ID'),all.x=TRUE)


gsel<-c('Phenylalanine, tyrosine and tryptophan biosynthesis',
        'Pyrimidine metabolism',
        'Purine ribonucleosides degradation',
        'One carbon pool by folate',
        'Salvage pathways of pyrimidine ribonucleotides',
        'NADH to cytochrome bd oxidase electron transport II',
        'Succinate to cytochrome bd oxidase electron transfer',
        'Succinate to cytochrome bo oxidase electron transfer')

#Minimum number of hits
gtres<-4

circms<-subset(circms,! Term %in% c('Metabolic pathways',
                                    'Microbial metabolism in diverse environments',
                                    NA,
                                    'Biosynthesis of secondary metabolites',
                                    '') &
                 (Size_MICo5>gtres | Term %in% gsel))


bacmicnoWT<-subset(bacmic,Gene!='WT')
tot<-length(bacmicnoWT$Gene)
hits<-length(subset(bacmicnoWT,MIC>5)$Gene)
tot
hits

#Hits by FDR
bacInt<-length(subset(bacmic,MIC>micthres & GT_FDR<0.05)$Gene)
bacSyn<-length(subset(bacmic,MIC>micthres & GT_Interaction < 0 &GT_FDR<0.05)$Gene)
bacAnt<-length(subset(bacmic,MIC>micthres & GT_Interaction > 0 &GT_FDR<0.05)$Gene)



circms$Coverage_EC<-circms$Size_MICo5/circms$Size_EC
circms$Coverage_KeioS<-circms$Size_MICo5/circms$Size_KeioS


circms$Coverage_BacInteracting<-circms$Size_BacInteracting/circms$Size_MICo5
circms$Coverage_BacSynergistic<-circms$Size_BacSynergistic/circms$Size_MICo5
circms$Coverage_BacAntagonistic<-circms$Size_BacAntagonistic/circms$Size_MICo5





circms$Enrichment_MIC_p.value<-phyper(circms$Size_MICo5-1,hits,tot-hits,circms$Size_KeioS,lower.tail = FALSE)



circms$Enrichment_MIC_FDR<-p.adjust(circms$Enrichment_MIC_p.value,method = 'fdr')
  
circms$Enrichment_BacInt_p.value<-phyper(circms$Size_BacInteracting-1,bacInt,hits-bacInt,circms$Size_MICo5,lower.tail = FALSE)
circms$Enrichment_BacSyn_p.value<-phyper(circms$Size_BacSynergistic-1,bacSyn,hits-bacSyn,circms$Size_MICo5,lower.tail = FALSE)
circms$Enrichment_BacAnt_p.value<-phyper(circms$Size_BacAntagonistic-1,bacAnt,hits-bacAnt,circms$Size_MICo5,lower.tail = FALSE)



subset(circms,Enrichment_BacInt_p.value<0.05)[,c('Size_EC','Size_KeioS','Size_MICo5','Size_BacInteracting','Enrichment_BacInt_p.value','Category','Term_ID','Term')]
subset(circms,Enrichment_BacSyn_p.value<0.05)[,c('Size_EC','Size_KeioS','Size_MICo5','Size_BacSynergistic','Enrichment_BacSyn_p.value','Category','Term_ID','Term')]
subset(circms,Enrichment_BacAnt_p.value<0.05)[,c('Size_EC','Size_KeioS','Size_MICo5','Size_BacAntagonistic','Enrichment_BacAnt_p.value','Category','Term_ID','Term')]

subset(circms,Enrichment_MIC_FDR<0.05 )[,c('Size_EC','Size_KeioS','Size_MICo5','Category','Term_ID','Term','Enrichment_MIC_FDR')]
subset(circms,Enrichment_MIC_p.value<0.05 )[,c('Size_EC','Size_KeioS','Size_MICo5','Category','Term_ID','Term','Enrichment_MIC_FDR','Enrichment_MIC_p.value')]



#PLP test
phyper(9,hits,tot-hits,43,lower.tail = FALSE)




#subset(enr,Enrichment_BacSyn_FDR<0.05)
subset(circms,Enrichment_BacSyn_p.value<0.05)
subset(circms,Enrichment_BacAnt_p.value<0.05)




#Redo for saving

circms_re<-subset(circms,!is.na(Category) & Size_KeioS>3)
circms_re<-circms_re[, -grep("logOD", colnames(circms_re))]
circms_re<-circms_re[, -grep("WTDiff", colnames(circms_re))]
circms_re<-circms_re[, -grep("OD_Mean", colnames(circms_re))]
circms_re<-circms_re[, -grep("LB_22hr", colnames(circms_re))]
circms_re<-circms_re[, -grep("GT_Interaction", colnames(circms_re))]
circms_re<-circms_re[, -grep("BacInt", colnames(circms_re))]

colnames(circms_re)

# 
# circms_rew<-circms_re[,c(1:6,
#                         match(c('Size_BacAntagonistic'),colnames(circms_re)):match(c('Enrichment_BacAnt_p.value'),colnames(circms_re)),
#                         7: (match(c('Size_BacAntagonistic'),colnames(circms_re))-1) )]


circms_rew<-circms_re[,c(1:6,16,17,20,21)]



colnames(circms_rew)

annotexpl<-read.table(paste(ddir,'/Enrichment_distribution_column_explanation.csv',sep=''),
                    sep=',',quote = '"',header = TRUE,stringsAsFactors=FALSE)

MICexpl<-read.table(paste(ddir,'/MICs_column_explanation.csv',sep=''),
                      sep=',',quote = '"',header = TRUE,stringsAsFactors=FALSE)



explrm<-subset(annotexpl,Column %in% colnames(circms_rew))

explrm<-explrm[match(colnames(circms_rew),explrm$Column),]

write.csv(circms_rew,
          paste(ddir,'/Enrichment_Distributions_for_terms_MICo',micthres,'.csv',sep=''),
          row.names = FALSE)

write.xlsx2(explrm, file=paste(ddir,'/Enrichment_Distributions_for_terms_MICo',micthres,'.xlsx',sep=''),
           sheetName="Readme",row.names = FALSE,showNA=FALSE)
write.xlsx2(circms_rew, file=paste(ddir,'/Enrichment_Distributions_for_terms_MICo',micthres,'.xlsx',sep=''),
           sheetName="Data", append=TRUE,row.names = FALSE,showNA=FALSE)


write.xlsx2(explrm, file='/Users/Povilas/Projects/B-D-H paper/Without Chemical screen/tables/Scott et al_Table S2.xlsx',
            sheetName="Readme_Enrichment",row.names = FALSE,showNA=FALSE,append = TRUE)
write.xlsx2(circms_rew, file='/Users/Povilas/Projects/B-D-H paper/Without Chemical screen/tables/Scott et al_Table S2.xlsx',
            sheetName="Enrichment", append=TRUE,row.names = FALSE,showNA=FALSE)


colnames(circm)
colnames(circms_rew)


#Circms instead of circm, filtered for meaningful terms
#Raw data not saved


enr<-merge(circm[,c('Category','Term','Term_ID','Gene','Plate','Well',
                    'MIC','MIC_SD',
                    'GT_Interaction','GT_p.value','GT_FDR','Size_EC','Size_KeioS')],
           circms_rew[,!colnames(circms_rew) %in% c('Size_EC','Size_KeioS')],
           by=c('Category','Term_ID','Term'),all.x=TRUE)

enr[enr=='NaN']<-NA
dim(enr)
head(enr)

annotexpl2<-merge(annotexpl,MICexpl, all=TRUE)


explrm<-subset(annotexpl2,Column %in% colnames(enr))

explrm<-explrm[match(colnames(enr),explrm$Column),]


write.csv(enr,paste(ddir,'/Enrichment_Complete_data_term_distribution_MICo',micthres,'.csv',sep=''))


write.xlsx2(explrm, file=paste(ddir,'/Enrichment_Complete_data',micthres,'.xlsx',sep=''),
            sheetName="Readme",row.names = FALSE,showNA=FALSE)
write.xlsx2(enr, file=paste(ddir,'/Enrichment_Complete_data',micthres,'.xlsx',sep=''),
            sheetName="Data", append=TRUE,row.names = FALSE,showNA=FALSE)


write.xlsx2(explrm, file='/Users/Povilas/Projects/B-D-H paper/Without Chemical screen/tables/Scott et al_Table S2.xlsx',
            sheetName="Readme_Enrichment_with_genes",row.names = FALSE,showNA=FALSE,append = TRUE)
write.xlsx2(enr, file='/Users/Povilas/Projects/B-D-H paper/Without Chemical screen/tables/Scott et al_Table S2.xlsx',
            sheetName="Enrichment_with_genes", append=TRUE,row.names = FALSE,showNA=FALSE)







#Create necessary directories
envar<-paste('Enrichment_MICo',micthres,sep='')
dir.create(paste(odir,'/',envar,sep=''), showWarnings = TRUE, recursive = FALSE, mode = "0777")

micname<-expression(paste('C. elegans MIC\n median [5FU], ',mu,'M',sep=''))
mic<-expression(paste('C. elegans MIC [5FU], ',mu,'M',sep=''))



#,'Pyridoxal 5-phosphate salvage I'
#'Succinate to cytochrome bo oxidase electron transfer',

sel<-subset(enr,(Enrichment_MIC_FDR < 0.05 | Term %in% gsel))

fname<-paste(odir,'/',envar,'/MICsignificance/EcoKEGG_MIC_all_MICo5.pdf',sep = '')
fname2<-paste(odir,'/',envar,'/MICsignificance/EcoKEGG_MIC_all_MICo5_large.pdf',sep = '')
print(fname)
sel$Term<-gsub('metabolism','metab.',sel$Term)
sel$Term<-gsub('biosynthesis','bs.',sel$Term)


customorder<-c('Pyridoxal 5-phosphate bs. I',
               'Vitamin B6 metab.',
               'N10-formyl-tetrahydrofolate bs.',
               'One carbon pool by folate',
               'Chorismate bs. from 3-dehydroquinate',
               'Phenylalanine, tyrosine and tryptophan bs.',
               'TCA cycle I (prokaryotic)',
               'Citrate cycle (TCA cycle)',
               'Oxidative phosphorylation',
               'Succinate to cytochrome bd oxidase electron transfer',
               'Succinate to cytochrome bo oxidase electron transfer',
               'NADH to cytochrome bd oxidase electron transport II',
               'Flagellar assembly',
               'Salvage pathways of pyrimidine ribonucleotides',
               'Pyrimidine metab.',
               'Purine ribonucleosides degradation',
               'Purine metab.')

length(unique(customorder))
length(unique(sel$Term))

sel$Term<-factor(sel$Term, levels=rev(customorder))

#sel<-sel[order(sel$Term),]

sizecol<-c('black','red','yellow')

termMIC<-ggplot(sel,
                aes(y=MIC,x=Term,fill=Coverage_KeioS,color=Size_MICo5))+
  geom_hline(yintercept=5,color='red',alpha=0.5)+
  geom_violin(position='identity',scale = "width")+
  geom_jitter(height = 0)+
  coord_flip()+
  scale_colour_gradientn(colours = sizecol,breaks=c(0,10,20),limits=c(0,20))+#breaks=c(5,10,15,20),limits=c(2,20) #,guide='legend',guide='legend'
  labs(fill='Coverage',
       color='Number of KOs with MIC>5')+
  scale_fill_gradientn(colours=c('grey70',"blue",'cyan'),limits=c(0,1),breaks=c(0,0.2,0.4,0.6,0.8,1))+
  geom_text(aes(y=0,label=as.character(stars.pval(Enrichment_MIC_FDR))),
            color='black',hjust=0, vjust=0,size=3)+
  geom_text(aes(y=100,label=ifelse(Category=='KEGG','K','E')),
            color='black',hjust=0, vjust=0,size=3)+
  ylab(mic)+
  scale_y_continuous(breaks=c(1,2.5,5,10,25,50,100),limits=c(1,100),trans='log')+
  xlab('EcoKEGG')+
  theme(legend.box='horizontal',
        legend.position='top',
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())
termMIC


dev.copy2pdf(device=cairo_pdf,
             file=fname,
             width=6,height=8,
             useDingbats=FALSE)
dev.copy2pdf(device=cairo_pdf,
             file=fname2,
             width=10,height=8,
             useDingbats=FALSE)

# 
# 
# 
# ##MIC distributions and coverage
# print('MIC distributions and coverage')
# dir.create(paste(odir,'/',envar,'/MIC_coverage_and_distribution',sep=''), showWarnings = TRUE, recursive = FALSE, mode = "0777")
# for(s in 1:nrow(Stat)) {
#   sr <- Stat[s,]
#   for(t in 1:nrow(Types)) {
#     tr <- Types[t,]
#     for(m in 1:nrow(Thresholds)) {
#       mr<-Thresholds[m,]
#       sel<-subset(enr,Category==as.character(tr$tname) & MIC>mtres & Size_MICo5>gtres&
#                     MIC_S_Med>mttres &!is.na(Term) &
#                     ! Term %in% c('Metabolic pathways',
#                                   'Biosynthesis of secondary metabolites',
#                                   'Microbial metabolism in diverse environments'))
#       selcov<-subset(circms,Category==as.character(tr$tname) & Size_MICo5>gtres&
#                     MIC_S_Med>mttres &!is.na(Term) &
#                     ! Term %in% c('Metabolic pathways',
#                                   'Biosynthesis of secondary metabolites',
#                                   'Microbial metabolism in diverse environments'))
#       
#       gtitle<-paste('Enriched E. Coli MG1655 ',tr$ttitle,' - MIC distribution, ',sr$sabr,' MIC>',micthres,sep = '')
#       fname<-paste(odir,'/',envar,'/MIC_coverage_and_distribution/',tr$tabr,'_MIC_',sr$sabr,'_MICo',micthres,'_',mr$abr,'.pdf',sep = '')
#       #print(gtitle)
#       
#       termMICdist<-ggplot(sel,
#                           aes(y=MIC,x=reorder(sel$Term,sel[,as.character(mr$cols)]),fill=Size_MICo5,color=MIC_S_Med))+
#         geom_boxplot(position='identity')+
#         coord_flip()+
#         labs(fill='Number of genes',color='MIC median')+
#         scale_fill_gradient( high="red",low='white')+
#         #scale_color_gradient(limits=c(0,100), high="blue",low='gray')+
#         scale_colour_gradientn(colours = c('black','orange','red'),
#                                breaks=c(5,25,50),limits=c(5,50),guide='legend',name=micname)+
#         ylab(mic)+ylim(0,100)+
#         xlab(tr$ttitle)+ggtitle(gtitle)+theme(legend.position="left")
#       gtitle2<-paste(mr$names,' MIC>',micthres,sep = '')
#       termMICcov<-ggplot(selcov,
#                          aes(y=selcov[,as.character(mr$cols)],x=reorder(selcov$Term,selcov[,as.character(mr$cols)])))+
#         geom_bar(stat = "identity",fill='red',alpha=0.8)+
#         coord_flip()+
#         labs(fill='Number of hits')+
#         scale_fill_gradient( high="red",low='gray')+
#         ylab(mr$axis)+
#         xlab(tr$ttitle)+ggtitle(gtitle2)+theme(axis.text.y=element_blank(),
#                                                axis.title.y=element_blank(),
#                                                axis.ticks.y=element_blank(),
#                                                legend.position="none")
#       #termMICcov
#       #fill=Size_MICo5
#       gt1 <- ggplot_gtable(ggplot_build(termMICdist))
#       gt2 <- ggplot_gtable(ggplot_build(termMICcov))
#       #newWidth = unit.pmax(gt1$widths[2:3], gt2$widths[2:3])
#       #gt1$widths[2:3] = as.list(newWidth)
#       #gt2$widths[2:3] = as.list(newWidth)
#       
#       
#       print(fname)
#       cairo_pdf(fname,width=12,height=9)
#       grid.arrange(gt1, gt2, ncol=2,widths=c(3,1))
#       #multiplot(termMICdist, termMICcov,  cols=2)
#       dev.off()
#     }
#   }
# }
# 
# # 
# # p + theme(axis.line=element_blank(),axis.text.x=element_blank(),
# #           axis.text.y=element_blank(),axis.ticks=element_blank(),
# #           axis.title.x=element_blank(),
# #           axis.title.y=element_blank(),legend.position="none",
# #           panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
# #           panel.grid.minor=element_blank(),plot.background=element_blank())
# 
# ##Coverage
# #                    !is.na(`mr$cols`) & `mr$cols`!=0 &
# print('Coverage')
# dir.create(paste(odir,'/',envar,'/Coverage',sep=''), showWarnings = TRUE, recursive = FALSE, mode = "0777")
# for(m in 1:nrow(Thresholds)) {
#   mr<-Thresholds[m,]
#   for(s in 1:nrow(Stat)) {
#     sr <- Stat[s,]
#     for(t in 1:nrow(Types)) {
#       print(mr$cols)
#       tr <- Types[t,]
#       sel<-subset(circms,Category==as.character(tr$tname)  &
#                     Size_MICo5>gtres & ! is.na(Term) &
# 
#                     ! Term %in% c('Metabolic pathways',
#                                   'Biosynthesis of secondary metabolites',
#                                   'Microbial metabolism in diverse environments'))
# 
#       gtitle<-paste('Enriched E. Coli MG1655 ',tr$ttitle,' - ',mr$names,' MIC>',micthres,sep = '')
#       fname<-paste(odir,'/',envar,'/Coverage/',tr$tabr,'_',mr$cols,'_MICo',micthres,'.pdf',sep = '')
#       print(gtitle)
#       print(fname)
#       termBG<-ggplot(sel,
#                      aes(y=sel[,as.character(mr$cols)],x=reorder(sel$Term,sel[,as.character(mr$cols)]),
#                          fill=Size_MICo5))+
#         geom_bar(stat = "identity")+
#         coord_flip()+
#         labs(fill='Number of hits')+
#         scale_fill_gradient( high="red",low='grey')+
#         ylab(mr$axis)+#ylim(0,mr$mscale)+
#         xlab(tr$ttitle)+ggtitle(gtitle)
#       #termBG
#       cairo_pdf(fname,width=9,height=9)
#       print(termBG)
#       dev.off()
#       #scale_color_gradient(limits=c(0,8), high="blue",low='gray')+
#       #,
#       #ggsave(plot=termBG,file=fname,width=9,height=9)
#       #dev.copy2pdf(device=cairo_pdf,file=fname,width=12,height=12)
#     }
#   }
# }
# 
# #mtres<-0
# ##MIC distributions
# print('MIC distributions')
# dir.create(paste(odir,'/',envar,'/MIC',sep=''), showWarnings = TRUE, recursive = FALSE, mode = "0777")
# for(s in 1:nrow(Stat)) {
#   sr <- Stat[s,]
#   for(t in 1:nrow(Types)) {
#     tr <- Types[t,]
#     # & Size_MICo5>gtres&MIC_S_Mean>mttres
#     sel<-subset(enr,Category==as.character(tr$tname) & MIC>mtres &!is.na(Term)& Size_MICo5>gtres &
#                   ! Term %in% c('Metabolic pathways',
#                                 'Biosynthesis of secondary metabolites',
#                                 'Microbial metabolism in diverse environments'))
#     #& logp_sum>2 & logp_sum>2  & logp_sum>2  & logp_sum>2
#     gtitle<-paste('Enriched E. Coli MG1655 ',tr$ttitle,' - MIC distribution, ',sr$sabr,' MIC>',micthres,sep = '')
#     fname<-paste(odir,'/',envar,'/MIC/',tr$tabr,'_MIC_',sr$sabr,'_MICo',micthres,'.pdf',sep = '')
#     print(gtitle)
#     print(fname)
#     sel$Term<-gsub('metabolism','met.',sel$Term)
#     sel$Term<-gsub('biosynthesis','bs.',sel$Term)
#     termMIC<-ggplot(sel,
#                     aes(y=MIC,x=reorder(Term,Coverage_KeioS),fill=Coverage_KeioS,color=Size_KeioS))+
#       geom_boxplot(position='identity')+
#       geom_point()+
#       coord_flip()+
#       scale_colour_gradientn(colours = c('black','yellow','red'),
#                              breaks=c(5,10,15,20),limits=c(2,20),guide='legend')+#breaks=c(5,10,15,20),limits=c(2,20)
#       labs(fill='Coverage',
#            color='Number of genes with MIC>5')+
#       scale_fill_gradientn(colours=c('grey70',"blue",'cyan'),limits=c(0,1),breaks=c(0,0.2,0.4,0.6,0.8,1)) +
#       #scale_color_gradient(limits=c(0,100), high="blue",low='gray')+
#       ylab(mic)+#ylim(1,100)+
#       scale_y_continuous(breaks=c(1,2.5,5,10,25,50,100),trans='log',limits=c(1,100))+
#       #annotation_logticks(base = 10,color='grey50')+#,sides=''
#       xlab(tr$ttitle)+ggtitle(gtitle)+
#       theme(legend.box='horizontal',legend.position='top')
#     cairo_pdf(fname,width=10,height=8)
#     print(termMIC)
#     dev.off()
#     #      
#   }
# }
# 
# 
# mtres<-0
# #mtres<-5
# ##MIC significance
# print('MIC distributions')
# dir.create(paste(odir,'/',envar,'/MICsignificance',sep=''), showWarnings = TRUE, recursive = FALSE, mode = "0777")
# for(s in 1:nrow(Stat)) {
#   sr <- Stat[s,]
#   for(t in 1:nrow(Types)) {
#     tr <- Types[t,]
#     # & Size_MICo5>gtres&MIC_S_Mean>mttres & Enrichmentp<0.1
#     sel<-subset(enr,Category==as.character(tr$tname) &
#                   MIC>mtres &
#                   #Coverage_KeioS>0 &
#                   (Enrichment_MIC_p.value < 0.1 | Term=='Pyrimidine metabolism') &
#                   !is.na(Term) &
#                   Size_MICo5>gtres  &
#                   ! Term %in% c('Metabolic pathways',
#                                 'Biosynthesis of secondary metabolites',
#                                 'Microbial metabolism in diverse environments'))
#     gtitle<-paste('Enriched E. Coli MG1655 ',tr$ttitle,' - MIC significance distribution, ',sr$sabr,' MIC>',micthres,sep = '')
#     fname<-paste(odir,'/',envar,'/MICsignificance/',tr$tabr,'_MIC_',sr$sabr,'_MICo',micthres,'.pdf',sep = '')
#     fname2<-paste(odir,'/',envar,'/MICsignificance/',tr$tabr,'_MIC_',sr$sabr,'_MICo',micthres,'_large.pdf',sep = '')
#     print(gtitle)
#     print(fname)
#     sel$Term<-gsub('metabolism','metab.',sel$Term)
#     sel$Term<-gsub('biosynthesis','bs.',sel$Term)
#     termMIC<-ggplot(sel,
#                     aes(y=MIC,x=reorder(Term,Coverage_KeioS),fill=Coverage_KeioS,color=Size_MICo5))+
#       #geom_boxplot(position='identity')+
#       geom_violin(position='identity')+
#       #geom_boxplot(stat='identity')+
#       #geom_point()+
#       coord_flip()+
#       geom_hline(yintercept=5,color='red',alpha=0.5)+
#       scale_colour_gradientn(colours = c('black','yellow','red'),breaks=c(5,10,20,30),limits=c(2,31),guide='legend')+#breaks=c(5,10,15,20),limits=c(2,20) #,guide='legend'
#       labs(fill='Coverage',
#            color='Number of genes with MIC>5')+
#       scale_fill_gradientn(colours=c('grey70',"blue",'cyan'),limits=c(0,1),breaks=c(0,0.2,0.4,0.6,0.8,1))+
#       geom_text(aes(y=0,label=as.character(stars.pval(Enrichment_MIC_FDR))),
#                 color='black',hjust=0, vjust=0,size=3)+#-Size_MICo5
#       ylab(mic)+
#       scale_y_continuous(breaks=c(1,2.5,5,10,25,50,100),limits=c(1,100),trans='log')+#,trans='log'
#       #annotation_logticks(base = 10,color='grey50')+#,sides=''
#       xlab(tr$ttitle)+ggtitle(gtitle)+
#       theme(legend.box='horizontal',
#             legend.position='top',
#             panel.grid.major=element_blank(),
#             panel.grid.minor=element_blank())
# 
#     cairo_pdf(fname,width=6,height=8)
#     print(termMIC)
#     #grid.arrange(termMIC, termMIC2, ncol=2,widths=c(2,2))
#     dev.off()
#     cairo_pdf(fname2,width=12,height=8)
#     print(termMIC)
#     #grid.arrange(termMIC, termMIC2, ncol=2,widths=c(2,2))
#     dev.off()
#     #      
#   }
# }
# 
# 


