library('ggplot2')
library('gplots')
library('plyr')
library('reshape2')
library('tidyr')
library('Vennerable')
library('quantreg')

library('splitstackshape')
library('VennDiagram')
library('splitstackshape')


odir<-'Figures_final'
ddir<-'Data_final'

bacmic<-read.table(paste(ddir,'/MICs_and_bacterial_growth-Complete.csv',sep=''),sep=',',header=TRUE)
unicom<-read.table(paste('Data/MICs_Unknown_function.csv',sep=''),sep=',',header=TRUE,stringsAsFactors = FALSE)
allwb<-read.table(paste(ddir,'/Biolog_Combined_Summary_Statistics.csv',sep=''),sep=',',header=TRUE)


keioinfo<-read.table('../Keio_library/Keio_library_fully_annotated.csv',sep=',',quote = '"',header = TRUE,stringsAsFactors=FALSE)


allk<-subset(keioinfo,!Gene %in% c('Blank','WT'))$Gene
PLP<-read.table('Data/Genes_using_PLP.csv',sep=',',header=FALSE)
all<-subset(bacmic,!is.na(Gene))$Gene
PLPg<-PLP$V1
MIC1<-subset(bacmic,!is.na(Gene) & MIC>1)$Gene
MIC25<-subset(bacmic,!is.na(Gene) & MIC>2.5)$Gene
MIC5<-subset(bacmic,!is.na(Gene) & MIC>5)$Gene


PLPl<-list('PLP using'=as.character(PLPg),'Screened Keio library'=all,'MIC>5'=MIC5)
plot(Venn(PLPl), doWeights = FALSE)#,type='ellipses',type='ellipses',type='ellipses'
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/PLP_using_enzymes_ScreenedKeio.pdf",sep=''),width=9,height=9)


PLPl<-list('PLP using'=as.character(PLPg),'Whole Keio library'=allk,'MIC>5'=MIC5)
plot(Venn(PLPl), doWeights = FALSE)#,type='ellipses',type='ellipses',type='ellipses'
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/PLP_using_enzymes_WholeKeio.pdf",sep=''),width=9,height=9)




PLPc<-PLPl[c('PLP using','Whole Keio library','MIC>5')]
PLP_df<-plyr::ldply(PLPc, cbind)
PLP_df<-rename(PLP_df,c('.id'='List','1'='Gene'))


plps<-subset(bacmic, Gene %in% unique(PLP_df$Gene))
plps$'PLP using'<-ifelse(plps$Gene %in% PLPc$'PLP using',TRUE,FALSE)
plps$'MIC>1'<-ifelse(plps$Gene %in% PLPc$'MIC>1',TRUE,FALSE)
plps$'MIC>5'<-ifelse(plps$Gene %in% PLPc$'MIC>5',TRUE,FALSE)
plps<-plps[,c(1,19:21,2:7,9:18,8)]
write.csv(plps,paste(ddir,'/PLP_use_MIC_overlap.csv',sep=''))



unw<-list('Whole Keio library'=allk,
          'Unknown function'=subset(unicom,Unknown.function=='TRUE' & Gene %in% all)$Gene,
          'MIC>5'=MIC5)
plot(Venn(unw), doWeights = FALSE)
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Unknown_enzymes_3sets_WholeKeio.pdf",sep=''),width=9,height=9)


unw<-list('Whole Keio library'=all,
          'Unknown function'=subset(unicom,Unknown.function=='TRUE' & Gene %in% all)$Gene,
          'MIC>5'=MIC5)
plot(Venn(unw), doWeights = FALSE)
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Unknown_enzymes_3sets_ScreenedKeio.pdf",sep=''),width=9,height=9)






fitqr<-rq(OD_T_Mean ~ OD_C_Mean,data=bacmic,tau=c(0.05,0.95))
bgli<-coefficients(fitqr)[1,][[1]]
bgui<-coefficients(fitqr)[1,][[2]]
bgls<-coefficients(fitqr)[2,][[1]]
bgus<-coefficients(fitqr)[2,][[2]]

fitqr2<-rq(OD_T_Mean ~ OD_C_Mean,data=bacmic,tau=c(0.10,0.90))
bgli2<-coefficients(fitqr2)[1,][[1]]
bgui2<-coefficients(fitqr2)[1,][[2]]
bgls2<-coefficients(fitqr2)[2,][[1]]
bgus2<-coefficients(fitqr2)[2,][[2]]


q05<-quantile(bacmic$MIC,0.05,na.rm=TRUE)[[1]]
q10<-quantile(bacmic$MIC,0.1,na.rm=TRUE)[[1]]
q90<-quantile(bacmic$MIC,0.9,na.rm=TRUE)[[1]]
q95<-quantile(bacmic$MIC,0.95,na.rm=TRUE)[[1]]
q99<-quantile(bacmic$MIC,0.99,na.rm=TRUE)[[1]]



bacODsigup<-subset(bacmic,CTODDiff_norm_Mean>0 &CTODDiff_norm_pval<0.05)$Gene
bacODsigdown<-subset(bacmic,CTODDiff_norm_Mean<0 &CTODDiff_norm_pval<0.05)$Gene


bacres95<-subset(bacmic,OD_T_Mean>OD_C_Mean*bgus+bgui)$Gene
bacsens05<-subset(bacmic,OD_T_Mean<OD_C_Mean*bgls+bgli)$Gene

bacres90<-subset(bacmic,OD_T_Mean>OD_C_Mean*bgus2+bgui2)$Gene
bacsens10<-subset(bacmic,OD_T_Mean<OD_C_Mean*bgls2+bgli2)$Gene

wormres95<-subset(bacmic,MIC>q95)$Gene
wormsens05<-subset(bacmic,MIC<q05)$Gene

wormres90<-subset(bacmic,MIC>q90)$Gene
wormres25<-subset(bacmic,MIC>2.5)$Gene
wormresall<-subset(bacmic,MIC>1)$Gene
wormresall5<-subset(bacmic,MIC>5)$Gene
wormsens10<-subset(bacmic,MIC<q10)$Gene

resmic5=list('Bacteria sensitive bottom 5%'=bacsens05,
            'Bacteria resistant top 5%'=bacres95,
            'Worms resistant top 5%'=wormres95,
            'Worms sensitive bottom 5%'=wormsens05)




resmic10=list('Worms resistant top 10%'=wormres90,
              'Bacteria sensitive bottom 10%'=bacsens10,
              'Bacteria resistant top 10%'=bacres90) # ,'Worms sensitive bottom 10%'=wormsens10 ,'Worms resistant all'=wormresall5


bacsigmictop10=list('Worms resistant top 10%'=wormres90,
              'Bacteria significant sensitive (Trend)'=bacODsigdown,
              'Bacteria significantly resistant (Trend)'=bacODsigup) # ,'Worms sensitive bottom 10%'=wormsens10 ,'Worms resistant all'=wormresall5


bacsigmico5=list('Worms MIC>5'=wormresall5,
                    'Bacteria significant sensitive (Trend)'=bacODsigdown,
                    'Bacteria significantly resistant (Trend)'=bacODsigup) # ,'Worms sensitive bottom 10%'=wormsens10 ,'Worms resistant all'=wormresall5




plot(Venn(resmic10), doWeights = FALSE)#,type='ellipses'  ,type='ellipses' ,type='ellipses'
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Venn_Wormtop10-Bacteria_ByBacPercentage.pdf",sep=''),
             width=9,height=9)


plot(Venn(bacsigmictop10), doWeights = FALSE)
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Venn_Wormtop10-Bacteria_BySignificance.pdf",sep=''),
             width=9,height=9)


plot(Venn(bacsigmico5), doWeights = FALSE)
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Venn_WormMICo5-Bacteria_BySignificance.pdf",sep=''),
             width=9,height=9)





resmic_df<-plyr::ldply(resmic10, cbind)
resmic_df<-rename(resmic_df,c('.id'='List','1'='Gene'))

#Summary
rmc<-subset(bacmic, Gene %in% unique(resmic_df$Gene))
rmc$'Bacteria resistant top 5%'<-ifelse(rmc$Gene %in% subset(resmic_df,List=="Bacteria resistant top 5%")$Gene,TRUE,FALSE)
rmc$'Bacteria sensitive bottom 5%'<-ifelse(rmc$Gene %in% subset(resmic_df,List=="Bacteria sensitive bottom 5%")$Gene,TRUE,FALSE)
rmc$'Worms resistant top 5%'<-ifelse(rmc$Gene %in% subset(resmic_df,List=="Worms resistant top 5%")$Gene,TRUE,FALSE)
rmc$'Worms sensitive bottom 5%'<-ifelse(rmc$Gene %in% subset(resmic_df,List=="Worms sensitive bottome 5%%")$Gene,TRUE,FALSE)
rmc<-rmc[,c(1,19:21,2:7,9:18,8)]
write.csv(rmc,paste(ddir,'/Venn_Worm-Bacteria_resistance_overlap.csv',sep=''))



#Biolog

#Venn
# 



w1<-subset(allwb,W_Median>1)$UniqueName
w2<-subset(allwb,W_Median>2)$UniqueName
w3<-subset(allwb,W_Median>3)$UniqueName

sig_C<-subset(allwb,NGMDiff_C_pval<0.05)$UniqueName
sig_T<-subset(allwb,NGMDiff_T_pval<0.05)$UniqueName
sig_CT<-subset(allwb,CTDiff_pval<0.05)$UniqueName
sig_CTnorm<-subset(allwb,CTDiff_norm_pval<0.05)$UniqueName

sig_Cp<-subset(allwb,NGMDiff_C_pval<0.05 & NGMDiff_C_Mean>0)$UniqueName
sig_Tp<-subset(allwb,NGMDiff_T_pval<0.05 & NGMDiff_T_Mean>0)$UniqueName
sig_CTp<-subset(allwb,CTDiff_pval<0.05 & CTDiff_Mean>0)$UniqueName
sig_CTnormp<-subset(allwb,CTDiff_norm_pval<0.05 & CTDiff_norm_Mean>0)$UniqueName

sig_Cm<-subset(allwb,NGMDiff_C_pval<0.05 & NGMDiff_C_Mean<0)$UniqueName
sig_Tm<-subset(allwb,NGMDiff_T_pval<0.05 & NGMDiff_T_Mean<0)$UniqueName
sig_CTm<-subset(allwb,CTDiff_pval<0.05 & CTDiff_Mean<0)$UniqueName
sig_CTnormm<-subset(allwb,CTDiff_norm_pval<0.05 & CTDiff_norm_Mean<0)$UniqueName

#Just sets
sigMet<-list('Treatment: Metabolite/NGM'=sig_T,'Control: Metabolite/NGM'=sig_C,'Treatment/Trend (NGM relative)'=sig_CTnorm,'Treatment/Control (NGM relative)'=sig_CT)
plot(Venn(sigMet), doWeights = FALSE,type='ellipses')
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"Bacteria_Venn_Significant-metabolites.pdf",sep=''),width=9,height=9)

#Positive negative Control and treatment
sigMetpm<-list('Control: Metabolite/NGM -'=sig_Cm,'Control: Metabolite/NGM +'=sig_Cp,'Treatment: Metabolite/NGM -'=sig_Tm,'Treatment: Metabolite/NGM +'=sig_Tp)
plot(Venn(sigMetpm), doWeights = FALSE,type='ellipses')
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"Bacteria_Venn_Significant-positive|negative-metabolites.pdf",sep=''),width=9,height=9)


#Worms>1 CT pos neg
sigW1CT<-list('Worm growth rescue>1'=w1,'Treatment/Control (NGM relative) +'=sig_CTp,'Treatment/Control (NGM relative) -'=sig_CTm)#
plot(Venn(sigW1CT),doWeights=FALSE,doEuler=TRUE)#,type='ellipses'
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"Bacteria_Venn_Worms1-CTboth.pdf",sep=''),width=9,height=9)


#Draw simple
dev.off()
drawEuler(sig_CTp,w1,sig_CTm,c('Treatment/Control (NGM relative) +','Worm growth rescue >1','Treatment/Control (NGM relative) -'))
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"Bacteria_Venn_Worms1-CTboth_simpleEuler.pdf",sep=''),width=9,height=9)
dev.off()


#Worms>1 CTnorm pos neg
sigW1CTnorm<-list('Worm growth rescue >1'=w1,'Treatment/Trend (NGM relative) -'=sig_CTnormm,'Treatment/Trend (NGM relative) +'=sig_CTnormp)
plot(Venn(sigW1CTnorm), doWeights = FALSE)#,type='ellipses'
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"Bacteria_Venn_Worms1-CTnormboth.pdf",sep=''),width=9,height=9)

#Draw simple
dev.off()
drawEuler(sig_CTnormp,w1,sig_CTnormm,c('Treatment/Trend (NGM relative) +','Worm growth rescue >1','Treatment/Trend (NGM relative) -'))
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"Bacteria_Venn_Worms1-CTnormboth_simpleEuler.pdf",sep=''),width=9,height=9)
dev.off()


#Worms>2 CTnorm pos neg
sigW2CTnorm<-list('Worm growth rescue>2'=w2,'Treatment/Trend (NGM relative) -'=sig_CTnormm,'Treatment/Trend (NGM relative) +'=sig_CTnormp)
plot(Venn(sigW2CTnorm), doWeights = FALSE)#,type='ellipses'
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"Bacteria_Venn_Worms2-CTnormboth.pdf",sep=''),width=9,height=9)

#Draw simple
dev.off()
drawEuler(sig_CTnormp,w2,sig_CTnormm,c('Treatment/Trend (NGM relative) +','Worm growth rescue >2','Treatment/Trend (NGM relative) -'))
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"Bacteria_Venn_Worms2-CTnormboth_simpleEuler.pdf",sep=''),width=9,height=9)
dev.off()


#Make significant hit comparison table
lists<-list('Control: Metabolite+NGM/NGM'=sig_C,'Control: Metabolite+NGM/NGM positive'=sig_Cp,'Control: Metabolite+NGM/NGM negative'=sig_Cm,
            'Treatment: Metabolite+NGM/NGM'=sig_T,'Treatment: Metabolite+NGM/NGM positive'=sig_Tp,'Treatment: Metabolite+NGM/NGM negative'=sig_Tm,
            'Treatment/Control'=sig_CT,'Treatment/Control positive'=sig_CTp,'Treatment/Control negative'=sig_CTm,
            'Treatment/Trend'=sig_CTnorm,'Treatment/Trend positive'=sig_CTnormp,'Treatment/Trend negative'=sig_CTnormm)

sigcols<-list('Control: Metabolite+NGM/NGM'='NGMDiff_C_pval','Control: Metabolite+NGM/NGM positive'='NGMDiff_C_pval','Control: Metabolite+NGM/NGM negative'='NGMDiff_C_pval',
              'Treatment: Metabolite+NGM/NGM'='NGMDiff_T_pval','Treatment: Metabolite+NGM/NGM positive'='NGMDiff_T_pval','Treatment: Metabolite+NGM/NGM negative'='NGMDiff_T_pval',
              'Treatment/Control'='CTDiff_pval','Treatment/Control positive'='CTDiff_pval','Treatment/Control negative'='CTDiff_pval',
              'Treatment/Trend'='CTDiff_norm_pval','Treatment/Trend positive'='CTDiff_norm_pval','Treatment/Trend negative'='CTDiff_norm_pval')

datacols<-list('Control: Metabolite+NGM/NGM'='NGMDiff_C_Mean','Control: Metabolite+NGM/NGM positive'='NGMDiff_C_Mean','Control: Metabolite+NGM/NGM negative'='NGMDiff_C_Mean',
               'Treatment: Metabolite+NGM/NGM'='NGMDiff_T_Mean','Treatment: Metabolite+NGM/NGM positive'='NGMDiff_T_Mean','Treatment: Metabolite+NGM/NGM negative'='NGMDiff_T_Mean',
               'Treatment/Control'='CTDiff_Mean','Treatment/Control positive'='CTDiff_Mean','Treatment/Control negative'='CTDiff_Mean',
               'Treatment/Trend'='CTDiff_norm_Mean','Treatment/Trend positive'='CTDiff_norm_Mean','Treatment/Trend negative'='CTDiff_norm_Mean')

allsigdf<-data.frame('SigGroup'=NA,'Plate'=NA,'Well'=NA,'UniqueName'=NA,'pvalue'=NA,'SigLabel'=NA)
for (nm in names(lists)){
  print(nm)
  lst<-lists[nm]
  print(length(lst[[1]]))
  sigcol<-as.character(sigcols[nm][[1]])
  sigdat<-as.character(datacols[nm][[1]])
  if (length(grep("positive",nm))>0){
    datasel<-allwb[allwb[,sigcol]<0.05,]
    datasel<-datasel[datasel[,sigdat]>0,]
  } else if (length(grep("negative",nm))>0) {
    datasel<-allwb[allwb[,sigcol]<0.05,]
    datasel<-datasel[datasel[,sigdat]<0,]
  } else {
    datasel<-allwb[allwb[,sigcol]<0.05,]
  }
  
  print(dim(datasel))
  datasel$SigGroup<-nm
  print(as.character(sigcols[nm][[1]]))
  
  datasel[,'pvalue']<-datasel[,sigcol]
  datasel$SigLabel<-ifelse(stars.pval(datasel$pvalue)!=' ',stars.pval(datasel$pvalue),'')
  allsigdf<-merge(allsigdf,datasel[,c('SigGroup','Plate','Well','UniqueName','pvalue','SigLabel')],all.x=TRUE,all.y=TRUE)
}

allsigdf<-subset(allsigdf,Plate!='NA')
significance<-dcast(allsigdf,Plate+Well+UniqueName~SigGroup,value.var = c('SigLabel'),fill=as.character(''),fun.aggregate = NULL)
significancecomp<-merge(allwb[,c('Plate','Well','UniqueName','W_Median')],significance,by=c('Plate','Well','UniqueName'),all.x=TRUE)
write.csv(significancecomp,paste(ddir,'/Biolog_Significant_hit_comparison.csv',sep=''))






