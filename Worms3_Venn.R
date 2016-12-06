library('ggplot2')
library('gplots')
library('plyr')
library('reshape2')
library('tidyr')
library('Vennerable')
library('quantreg')

library('VennDiagram')
library('splitstackshape')
library(xlsx)


drawEuler<-function(set1,set2,set3,labels){
  draw.triple.venn(length(set1),length(set2),length(set3),
                   length(intersect(set1,set2)),length(intersect(set2,set3)),length(intersect(set3,set1)),length(intersect(intersect(set3,set1),set2)),
                   category=labels,
                   col = c("red",'blue', 'green'),
                   cat.pos=c(0,0,0))
}



odir<-'Figures_final'
ddir<-'Data_final'

bacmic<-read.table(paste(ddir,'/MICs_and_bacterial_growth-Complete.csv',sep=''),sep=',',header=TRUE)
unicom<-read.table(paste('Data/MICs_Unknown_function.csv',sep=''),sep=',',header=TRUE,stringsAsFactors = FALSE)
keioinfo<-read.table('../Keio_library/Keio_library_fully_annotated.csv',sep=',',quote = '"',header = TRUE,stringsAsFactors=FALSE)
PLP<-read.table(paste(ddir,'/Genes_using_PLP.csv',sep=''),sep=',',header=FALSE)


unknown<-read.table('../EColi_annotation/Unknown/Unknown.tab',sep='\t',quote = '"',header = TRUE,stringsAsFactors=FALSE)
phantom<-read.table('../EColi_annotation/Unknown/Phantom.tab',sep='\t',quote = '"',header = TRUE,stringsAsFactors=FALSE)
pseudo<-read.table('../EColi_annotation/Unknown/Pseudo.tab',sep='\t',quote = '"',header = TRUE,stringsAsFactors=FALSE)
orf<-read.table('../EColi_annotation/Unknown/ORF.tab',sep='\t',quote = '"',header = TRUE,stringsAsFactors=FALSE)



bacmic$Index<-paste(bacmic$Gene,bacmic$Plate,bacmic$Well,sep='')
keioinfo$Index<-paste(keioinfo$Gene,keioinfo$Plate,keioinfo$Well,sep='')

keioclean<-subset(keioinfo,!Gene %in% c('Blank','WT','present',NA) &! Plate %in% c('91','93','95'))
allk<-keioclean$Index

all<-subset(bacmic,!Gene %in% c('WT',NA))$Index

excluded<-setdiff(allk,all)
PLPg<-subset(keioclean,Gene %in% PLP$V1)$Index
MIC1<-subset(bacmic,!is.na(Gene) & MIC>1)$Index
MIC25<-subset(bacmic,!is.na(Gene) & MIC>2.5)$Index
MIC5<-subset(bacmic,!is.na(Gene) & MIC>5)$Index
unkn<-subset(keioclean,Gene %in% unknown$Name)$Index
psd<-subset(keioclean,Gene %in% pseudo$Name)$Index
phntm<-subset(keioclean,Gene %in% phantom$Name)$Index
orfl<-subset(keioclean,Gene %in% orf$Name)$Index

allu<-union(union(unkn,psd),union(phntm,orfl))


PLPla<-list('Screened Keio library'=all,'PLP using'=PLPg,'MIC>5'=MIC5)
plot(Venn(PLPla),show = list(Faces = FALSE),doWeights = FALSE)
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Venn_PLP_using_enzymes_ScreenedKeio.pdf",sep=''),width=6,height=6)



#Draw simple
dev.off()
drawEuler(all,
          PLPg,
          MIC5,
          c('Screened Keio library',
            'PLP using',
            'MIC>5'))
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Venn_PLP_using_enzymes_ScreenedKeio_Euler.pdf",sep=''),width=4,height=4)
dev.off()



PLPlk<-list('PLP using'=PLPg,'Whole Keio library'=allk,'MIC>5'=MIC5)
plot(Venn(PLPlk),show = list(Faces = FALSE), doWeights = FALSE)
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Venn_PLP_using_enzymes_WholeKeio.pdf",sep=''),width=6,height=6)



#Unknowns

unwk<-list('Whole Keio library'=setdiff(allk,'WT'),
          'Unknown function'=unkn,
          'MIC>5'=MIC5)
plot(Venn(unwk),show = list(Faces = FALSE), doWeights = FALSE)
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Venn_Unknown_enzymes_WholeKeio.pdf",sep=''),width=6,height=6)


unw<-list('Screened Keio library'=setdiff(all,'WT'),
          'Unknown function'=unkn,
          'MIC>5'=MIC5)
plot(Venn(unw),show = list(Faces = FALSE), doWeights = FALSE)
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Venn_Unknown_enzymes_ScreenedKeio.pdf",sep=''),width=6,height=6)



#Draw simple
dev.off()
drawEuler(setdiff(all,'WT'),
          unkn,
          MIC5,
          c('Screened Keio library',
            'Unknown function',
            'MIC>5'))
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Venn_Unknown_enzymes_ScreenedKeio_Euler.pdf",sep=''),width=4,height=4)
dev.off()



proba<-list('Whole Keio library'=allk,
           'Problematic'=allu,
           'MIC>5'=MIC5)
plot(Venn(proba),show = list(Faces = FALSE), doWeights = FALSE)
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Venn_Unknown-pseudo-orf-phantom_WholeKeio.pdf",sep=''),width=6,height=6)


#Draw simple
dev.off()
drawEuler(allk,
          allu,
          MIC5,
          c('Whole Keio library',
            'Problematic',
            'MIC>5'))
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Venn_Unknown-pseudo-orf-phantom_WholeKeio_Euler.pdf",sep=''),width=4,height=4)
dev.off()



prob<-list('Screened Keio library'=all,
          'Problematic'=intersect(allu,all),
          'MIC>5'=MIC5)
plot(Venn(prob),show = list(Faces = FALSE), doWeights = FALSE)
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Venn_Unknown-pseudo-orf-phantom_ScreenedKeio.pdf",sep=''),width=6,height=6)


#Draw simple
dev.off()
drawEuler(all,
          allu,
          MIC5,
          c('Screened Keio library',
            'Problematic',
            'MIC>5'))
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Venn_Unknown-pseudo-orf-phantom_ScreenedKeio_Euler.pdf",sep=''),width=4,height=4)
dev.off()





q05<-quantile(bacmic$MIC,0.05,na.rm=TRUE)[[1]]
q10<-quantile(bacmic$MIC,0.1,na.rm=TRUE)[[1]]
q90<-quantile(bacmic$MIC,0.9,na.rm=TRUE)[[1]]
q95<-quantile(bacmic$MIC,0.95,na.rm=TRUE)[[1]]
q99<-quantile(bacmic$MIC,0.99,na.rm=TRUE)[[1]]



bacGTsigup<-subset(bacmic,GT_Interaction>0 & GT_p.value<0.05)$Index
bacGTsigdown<-subset(bacmic,GT_Interaction<0 & GT_p.value<0.05)$Index

bacGTFDRup<-subset(bacmic,GT_Interaction>0 & GT_FDR<0.05)$Index
bacGTFDRdown<-subset(bacmic,GT_Interaction<0 & GT_FDR<0.05)$Index
# 
# bacODsigup<-subset(bacmic,CTODDiff_Mean>0 &CTODDiff_pval<0.05)$Index
# bacODsigdown<-subset(bacmic,CTODDiff_Mean<0 &CTODDiff_pval<0.05)$Index



wormres95<-subset(bacmic,MIC>q95)$Index
wormsens05<-subset(bacmic,MIC<q05)$Index

wormres90<-subset(bacmic,MIC>q90)$Index
wormres25<-subset(bacmic,MIC>2.5)$Index
wormresall<-subset(bacmic,MIC>1)$Index
wormresall5<-subset(bacmic,MIC>5)$Index
wormsens10<-subset(bacmic,MIC<q10)$Index


resmic5=list('Bacteria sensitive\nbottom 5%'=bacsens05,
            'Bacteria resistant\ntop 5%'=bacres95,
            'Worms resistant\ntop 5%'=wormres95,
            'Worms sensitive\nbottom 5%'=wormsens05)

bacFDRmictop10=list('Worms resistant top 10%'=wormres90,
                    'KO synergistic'=bacGTFDRdown,
                    'KO antagonistic'=bacGTFDRup) # ,'Worms sensitive bottom 10%'=wormsens10 ,'Worms resistant all'=wormresall5


bacFDRmictop5=list('Worms resistant top 5%'=wormres95,
                   'KO synergistic'=bacGTFDRdown,
                   'KO antagonistic'=bacGTFDRup) # ,'Worms sensitive bottom 10%'=wormsens10 ,'Worms resistant all'=wormresall5



bacsigmictop10=list('Worms resistant top 10%'=wormres90,
              'KO synergistic'=bacGTsigdown,
              'KO antagonistic'=bacGTsigup) # ,'Worms sensitive bottom 10%'=wormsens10 ,'Worms resistant all'=wormresall5


bacsigmictop5=list('Worms resistant top 5%'=wormres95,
                    'KO synergistic'=bacGTsigdown,
                    'KO antagonistic'=bacGTsigup) # ,'Worms sensitive bottom 10%'=wormsens10 ,'Worms resistant all'=wormresall5


bacsigmico5=list('Worms MIC>5'=wormresall5,
                    'KO synergistic'=bacGTsigdown,
                    'KO antagonistic'=bacGTsigup) # ,'Worms sensitive bottom 10%'=wormsens10 ,'Worms resistant all'=wormresall5


plot(Venn(bacFDRmictop10),show = list(Faces = FALSE), doWeights = FALSE)
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Venn_Wormtop10-Bacteria_ByFDR.pdf",sep=''),
             width=6,height=6)


plot(Venn(bacFDRmictop5),show = list(Faces = FALSE), doWeights = FALSE)
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Venn_Wormtop5-Bacteria_ByFDR.pdf",sep=''),
             width=6,height=6)



plot(Venn(bacsigmictop10),show = list(Faces = FALSE), doWeights = FALSE)
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Venn_Wormtop10-Bacteria_BySignificance.pdf",sep=''),
             width=6,height=6)


plot(Venn(bacsigmictop5),show = list(Faces = FALSE), doWeights = FALSE)
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Venn_Wormtop5-Bacteria_BySignificance.pdf",sep=''),
             width=6,height=6)

plot(Venn(bacsigmico5),show = list(Faces = FALSE), doWeights = FALSE)
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Venn_WormMICo5-Bacteria_BySignificance.pdf",sep=''),
             width=6,height=6)
# 
# #Draw simple
# dev.off()
# drawEuler(wormres90,
#           bacODnormsigdown,
#           bacODnormsigup,
#           c('Worms resistant top 10%',
#             'Bacteria significantly\nsensitive (Trend)',
#             'Bacteria significantly\nresistant (Trend)'))
# dev.copy2pdf(device=cairo_pdf,
#              file=paste(odir,"/Venn_Wormtop10-Bacteria_BySignificance_Euler.pdf",sep=''),width=4,height=4)
# dev.off()



unknowns<-list('Screened'=all,
               'Excluded'=excluded,
               'Unknown function'=unkn,
               'Pseudo gene'=psd,
               'Phantom gene'=phntm,
               'ORF'=orfl)


res<-list('Worms resistant top 5%'=wormres95,
          'Worms resistant top 10%'=wormres90,
          'Bacteria KO synergistic to 5FU'=bacGTsigdown,
          'Bacteria KO antagonistic to 5FU'=bacGTsigup,
          'Bacteria KO synergistic to 5FU (FDR)'=bacGTFDRdown,
          'Bacteria KO antagonistic to 5FU (FDR)'=bacGTFDRup)


lists<-c(unknowns,PLPla['PLP using'],res)

inclusion<-function(data,ilist){
  for (n in names(ilist)){
    print(n)
    data[data$Index %in% ilist[n][[1]],as.character(n)]<-'TRUE'
    data[!data$Index %in% ilist[n][[1]],as.character(n)]<-''
  }
  return(data)
}

keioMIC<-merge(keioclean[,c('Gene','Plate','Well','Index')],bacmic,by=c('Gene','Plate','Well','Index'),all.x=TRUE,all.y=TRUE)

bm<-inclusion(keioMIC,lists)

bmv<-bm[,colnames(bm) %in% c('Gene','Plate','Well','MIC',names(lists))]



vennexpl<-read.table(paste(ddir,'/Venn_Worm_Bacteria_column_explanation.csv',sep=''),
                      sep=',',quote = '"',header = TRUE,stringsAsFactors=FALSE)




explrm<-subset(vennexpl,Column %in% colnames(bmv))
explrm<-explrm[match(colnames(bmv),explrm$Column),]

write.csv(bmv,paste(ddir,'/Venn_Worm_Bacteria.csv',sep=''),row.names = FALSE)

write.xlsx2(explrm, file=paste(ddir,'/Venn_Worm_Bacteria.xlsx',sep=''),
           sheetName="Readme",row.names = FALSE,showNA=FALSE)
write.xlsx2(bmv, file=paste(ddir,'/Venn_Worm_Bacteria.xlsx',sep=''),
           sheetName="Data", append=TRUE,row.names = FALSE,showNA=FALSE)






#Biolog update later
allwb<-read.table(paste(ddir,'/Biolog_Combined_Summary_Statistics.csv',sep=''),sep=',',header=TRUE)

odir<-'Figures_final/Biolog'

#Venn
# 

w1<-subset(allwb,W_Median>1)$UniqueName
w2<-subset(allwb,W_Median>2)$UniqueName
w3<-subset(allwb,W_Median>3)$UniqueName


sig_A_MT_p<-subset(allwb,MT_p.value < 0.05 & MT_Interaction > 0)$UniqueName
sig_S_MT_p<-subset(allwb,MT_p.value < 0.05 & MT_Interaction < 0)$UniqueName
sig_A_MT_FDR<-subset(allwb,MT_FDR < 0.05 & MT_Interaction > 0)$UniqueName
sig_S_MT_FDR<-subset(allwb,MT_FDR < 0.05 & MT_Interaction < 0)$UniqueName



#Just sets
sigMet<-list('Worms > 1'=w1,
              'Synergistic'=sig_S_MT_FDR,
              'Antagonistic'=sig_A_MT_FDR)
plot(Venn(sigMet),show = list(Faces = FALSE), doWeights = FALSE)#
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Biolog_Venn_FDR.pdf",sep=''),width=6,height=6)


#Just sets
sigMetp<-list('Worms > 1'=w1,
             'Synergistic'=sig_S_MT_p,
             'Antagonistic'=sig_A_MT_p)
plot(Venn(sigMetp),show = list(Faces = FALSE), doWeights = FALSE)#
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Biolog_Venn_pval.pdf",sep=''),width=6,height=6)

# 
# 
# #Make significant hit comparison table
# lists<-list('Control: Metabolite+NGM/NGM'=sig_C,'Control: Metabolite+NGM/NGM positive'=sig_Cp,'Control: Metabolite+NGM/NGM negative'=sig_Cm,
#             'Treatment: Metabolite+NGM/NGM'=sig_T,'Treatment: Metabolite+NGM/NGM positive'=sig_Tp,'Treatment: Metabolite+NGM/NGM negative'=sig_Tm,
#             'Treatment/Control'=sig_CT,'Treatment/Control positive'=sig_CTp,'Treatment/Control negative'=sig_CTm,
#             'Treatment/Trend'=sig_CTnorm,'Treatment/Trend positive'=sig_CTnormp,'Treatment/Trend negative'=sig_CTnormm)
# 
# sigcols<-list('Control: Metabolite+NGM/NGM'='NGMDiff_C_pval','Control: Metabolite+NGM/NGM positive'='NGMDiff_C_pval','Control: Metabolite+NGM/NGM negative'='NGMDiff_C_pval',
#               'Treatment: Metabolite+NGM/NGM'='NGMDiff_T_pval','Treatment: Metabolite+NGM/NGM positive'='NGMDiff_T_pval','Treatment: Metabolite+NGM/NGM negative'='NGMDiff_T_pval',
#               'Treatment/Control'='CTDiff_pval','Treatment/Control positive'='CTDiff_pval','Treatment/Control negative'='CTDiff_pval',
#               'Treatment/Trend'='CTDiff_norm_pval','Treatment/Trend positive'='CTDiff_norm_pval','Treatment/Trend negative'='CTDiff_norm_pval')
# 
# datacols<-list('Control: Metabolite+NGM/NGM'='NGMDiff_C_Mean','Control: Metabolite+NGM/NGM positive'='NGMDiff_C_Mean','Control: Metabolite+NGM/NGM negative'='NGMDiff_C_Mean',
#                'Treatment: Metabolite+NGM/NGM'='NGMDiff_T_Mean','Treatment: Metabolite+NGM/NGM positive'='NGMDiff_T_Mean','Treatment: Metabolite+NGM/NGM negative'='NGMDiff_T_Mean',
#                'Treatment/Control'='CTDiff_Mean','Treatment/Control positive'='CTDiff_Mean','Treatment/Control negative'='CTDiff_Mean',
#                'Treatment/Trend'='CTDiff_norm_Mean','Treatment/Trend positive'='CTDiff_norm_Mean','Treatment/Trend negative'='CTDiff_norm_Mean')
# 
# allsigdf<-data.frame('SigGroup'=NA,'Plate'=NA,'Well'=NA,'UniqueName'=NA,'pvalue'=NA,'SigLabel'=NA)
# for (nm in names(lists)){
#   print(nm)
#   lst<-lists[nm]
#   print(length(lst[[1]]))
#   sigcol<-as.character(sigcols[nm][[1]])
#   sigdat<-as.character(datacols[nm][[1]])
#   if (length(grep("positive",nm))>0){
#     datasel<-allwb[allwb[,sigcol]<0.05,]
#     datasel<-datasel[datasel[,sigdat]>0,]
#   } else if (length(grep("negative",nm))>0) {
#     datasel<-allwb[allwb[,sigcol]<0.05,]
#     datasel<-datasel[datasel[,sigdat]<0,]
#   } else {
#     datasel<-allwb[allwb[,sigcol]<0.05,]
#   }
#   
#   print(dim(datasel))
#   datasel$SigGroup<-nm
#   print(as.character(sigcols[nm][[1]]))
#   
#   datasel[,'pvalue']<-datasel[,sigcol]
#   datasel$SigLabel<-ifelse(stars.pval(datasel$pvalue)!=' ',stars.pval(datasel$pvalue),'')
#   allsigdf<-merge(allsigdf,datasel[,c('SigGroup','Plate','Well','UniqueName','pvalue','SigLabel')],all.x=TRUE,all.y=TRUE)
# }
# 
# allsigdf<-subset(allsigdf,Plate!='NA')
# significance<-dcast(allsigdf,Plate+Well+UniqueName~SigGroup,value.var = c('SigLabel'),fill=as.character(''),fun.aggregate = NULL)
# significancecomp<-merge(allwb[,c('Plate','Well','UniqueName','W_Median')],significance,by=c('Plate','Well','UniqueName'),all.x=TRUE)
# 
# 
# 
# 
# blexpl<-read.table(paste(ddir,'/Venn_Biolog_column_explanation.csv',sep=''),
#                      sep=',',quote = '"',header = TRUE,stringsAsFactors=FALSE)
# 
# 
# 
# 
# explrm<-subset(blexpl,Column %in% colnames(significancecomp))
# explrm<-explrm[match(colnames(significancecomp),explrm$Column),]
# 
# 
# write.csv(significancecomp,paste(ddir,'/Biolog_Significant_hit_comparison.csv',sep=''),row.names = FALSE,na = "")
# 
# 
# write.xlsx2(explrm, file=paste(ddir,'/Biolog_Significant_hit_comparison.xlsx',sep=''),
#            sheetName="Readme",row.names = FALSE,showNA=FALSE)
# write.xlsx2(significancecomp, file=paste(ddir,'/Biolog_Significant_hit_comparison.xlsx',sep=''),
#            sheetName="Data", append=TRUE,row.names = FALSE,showNA=FALSE)
# 
# 
# 
# write.xlsx2(explrm, file='/Users/Povilas/Projects/B-D-H paper/figures and data/figure 5/final files/table S5.xlsx',
#            sheetName="Biolog_Hits_Readme", append=TRUE,row.names = FALSE,showNA=FALSE)
# write.xlsx2(significancecomp, file='/Users/Povilas/Projects/B-D-H paper/figures and data/figure 5/final files/table S5.xlsx',
#            sheetName="Biolog_Hits_Data", append=TRUE,row.names = FALSE,showNA=FALSE)
# 
# 

