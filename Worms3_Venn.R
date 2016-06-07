library('ggplot2')
library('gplots')
library('plyr')
library('reshape2')
library('tidyr')
library('Vennerable')
#library('GOplot')
library('quantreg')
#library('KEGGREST')
#library('devtools')
library('splitstackshape')


odir<-'Figures_final'
ddir<-'Data_final'

qbacmicq<-read.table(paste(ddir,'/MICs_and_bacterial_growth-Clean.csv',sep=''),sep=',',header=TRUE)
unicom<-read.table(paste('Data/MICs_Unknown_function.csv',sep=''),sep=',',header=TRUE,stringsAsFactors = FALSE)




PLP<-read.table('Data/Genes_using_PLP.csv',sep=',',header=FALSE)
all<-subset(qbacmicq,!is.na(Gene))$Gene
PLPg<-PLP$V1
MIC1<-subset(qbacmicq,!is.na(Gene) & MIC>1)$Gene
MIC25<-subset(qbacmicq,!is.na(Gene) & MIC>2.5)$Gene
MIC5<-subset(qbacmicq,!is.na(Gene) & MIC>5)$Gene


PLPl<-list('PLP using'=as.character(PLPg),'Whole Keio library'=all,'MIC>5'=MIC5)

plot(Venn(PLPl), doWeights = FALSE)#,type='ellipses',type='ellipses',type='ellipses'
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/PLP_using_enzymes.pdf",sep=''),width=9,height=9)


PLPc<-PLPl[c('PLP using','Whole Keio library','MIC>5')]
PLP_df<-plyr::ldply(PLPc, cbind)
PLP_df<-rename(PLP_df,c('.id'='List','1'='Gene'))


plps<-subset(qbacmicq, Gene %in% unique(PLP_df$Gene))
plps$'PLP using'<-ifelse(plps$Gene %in% PLPc$'PLP using',TRUE,FALSE)
plps$'MIC>1'<-ifelse(plps$Gene %in% PLPc$'MIC>1',TRUE,FALSE)
plps$'MIC>5'<-ifelse(plps$Gene %in% PLPc$'MIC>5',TRUE,FALSE)
plps<-plps[,c(1,19:21,2:7,9:18,8)]
write.csv(plps,paste(ddir,'/PLP_use_MIC_overlap.csv',sep=''))



unw<-list('Whole Keio library'=all,'Unknown function'=subset(unicom,Unknown.function=='TRUE' & Gene %in% all)$Gene,'MIC>5'=MIC25)
plot(Venn(unw), doWeights = FALSE)#,type='ellipses',type='ellipses' ,type='ellipses'  ,doEuler=TRUE,
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Unknown_enzymes_3sets.pdf",sep=''),width=9,height=9)





fitqr<-rq(NGM_D ~ NGM_C,data=qbacmicq,tau=c(0.05,0.95))
bgli<-coefficients(fitqr)[1,][[1]]
bgui<-coefficients(fitqr)[1,][[2]]
bgls<-coefficients(fitqr)[2,][[1]]
bgus<-coefficients(fitqr)[2,][[2]]

fitqr2<-rq(NGM_D ~ NGM_C,data=qbacmicq,tau=c(0.10,0.90))
bgli2<-coefficients(fitqr2)[1,][[1]]
bgui2<-coefficients(fitqr2)[1,][[2]]
bgls2<-coefficients(fitqr2)[2,][[1]]
bgus2<-coefficients(fitqr2)[2,][[2]]


q05<-quantile(qbacmicq$MIC,0.05,na.rm=TRUE)[[1]]
q10<-quantile(qbacmicq$MIC,0.1,na.rm=TRUE)[[1]]
q90<-quantile(qbacmicq$MIC,0.9,na.rm=TRUE)[[1]]
q95<-quantile(qbacmicq$MIC,0.95,na.rm=TRUE)[[1]]
q99<-quantile(qbacmicq$MIC,0.99,na.rm=TRUE)[[1]]



bacres95<-subset(qbacmicq,NGM_D>NGM_C*bgus+bgui)$Gene
bacsens05<-subset(qbacmicq,NGM_D<NGM_C*bgls+bgli)$Gene

bacres90<-subset(qbacmicq,NGM_D>NGM_C*bgus2+bgui2)$Gene
bacsens10<-subset(qbacmicq,NGM_D<NGM_C*bgls2+bgli2)$Gene

wormres95<-subset(qbacmicq,MIC>q95)$Gene
wormsens05<-subset(qbacmicq,MIC<q05)$Gene

wormres90<-subset(qbacmicq,MIC>q90)$Gene
wormres25<-subset(qbacmicq,MIC>2.5)$Gene
wormresall<-subset(qbacmicq,MIC>1)$Gene
wormresall5<-subset(qbacmicq,MIC>5)$Gene
wormsens10<-subset(qbacmicq,MIC<q10)$Gene

resmic5=list('Bacteria sensitive bottom 5%'=bacsens05,
            'Bacteria resistant top 5%'=bacres95,
            'Worms resistant top 5%'=wormres95,
            'Worms sensitive bottom 5%'=wormsens05)

resmic10=list('Worms resistant top 10%'=wormres90,
              'Bacteria sensitive bottom 10%'=bacsens10,
              'Bacteria resistant top 10%'=bacres90) # ,'Worms sensitive bottom 10%'=wormsens10 ,'Worms resistant all'=wormresall5

resmic_df<-plyr::ldply(resmic10, cbind)
resmic_df<-rename(resmic_df,c('.id'='List','1'='Gene'))

#Summary
rmc<-subset(qbacmicq, Gene %in% unique(resmic_df$Gene))
rmc$'Bacteria resistant top 5%'<-ifelse(rmc$Gene %in% subset(resmic_df,List=="Bacteria resistant top 5%")$Gene,TRUE,FALSE)
rmc$'Bacteria sensitive bottom 5%'<-ifelse(rmc$Gene %in% subset(resmic_df,List=="Bacteria sensitive bottom 5%")$Gene,TRUE,FALSE)
rmc$'Worms resistant top 5%'<-ifelse(rmc$Gene %in% subset(resmic_df,List=="Worms resistant top 5%")$Gene,TRUE,FALSE)
rmc$'Worms sensitive bottom 5%'<-ifelse(rmc$Gene %in% subset(resmic_df,List=="Worms sensitive bottome 5%%")$Gene,TRUE,FALSE)
rmc<-rmc[,c(1,19:21,2:7,9:18,8)]
write.csv(rmc,paste(ddir,'/Venn_Worm-Bacteria_resistance_overlap.csv',sep=''))



plot(Venn(resmic10), doWeights = FALSE)#,type='ellipses'  ,type='ellipses' ,type='ellipses'
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Venn_Worm-Bacteria_resistance_overlap.pdf",sep=''),
             width=9,height=9)

