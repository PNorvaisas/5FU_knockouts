library('ggplot2')

ddir<-'Data_v2'
data<-read.table(paste(ddir,'/MICs_and_bacterial_growth-Clean.csv',sep=''),sep=',',header = TRUE,stringsAsFactors = FALSE)

sel<-data[,c('Gene','MIC')]
sel$Gene<-toupper(sel$Gene)
sel<-subset(sel,!is.na(MIC))

sel$MIC<-jitter(sel$MIC,amount=0.01)

sel$lMIC<-log(sel$MIC,2)
sel$MICs<-sel$MIC/100
sel$pMIC<-sel$MIC
sel[sel$pMIC==1,'pMIC']<--1



ggplot(sel,aes(x=MIC))+geom_histogram(binwidth=1)

seld<-sel

seld$Gene<-paste(seld$Gene,'D',sep='')
seld$MIC<--seld$MIC
seld$lMIC<--seld$lMIC

doub<-merge(sel,seld,all.x = TRUE,all.y='TRUE')



#sel2<-subset(sel,MIC>1)
write.table(sel[,c('Gene','MIC')],paste(ddir,'/Gene_sets/MICS_GSE.rnk',sep=''),sep='\t',row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(sel[,c('Gene','MICs')],paste(ddir,'/Gene_sets/MICS_GSE_scaled.rnk',sep=''),sep='\t',row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(sel[,c('Gene','lMIC')],paste(ddir,'/Gene_sets/MICS_GSE-l.rnk',sep=''),sep='\t',row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(sel[,c('Gene','pMIC')],paste(ddir,'/Gene_sets/MICS_GSE_phenotype.rnk',sep=''),sep='\t',row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(subset(sel,lMIC!=0)[,c('Gene','lMIC')],paste(ddir,'/Gene_sets/MICS_GSE-l_over1.rnk',sep=''),sep='\t',row.names = FALSE,col.names = FALSE,quote = FALSE)

#write.table(sel2,'Data/MICS_GSE2.rnk',sep='\t',row.names = FALSE,col.names = FALSE,quote = FALSE)

write.table(doub[,c('Gene','MIC')],paste(ddir,'/Gene_sets/MICS_GSE_dbl.rnk',sep=''),sep='\t',row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(doub[,c('Gene','lMIC')],paste(ddir,'/Gene_sets/MICS_GSE-l_dbl.rnk',sep=''),sep='\t',row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(subset(doub,lMIC!=0)[,c('Gene','lMIC')],paste(ddir,'/Gene_sets/MICS_GSE-l_dbl_over1.rnk',sep=''),sep='\t',row.names = FALSE,col.names = FALSE,quote = FALSE)


