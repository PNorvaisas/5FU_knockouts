library('ggplot2')

data<-read.table('Data/MICs_and_bacterial_growth-Unique_clean.csv',sep=',',header = TRUE,stringsAsFactors = FALSE)

sel<-data[,c('Gene','MIC')]
sel$Gene<-toupper(sel$Gene)
sel<-subset(sel,!is.na(MIC))

sel$lMIC<-log(sel$MIC,2)



ggplot(sel,aes(x=lMIC))+geom_histogram(binwidth=1)

seld<-sel

seld$Gene<-paste(seld$Gene,'D',sep='')
seld$MIC<--seld$MIC
seld$lMIC<--seld$lMIC

doub<-merge(sel,seld,all.x = TRUE,all.y='TRUE')



#sel2<-subset(sel,MIC>1)
write.table(sel[,c('Gene','MIC')],'Data/MICS_GSE.rnk',sep='\t',row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(sel[,c('Gene','lMIC')],'Data/MICS_GSE-l.rnk',sep='\t',row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(subset(sel,lMIC!=0)[,c('Gene','lMIC')],'Data/MICS_GSE-l_over1.rnk',sep='\t',row.names = FALSE,col.names = FALSE,quote = FALSE)

#write.table(sel2,'Data/MICS_GSE2.rnk',sep='\t',row.names = FALSE,col.names = FALSE,quote = FALSE)

write.table(doub[,c('Gene','MIC')],'Data/MICS_GSE_dbl.rnk',sep='\t',row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(doub[,c('Gene','lMIC')],'Data/MICS_GSE-l_dbl.rnk',sep='\t',row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(subset(doub,lMIC!=0)[,c('Gene','lMIC')],'Data/MICS_GSE-l_dbl_over1.rnk',sep='\t',row.names = FALSE,col.names = FALSE,quote = FALSE)
