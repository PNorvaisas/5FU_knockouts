library('ggplot2')
library('gplots')
library('plyr')
library('reshape2')
library(tidyr)
library(BSDA)
library(xlsx)
#library(quantreg)
#library(ellipse)


#Output folder:
odir<-'Figures_final'
ddir<-'Data_final'


#Read Keio info table with identifiers
keioinfo<-read.table('../Keio_library/Keio_library_fully_annotated.csv',sep=',',quote = '"',header = TRUE,stringsAsFactors=FALSE)
keioinfo$X<-NULL
keioinfo<-subset(keioinfo,!Plate %in% c('91','93','95'))
rdupl<-as.factor(unique(keioinfo[which(duplicated(keioinfo$Gene)),]$Gene))
allgenes<-subset(keioinfo,!is.na(Gene) & !Gene %in% c('present','WT'))$Gene

nodata<-read.table(paste(ddir,'/Worms_nodata.csv',sep=''),sep=',',quote = '"',header = TRUE,stringsAsFactors=FALSE)



mics<-read.table(paste(ddir,'/Worm_MICs.csv',sep=''),sep=',',quote = '"',header = TRUE,stringsAsFactors=FALSE)
colnames(mics)<-gsub('X','',colnames(mics))


bacall<-read.table(paste(ddir,'/Bacterial_growth_summary_with_yjjG.csv',sep=''),sep=',',quote = '"',header = TRUE,stringsAsFactors=FALSE)

bacall<-read.table(paste(ddir,'/Bacterial_growth_summary.csv',sep=''),sep=',',quote = '"',header = TRUE,stringsAsFactors=FALSE)



#poorgrowth<-unique(subset(bacall,C_OD_Mean<0.05)$Gene)
poorgrowth<-c('atpB','atpE','atpF','atpG','atpE','lpd','sucA','atpE','ydcT','glnA')


#lowmicbac<-c('dcp','fre','glgC','oxyR','ybdR','yecE','yehA','yghT','ytjB')

bacmict<-merge(mics,bacall,id=c('Gene','Plate','Well'),all.x = TRUE)



#Remove
keioundis<-read.table('../Keio_library/Keio_undisrupted.csv',sep=',',quote = '"',header = TRUE,stringsAsFactors=FALSE)
rmdev<-read.table('Keio_development_delays.csv',sep=',',quote = '"',header = TRUE,stringsAsFactors=FALSE)




#Remove undisrupted and poor growing strains
#Make sure you can use Gene as unique identifier

#Remove knockouts that cause developmental delays


#Need rdupl
#Need nodata
#Need rdupl

intersect(rdupl,rmdev$Gene)
intersect(rdupl,nodata$Gene)
intersect(rdupl,keioundis$Gene)
intersect(rdupl,poorgrowth)

bacmic<-subset(bacmict,!Gene %in% keioundis$Gene & !Gene %in% poorgrowth )
bacmic<-subset(bacmic,!Gene %in% setdiff(rmdev$Gene,'wbbL'))
bacmic<-bacmic[!(bacmic$Gene=='wbbL' & bacmic$Plate=='21' & bacmic$Well=='A7'), ]



#Data fully merged!!

poorgr<-subset(bacmict,Gene %in% setdiff(poorgrowth,'Blank'))
poorgr$NComment<-'Eliminated; no growth'

devdel<-subset(bacmict,Gene %in% setdiff(rmdev$Gene,'wbbL') | (Gene=='wbbL' & Plate=='21' & Well=='A7'))
devdel$NComment<-'Causes developmental delay'


undisrupted<-subset(bacmict,Gene %in% keioundis$Gene)
undisrupted$NComment<-'Undisrupted gene, possibly essential (Yamamoto2009)'

#Undisrupted and not in data
unmis<-setdiff(keioundis$Gene,undisrupted$Gene)

unno<-subset(nodata,Gene %in% keioundis$Gene)
unno$NComment<-'Undisrupted gene'

otherno<-subset(nodata,!Gene %in% keioundis$Gene)
otherno$NComment<-'Eliminated; no growth'

remsett<-merge(undisrupted,unno,all.x=TRUE,all.y=TRUE)
remsetf<-merge(remsett,otherno,all.x=TRUE,all.y=TRUE)
remsetff<-merge(remsetf,devdel,all.x=TRUE,all.y=TRUE)
remsetfff<-merge(remsetff,poorgr,all.x=TRUE,all.y=TRUE)
remsetf<-remsetfff[,!colnames(remsetfff) %in% c('Row','Column','Comment')]
remsetf<-rename(remsetf,c("NComment"="Comment"))

remset<-remsetf[,c('Gene','Plate','Well','JW_id','ECK','bno',
                   'LB_22hr','MOPS_24hr','MOPS_48hr',
                   'C_OD_Mean','C_OD_SD','T_OD_Mean','T_OD_SD','Comment')]


dim(remset)


#All unique gens and one instance of WT
length(allgenes)+1
#Genes in screen
length(bacmic$Gene)
#Genes removed
length(remset$Gene)



#Do we have everything?
length(allgenes)+1==length(bacmic$Gene)+length(remset$Gene)




#Sanity checks:
#Check whether the removal is OK
poorgrowth



micexpl<-read.table(paste(ddir,'/MICs_column_explanation.csv',sep=''),
                    sep=',',quote = '"',header = TRUE,stringsAsFactors=FALSE)


bacmic<-bacmic[,setdiff(colnames(bacmic),c('0','1','2.5','5','UniACC'))] #,'EG','GI'
colnames(bacmic)
bacmicw<-bacmic[,c(1:5,8:10,6:7,11:length(colnames(bacmic)))]

colnames(bacmicw)

explrd<-subset(micexpl,Column %in% colnames(bacmicw))
explrd<-explrd[match(colnames(bacmicw),explrd$Column),]


write.csv(bacmicw,paste(ddir,'/MICs_and_bacterial_growth-Complete.csv',sep=''),row.names = FALSE)
write.xlsx2(explrd, file=paste(ddir,'/MICs_and_bacterial_growth-Complete.xlsx',sep=''), sheetName="Readme",row.names = FALSE,showNA=FALSE)
write.xlsx2(bacmicw, file=paste(ddir,'/MICs_and_bacterial_growth-Complete.xlsx',sep=''), sheetName="Data", append=TRUE,row.names = FALSE,showNA=FALSE)

write.xlsx2(explrd, file='/Users/Povilas/Projects/B-D-H paper/figures and data/figure 2/final files/Table S2.xlsx',
            sheetName="Readme_All-data",row.names = FALSE,showNA=FALSE,append = TRUE)
write.xlsx2(bacmicw, file='/Users/Povilas/Projects/B-D-H paper/figures and data/figure 2/final files/Table S2.xlsx',
            sheetName="All-data", append=TRUE,row.names = FALSE,showNA=FALSE)


#Knockouts removed from screen

remsetw<-remset[,! colnames(remset) %in% c('C_OD_Mean','C_OD_SD','T_OD_Mean','T_OD_SD')]

explrm<-subset(micexpl,Column %in% colnames(remsetw))
explrm<-explrm[match(colnames(remsetw),explrm$Column),]
write.csv(remsetw,paste(ddir,'/All_removed_from_screen.csv',sep=''),row.names = FALSE)
write.xlsx2(explrm, file=paste(ddir,'/All_removed_from_screen.xlsx',sep=''), sheetName="Readme",row.names = FALSE,showNA=FALSE)
write.xlsx2(remsetw, file=paste(ddir,'/All_removed_from_screen.xlsx',sep=''), sheetName="Data", append=TRUE,row.names = FALSE,showNA=FALSE)


write.xlsx2(explrm, file='/Users/Povilas/Projects/B-D-H paper/figures and data/figure 2/final files/Table S2.xlsx',
            sheetName="Readme_Removed",row.names = FALSE,showNA=FALSE,append = TRUE)
write.xlsx2(remsetw, file='/Users/Povilas/Projects/B-D-H paper/figures and data/figure 2/final files/Table S2.xlsx',
            sheetName="Removed", append=TRUE,row.names = FALSE,showNA=FALSE)



#Duplicates
duplkeio<-subset(keioinfo,Gene %in% rdupl & !Gene %in% c('WT',NA,'present'))
duplset<-duplkeio[order(duplkeio$Gene),]
write.csv(duplset,paste(ddir,'/Gene_duplicates_in_library.csv',sep=''),row.names = FALSE)
duplscreen<-subset(bacmic,Gene %in% rdupl & !Gene %in% c('WT',NA,'present'))
duplscreen<-duplscreen[order(duplscreen$Gene),]
write.csv(duplscreen,paste(ddir,'/Gene_duplicates_in_screen.csv',sep=''),row.names = FALSE)
#write.csv(duplset,paste(ddir,'/Real_duplicates.csv',sep=''))

