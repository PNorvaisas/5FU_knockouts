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


lmsum<-function(m){fres<-summary(m)
l <- list(b = as.double(coef(m)[1], digits = 2),
          a = as.double(coef(m)[2], digits = 2),
          r2 = as.double(summary(m)$r.squared, digits = 3),
          p2 = as.double(pf(fres$fstatistic[1], fres$fstatistic[2], fres$fstatistic[3],lower.tail = FALSE)[[1]], digits = 5))
return(l)
}

trim <- function (x) gsub("^\\s+|\\s+$", "", x)



evalmic=function(m){
  print(m[1][1])
  m<-as.numeric(m)
  mic=1
  if (!is.na(m[[3]]) & m[[3]]>0) {
    mic=2.5
  }
  if (!is.na(m[[4]]) & m[[4]]>0) {
    mic=5
  }
  if (!is.na(m[[5]]) & m[[5]]>0) {
    mic=10
  }
  return(mic)
}

evalmic2=function(all){
  for (i in rownames(all)){
    mic=1
    if (!is.na(all[i,'1']) & all[i,'1']>0) {
      mic=2.5
    }
    if (!is.na(all[i,'2.5']) & all[i,'2.5']>0) {
      mic=5
    }
    if (!is.na(all[i,'5']) & all[i,'5']>0) {
      mic=10
    }
    all[i,'MIC']<-mic
  }
  return(all)
}



#Read Keio info table with identifiers
keioinfo<-read.table('../Keio_library/Keio_library_fully_annotated.csv',sep=',',quote = '"',header = TRUE,stringsAsFactors=FALSE)
keioinfo$X<-NULL
keioinfo<-subset(keioinfo,!Plate %in% c('91','93','95'))





#Get data that's already scored
scr1r<-read.table('Primary_screen_PN_clean_fixed.csv',sep=',',quote = '"',header = TRUE,stringsAsFactors=FALSE)

#scr1r<-subset(scr1r, ! Gene %in% c('XXXXXXX','no bact','empty','',NA))
scr1r$Plate<-scr1r$Keio.Plate.no.
scr1r$Well<-scr1r$Position
scr1r$Keio.Plate.no.<-NULL
scr1r$Position<-NULL
scr1r$Faults<-NULL
scr1r$Details<-NULL
scr1r[scr1r$Plate %in% c('-','',' '),'Plate']<-NA
scr1r[scr1r$Well %in% c('-','',' '),'Well']<-NA
scr1r[scr1r$Gene %in% c('WT?','WT control', 'dodgy "WT"'),'Gene']<-'WT'
scr1r[scr1r$Gene %in% c('XXXXXXX','no bact','empty',''),'Gene']<-NA
#scr1r[scr1r$Gene=='WT',c('Keio.Plate.no.','Position')]<-NA



##Apply bandage
timfix<-read.csv('Tim_fixed_scr1.csv',sep=',',quote = '"',header = TRUE,stringsAsFactors=FALSE)
timfix$Action<-NULL
timfix$Tim.s.comments<-NULL
timfix$Index<-apply(timfix[,c('Plate','Well')],1,paste,collapse='-')
scr1r$Index<-apply(scr1r[,c('Plate','Well')],1,paste,collapse='-')

scr1r[match(timfix$Index,scr1r$Index),'Gene']<-timfix$Gene_fix
scr1r$Index<-NULL

scr1<-rename(scr1r, c("X0"="0", "X1"="1", "X2.5"="2.5", "X5"="5"))
scr1m<-melt(scr1[,c('Gene','Plate','Well','0','1','2.5','5')],id=c('Gene','Plate','Well'),variable.name = 'Measure',value.name='Score')
scr1m$Measure<-as.character(scr1m$Measure)
scr1m$Score<-as.character(scr1m$Score)

#Secondary screen
scr2<-read.table('Secondary_screen_PN_new.csv',sep=',',quote = '"',header = TRUE,stringsAsFactors = FALSE)
scr2$Starving<-as.factor(scr2$Starving)
#Fix Gene-location relationships. Well must come first!
scr2[scr2$Gene=='yedN','Plate']<-'31'
scr2[scr2$Gene=='yedN','Well']<-'G12'
scr2[scr2$Gene=='dcuC' & scr2$Plate=='4','Well']<-'H3'
scr2[scr2$Gene=='dcuC' & scr2$Plate=='4','Plate']<-'89'
scr2[scr2$Gene=='dcuC' & scr2$Plate=='2','Well']<-'A4'
scr2[scr2$Gene=='dcuC' & scr2$Plate=='2','Plate']<-'89'

#scr2<-scr2[,c('Gene','Plate','Well','MIC1','MIC2','MIC3')]
scr2m<-melt(scr2[,c('Gene','Plate','Well','MIC1','MIC2','MIC3')],id=c('Gene','Plate','Well'),variable.name = 'Replicate',value.name='MIC')
scr2m$Replicate<-as.character(scr2m$Replicate)
scr2m$MIC<-as.character(scr2m$MIC)

scr2m<-subset(scr2m,!Gene %in% c('no bacteria','empty'))


#Supplemented screen 3
scr3fix<-read.table('3rd_screen_location_fix.csv',sep=',',quote = '"',header = TRUE,stringsAsFactors = FALSE)
scr3Scores<-read.table('3rd_screen_Scores.csv',sep=',',quote = '"',header = TRUE,stringsAsFactors = FALSE)
scr3Scores<-rename(scr3Scores, c("Keio.Plate.no."="Plate", "Position"="Well", "Well"="NWell"))

#Fix mistakes with gene-location mismatches
#Tim picked knockouts by provided plate-well locations, therefore, names have to be updated!
scr3Scores<-merge(scr3Scores,scr3fix[,c('Plate','Well','Right_gene')],by=c('Plate','Well'),all.x=TRUE)
scr3Scores[!is.na(scr3Scores$Right_gene),c('Gene')]<-scr3Scores[!is.na(scr3Scores$Right_gene),c('Right_gene')]
scr3Scores$Right_gene<-NULL


scr3Scores[scr3Scores$Gene=='upp cont','Plate']<-'61'
scr3Scores[scr3Scores$Gene=='upp cont','Well']<-'B1'
scr3Sm<-melt(scr3Scores[,c('Gene','Plate','Well',
                           "X0.1","X0.2","X0.3","X1.1","X1.2","X1.3",
                           "X2.5.1","X2.5.2","X2.5.3","X5.1","X5.2","X5.3")],
             id=c('Gene','Plate','Well'),variable.name = 'Measure',value.name='Score')
scr3Sm$Measure<-as.character(scr3Sm$Measure)
scr3Sm$Score<-as.character(scr3Sm$Score)
scr3Sm[scr3Sm$Measure %in% c('X0.1','X0.2','X0.3'),'Measure']<-'0'
scr3Sm[scr3Sm$Measure %in% c('X1.1','X1.2','X1.3'),'Measure']<-'1'
scr3Sm[scr3Sm$Measure %in% c('X2.5.1','X2.5.2','X2.5.3'),'Measure']<-'2.5'
scr3Sm[scr3Sm$Measure %in% c('X5.1','X5.2','X5.3'),'Measure']<-'5'
scr3Sm[scr3Sm$Score %in% c('3'),'Score']<-'9'
scr3Sm[scr3Sm$Score %in% c('2'),'Score']<-'6'
scr3Sm[scr3Sm$Score %in% c('1'),'Score']<-'3'



#3rd screen MIC values all
scr3Mics<-read.table('3rd_screen_MICs.csv',sep=',',quote = '"',header = TRUE,stringsAsFactors = FALSE)
scr3Mics<-rename(scr3Mics, c("Keio.Plate.no."="Plate", "Position"="Well", "Well"="NWell"))
scr3Mics<-merge(scr3Mics,scr3fix[,c('Plate','Well','Right_gene')],by=c('Plate','Well'),all.x=TRUE)
#Fix mistakes with gene-location missmatches
#Tim picked knockouts by provided plate-well locations, therefore, names have to be updated!
scr3Mics[!is.na(scr3Mics$Right_gene),c('Gene')]<-scr3Mics[!is.na(scr3Mics$Right_gene),c('Right_gene')]
scr3Mics$Right_gene<-NULL
scr3Micsm<-melt(scr3Mics[,c('Gene','Plate','Well','MIC.1','MIC.2','MIC.3')],id=c('Gene','Plate','Well'),variable.name = 'Replicate',value.name='MIC')
scr3Micsm$Replicate<-as.character(scr3Micsm$Replicate)
scr3Micsm$MIC<-as.character(scr3Micsm$MIC)

#Part not really used
# scr3<-subset(scr3Micsm,!(is.na(as.numeric(scr3Micsm$MIC))))
# scr3[scr3$Gene=='WT cont','Gene']<-'WT'
# scr3[scr3$Gene=='upp cont','Gene']<-'upp'
# scr3$MIC<-as.numeric(scr3$MIC)
# scr3avg<-ddply(scr3, .(Gene,Plate,Well), summarise, MIC_avg=mean(MIC,na.rm = TRUE),MIC_SD=sd(MIC,na.rm = TRUE))


sc1g<-unique(scr1$Gene)
sc2g<-unique(scr2$Gene)
sc3Sg<-unique(scr3Scores$Gene)
sc3Mg<-unique(scr3Mics$Gene)


allscores<-merge(scr1m,scr3Sm,all.x=TRUE,all.y=TRUE)
allmics<-merge(scr2m,scr3Micsm,all.x=TRUE,all.y=TRUE)

#Clean values
allscores<-subset(allscores,! is.na(Gene) & !Gene %in% c('','XXXXX') & ! Score %in% c('','no worms','2, contamination','XXXX'))
allmics<-subset(allmics,! is.na(Gene) & !Gene %in% c('','XXXXX','empty','no bacteria') & !MIC %in% c('empty','no growth','XXXX'))


#Name fixes
allscores[allscores$Gene=='WT cont','Gene']<-'WT'
allscores[allscores$Gene=='upp cont','Gene']<-'upp'
allmics[allmics$Gene=='WT cont','Gene']<-'WT'
allmics[allmics$Gene=='upp cont','Gene']<-'upp'

#Location fixes
allscores[allscores$Gene=='WT','Plate']<-''
allscores[allscores$Gene=='WT','Well']<-''
allscores[allscores$Gene=='ubiX','Plate']<-'91,93,95 comb'
allscores[allscores$Gene=='ubiX','Well']<-'D5'
allmics[allmics$Gene=='WT','Plate']<-''
allmics[allmics$Gene=='WT','Well']<-''




#Real duplicates!
rdupl<-as.factor(unique(keioinfo[which(duplicated(keioinfo$Gene)),]$Gene))

#Get Gene names with an exception of duplicates!
allmics<-merge(allmics,subset(keioinfo,!Gene %in% rdupl)[,c('Gene','Plate','Well')],by='Gene',all.x=TRUE)

allmics[!is.na(allmics$Plate.y) &!is.na(allmics$Well.y),c('Plate.x','Well.x')]<-allmics[!is.na(allmics$Plate.y) &!is.na(allmics$Well.y),c('Plate.y','Well.y')]
allmics<-rename(allmics,c('Plate.x'='Plate','Well.x'='Well'))
allmics$Plate.y<-NULL
allmics$Well.y<-NULL
#Real duplicates in mics
rduplm<-c('dcuC','yedN','yhcE')

#rdupl<-as.factor(unique(allmics[which(duplicated(allmics[,c('')])),]$Gene))


scores<-subset(allscores,! is.na(Score))
scores$Score<-as.numeric(scores$Score)
unique(scores$Score)

micss<-subset(allmics,! is.na(MIC))
micss$MIC<-as.numeric(micss$MIC)
unique(micss$MIC)

#Calculate averages for Scores and MICs over replicates
scores_avg<-ddply(scores, .(Gene,Plate,Well,Measure), summarise, Score_avg=mean(Score,na.rm = TRUE),Score_sd=sd(Score,na.rm = TRUE))#,Score_N=length(Well)

#Averaging by gene names and positions
mics_avg<-ddply(micss, .(Gene,Plate,Well), summarise, MIC_avg=mean(MIC,na.rm = TRUE),MIC_SD=sd(MIC,na.rm = TRUE))#,MIC_N=length(Well)


alls<-dcast(scores_avg,Gene+Plate+Well ~Measure,mean,value.var = c('Score_avg'),fill = as.numeric(NA))

#Evaluate MICS
alls$MIC<-evalmic2(alls)$MIC


#Consistency checks
#duplm should be dcuC and yhcE

#Scr1 length+21 duplicates, +3 triplicates
gfreq<-sort(table(keioinfo$Gene))
length(unique(subset(scr1,!Gene %in% c(NA,'WT'))$Gene))+21+3

#Scr2 length + duplicate
length(unique(scr2m$Gene))+1


#Scr3 Scores length
scr3genes<-subset(scr3Scores,!Gene %in% c('XXXXX','WT cont'))$Gene
intersect(scr3genes,rdupl)
length(unique(scr3genes))+6

#Scr3 MICs length
scr3mgenes<-subset(scr3Mics,!Gene %in% c('XXXXX','WT cont','upp cont',''))$Gene
intersect(scr3mgenes,rdupl)
length(unique(scr3mgenes))+6


duplm<-as.factor(unique(mics_avg[which(duplicated(mics_avg$Gene)),]$Gene))
dupls<-as.factor(unique(alls[which(duplicated(alls$Gene)),]$Gene))

length(subset(alls,Gene!='WT')$Gene)
length(subset(mics_avg,Gene!='WT')$Gene)
length(subset(keioinfo,!is.na(Gene) & !Gene %in% c('present','WT'))$Gene)
allgenes<-subset(keioinfo,!is.na(Gene) & !Gene %in% c('present','WT'))$Gene

missg<-setdiff(allgenes,alls$Gene)
conc<-merge(alls,keioinfo,by=c('Gene','Plate','Well'),all.x=TRUE,all.y=TRUE)

#15 excluded due to prior knowledge of poor growth
nodata<-subset(conc,is.na(MIC) & !Gene %in% c('present','WT',NA))
write.csv(nodata,paste(ddir,'/Worms_nodata.csv',sep=''),row.names = FALSE)


#Needs Keio growth data
#nogrowth<-merge(nodata, keio,by=c('Gene'))
#write.csv(nogrowth,'Data_final/Excluded.csv')

#Merging by 'Plate','Well' and Gene
allfull<-merge(alls,mics_avg,by=c('Gene','Plate','Well'),all.x=TRUE,all.y=TRUE)

allfull$MIC<-ifelse(!is.na(allfull$MIC_avg) & (allfull$MIC_avg > allfull$MIC & allfull$MIC %in% c(5,10)) ,allfull$MIC_avg,allfull$MIC)
allfull$MIC_avg<-NULL

#Manual fixes
allfull[allfull$Gene=='WT','MIC']<-1
allfull[allfull$Gene=='WT','MIC_SD']<-0

allfull[allfull$Gene=='yjjG','MIC']<-5
allfull[allfull$Gene=='yjjG','5']<-0
allfull[allfull$Gene=='yjjG','MIC_SD']<-NA

# allfull[allfull$Gene=='upp','MIC']<-15
# allfull[allfull$Gene=='upp','MIC_SD']<-0


#Find duplicates
dupl<-as.factor(unique(allfull[which(duplicated(allfull$Gene)),]$Gene))


#We have all necessary MIC values!
allfull$HasNA<-ifelse(is.na(allfull$`0`)|is.na(allfull$`1`)|is.na(allfull$`2.5`)|is.na(allfull$`5`),TRUE,FALSE)

mismic<-subset(allfull,HasNA & is.na(MIC_SD))
length(subset(mismic,MIC>2.5)$Gene)



# 
# #Find duplication problems
#Check for duplicates
# dad<-table(allfull$Gene)
# ind<-table(keioinfo$Gene)
# match<-merge(data.frame(dad),data.frame(ind),by='Var1')
# match<-rename(match, c("Var1"="Gene", "Freq.x"="Data_freq", "Freq.y"="Info_freq"))
# strange<-match[match$Data_freq!=match$Info_freq & match$Gene!='WT',]$Gene
# 
# problems<-merge(subset(allfull,Gene %in% strange),
#                 subset(keioinfo,Gene %in% strange),
#                 by=c('Gene'),all.x=TRUE,all.y=TRUE)
# 
# ps<-subset(problems,Gene %in% strange)
# ps<-merge(ps,match,by='Gene',all.x=TRUE)
# ps<-merge(ps,keioinfo[,c('Gene','Plate','Well')],by.x=c('Plate.x','Well.x'),by.y=c('Plate','Well'),all.x=TRUE)
# write.csv(ps,'Data/Mismatch_check.csv')
# #All missmatches between data and Keio info must come from plates 91,93,95


allinfr<-merge(allfull,keioinfo,by=c('Gene','Plate','Well'),all.x=TRUE)
allinf<-allinfr[,! colnames(allinfr) %in% c('Gene.y','Row','Column','Comment','HasNA')]
#allinf<-rename(allinfr, c("Gene.x"="Gene"))


#Merge with Keio reference growth

keio<-read.table('Keio_growth.csv',sep=',',quote = '"',header = TRUE,stringsAsFactors=FALSE)
#keio<-keio[,colnames(keio) %in% c('Gene','LB_22hr','MOPS_24hr','MOPS_48hr')]
keio[is.na(keio)]<-NA
#keio<-subset(keio,!LB_22hr=='N.A.' & !MOPS_24hr=='N.A.' & !MOPS_48hr=='N.A.')
keio$LB_22hr<-as.numeric(as.character(keio$LB_22hr))
keio$MOPS_24hr<-as.numeric(as.character(keio$MOPS_24hr))
keio$MOPS_48hr<-as.numeric(as.character(keio$MOPS_48hr))

kdupl<-as.factor(unique(keio[which(duplicated(keio$Gene)),]$Gene))
kduplicates<-subset(keio,Gene %in% kdupl & Gene!='none')

#Merge with Keio library data
mics<-merge(allinf,keio,by.x=c('JW_id','ECK','Gene'),by.y=c('JW.id','ECK.number','Gene'),all.x=TRUE)
mics[is.na(mics$LB_22hr),c('LB_22hr','MOPS_24hr','MOPS_48hr')]<-keio[match(subset(mics,is.na(LB_22hr))$Gene,keio$Gene),c('LB_22hr','MOPS_24hr','MOPS_48hr')]



write.csv(mics,paste(ddir,'/Worm_MICs.csv',sep=''),row.names = FALSE)


