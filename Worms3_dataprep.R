library('ggplot2')
library('gplots')
library('plyr')
library('reshape2')
library(tidyr)
#library(quantreg)
#library(ellipse)


#Output folder:
odir<-'Figures_final'
ddir<-'Data_final'





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

keioundis<-read.table('../Keio_library/Keio_undisrupted.csv',sep=',',quote = '"',header = TRUE,stringsAsFactors=FALSE)






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
# scr3avg<-ddply(scr3, .(Gene,Plate,Well), summarise, MIC_avg=mean(MIC,na.rm = TRUE),MIC_sd=sd(MIC,na.rm = TRUE))


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
mics_avg<-ddply(micss, .(Gene,Plate,Well), summarise, MIC_avg=mean(MIC,na.rm = TRUE),MIC_sd=sd(MIC,na.rm = TRUE))#,MIC_N=length(Well)


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

#Needs Keio growth data
#nogrowth<-merge(nodata, keio,by=c('Gene'))
#write.csv(nogrowth,'Data_final/Excluded.csv')

#Merging by 'Plate','Well' and Gene
allfull<-merge(alls,mics_avg,by=c('Gene','Plate','Well'),all.x=TRUE,all.y=TRUE)

allfull$MIC<-ifelse(!is.na(allfull$MIC_avg) & (allfull$MIC_avg > allfull$MIC & allfull$MIC %in% c(5,10)) ,allfull$MIC_avg,allfull$MIC)
allfull$MIC_avg<-NULL

#Manual fixes
allfull[allfull$Gene=='WT','MIC']<-1
allfull[allfull$Gene=='WT','MIC_sd']<-0

# allfull[allfull$Gene=='upp','MIC']<-15
# allfull[allfull$Gene=='upp','MIC_sd']<-0


#Find duplicates
dupl<-as.factor(unique(allfull[which(duplicated(allfull$Gene)),]$Gene))


#We have all necessary MIC values!
allfull$HasNA<-ifelse(is.na(allfull$`0`)|is.na(allfull$`1`)|is.na(allfull$`2.5`)|is.na(allfull$`5`),TRUE,FALSE)

mismic<-subset(allfull,HasNA & is.na(MIC_sd))
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



#Bacterial growth OLD
#NGM media start OD=0.057392708
#Data from 24hr growth
#In ideal case get location map from 2nd screen and merge by location-gene
# bac<-read.table('Bacteria.csv',sep=',',quote = '"',header = TRUE,stringsAsFactors=FALSE)
# bac[bac$Gene=='BW',]$Gene<-'WT'
# bac<-rename(bac, c("Plate"="Nplate", "Well"="NWell"))
# bac<-subset(bac,!Nplate=='')
# 
# #Get Gene positions for not duplicated genes
# bact<-merge(bac,subset(keioinfo,!Gene %in% c('dcuC','yedN','WT'))[,c('Gene','Plate','Well')],by='Gene',all.x=TRUE)
# 
# #Fix values for duplicates
# bact[bact$Gene=='yedN','Plate']<-'31'
# bact[bact$Gene=='yedN','Well']<-'G12'
# bact[bact$Gene=='dcuC' & bact$Nplate=='4','Well']<-'H3'
# bact[bact$Gene=='dcuC' & bact$Nplate=='4','Plate']<-'89'
# bact[bact$Gene=='dcuC' & bact$Nplate=='2','Well']<-'A4'
# bact[bact$Gene=='dcuC' & bact$Nplate=='2','Plate']<-'89'
# bact[bact$Gene=='WT' & bact$Nplate=='C','Well']<-''
# bact[bact$Gene=='WT' & bact$Nplate=='C','Plate']<-''
# 




###Get annotations for bacterial plates
# 
# pl14<-bact[,c('Gene','Nplate','NWell','Plate','Well')]
# pl14<-pl14[!duplicated(pl14), ]
# pl58<-scr3Mics[,c('Gene','Nplate','NWell','Plate','Well')]
# pl58<-subset(pl58,!is.na(Nplate))
# pl58[pl58$Gene=='WT cont','Gene']<-'WT'
# pl58[pl58$Gene=='upp cont','Gene']<-'upp'
# pl58[pl58$Gene %in% c('upp','hyfF','hyfB','hyfR','hyfH','purU','ubiX'),c('Plate','Well')]<-keioinfo[match(pl58[pl58$Gene %in% c('upp','hyfF','hyfB','hyfR','hyfH','purU','ubiX'),c('Gene')],keioinfo$Gene),c('Plate','Well')]
# #Yet another Timfix
# pl58[pl58$Nplate==8 & pl58$NWell %in% c('D1','E1','F1','D12','E12','G1','H1','A2'),'Gene']<-'Blank'
# 
# bac18<-merge(pl14,pl58,all.x=TRUE,all.y=TRUE)
# bac18[bac18$Gene %in% c('XXXXX','XXXXXXX','empty','no bacteria'),'Gene']<-'Blank'
# bac18[bac18$Gene=='Blank','Well']<-''
# bac18[bac18$Gene=='Blank','Plate']<-''
# #bactest<-merge(bac18,keioinfo[,c('Gene','Plate','Well')],by.x=c('Plate','Well'),by.y=c('Plate','Well'))
# write.csv(bac18,paste(ddir,'/Plate_annotations.csv',sep=''))
# 


# 
# #Get bacterial growth from the 3rd screen
# bac3<-read.table('3rd_screen_bacterial_growth/3rd_screen_bacterial_linear.csv',sep=',',quote = '"',header = TRUE,stringsAsFactors=FALSE)
# bac3<-rename(bac3, c("Plate"="Nplate", "Well"="NWell",'OD'='OD_raw'))
# #Get gene names and positions
# bac3a<-merge(bac3,scr3Mics[,c('Plate','Well','Gene','Nplate','NWell')],by=c('Nplate','NWell'),all.x = TRUE)
# bac3a[bac3a$Gene=='WT cont','Gene']<-'WT'
# bac3a[bac3a$Gene=='upp cont','Gene']<-'upp'
# 
# bac3a[bac3a$Gene=='purU' ,'Well']<-'F7'
# bac3a[bac3a$Gene=='purU' ,'Plate']<-'49'
# bac3a[bac3a$Gene=='hyfH' ,'Well']<-'A7'
# bac3a[bac3a$Gene=='hyfH' ,'Plate']<-'57'
# 
# bac3a[bac3a$Drug=='0uM','Drug']<-'0'
# bac3a[bac3a$Drug=='100uM','Drug']<-'100'
# 
# bac3m<-subset(bac3a,Time=='24hrs')
# 
# bac3m$OD<-bac3m$OD_raw-0.057392708
# write.csv(bac3m[,c('Plate','Well','Gene','Nplate','NWell','Replicate','Drug','Time','OD_raw','OD')],paste(ddir,'/3rd_screen_bacterial_annotated.csv',sep=''))

#Get bacterial annotation
bacannot<-read.table('Data_v4/Plate_annotations.csv',sep=',',quote = '"',header = TRUE,stringsAsFactors=FALSE)
bacannot$X<-NULL

#Get 4th screen
bac41<-read.table('4th_bacterial_screen_batch-1/4th_screen_bacterial_linear_1.csv',sep=',',quote = '"',header = TRUE,stringsAsFactors=FALSE)
bac41<-rename(bac41, c("Plate"="Nplate", "Well"="NWell",'OD'='OD_raw'))

bac42<-read.table('4th_bacterial_screen_batch-2/4th_screen_bacterial_linear_2.csv',sep=',',quote = '"',header = TRUE,stringsAsFactors=FALSE)
bac42<-rename(bac42, c("Plate"="Nplate", "Well"="NWell",'OD'='OD_raw'))

bac42[bac42$Replicate==3,'Replicate']<-5
bac42[bac42$Replicate==1,'Replicate']<-3
bac42[bac42$Replicate==2,'Replicate']<-4
bac42[bac42$Replicate==5,'Replicate']<-2

bac4f<-read.table('4th_bacterial_screen_fill/4th_screen_bacterial_linear_fill.csv',sep=',',quote = '"',header = TRUE,stringsAsFactors=FALSE)
bac4f<-rename(bac4f, c("Plate"="Nplate", "Well"="NWell",'OD'='OD_raw'))


#Needs separate annotation because of changed layout
bac4fa<-merge(bac4f,bacannot,by=c('Nplate','NWell'),all.x = TRUE)
#
#G1,H1,A2,Eq  fecI, dacC, idcA, purA

bac4fix<-data.frame(Nplate=c(8,8,8,8),
                    NWell=c('G1','H1','A2','E1'),
                    Gene=c('fecI','dacC','ldcA','purA'))
bac4fix<-merge(bac4fix,keioinfo[,c('Gene','Plate','Well')],by='Gene')
bac4fix$Gene<-as.character(bac4fix$Gene)

bac4fa[bac4fa$NWell %in% c('G1','H1','A2','E1'),c('Gene','Plate','Well')]<-bac4fix[match(bac4fa[bac4fa$NWell %in% c('G1','H1','A2','E1'),]$NWell,bac4fix$NWell),c('Gene','Plate','Well')]
bac4fa[!bac4fa$NWell %in% c('G1','H1','A2','E1'),'Replicate']<-bac4fa[!bac4fa$NWell %in% c('G1','H1','A2','E1'),'Replicate']+4

bac4a<-merge(bac41,bac42,all.x=TRUE,all.y=TRUE)
bac4al<-merge(bac4a,bacannot,by=c('Nplate','NWell'),all.x = TRUE)

bac4all<-merge(bac4al,bac4fa,all.x=TRUE,all.y=TRUE)

bac4all[bac4all$Drug=='0uM','Drug']<-'0'
bac4all[bac4all$Drug=='50uM','Drug']<-'50'

#Remove outliers

bac4all<-subset(bac4all,!(Gene=='ydcT' & Replicate==1 & Drug==0))
bac4all<-subset(bac4all,!(Gene=='lpd' & Replicate==2 & Drug==0))
bac4all<-subset(bac4all,!(Gene=='glyA' & Replicate==1 & Drug==50))
bac4all<-subset(bac4all,!(Gene=='glyA' & Replicate==1 & Drug==50))
bac4all<-subset(bac4all,!(Gene=='oppB' & Replicate==1 & Drug==50))
bac4all<-subset(bac4all,!(Gene=='oppB' & Replicate==2 & Drug==50))
bac4all<-subset(bac4all,!(Gene=='oppB' & Replicate==1 & Drug==0))


# 
# 
# ggplot(subset(bac4all,Gene=='Blank'),aes(x=OD_raw))+geom_histogram()+
#   facet_grid(Replicate+'Replicate'~'Drug, uM'+Drug)


back<-mean(subset(bac4all,Gene=='Blank' & OD_raw<0.075)$OD_raw)
back_sd<-sd(subset(bac4all,Gene=='Blank' & OD_raw<0.075)$OD_raw)
bac4all$OD<-bac4all$OD_raw-back

bac4sumN<-ddply(bac4all,.(Gene,Plate,Well,Drug),summarise,
      Replicates=length(Gene))

bac4sum<-dcast(bac4all,Gene+Plate+Well+Nplate+NWell+Drug~Replicate,mean,value.var = c('OD'))
bac4sum$OD_avg<-apply(bac4sum[,c('1','2','3','4')],1,mean,na.rm=TRUE)
bac4sum$OD_sd<-apply(bac4sum[,c('1','2','3','4')],1,sd,na.rm=TRUE)

#dir.create(ddir, showWarnings = TRUE, recursive = FALSE, mode = "0777")
write.csv(bac4sum,paste(ddir,'/Bacterial_growth_summary.csv',sep=''))


head(bac4sum)

# 
# bacm<-melt(subset(bact,!Gene %in% c('XXXXXXX','no bacteria','empty',''))[,c('Plate','Well','Gene','Drug','X1','X2','X3','X4')],
#            id=c('Plate','Well','Gene','Drug'),variable.name = 'Replicate',value.name='OD_raw')
# bacm$Replicate<-as.character(bacm$Replicate)
# 
# bacm$Screen<-'2'
# bac3m$Screen<-'3'
# 
# bac4a$Screen<-'4'

#Merge datasets
# bacs<-merge(bacm[,c('Plate','Well','Gene','Drug','Replicate','OD_raw','Screen')],
#             bac3m[,c('Plate','Well','Gene','Drug','Replicate','OD_raw','Screen')],
#             all.x=TRUE,all.y=TRUE)
# bacs<-merge(bacs[,c('Plate','Well','Gene','Drug','Replicate','OD_raw','Screen')],
#             bac4a[,c('Plate','Well','Gene','Drug','Replicate','OD_raw','Screen')],
#             all.x=TRUE,all.y=TRUE)

bacs<-bac4all[,!colnames(bac4all) %in% c('Nplate','NWell')]

bacs[bacs$Gene=='WT' ,'Well']<-''
bacs[bacs$Gene=='WT' ,'Plate']<-''
bacs[bacs$Gene=='upp','Plate']<-'61'
bacs[bacs$Gene=='upp','Well']<-'B1'

# bacs[bacs$Screen %in% c('2','3'),'OD']<-bacs[bacs$Screen %in% c('2','3'),'OD_raw']-0.057392708
# bacs[bacs$Screen=='4','OD']<-bacs[bacs$Screen=='4','OD_raw']-back

#Check for genes not present in the library
unique(subset(bacs,!Gene %in% keioinfo$Gene)$Gene)

bacavg<-ddply(subset(bacs,Gene!='Blank'), .(Plate,Well,Gene,Drug), summarise, NGM=mean(OD,na.rm=TRUE),NGM_sd=sd(OD,na.rm=TRUE)) #Plate,Well,Drug,


# bothscr<-subset(bacavg,Gene %in% c('upp','hyfF','WT','hyfB'))
bothscr<-bacavg
bothscrall<-merge(subset(bothscr,Drug==0),subset(bothscr,Drug==50),by=c('Plate','Well','Gene'),suffixes = c("_C","_D"))
bothscrall$Drug_C<-NULL
bothscrall$Drug_D<-NULL
#write.csv(bothscrall,paste(ddir,'/2nd-3rd_comparison.csv',sep=''))


scrcomp<-ggplot(bothscrall,aes(x=NGM_C,y=NGM_D))+
  #stat_smooth(aes(group = 1),method = "lm")+
  geom_abline(intercept=0,slope=1,alpha=0.5,aes(color='grey'),linetype='longdash')+
  geom_errorbarh(aes(xmax=NGM_C+NGM_sd_C,xmin=NGM_C-NGM_sd_C),height=.001,alpha=0.2)+
  geom_errorbar(aes(ymax=NGM_D+NGM_sd_D,ymin=NGM_D-NGM_sd_D),width=0.001,alpha=0.2)+
  geom_point(size=1)+ylim(-0.1, .35)+xlim(-0.1,.35)+
  geom_text(aes(label=as.character(Gene)))+
  ylab(expression(paste('Knockout strain growth OD - 5FU treatment')))+
  xlab('Knockout strain growth OD - Control')+
  ggtitle(expression(paste('Growth of knockout strains in control and 5FU treatment')))
  #scale_x_continuous(breaks=seq(0,.3,by=.05))+
scrcomp
# dev.copy2pdf(device=cairo_pdf,
#              file=paste(odir,"/Control-Treatment_NGM_growth_KO-Screen_comparison.pdf",sep = ''),
#              width=9,height=9)

# 
# bothscrside<-merge(subset(bothscr,Screen=='2'),subset(bothscr,Screen=='4'),by=c('Plate','Well','Gene','Drug'),suffixes = c("_2","_4"))
# bothscrside$Drug<-factor(bothscrside$Drug,levels=c('0','100'),labels=c('Control','Treatment'))
# 
# scrcomps<-ggplot(bothscrside,aes(x=NGM_2,y=NGM_4))+
#   geom_abline(intercept=0,slope=1,alpha=0.5,aes(color='grey'),linetype='longdash')+
#   geom_errorbarh(aes(xmax=NGM_2+NGM_sd_2,xmin=NGM_2-NGM_sd_2),height=.001,alpha=0.2)+
#   geom_errorbar(aes(ymax=NGM_4+NGM_sd_4,ymin=NGM_4-NGM_sd_4),width=0.001,alpha=0.2)+
#   geom_point(size=1)+ylim(0, .35)+xlim(0,.35)+
#   geom_text(aes(label=as.character(Gene)))+
#   stat_smooth(aes(group = 1),method = "lm")+
#   #scale_x_continuous(breaks=seq(0,.3,by=.05))+
#   ylab(expression(paste('Knockout strain growth OD - Screen 4')))+
#   xlab('Knockout strain growth OD - Screen 2')+
#   ggtitle(expression(paste('Growth of knockout strains in control and 5FU treatment')))+
#   labs(color='Screen')+
#   facet_grid(.~Drug)
# scrcomps
# dev.copy2pdf(device=cairo_pdf,
#              file=paste(odir,"/Control-Treatment_NGM_growth_KO-Screen_sided.pdf",sep = ''),
#              width=9,height=9)
# 


#Chance to select from several screens
#bacavgsel<-subset(bacavg,Screen=='4')

bacavgsel<-bacavg

bacall<-merge(subset(bacavgsel,Drug==0),subset(bacavgsel,Drug==50),by=c('Plate','Well','Gene'),suffixes = c("_C","_D"))
bacall$Drug_C<-NULL
bacall$Drug_D<-NULL
bacall$Screen_C<-NULL
bacall$Screen_D<-NULL


bacmict<-merge(mics,bacall,id=c('Plate','Well','Gene'),all.x = TRUE)


#Remove undisrupted and poor growing strains
#Make sure you can use Gene as unique identifier
intersect(rdupl,nodata$Gene)
intersect(rdupl,keioundis$Gene)

bacmic<-subset(bacmict,!Gene %in% keioundis$Gene & !Gene=='ydcT' )
#Data fully merged!!

undisrupted<-subset(bacmict,Gene %in% keioundis$Gene)
undisrupted$NComment<-'Undisrupted gene'

#Undisrupted and not in data
unmis<-setdiff(keioundis$Gene,undisrupted$Gene)

unno<-subset(nodata,Gene %in% keioundis$Gene)
unno$NComment<-'Undisrupted gene'

#Also remove ydcT
othernot<-subset(nodata,!Gene %in% keioundis$Gene)
otherno<-merge(othernot,subset(bacmict,Gene=='ydcT'),all.x=TRUE,all.y=TRUE)
otherno$NComment<-'Eliminated; no growth'

remsett<-merge(undisrupted,unno,all.x=TRUE,all.y=TRUE)
remsetf<-merge(remsett,otherno,all.x=TRUE,all.y=TRUE)

remsetf<-remsetf[,!colnames(remsetf) %in% c('Row','Column','Comment')]
remsetf<-rename(remsetf,c("NComment"="Comment"))

remset<-remsetf[,c('Gene','Plate','Well','JW_id','ECK','bno','0','1','2.5','5',
                   'MIC','MIC_sd','NGM_C','NGM_sd_C','NGM_D','NGM_sd_D',
                   'LB_22hr','MOPS_24hr','MOPS_48hr','Comment')]

#All unique gens and one instance of WT
length(allgenes)+1
#Genes in screen
length(bacmic$Gene)
#Genes removed
length(remset$Gene)

#Do we have everything?
length(allgenes)+1==length(bacmic$Gene)+length(remset$Gene)


#Remove knockouts that cause developmental delays
rmdev<-read.table('Keio_development_delays.csv',sep=',',quote = '"',header = TRUE,stringsAsFactors=FALSE)




#qbacmicq - no outliers
qbacmicq<-subset(bacmic,NGM_C>0.05 | is.na(NGM_C))

qbacmicqdt<-subset(qbacmicq,!Gene %in% setdiff(rmdev$Gene,'wbbL'))
qbacmicqd<-qbacmicqdt[!(qbacmicqdt$Gene=='wbbL' & qbacmicqdt$Plate=='21' & qbacmicqdt$Well=='A7'), ]

#Sanity checks:
poor<-subset(bacmic,NGM_C<0.05)$Gene
removed<-setdiff(bacmic$Gene,qbacmicq$Gene)
#Check whether the removal is OK
intersect(removed,poor)
poor




#Missing bacterial growth for MIC>=5
subset(bacmic,MIC>2.5 & is.na(NGM_C))
#write.csv(subset(bacmic,MIC>2.5 & is.na(NGM_C)),paste(ddir,'/MICs_and_bacterial_growth_Missing-growth_MICover2.5.csv',sep=''))




head(subset(bac4sum[order(bac4sum$OD_sd,decreasing = TRUE),],Gene!='Blank'),n=20)

write.csv(bacmic,paste(ddir,'/MICs_and_bacterial_growth-All.csv',sep=''))
write.csv(qbacmicq,paste(ddir,'/MICs_and_bacterial_growth-Clean.csv',sep=''))
write.csv(qbacmicqd,paste(ddir,'/MICs_and_bacterial_growth-Clean-NoDevDelays.csv',sep=''))

write.csv(remset,'Data_final/All_removed_from_screen.csv')



bacmic_shortt<-bacmic[,!colnames(bacmic) %in% c('0','1','2.5','5','EG','GI')]
bacmic_short<-bacmic_shortt[,c(1:5,8:9,6:7,10:16)]
write.csv(bacmic_short,paste(ddir,'/MICs_and_bacterial_growth-All_simple.csv',sep=''))


duplset<-subset(keioinfo,Gene %in% rdupl & !Gene %in% c('WT',NA,'present'))
duplset<-duplset[order(duplset$Gene),]

#write.csv(duplset,paste(ddir,'/Real_duplicates.csv',sep=''))

# 
# #PT lists
# PT_clean<-subset(qbacmicq,MIC>2.5 & !is.na(Gene))[,c('Gene','MIC')]
# PT_all<-subset(qbacmicq, !is.na(Gene) & !is.na(MIC) )[,c('Gene','MIC')]
# write.table(PT_clean,file='Data/Gene_MIC_above2.5_ForPT.tab',sep='\t',row.names = FALSE,quote=FALSE)
# write.table(PT_all,file='Data_Gene_MIC_all_ForPT.tab',sep='\t',row.names = FALSE,quote=FALSE)


