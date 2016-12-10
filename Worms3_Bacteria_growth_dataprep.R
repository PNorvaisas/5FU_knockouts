library('ggplot2')
library('gplots')
library('plyr')
library('reshape2')
library(tidyr)
library(BSDA)
library(xlsx)

library(sp)
library(car)
library('rafalib')
library(multcomp)
library('contrast')

theme_set(theme_light())

theme_update(panel.background = element_rect(colour = "black"),
             axis.text = element_text(colour = "black"))



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

#Read Keio info table with identifiers
keioinfo<-read.table('../Keio_library/Keio_library_fully_annotated.csv',sep=',',quote = '"',header = TRUE,stringsAsFactors=FALSE)
keioinfo$X<-NULL
keioinfo<-subset(keioinfo,!Plate %in% c('91','93','95'))


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
#bac4all<-subset(bac4all,!(Gene=='glyA' & Replicate==1 & Drug==50))
bac4all<-subset(bac4all,!(Gene=='oppB' & Replicate==1 & Drug==50))
bac4all<-subset(bac4all,!(Gene=='oppB' & Replicate==2 & Drug==50))
bac4all<-subset(bac4all,!(Gene=='oppB' & Replicate==1 & Drug==0))


# ggplot(subset(bac4all,Gene=='Blank'),aes(x=OD_raw))+geom_histogram()+
#   facet_grid(Replicate+'Replicate'~'Drug, uM'+Drug)


back<-mean(subset(bac4all,Gene=='Blank' & OD_raw<0.075)$OD_raw)
back_sd<-sd(subset(bac4all,Gene=='Blank' & OD_raw<0.075)$OD_raw)
bac4all$OD<-bac4all$OD_raw-back
bac4all$logOD<-log2(bac4all$OD)
bac4all$OD_raw<-NULL

bac4all[bac4all$Gene=='WT' ,'Well']<-''
bac4all[bac4all$Gene=='WT' ,'Plate']<-''
bac4all[bac4all$Gene=='upp','Plate']<-'61'
bac4all[bac4all$Gene=='upp','Well']<-'B1'

bac4all$Drug<-revalue(bac4all$Drug,c('0'='C','50'='T'))

bac4all$Index<-paste(bac4all$Gene,bac4all$Plate,bac4all$Well,sep='_')
bac4all$Gene<-as.factor(bac4all$Gene)
bac4all$Drug<-as.factor(bac4all$Drug)
bac4all$Index<-as.factor(bac4all$Index)

bac4all$Gene<-relevel(bac4all$Gene,ref = 'WT')
bac4all$Index<-relevel(bac4all$Index,ref = 'WT__')

#Is yjjG present
yjjG<-FALSE

#No duplicates, save to remove by gene name
if (yjjG==FALSE) {
  lowmicbac<-c('dcp','fre','glgC','oxyR','ybdR','yecE','yehA','yghT','ytjB','yjjG')
} else {
  lowmicbac<-c('dcp','fre','glgC','oxyR','ybdR','yecE','yehA','yghT','ytjB')
}
#

table(as.character(subset(bac4all,Gene %in% lowmicbac)$Gene))

poorgrowth<-c('atpB','atpE','atpF','atpG','atpE','lpd','sucA','atpE','ydcT','glnA')



dmw<-subset(bac4all,Gene!='Blank' & ! Gene %in% lowmicbac & ! Gene %in% poorgrowth)






WTref<-subset(dmw,Gene=='WT')
WTm<-melt(WTref,measure.vars = c('OD','logOD'),
          variable.name = 'Measure',value.name='Value')

WT_sum<-ddply(WTm,.(Drug,Measure),
              summarise,Mean=mean(Value,na.rm=TRUE),
              SD=sd(Value,na.rm=TRUE),
              N=length(Value))

WTsm<-melt(WT_sum,measure.vars = c('Mean','SD','N'),
           variable.name = 'Stat',value.name='Value')

WT<-dcast(WTsm,Drug~ Measure+Stat,value.var = 'Value')


# 
# dmwt<-merge(dmw,WT,by='Drug',all.x=TRUE,all.y=TRUE)

#Substract WT effect
dmw$WTDiff<-ifelse(dmw$Drug=='C',
                   dmw$logOD-WT[WT$Drug=='C',c("logOD_Mean")],
                   dmw$logOD-WT[WT$Drug=='T',c("logOD_Mean")])



#dmwt<-dmw[,c('Gene','Plate','Well','Drug','Replicate','WTDiff')]
#Change to dcast


#It's a bad idea to try and do subtractions one by one

ggplot(dmw,aes(x=WTDiff,fill=Drug))+geom_histogram(position='identity',alpha=0.5)

dmwm<-melt(dmw,id.vars = c('Gene','Index','Plate','Well','Drug','Replicate','Nplate','NWell'),variable.name = 'Measure',value.name = 'Value')
dmws<-dcast(dmwm,Gene+Index+Plate+Well+Nplate+NWell+Replicate~Drug+Measure,
                  fun.aggregate = NULL,value.var = c('Value'),fill = as.numeric(NA))

if (yjjG==TRUE) {
  write.csv(dmws,paste(ddir,'/Bacterial_growth_with_yjjG.csv',sep=''),row.names = FALSE)
} else {
  write.csv(dmws,paste(ddir,'/Bacterial_growth.csv',sep=''),row.names = FALSE)
}



# ggplot(dmws,aes(x=WTDiff_C))

# dmws<-merge(subset(dmw,Drug=='0'),subset(dmw,Drug=='50'),
#             by=c('Gene','Plate','Well','Nplate','NWell','Replicate'),all.x = TRUE,all.y = TRUE,suffixes=c('_C','_T'))


dmws$GT_Interaction_temp<-dmws$T_WTDiff-dmws$C_WTDiff
dmws$ODDiff_temp<-dmws$T_OD-dmws$C_OD

#poorgrowth<-unique(subset(dmws,C_OD<0.05)$Gene)

#write.csv(dmws,paste(ddir,'/Bacterial_growth.csv',sep=''),row.names = FALSE)


#Poor growers are not included in the fit
#dmwsc<-subset(dmws,Gene %in% poorgrowth)


dmwscstack<-melt(dmws,id=c('Gene','Index','Plate','Well','Nplate','NWell','Replicate'),
                 variable.name = 'Measure',value.name='Value')


dmwavg<-ddply(dmwscstack, .(Gene,Index,Plate,Well,Measure),summarise,
               Mean=mean(Value,na.rm = TRUE),SD=sd(Value,na.rm = TRUE),
               N=length(Value)) #Plate,Well,Drug,

dmwavgm<-melt(dmwavg,measure.vars = c('Mean','SD','N'),
               variable.name = 'Stat',value.name='Value')

dmwallsum<-dcast(dmwavgm,Gene+Index+Plate+Well~ Measure+Stat,value.var = 'Value')


#dmwallsumc<-subset(dmwallsum,!Gene %in% poorgrowth)



fitGTt<-lm(GT_Interaction_temp_Mean ~ C_WTDiff_Mean,data=dmwallsum)
GTrest<-summary(fitGTt)
GTrest
outlierTest(fitGTt)
otGT<-outlierTest(fitGTt)
outlistGT<-c(names(otGT$rstudent))

if (yjjG==FALSE) {
  outlistGT<-c()
}

subset(dmwallsum,rownames(dmwallsum) %in% outlistGT)




fitGT<-lm(GT_Interaction_temp_Mean ~ C_WTDiff_Mean,
           data=subset(dmwallsum,!rownames(dmwallsum) %in% outlistGT))
GTres<-summary(fitGT)
GTres

bGT<-GTres$coefficients[[1]]
aGT<-GTres$coefficients[[2]]

fitODt<-lm(ODDiff_temp_Mean ~ C_OD_Mean,data=dmwallsum)
ODrest<-summary(fitODt)
ODrest
outlierTest(fitODt)

otOD<-outlierTest(fitODt)
outlistOD<-c(names(otOD$rstudent))

if (yjjG==FALSE) {
  outlistOD<-c()
}

subset(dmwallsum,rownames(dmwallsum) %in% outlistOD)



fitOD<-lm(ODDiff_temp_Mean ~ C_OD_Mean,
           data=subset(dmwallsum,!rownames(dmwallsum) %in% outlistOD))
ODres<-summary(fitOD)
ODres

bOD<-ODres$coefficients[[1]]
aOD<-ODres$coefficients[[2]]


#Adjustment
dmwallsum$GT_Interaction<-dmwallsum$GT_Interaction_temp_Mean-(dmwallsum$C_WTDiff_Mean*aGT+bGT)
dmwallsum$ODDiff<-dmwallsum$ODDiff_temp_Mean-(dmwallsum$C_OD_Mean*aOD+bOD)


#Check for genes not present in the library
unique(subset(dmwallsum,!Gene %in% keioinfo$Gene)$Gene)

dmwallsum$N<-dmwallsum$C_OD_N
dmwallsum$ODDiff_SD<-dmwallsum$ODDiff_temp_SD
dmwallsum$GT_Interaction_SD<-dmwallsum$GT_Interaction_temp_SD

dmwallsum<-dmwallsum[, -grep("_N", colnames(dmwallsum))]
#dmwallsum<-dmwallsum[, -grep("OD_._pval", colnames(dmwallsum))]

ggplot(dmwallsum,aes(x=GT_Interaction))+geom_histogram()




#Linear modelling

bacfitWT<-lm(logOD~Drug+Drug:Index,data=dmw)
resbacWT<-summary(bacfitWT)

make.bacCT.frame<-function(res,len){
  df<-data.frame(res$coefficients)
  df<-df[(len+1):dim(df)[1],]
  names<-gsub('Index','',rownames(df))
  df$ID<-gsub('Drug','',names)
  df = transform(df, ID=colsplit(df$ID, ":", c('Drug', 'Index')) )
  df$Index<-df$ID$Index
  df$Drug<-df$ID$Drug
  df$ID<-NULL
  rownames(df)<-NULL
  return(df)
}

#Make WTDiff comparisons


bacCT<-make.bacCT.frame(resbacWT,2)
bacCT<-rename(bacCT,c('Estimate'='WTDiff',
                      'Std..Error'='WTDiff_SE',
                      't.value'='WTDiff_t.value',
                      'Pr...t..'='WTDiff_p.value'))
bacCT$WTDiff_FDR<-p.adjust(bacCT$WTDiff_p.value,method = 'fdr')
bacCTm<-melt(bacCT,id.vars = c('Index','Drug'),variable.name = 'Stats',value.name = 'Value')


bacCTf<-dcast(bacCTm,Index~Drug+Stats,value.var = 'Value')

fitbCTt<-lm(T_WTDiff~C_WTDiff,data=bacCTf)
outlierTest(fitbCTt)
otbCT<-outlierTest(fitbCTt)
outlistbCT<-c(names(otbCT$rstudent))

if (yjjG==FALSE) {
  outlistbCT<-c()
}

bacCTf[outlistbCT,]





bacCTfclean<-subset(bacCTf,!rownames(bacCTf) %in% outlistbCT)

fitbCT<-lm(T_WTDiff~C_WTDiff,data=bacCTfclean)
summary(fitbCT)


ggplot(bacCTf,aes(x=C_WTDiff,y=T_WTDiff))+
  geom_abline(slope = 1,intercept = 0,color='grey50',alpha=0.5)+
  geom_abline(slope = fitbCT$coefficients[[2]],intercept = fitbCT$coefficients[[1]],color='red',alpha=0.5)+
  geom_errorbar(aes(ymin=T_WTDiff-T_WTDiff_SE,ymax=T_WTDiff+T_WTDiff_SE),color='grey70')+
  geom_errorbarh(aes(xmin=C_WTDiff-C_WTDiff_SE,xmax=C_WTDiff+C_WTDiff_SE),color='grey70')+
  geom_point()

bacCTf<-transform(bacCTf, Index=colsplit(bacCTf$Index, "_", c('Gene','Plate','Well')))
bacCTf$Gene<-bacCTf$Index$Gene
bacCTf$Plate<-bacCTf$Index$Plate
bacCTf$Well<-bacCTf$Index$Well
bacCTf$Index<-NULL


#Get GT interaction

make.coef.frame<-function(res,len){
  tot<-dim(res$coefficients)[1]
  allcof<-res$coefficients[(len+1):tot]
  #print(allcof)
  #print(length(allcof))
  coflen<-length(allcof)
  w1<-allcof[1:(coflen/2)]
  #print(w1)
  w3<-allcof[(coflen/2+1):coflen]
  #print(w3)
  names<-gsub('Index','',rownames(res$coefficients)[(len+1):((coflen/2)+len)])
  #print(names)
  cf<-data.frame(Index=names,w1=w1,w3=w3)
  return(cf)
}

fitbaclog<-lm(logOD~Drug*Index,data=dmw)
resbaclog<-summary(fitbaclog)


cfbac<-make.coef.frame(resbaclog,2)

#Get fit to find outliers
genfitbact=lm(w3~w1,data=cfbac)
genresbact=summary(genfitbact)

outlierTest(genfitbact)
otbl<-outlierTest(genfitbact)
qqPlot(genfitbact, main="QQ Plot")

outlistbaclog<-c(names(otbl$rstudent))
#Outliers:
if (yjjG==FALSE) {
  outlistbaclog<-c()
}

print('Outliers in knockouts')
cfbac[outlistbaclog,]

#No outliers


#Get final fit
genfitbac=lm(w3~w1,data=subset(cfbac,!rownames(cfbac) %in% outlistbaclog))
genresbac=summary(genfitbac)
genresbac
bbl<-genresbac$coefficients[[1]]
abl<-genresbac$coefficients[[2]]


ggplot(cfbac,aes(x=w3))+geom_histogram(position='identity',alpha=0.5)


cfbac$w3a<-cfbac$w3-(abl*cfbac$w1+bbl)
remgens<-as.character(cfbac[outlistbaclog,'Index'])
print('Outliers')
remgens

cleancfbac<-subset(cfbac,!Index %in% remgens)
genfitbaca=lm(w3a~w1,data=cleancfbac)
genresbaca=summary(genfitbaca)
bbla<-genresbaca$coefficients[[1]]
abla<-genresbaca$coefficients[[2]]

ggplot(cfbac,aes(x=w1,y=w3a))+geom_point()+
  geom_abline(slope=abla,intercept=bbla,color='maroon')


ggplot(cfbac,aes(x=w3a))+geom_histogram(position='identity',alpha=0.5)

#Adjust coefficients
genesbac<-subset(cfbac)$Index
allcontrbac<-paste("`DrugT:Index",genesbac, "`-(", abl,")*`Index",genesbac,"`=",bbl,sep='')#,b12
adjconfbac<-glht(fitbaclog,linfct=allcontrbac)
adjresbac<-summary(adjconfbac,test = adjusted("none"))
#adjres122<-summary(adjconf12,test = adjusted("fdr"))


finalresbac<-data.frame('GT_Interaction'=adjresbac$test$coefficients,
                        'GT_SE'=adjresbac$test$sigma,
                        'GT_p.value'=adjresbac$test$pvalues)#,'FDR'=adjres122$test$pvalues
finalresbac$GT_Interaction<-finalresbac$GT_Interaction-bbl


finalresbac$Index<-genesbac
finalresbac<-transform(finalresbac, Index=colsplit(finalresbac$Index, "_", c('Gene','Plate','Well')))
finalresbac$Gene<-finalresbac$Index$Gene
finalresbac$Plate<-finalresbac$Index$Plate
finalresbac$Well<-finalresbac$Index$Well
finalresbac$Index<-NULL



rownames(finalresbac)<-NULL
finalresbac$GT_FDR<-p.adjust(finalresbac$GT_p.value,method = 'fdr')

finalresbac<-finalresbac[,c('Gene','Plate','Well','GT_Interaction','GT_SE','GT_p.value','GT_FDR')]

ggplot(finalresbac,aes(x=GT_Interaction,y=-log10(GT_p.value)))+geom_point()+
  geom_hline(yintercept = -log10(0.05),color='red',alpha=0.5)+
  geom_text(aes(label=Gene))

ggplot(finalresbac,aes(x=GT_Interaction,y=-log10(GT_FDR)))+geom_point()+
  geom_hline(yintercept = -log10(0.05),color='red',alpha=0.5)+
  geom_text(aes(label=Gene))


#fitp<-lm(-log10(CT_FDR)~-log10(CT_Pval),data<-subset(finalres,CT_Pval<0.05))

ggplot(finalresbac,aes(x=GT_p.value,y=GT_FDR))+geom_point()+
  geom_hline(yintercept = 0.05,color='red',alpha=0.5)+
  geom_vline(xintercept = 0.05,color='red',alpha=0.5)


#All stats results



head(finalresbac)
head(bacCTf)

head(dmwallsum)

rawsum<-c('Gene','Plate','Well',
          'C_OD_Mean','C_OD_SD',
          'T_OD_Mean','T_OD_SD',
          'C_logOD_Mean','C_logOD_SD',
          'T_logOD_Mean','T_logOD_SD',
          'C_WTDiff_SD','T_WTDiff_SD',
          'GT_Interaction_SD')
dmwallsum[dmwallsum$Gene=='WT',c('C_WTDiff_SD','T_WTDiff_SD','GT_Interaction_SD')]<-NA


bacalltt<-merge(dmwallsum[,colnames(dmwallsum) %in% c(rawsum,'N')],bacCTf,by=c('Gene','Plate','Well'),all.x=TRUE)
bacallt<-merge(bacalltt,finalresbac,by=c('Gene','Plate','Well'),all.x=TRUE)
bacallt<-rename(bacallt,c('GT_Interaction_SD'='GT_SD'))

colnames(bacallt)
bacall<-bacallt[,c(1:3,14,4:7,9:12,
                   16,8,17:20,
                   21,13,22:25,
                   26,15,27:29)]
colnames(bacall)


if (yjjG==TRUE) {
  write.csv(bacall,paste(ddir,'/Bacterial_growth_summary_with_yjjG.csv',sep=''),row.names = FALSE)
} else {
  write.csv(bacall,paste(ddir,'/Bacterial_growth_summary.csv',sep=''),row.names = FALSE)
  
}


