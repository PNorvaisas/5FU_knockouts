library('ggplot2')
library('gplots')
library('plyr')
library('reshape2')
library('tidyr')
library('gtools')

library(sp)
library(xlsx)
library(car)
library('rafalib')
library(multcomp)
library('contrast')

lmsum<-function(m){fres<-summary(m)
l <- list(b = as.double(coef(m)[1], digits = 2),
          a = as.double(coef(m)[2], digits = 2),
          r2 = as.double(summary(m)$r.squared, digits = 3),
          p2 = as.double(pf(fres$fstatistic[1], fres$fstatistic[2], fres$fstatistic[3],lower.tail = FALSE)[[1]], digits = 5))
return(l)
}

make.CT.frame<-function(res,len){
  df<-data.frame(res$coefficients)
  df<-df[(len+1):dim(df)[1],]
  names<-gsub('UniqueName','',rownames(df))
  df$ID<-gsub('Type','',names)
  df = transform(df, ID=colsplit(df$ID, ":", c('Type', 'UniqueName')) )
  df$UniqueName<-df$ID$UniqueName
  df$Type<-df$ID$Type
  df$ID<-NULL
  rownames(df)<-NULL
  return(df)
}


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
  names<-gsub('UniqueName','',rownames(res$coefficients)[(len+1):((coflen/2)+len)])
  #print(names)
  cf<-data.frame(UniqueName=names,w1=w1,w3=w3)
  return(cf)
}


trim <- function (x) gsub("^\\s+|\\s+$", "", x)

setwd("~/Projects/2015-Metformin/Worms")

theme_set(theme_light())


ddir<-'Data_final'
odir<-'Figures_final/Biolog/'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")

#Please add descriptors used in Design file to this list!
#Plase add all non-numeric or descriptory columns in Summary file to this list!
design<-c('File','Plate','Strain','Type','Media','Replicate','Well','Index',
          'Data','Name','EcoCycID','KEGG_ID','CAS_ID','Group','Description','Metformin_mM','Sugar_mM','Descriptor')

data<-read.table('./Biolog/Bacteria_data/Bacteria_collected/Summary.csv',sep=',',quote = "\"",header=TRUE)


cols<-colnames(data)
realanot<-intersect(design,cols)

wellspec<-c('Well','Index','Name','EcoCycID','KEGG_ID','CAS_ID','Group','Description')


#File	Plate	Strain	Type	Media	Replicate	Well	Index	Data	Name	EcoCycID	Group	Description

dm<-melt(data,id=realanot,
         variable.name='Descriptor',value.name='Value')

dm$ReplicateB<-paste('B',as.character(dm$Replicate),sep='')

dm$Type<-ifelse(dm$Type=='Control','C','T')

byrepsided<-dcast(dm,Plate+Strain+Well+Index+Name+EcoCycID+KEGG_ID+CAS_ID+Descriptor~Type+ReplicateB,
             fun.aggregate = NULL,value.var = c('Value'),fill = as.numeric(NA))


#Summary
dm_sum<-ddply(dm,.(Plate,Type,Strain,Well,Index,Name,EcoCycID,KEGG_ID,CAS_ID,Descriptor),
              summarise,Mean=mean(Value,na.rm=TRUE),SD=sd(Value,na.rm=TRUE),N=length(Value))

#Give quick summary
ref<-subset(dm_sum,Well=='A1')
ref<-rename(ref, c("Mean"="A1_Mean", "SD"="A1_SD", "N"="A1_N"))
ref<-ref[,!colnames(ref) %in% c('Well','Index','Name','EcoCycID','KEGG_ID','CAS_ID')]

dm_ref<-merge(dm_sum,ref,by=c('Plate','Type','Strain','Descriptor'),all.x=TRUE)



dmrefs<-melt(dm_ref,measure.vars = c('Mean','SD','N','A1_Mean','A1_SD','A1_N'),
         variable.name='Stats',value.name='Value')

sumsided<-dcast(dmrefs,Plate+Strain+Well+Index+Name+EcoCycID+KEGG_ID+CAS_ID+Descriptor~Type+Stats,
                  fun.aggregate = NULL,value.var = c('Value'),fill = as.numeric(NA))

sumsidedw<-subset(sumsided,Descriptor %in% c('Int_750nm','Int_750nm_log'))[, -grep("_N", colnames(sumsided))]
byrepsidedw<-subset(byrepsided,Descriptor %in% c('Int_750nm','Int_750nm_log'))

write.csv(byrepsidedw,paste(ddir,'/Biolog_Bacteria_Summary_replicates_by-treatment.csv',sep=''))
write.csv(sumsidedw,paste(ddir,'/Biolog_Bacteria_Summary_averages_by-treatment.csv',sep=''))


#Worm data
wormd<-read.table('./Biolog/Worm_data/Biolog_Worm_linearised.csv',sep=',',quote = "\"",header=TRUE)
bioinfo<-read.table('../Biolog_results/Biolog_metabolites_Ecocyc.csv',sep=',',quote = "\"",header=TRUE)
bioinfo$Metabolite<-trim(bioinfo$Metabolite)

wormdata<-merge(wormd,bioinfo,by=c('Plate','Well'),all.x=TRUE)
wormdata$Name<-NULL
wormdata<-rename(wormdata,c('Metabolite'='Name'))

wormbyrep<-dcast(wormdata,Plate+Well+Index+Name+EcoCycID+KEGG_ID+CAS_ID~Replicate,
                 fun.aggregate = NULL,value.var = c('Value'),fill = as.numeric(NA))


wormsum<-ddply(wormdata,.(Plate,Well),
               summarise,W_Mean=mean(Value,na.rm=TRUE),W_SD=sd(Value,na.rm=TRUE),
               W_Median=median(Value,na.rm=TRUE))

wormall<-merge(wormbyrep,wormsum,by=c('Plate','Well'))
wormall<-rename(wormall,c('1'='W1','2'='W2','3'='W3','4'='W4'))
write.csv(wormall,paste(ddir,'/Biolog_Worms_Summary.csv',sep=''))


intsa<-subset(dm,Descriptor %in% c('Int_750nm_log','Int_750nm'))
#,'Int_750nm','Max_750nm','Max_750nm_log','X24h_750nm','X24h_750nm_log'
intsa$Type<-as.factor(ints$Type)

#Linear modelling
intsa$Name <- relevel(ints$Name, ref = "Negative Control")
#ints$UniqueName<-ints$Name)
intsa$UniqueName<-as.factor(gsub('`','',ifelse(ints$Plate=='PM5',
                                  paste(as.character(ints$Name),'-PM5',sep=''),
                                  as.character(ints$Name))))
intsa$UniqueName<-relevel(intsa$UniqueName, ref = "Negative Control")
intsa$PlateGroup<-ifelse(intsa$Plate=='PM5','PM5','PM1&PM2A')

#Calculate ODInt stats



intsr<-subset(intsa,Descriptor %in% c('Int_750nm'))





#Check Negative control in plates
ints<-subset(intsa,Descriptor %in% c('Int_750nm_log'))

back=subset(ints,Name=='Negative Control')

ggplot(back,aes(y=Value,x=Plate,colour=Type))+
  geom_boxplot()+geom_point()+
  ylab('Growth')

fitneg<-lm(Value~Plate+Type:Plate,data=subset(ints,Name=='Negative Control'))#Type+Plate+
summary(fitneg)



#Get separate datasets
ints12<-subset(ints,PlateGroup=='PM1&PM2A' )
ints12$UniqueName <- relevel(ints12$UniqueName, ref = "Negative Control")

ints5<-subset(ints,PlateGroup=='PM5')
ints5$UniqueName <- relevel(ints5$UniqueName, ref = "Negative Control-PM5")


#Make NGMDiff comparisons
fitCT12<-lm(Value~Type+Type:UniqueName,data=ints12)
resCT12<-summary(fitCT12)


fitCT5<-lm(Value~Type+Type:UniqueName,data=ints5)
resCT5<-summary(fitCT5)



resCT5df<-make.CT.frame(resCT5,2)
resCT5df$PlateGroup<-'PM5'
resCT12df<-make.CT.frame(resCT12,2)
resCT12df$PlateGroup<-'PM1&PM2A'

 
resCT<-merge(resCT5df,resCT12df,all=TRUE)

resCT<-rename(resCT,c('Estimate'='NGMDiff','Std..Error'='SE','Pr...t..'='p.value'))
resCT$FDR<-p.adjust(resCT$p.value,method = 'fdr')
resCTm<-melt(resCT,id.vars = c('UniqueName','PlateGroup','Type'),variable.name = 'Stats',value.name = 'Value')


resCTf<-dcast(resCTm,UniqueName+PlateGroup~Type+Stats,value.var = 'Value')

fitCTt<-lm(T_NGMDiff~C_NGMDiff,data=resCTf)
outlierTest(fitCTt)
otCT<-outlierTest(fitCTt)
outlistCT<-c(names(otCT$rstudent))

subset(resCTf,rownames(resCTf) %in% outlistCT)

resCTfclean<-subset(resCTf,!rownames(resCTf) %in% outlistCT & !UniqueName %in% c('2-Hydroxy Benzoic Acid','L-Leucine'))

fitCT<-lm(T_NGMDiff~C_NGMDiff,data=resCTfclean)
summary(fitCT)
#plot(fitCT)


ggplot(resCTf,aes(x=C_NGMDiff,y=T_NGMDiff))+
  geom_abline(slope = 1,intercept = 0,color='grey50',alpha=0.5)+
  geom_abline(slope = fitCT$coefficients[[2]],intercept = fitCT$coefficients[[1]],color='red',alpha=0.5)+
  geom_errorbar(aes(ymin=T_NGMDiff-T_SE,ymax=T_NGMDiff+T_SE),color='grey70')+
  geom_errorbarh(aes(xmin=C_NGMDiff-C_SE,xmax=C_NGMDiff+C_SE),color='grey70')+
  geom_point()


#Get Treatment vs Metabolite interaction

fit12<-lm(Value~Type*UniqueName,data=ints12)
result12<-summary(fit12)

fit5<-lm(Value~Type*UniqueName,data=ints5)
result5<-summary(fit5)


cf12<-make.coef.frame(result12,2)
cf12$PlateGroup<-'PM1&PM2A'
cf5<-make.coef.frame(result5,2)
cf5$PlateGroup<-'PM5'


cf<-merge(cf12,cf5,all=TRUE)

#Get fit to find outliers
genfit12t=lm(w3~w1,data=cf12)
genres12t=summary(genfit12t)
genres12t

#plot(genfit12t)
#influencePlot(genfit12t,	id.method="identify", main="Influence Plot", sub="Circle size is proportial to Cook's Distance" )
outlierTest(genfit12t)
ot12<-outlierTest(genfit12t)
qqPlot(genfit12t, main="QQ Plot")

outlist12<-c(names(ot12$rstudent))
#Outliers:
print('Outliers in PM1&PM2A')
cf12[outlist12,]


#Get final fit
genfit12=lm(w3~w1,data=subset(cf12,!rownames(cf12) %in% outlist12 & !UniqueName %in% c('2-Hydroxy Benzoic Acid','L-Leucine')))
genres12=summary(genfit12)
genres12
b12<-genres12$coefficients[[1]]
a12<-genres12$coefficients[[2]]


#Get fit to find outliers
genfit5t=lm(w3~w1,data=cf5)
genres5t=summary(genfit5t)
genres5t
#plot(genfit5t)
#influencePlot(genfit12t,	id.method="identify", main="Influence Plot", sub="Circle size is proportial to Cook's Distance" )
#Bonferroni outlier test

ot5<-outlierTest(genfit5t)
ot5
qqPlot(genfit5t, main="QQ Plot")

outlist5<-c(names(ot5$rstudent))
#Outliers:
print('Outliers in PM5')
cf5[outlist5,]

#Get final fit
genfit5=lm(w3~w1,data=subset(cf5,!rownames(cf5) %in% outlist5))
genres5=summary(genfit5)
genres5
b5<-genres5$coefficients[[1]]
a5<-genres5$coefficients[[2]]

#Compare w3~w1 between plates
fitw13<-lm(w3~w1*PlateGroup,data=cf)
resw13<-summary(fitw13)
resw13

ggplot(cf,aes(x=w1,y=w3,color=PlateGroup))+geom_point()+
  geom_abline(slope=a12,intercept=b12,color='red')+
  geom_abline(slope=a5,intercept=b5,color='blue')


cf$w3a<-ifelse(cf$PlateGroup=='PM5',cf$w3-(a5*cf$w1+b5),cf$w3-(a12*cf$w1+b12))
remmets<-c('2-Hydroxy Benzoic Acid','L-Leucine',as.character(cf5[outlist5,'UniqueName']),as.character(cf12[outlist12,'UniqueName']))
print('Outliers')
remmets
cleancf<-subset(cf,!UniqueName %in% remmets)
genfit125=lm(w3a~w1,data=cleancf)
genres125=summary(genfit125)
b125<-genres125$coefficients[[1]]
a125<-genres125$coefficients[[2]]

ggplot(cf,aes(x=w1,y=w3a,color=PlateGroup))+geom_point()+
  geom_abline(slope=a125,intercept=b125,color='maroon')


ggplot(cf,aes(x=w3a,fill=PlateGroup))+geom_histogram(position='identity',alpha=0.5)

#AdjustPM1&PM2A
names12<-subset(cf,PlateGroup=='PM1&PM2A')$UniqueName
allcontr12<-paste("`TypeT:UniqueName",names12, "`-(", a12,")*`UniqueName",names12,"`=",b12,sep='')#,b12
adjconf12<-glht(fit12,linfct=allcontr12)
adjres12<-summary(adjconf12,test = adjusted("none"))
#adjres122<-summary(adjconf12,test = adjusted("fdr"))

#Adjust PM5
names5<-subset(cf,PlateGroup=='PM5')$UniqueName
allcontr5<-paste("`TypeT:UniqueName",names5, "`-(", a5,")*`UniqueName",names5,"`=",b5,sep='')#,b5
adjconf5<-glht(fit5,linfct=allcontr5)
adjres5<-summary(adjconf5,test = adjusted("none"))
#adjres52<-summary(adjconf5,test = adjusted("fdr"))


finalres12<-data.frame('MT_Interaction'=adjres12$test$coefficients,
                       'MT_SE'=adjres12$test$sigma,
                       'MT_p.value'=adjres12$test$pvalues)#,'FDR'=adjres122$test$pvalues
finalres12$UniqueName<-names12
finalres12$PlateGroup<-'PM1&PM2A'


finalres5<-data.frame('MT_Interaction'=adjres5$test$coefficients,
                      'MT_SE'=adjres5$test$sigma,
                      'MT_p.value'=adjres5$test$pvalues)#,'FDR'=adjres52$test$pvalues
finalres5$UniqueName<-names5
finalres5$PlateGroup<-'PM5'

finalres12$MT_Interaction<-finalres12$MT_Interaction-b12
finalres5$MT_Interaction<-finalres5$MT_Interaction-b5
finalres<-merge(finalres12,finalres5,all=TRUE)

rownames(finalres)<-NULL
finalres$MT_FDR<-p.adjust(finalres$MT_p.value,method = 'fdr')

finalres<-finalres[,c('UniqueName','PlateGroup','MT_Interaction','MT_SE','MT_p.value','MT_FDR')]

ggplot(finalres,aes(x=MT_Interaction,y=-log10(MT_p.value),color=PlateGroup))+geom_point()+
  geom_hline(yintercept = -log10(0.05),color='red',alpha=0.5)+
  geom_text(aes(label=UniqueName))

#fitp<-lm(-log10(CT_FDR)~-log10(CT_Pval),data<-subset(finalres,CT_Pval<0.05))

#-log10(CT_FDR)
ggplot(finalres,aes(x=MT_p.value,y=MT_FDR,color=PlateGroup))+geom_point()+
  geom_hline(yintercept = 0.05,color='red',alpha=0.5)+
  geom_vline(xintercept = 0.05,color='red',alpha=0.5)
  #xlim(0,15)+ylim(0,15)
  #+
  #geom_text(aes(label=UniqueName))

#Combined stats from linear model
allstat<-merge(resCTf,finalres[,! colnames(finalres) %in% c("PlateGroup") ],by='UniqueName')






#Calculate stats manually

intrefs<-subset(intsr,Well=='A1')
intrefs<-rename(intrefs, c("Value"="A1"))
intrefs<-intrefs[,union(setdiff(realanot,wellspec),c('A1'))]


in_ref<-merge(intsr,intrefs,by=c('File','Strain','Plate','Type','Data','Replicate'),
              all.x = TRUE,all.y=TRUE)

in_ref$logODInt<-log(in_ref$Value,2)
in_ref$NGMDiff<-in_ref$logODInt-log(in_ref$A1,2)



in_ref<-rename(in_ref,c('Value'='ODInt'))

inrefm<-melt(in_ref[,!colnames(in_ref) %in% c('A1')],measure.vars = c('ODInt','logODInt','NGMDiff'),
             variable.name='Measure',value.name='Value')

insided<-dcast(inrefm,PlateGroup+Plate+Well+Replicate+ReplicateB+Name+UniqueName+EcoCycID+KEGG_ID+CAS_ID+Group+Description~Measure+Type,
                fun.aggregate = NULL,value.var = c('Value'),fill = as.numeric(NA))


PM1Names<-unique(subset(insided,Plate %in% c('PM1'))$Name)
PM2Names<-unique(subset(insided,Plate %in% c('PM2A'))$Name)
PM5Names<-unique(subset(insided,Plate=='PM5')$Name)
PM12Names<-unique(subset(insided,Plate %in% c('PM1','PM2A'))$Name)


dubl<-setdiff(intersect(PM5Names,PM12Names),c('Negative Control','Positive Control'))


exclude<-c('2-Hydroxy Benzoic Acid','L-Leucine')


insided$MT_Interaction<-insided$NGMDiff_T-insided$NGMDiff_C


insidedm<-melt(insided,id.vars = c('PlateGroup','Plate','Well','Replicate','ReplicateB',
                                     'Name','UniqueName','EcoCycID','KEGG_ID','CAS_ID','Group','Description'),
               variable.name='Measure',value.name='Value')

insidstat<-ddply(insidedm,.(PlateGroup,Plate,Well,Name,UniqueName,EcoCycID,KEGG_ID,CAS_ID,Group,Description,Measure),
                 summarise,
                 Mean=mean(Value,na.rm=TRUE),
                 SD=sd(Value,na.rm=TRUE),
                 p.value=ifelse(sd(Value,na.rm=TRUE)==0,NA,t.test(Value,mu=0)$p.value))

insidstatm<-melt(insidstat,measure.vars = c('Mean','SD','p.value'),
                 variable.name='Stat',value.name='Value')


insidsum<-dcast(insidstatm,PlateGroup+Plate+Well+Name+UniqueName+EcoCycID+KEGG_ID+CAS_ID+Group+Description~Measure+Stat,
                fun.aggregate = NULL,value.var = c('Value'),fill = as.numeric(NA))

insidsumA<-insidsum
insidsumA$PlateGroup<-'All'

insidsumd<-merge(insidsum,insidsumA,all.x=TRUE,all.y=TRUE)


sumfitsjoined<-ddply(subset(insidsumd,!UniqueName %in% exclude &
                              !Name %in% c('Positive Control','Negative Control')), .(PlateGroup), summarise,
                     logODInt_a=lmsum(lm(`logODInt_T_Mean` ~ `logODInt_C_Mean`))$a,
                     logODInt_b=lmsum(lm(`logODInt_T_Mean` ~ `logODInt_C_Mean`))$b,
                     logODInt_r2=lmsum(lm(`logODInt_T_Mean` ~ `logODInt_C_Mean`))$r2,
                     logODInt_p=lmsum(lm(`logODInt_T_Mean` ~ `logODInt_C_Mean`))$p2,
                     NGMDiff_a=lmsum(lm(`NGMDiff_T_Mean` ~ `NGMDiff_C_Mean`))$a,
                     NGMDiff_b=lmsum(lm(`NGMDiff_T_Mean` ~ `NGMDiff_C_Mean`))$b,
                     NGMDiff_r2=lmsum(lm(`NGMDiff_T_Mean` ~ `NGMDiff_C_Mean`))$r2,
                     NGMDiff_p=lmsum(lm(`NGMDiff_T_Mean` ~ `NGMDiff_C_Mean`))$p2,
                     MT_Interaction_a=lmsum(lm(`MT_Interaction_Mean` ~ `NGMDiff_C_Mean`))$a,
                     MT_Interaction_b=lmsum(lm(`MT_Interaction_Mean` ~ `NGMDiff_C_Mean`))$b,
                     MT_Interaction_r2=lmsum(lm(`MT_Interaction_Mean` ~ `NGMDiff_C_Mean`))$r2,
                     MT_Interaction_p=lmsum(lm(`MT_Interaction_Mean` ~ `NGMDiff_C_Mean`))$p2)
sumfitsjoined$Type<-'Summary'

ggplot(insidsum,aes(x=NGMDiff_C_Mean,y=MT_Interaction_Mean,color=PlateGroup))+geom_point()

insidsumf<-merge(insidsum,sumfitsjoined[sumfitsjoined$PlateGroup %in% c('PM1&PM2A','PM5'),c('PlateGroup','MT_Interaction_a','MT_Interaction_b')],by='PlateGroup',all.x=TRUE,all.y=TRUE)

insidsumf$MT_Interaction_norm_Mean<-insidsumf$MT_Interaction_Mean-(insidsumf$NGMDiff_C_Mean*insidsumf$MT_Interaction_a+insidsumf$MT_Interaction_b)
insidsumf$MT_Interaction_norm_SD<-insidsumf$MT_Interaction_SD


#insidsumf<-insidsumf[, -grep("logODInt_._p.value", colnames(insidsumf))]
insidsumf<-insidsumf[,!colnames(insidsumf) %in% c('MT_Interaction_a','MT_Interaction_b')]


simpleinc<-c('PlateGroup','Plate','Well','Name','UniqueName',
             'EcoCycID','KEGG_ID','CAS_ID','Group','Description',
             'ODInt_C_Mean','ODInt_C_SD','ODInt_T_Mean','ODInt_T_SD',
             'logODInt_C_Mean','logODInt_C_SD','logODInt_T_Mean','logODInt_T_SD',
             'NGMDiff_C_SD','NGMDiff_T_SD','MT_Interaction_SD')
insidsumex<-insidsumf[,simpleinc]
#Join raw data and fits for normalisation
#insided_f<-merge(insided,repfitsjoined[,c('PlateGroup','NGMDiff_a','NGMDiff_b')],by='PlateGroup',all.x=TRUE,all.y=TRUE)



#insided_f$NGMDiff_T_norm<-(insided_f$NGMDiff_T-insided_f$NGMDiff_b)/insided_f$NGMDiff_a



completestats<-merge(insidsumex,allstat,by=c('PlateGroup','UniqueName'),all.x=TRUE,all.y=TRUE)
completestats<-rename(completestats,c('C_NGMDiff'='NGMDiff_C_Mean',
                                      'C_SE'='NGMDiff_C_SE',
                                      'C_t.value'='NGMDiff_C_t.value',
                                      'C_p.value'='NGMDiff_C_p.value',
                                      'C_FDR'='NGMDiff_C_FDR',
                                      'T_NGMDiff'='NGMDiff_T_Mean',
                                      'T_SE'='NGMDiff_T_SE',
                                      'T_t.value'='NGMDiff_T_t.value',
                                      'T_p.value'='NGMDiff_T_p.value',
                                      'T_FDR'='NGMDiff_T_FDR',
                                      'MT_Interaction_SD'='MT_SD'))

allfits<-sumfitsjoined


#Join stats


allrep<-merge(wormall[,c('Plate','Well','W1','W2','W3','W4','W_Mean','W_SD','W_Median')],
              subset(byrepsided,Descriptor=='Int_750nm')[,c('Plate','Well',
                                                                'C_B1','C_B2','C_B3','C_B4',
                                                                'T_B1','T_B2','T_B3','T_B4')],
              by=c('Plate','Well'))

  
allwbr<-merge(allrep,completestats,by=c('Plate','Well'))
colnames(allwbr)
allwb<-allwbr[,c(18,1:2,19:25,3:17,26:33,37,34,38:41,42,35,43:46,47,36,48:50)]
colnames(allwb)

exclude<-c('2-Hydroxy Benzoic Acid','L-Leucine')
allwb<-subset(allwb,!UniqueName %in% exclude)

subset(allwb,UniqueName %in% exclude)


allwbw<-allwb[, -grep("_t.value", colnames(allwb))]


PM12names<-subset(allwb,PlateGroup=='PM1&PM2A' & Name!='Negative Control')$Name
PM5names<-subset(allwb,PlateGroup=='PM5'& Name!='Negative Control')$Name


intersect(PM12names,PM5names)
subset(allwb,Name=='Negative Control')
#Unique metabolites
length(unique(allwb$Name))-1
length(subset(allwb,W_Median>1)$Name)
length(unique(subset(allwb,W_Median>1)$Name))





write.csv(allwb,paste(ddir,'/Biolog_Combined_Summary_Statistics.csv',sep=''),row.names = FALSE)



#3 negative controls, 1 positive control, 




colnames(allwb)

allwbw<-allwb[, -grep("_t.value", colnames(allwb))]
allwbw<-allwbw[, -grep("C_B", colnames(allwbw))]
allwbw<-allwbw[, -grep("T_B", colnames(allwbw))]
colnames(allwbw)

colnames(allwbw)<-gsub('_Mean','',colnames(allwbw))
colnames(allwbw)<-gsub('logODInt_C','C_logODInt',colnames(allwbw))
colnames(allwbw)<-gsub('logODInt_T','T_logODInt',colnames(allwbw))
colnames(allwbw)<-gsub('ODInt_C','C_ODInt',colnames(allwbw))
colnames(allwbw)<-gsub('ODInt_T','T_ODInt',colnames(allwbw))
colnames(allwbw)<-gsub('NGMDiff_C','C_NGMDiff',colnames(allwbw))
colnames(allwbw)<-gsub('NGMDiff_T','T_NGMDiff',colnames(allwbw))

colnames(allwbw)





bioexpl<-read.table(paste(ddir,'/Biolog_column_explanation.csv',sep=''),
                    sep=',',quote = '"',header = TRUE,stringsAsFactors=FALSE)


#Officia file with explanations
explrd<-subset(bioexpl,Column %in% colnames(allwbw))
explrd<-explrd[match(colnames(allwbw),explrd$Column),]

write.xlsx2(explrd, file=paste(ddir,'/Biolog_Combined_Summary_Statistics.xlsx',sep=''), sheetName="Readme",row.names = FALSE,showNA=FALSE)
write.xlsx2(allwbw, file=paste(ddir,'/Biolog_Combined_Summary_Statistics.xlsx',sep=''), sheetName="Data", append=TRUE,row.names = FALSE)#showNA=FALSE


#Update manually
write.xlsx2(explrd, file='/Users/Povilas/Projects/B-D-H paper/First submission MS files/tables/Scott et al_Table S5.xlsx',
            sheetName="Biolog_Readme", append=TRUE,row.names = FALSE,showNA=FALSE)
write.xlsx2(allwbw, file='/Users/Povilas/Projects/B-D-H paper/First submission MS files/tables/Scott et al_Table S5.xlsx',
            sheetName="Biolog_Data", append=TRUE,row.names = FALSE)#showNA=FALSE




