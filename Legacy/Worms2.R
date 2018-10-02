library('ggplot2')
library('gplots')
library('plyr')
library('reshape2')
library(tidyr)
library(quantreg)
#Vennerable installation: install.packages("Vennerable", repos="http://R-Forge.R-project.org")

#library(quantreg)

lm_eqn = function(m) {
  fres<-summary(m)
  l <- list(a = format(coef(m)[1], digits = 2),
            b = format(abs(coef(m)[2]), digits = 2),
            r2 = format(summary(m)$r.squared, digits = 3),
            p2 = format(pf(fres$fstatistic[1], fres$fstatistic[2], fres$fstatistic[3],lower.tail = FALSE)[[1]], digits = 3));
  
  if (coef(m)[2] >= 0)  {
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2*","~~italic(p)~"="~p2,l)
  } else {
    eq <- substitute(italic(y) == a - b %.% italic(x)*","~~italic(r)^2~"="~r2*","~~italic(p)~"="~p2,l)    
  }
  
  as.character(as.expression(eq));                 
}


evalmic=function(m){
  m<-as.numeric(m)
  mic=1
  if (!is.na(m[[3]]) & m[[3]]>0) {
    mic=2.5
  }
  if (!is.na(m[[4]]) & m[[4]]>0) {
    mic=5
  }
  if (!is.na(m[[5]]) & m[[5]]>0) {
    mic=5
  }
  return(mic)
}

#Primary screen analysis
#scr1r<-read.table('Primary_screen_raw_fixed.csv',sep=',',quote = '"',header = TRUE,stringsAsFactors=FALSE)
#scr1r$Details<-NULL

#scr1r[!scr1r$X0 %in% c('+++','++','+',''),'X0' ]<-NA
#scr1r[!scr1r$X1 %in% c('+++','++','+',''),'X1' ]<-NA
#scr1r[!scr1r$X2.5 %in% c('+++','++','+',''),'X2.5' ]<-NA
#scr1r[!scr1r$X5 %in% c('+++','++','+',''),'X5' ]<-NA
#scr1r[scr1r[,c('X0','X1','X2.5','X5')]=='+++',]
#scr1r$X0<-factor(scr1$X0,levels=c('+++','++','+',''),labels=c('9','6','3','0'))

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


allg<-unique(scr1$Gene)
concs<-c(0,1,2.5,5)

dupl<-as.factor(unique(scr1[which(duplicated(scr1$Gene)),]$Gene))
dupl<-dupl[!dupl%in%c('WT',NA)]
duplicates<-subset(subset(scr1[order(scr1$Gene),], Gene %in% dupl ),Gene!='WT')

scr1m<-melt(scr1r[,c('Gene','0','1','2.5','5')],id='Gene',variable.name = 'Measure',value.name='Score')
scr1m$Measure<-as.character(scr1m$Measure)
scr1m$Score<-as.character(scr1m$Score)



#Secondary screen
scr2<-read.table('Secondary_screen_PN_new.csv',sep=',',quote = '"',header = TRUE,stringsAsFactors = FALSE)
scr2$Starving<-as.factor(scr2$Starving)
scr2m<-melt(scr2[,c('Gene','MIC1','MIC2','MIC3')],id='Gene',variable.name = 'Replicate',value.name='MIC')
scr2m$Replicate<-as.character(scr2m$Replicate)
scr2m$MIC<-as.character(scr2m$MIC)



# scr2avg<-ddply(scr2m, .(Gene), summarise, MIC_avg=mean(MIC,na.rm = TRUE),MIC_sd=sd(MIC,na.rm = TRUE))
# scr2avg[match(scr2$Gene,scr2avg$Gene),'Starving']<-scr2$Starving



#Supplemented screen 3

scr3Scores<-read.table('3rd_screen_Scores.csv',sep=',',quote = '"',header = TRUE,stringsAsFactors = FALSE)
scr3Scores<-scr3Scores[,c('Gene',"X0.1","X0.2","X0.3","X1.1","X1.2","X1.3","X2.5.1","X2.5.2","X2.5.3","X5.1","X5.2","X5.3")]
scr3Sm<-melt(scr3MScores,id='Gene',variable.name = 'Measure',value.name='Score')
scr3Sm$Measure<-as.character(scr3Sm$Measure)
scr3Sm$Score<-as.character(scr3Sm$Score)
scr3Sm[scr3Sm$Measure %in% c('X0.1','X0.2','X0.3'),'Measure']<-'0'
scr3Sm[scr3Sm$Measure %in% c('X1.1','X1.2','X1.3'),'Measure']<-'1'
scr3Sm[scr3Sm$Measure %in% c('X2.5.1','X2.5.2','X2.5.3'),'Measure']<-'2.5'
scr3Sm[scr3Sm$Measure %in% c('X5.1','X5.2','X5.3'),'Measure']<-'5'
scr3Sm[scr3Sm$Score %in% c('3'),'Score']<-'9'
scr3Sm[scr3Sm$Score %in% c('2'),'Score']<-'6'
scr3Sm[scr3Sm$Score %in% c('1'),'Score']<-'3'


scr3Mics<-read.table('3rd_screen_MICs.csv',sep=',',quote = '"',header = TRUE,stringsAsFactors = FALSE)
scr3Mics<-scr3Mics[,c('Gene','MIC.1','MIC.2','MIC.3')]
scr3Micsm<-melt(scr3Mics,id='Gene',variable.name = 'Replicate',value.name='MIC')
scr3Micsm$Replicate<-as.character(scr3Micsm$Replicate)
scr3Micsm$MIC<-as.character(scr3Micsm$MIC)


scr3<-subset(scr3Micsm,!(is.na(as.numeric(scr3Micsm$MIC))))
scr3[scr3$Gene=='WT cont','Gene']<-'WT'
scr3[scr3$Gene=='upp cont','Gene']<-'upp'
scr3$MIC<-as.numeric(scr3$MIC)
scr3avg<-ddply(scr3, .(Gene), summarise, MIC_avg=mean(MIC,na.rm = TRUE),MIC_sd=sd(MIC,na.rm = TRUE))



allscores<-merge(scr1m,scr3Sm,all.x=TRUE,all.y=TRUE)
allmics<-merge(scr2m,scr3Micsm,all.x=TRUE,all.y=TRUE)

allscores<-subset(allscores,! is.na(Gene) & Gene!='' & Score!='')
allmics<-subset(allmics,! is.na(Gene) & Gene!='' & Score!='')

allscores[allscores$Gene=='WT cont','Gene']<-'WT'
allscores[allscores$Gene=='upp cont','Gene']<-'upp'

allmics[allmics$Gene=='WT cont','Gene']<-'WT'
allmics[allmics$Gene=='upp cont','Gene']<-'upp'

scores<-subset(allscores,!(is.na(as.numeric(allscores$Score))))
scores$Score<-as.numeric(scores$Score)

micss<-subset(allmics,!(is.na(as.numeric(allmics$MIC))))
micss$MIC<-as.numeric(micss$MIC)


scores_avg<-ddply(scores, .(Gene,Measure), summarise, Score_avg=mean(Score,na.rm = TRUE),Score_sd=sd(Score,na.rm = TRUE))

mics_avg<-ddply(micss, .(Gene), summarise, MIC_avg=mean(MIC,na.rm = TRUE),MIC_sd=sd(MIC,na.rm = TRUE))

all<-dcast(scores_avg,Gene ~Measure,mean,value.var = c('Score_avg'))


all$MIC<-apply(scoresc,1,evalmic)

allfull<-merge(all,mics_avg,by='Gene',all.x=TRUE)
allfull$MIC<-ifelse(allfull$MIC_avg>allfull$MIC,allfull$MIC_avg,allfull$MIC)




sc1g<-unique(scr1$Gene)
sc2g<-unique(scr2$Gene)
sc3Mg<-unique(scr3Mix$Gene)
sc3JMg<-unique(scr3JMics$Gene)




ones<-subset(scr1,(scr1$'X0'>0 | scr1$'X0'==0 | is.na(scr1$'X0')) & (scr1$'X1'==0 | is.na(scr1$'X1')) & (scr1$'X2.5'==0 | is.na(scr1$'X2.5')) & (scr1$'X5'==0 | is.na(scr1$'X5')) )$Gene #scr1$'X2.5'==0 & (scr1$'X2.5'==0 | is.na(scr1$'X2.5')) & (scr1$'X5'==0 | is.na(scr1$'X5')) 
ones<-unique(ones)

twofives<-subset(scr1,(scr1$'X0'>0 | is.na(scr1$'X0')) & (scr1$'X1'>0 | is.na(scr1$'X1')) & scr1$'X2.5'==0 &  scr1$'X5'==0)$Gene
twofives<-setdiff(unique(twofives),ones)#unique(twofives)#
#Check for duplicates

fivers<-subset(scr1,scr1$'X0'>0 & scr1$'X1'>0 & scr1$'X2.5'>0 &  scr1$'X5'==0)$Gene
fivers<-unique(fivers)#setdiff(unique(fivers),twofives)
#Check for duplicates


intersect(ones,twofives)
intersect(twofives,fivers)
intersect(ones,fivers)


scr1done<-union(ones,twofives)
scr2done<-unique(scr2$Gene)

missing<-setdiff(allg,scr1done)
missing<-setdiff(missing,scr2done)
missing<-setdiff(missing,fivers)


## 9 0 0 3 pattern
## Needs discussion
dubious1<-subset(scr1,Gene %in% missing & (scr1$'X0'==9 | is.na(scr1$'X0')) &scr1$'X1'==0 & scr1$'X2.5'==0 & scr1$'X5'==3)$Gene
## 9 0 3 0 pattern
## Needs discussion
dubious2<-subset(scr1,Gene %in% missing & (scr1$'X0'==9 | is.na(scr1$'X0')) &scr1$'X1'==0 & scr1$'X2.5'==3 & scr1$'X5'==0)$Gene
dubiousall<-union(dubious1,dubious2)
emissing<-setdiff(missing,dubiousall)
#Dubious retested at 1uM - 5uM and greater depending on the results
#Greater tested at greater than 5uM


greater<-subset(scr1,Gene %in% emissing & (scr1$'X0'>0 | is.na(scr1$'X0')) & scr1$'X2.5'>0 & ( scr1$'X5'>0 | is.na(scr1$'X5')))$Gene
mmissing<-setdiff(emissing,greater)

almost5<-subset(scr1,Gene %in% mmissing & scr1$'X0'>0 & scr1$'X1'>0 & scr1$'X2.5'==0 & scr1$'X5'>0)$Gene
mmissing<-setdiff(mmissing,almost5)



allredo<-subset(scr1,Gene %in% missing)
allredo$Details<-NULL
allredo$Faults<-NULL

allredo[allredo$Gene %in% almost5,'Reason' ]<-'Inconsistent for replicates'
allredo[allredo$Gene %in% mmissing,'Reason' ]<-'Inconsistent or missing'
allredo[allredo$Gene %in% dubiousall,'Reason' ]<-'Sporadic development'
allredo[allredo$Gene %in% fivers,'Reason' ]<-'Assumed MIC=5'
allredo[allredo$Gene %in% greater,'Reason' ]<-'Needs testing at [5FU]>5uM'
allredo[allredo$Gene=='yjgX','Reason' ]<-'Inconsistent results in replicates'

allredo[allredo$Gene %in% almost5,'Suggested concentration range' ]<-'1-5uM or 5-100uM'
allredo[allredo$Gene %in% mmissing,'Suggested concentration range' ]<-'1-5uM and possibly more'
allredo[allredo$Gene %in% dubiousall,'Suggested concentration range' ]<-'1-5uM and possibly more'
allredo[allredo$Gene %in% fivers,'Suggested concentration range' ]<-'>5uM'
allredo[allredo$Gene %in% greater,'Suggested concentration range' ]<-'5-100uM'
allredo[allredo$Gene=='yjgX','Suggested concentration range' ]<-'1-5uM and possibly more'

#What needs to be redone
#write.csv(subset(scr1,Gene %in% mmissing),'Missing_0-5.csv')
#write.csv(subset(scr1,Gene %in% dubiousall),'Sporadic_2.5-5.csv')
#write.csv(subset(scr1,Gene %in% greater),'Greater_5-100.csv')
#write.csv(subset(scr1,Gene %in% fivers),'Fivers_5-100.csv')


write.csv(allredo,'Data/All_redo_new.csv')


#mc<-data.frame('Gene'=unique(scr1$Gene))
#mc<-merge(mc,scr1[,c('Gene','Plate','Well','X0','X1','X2.5','X5')],by=c('Gene'),all.x=TRUE)


#Wrong gene indices!!
mics<-merge(scr1[,colnames(scr1)[c(1:6,12:14,11)]],scr2avg[,c('Gene','MIC_avg','MIC_sd','Starving')],by=c('Gene'),all.x=TRUE)
colnames(mics)<-c('Gene','Plate','Well','JW_id','ECK','bno','EG','GI','UniACC','Comment','MIC','MIC_sd','Starving')

subset(scr1,Gene %in% c('yedN','ubiX'))

mics[mics$Gene %in% setdiff(ones,c('yedN')),'MIC' ]<-1
mics[mics$Gene %in% twofives,'MIC' ]<-2.5
mics[mics$Gene %in% fivers,'MIC' ]<-5


#Uncomment if you want to avoid bias
#Will have to be changd with values from supplementary screen
mics[mics$Gene %in% greater,'MIC' ]<-10
mics$Starving<-as.character(mics$Starving)
mics$Starving<-ifelse(is.na(mics$Starving),'No',mics$Starving)
mics$Starving<-factor(mics$Starving,levels=c('Yes','No'),labels=c('Yes','No'))

#Add new data! Merge by gene name




#Bacterial growth
#NGM media start OD=0.057392708
bac<-read.table('Bacteria.csv',sep=',',quote = '"',header = TRUE,stringsAsFactors=FALSE)

keio<-read.table('Keio_growth.csv',sep=',',quote = '"',header = TRUE,stringsAsFactors=FALSE)
#keio<-keio[,colnames(keio) %in% c('Gene','LB_22hr','MOPS_24hr','MOPS_48hr')]
keio[is.na(keio)]<-NA
#keio<-subset(keio,!LB_22hr=='N.A.' & !MOPS_24hr=='N.A.' & !MOPS_48hr=='N.A.')
keio$LB_22hr<-as.numeric(as.character(keio$LB_22hr))
keio$MOPS_24hr<-as.numeric(as.character(keio$MOPS_24hr))
keio$MOPS_48hr<-as.numeric(as.character(keio$MOPS_48hr))

mics<-merge(mics,keio,by.x=c('JW_id','ECK','Gene'),by.y=c('JW.id','ECK.number','Gene'),all.x=TRUE)
mics[is.na(mics$LB_22hr),c('LB_22hr','MOPS_24hr','MOPS_48hr')]<-keio[match(subset(mics,is.na(LB_22hr))$Gene,keio$Gene),c('LB_22hr','MOPS_24hr','MOPS_48hr')]



bac$Lookup<-NULL
bac$Row<-NULL
bac$Col<-NULL

bac[bac$Gene=='BW',]$Gene<-'WT'

bacm<-melt(bac,id=colnames(bac)[1:4],variable.name = 'Replicate',value.name='OD')
bacm$OD<-bacm$OD-0.057392708
bacavg<-ddply(bacm, .(Gene,Drug), summarise, NGM=mean(OD,na.rm=TRUE),NGM_sd=sd(OD,na.rm=TRUE)) #Plate,Well,Drug,
bacavg<-subset(bacavg,! Gene %in% c('XXXXXXX','no bacteria','empty'))

bacall<-merge(subset(bacavg,Drug==0),subset(bacavg,Drug==100),by=c('Gene'),suffixes = c("_C","_D"))
bacall$Drug_C<-NULL
bacall$Drug_D<-NULL

bacall<-subset(bacall,! Gene %in% c('XXXXXXX','no bacteria','empty'))

bacmic<-merge(mics,bacall,id=c('Gene'),all.x = TRUE)


#Get only unique instances of knockouts
ubacmic<-bacmic[!rev(duplicated(rev(bacmic$Gene))),]
ubacmic<-ubacmic[,colnames(ubacmic)[c(1:3,6:20)]]
qbacmic<-subset(ubacmic,!is.na(Gene))

bcq05<-quantile(qbacmic$NGM_C,0.05,na.rm=TRUE)[[1]]
bcq95<-quantile(qbacmic$NGM_C,0.95,na.rm=TRUE)[[1]]
bdq05<-quantile(qbacmic$NGM_D,0.05,na.rm=TRUE)[[1]]
bdq95<-quantile(qbacmic$NGM_D,0.95,na.rm=TRUE)[[1]]
blq05<-quantile(qbacmic$LB_22hr,0.05,na.rm=TRUE)[[1]]
blq95<-quantile(qbacmic$LB_22hr,0.95,na.rm=TRUE)[[1]]

fitbac<-lm(NGM_D ~ NGM_C,data=qbacmic)
#confint(fitbac,'(Intercept)',level=0.95)[[2]]

#coefficients(fitbac)[[2]]

fitqr<-rq(NGM_D ~ NGM_C,data=qbacmic,tau=c(0.05,0.95))
bgli<-coefficients(fitqr)[1,][[1]]
bgui<-coefficients(fitqr)[1,][[2]]
bgls<-coefficients(fitqr)[2,][[1]]
bgus<-coefficients(fitqr)[2,][[2]]

#qbacmicq - no duplicates
qbacmicq<-subset(qbacmic,! Gene %in% c('atpB','atpE','atpF','atpG','aceE','lpd','glnA'))

bacres<-subset(qbacmicq,NGM_D>NGM_C*bgus+bgui)
bacsens<-subset(qbacmicq,NGM_D<NGM_C*bgls+bgli)


baccor<-ggplot(qbacmic,aes(x=NGM_C,y=NGM_D,color=MIC))+geom_point(size=1)+ylim(0, .25)#+xlim(0, .3)
baccor<-baccor+ylab('Knockout strain growth OD - 100uM 5FU')+xlab('Knockout strain growth OD - Control')+
  ggtitle('Growth of knockout strains in control and 100uM 5FU treatment')+
  stat_smooth(aes(group = 1),method = "lm")+
  geom_abline(intercept=0,slope=1,alpha=0.5,aes(color='grey'),linetype='longdash')+
  geom_text(aes(label=ifelse(NGM_D>NGM_C*bgus+bgui | NGM_D < NGM_C*bgls+bgli | NGM_C<0.03 ,Gene,'')),
            hjust=-0.1, vjust=-0.1,size=3)+
  geom_errorbarh(aes(xmax=NGM_C+NGM_sd_C,xmin=NGM_C-NGM_sd_C),height=.001,alpha=0.2)+
  geom_errorbar(aes(ymax=NGM_D+NGM_sd_D,ymin=NGM_D-NGM_sd_D),width=0.001,alpha=0.2)+
  geom_abline(intercept=bgli,slope=bgls,alpha=0.5,color='red')+
  geom_abline(intercept=bgui,slope=bgus,alpha=0.5,color='red')+
  annotate("text", 0.25,0.25*bgls+bgli+0.005, label = "5%",color='red')+
  annotate("text", 0.25,0.25*bgus+bgui+0.005, label = "95%",color='red')+
  scale_x_continuous(breaks=seq(0,.3,by=.05))+
  labs(color='MIC [5FU], uM')+
  annotate('text',x = 0.125, y = 0.25, label = lm_eqn(fitbac), parse = TRUE)
baccor
dev.copy2pdf(device=cairo_pdf,file="Figures/Control-Treatment_NGM_growth.pdf",width=9,height=9)



#baccor<-baccor+stat_quantile(quantiles=c(.05,.95),alpha=0.5)
#baccor<-baccor+annotate("text", 0.3, bdq05+0.005, label = "5%",color='blue')+annotate("text", 0.3, bdq95+0.005, label = "95%",color='blue')
#baccor<-baccor+annotate("text", bcq05+0.01, 0.25, label = "5%",color='red')+annotate("text", bcq95+0.01, 0.25, label = "95%",color='red')
#baccor<-baccor+geom_hline(yintercept=bdq05,color='blue',alpha=0.5,linetype='longdash')+geom_hline(yintercept=bdq95,color='blue',alpha=0.5,linetype='longdash')
#baccor<-baccor+geom_vline(xintercept=bcq05,color='red',alpha=0.5,linetype='longdash')+geom_vline(xintercept=bcq95,color='red',alpha=0.5,linetype='longdash')


#baccor+stat_ellipse(level=0.95)



bacmed<-melt(qbacmic[,colnames(qbacmic)[c(1,12:15,17)]],id=c('Gene'),variable.name = 'Media',value.name='OD')
bacmed$Media<-factor(bacmed$Media,levels = c('NGM_C','NGM_D','LB_22hr','MOPS_24hr','MOPS_48hr'),labels=c('NGM - 24h','NGM + 100uM 5FU','LB - 22hr','MOPS - 24hr','MOPS - 48hr'))

bacmed<-subset(bacmed,!Media %in% c('MOPS - 48hr')) #, 'MOPS - 48hr'

bachist<-ggplot(bacmed,aes(x=OD,fill=Media))+geom_histogram(aes(y=0.01*..density..),position='identity',alpha=0.5,binwidth = 0.01)
bachist<-bachist+labs(fill='Media')+xlab('OD')+ylab('')+scale_y_continuous(limits=c(0,0.10), labels = scales::percent)
bachist+ggtitle('Distribution of strain growth')
dev.copy2pdf(device=cairo_pdf,file="Figures/Bac_growth_disribution.pdf",width=9,height=9)

#bacmic$MIC_avg<-bacmic$MIC
#bacmic<-subset(bacmic,! Gene %in% c('XXXXXXX','no bacteria','empty'))




#qbacmicq<-qbacmic




bdqq05<-quantile(qbacmicq$NGM_D,0.05,na.rm=TRUE)[[1]]
bdqq95<-quantile(qbacmicq$NGM_D,0.95,na.rm=TRUE)[[1]]
bcqq05<-quantile(qbacmicq$NGM_C,0.05,na.rm=TRUE)[[1]]
bcqq95<-quantile(qbacmicq$NGM_C,0.95,na.rm=TRUE)[[1]]
blqq05<-quantile(qbacmicq$LB_22hr,0.05,na.rm=TRUE)[[1]]
blqq95<-quantile(qbacmicq$LB_22hr,0.95,na.rm=TRUE)[[1]]

fitD<-lm(NGM_D ~ MIC,qbacmicq)
fitNDqr<-rq(NGM_D ~ MIC,data=qbacmicq,tau=c(0.05,0.95))
mndli<-coefficients(fitNDqr)[1,][[1]]
mndui<-coefficients(fitNDqr)[1,][[2]]
mndls<-coefficients(fitNDqr)[2,][[1]]
mndus<-coefficients(fitNDqr)[2,][[2]]


mbcD<-ggplot(qbacmicq,aes(x=MIC,y=NGM_D))+geom_point(size=1)+
  stat_smooth(aes(group = 1),method = "lm")+
  geom_errorbarh(aes(xmax=MIC+MIC_sd,xmin=MIC-MIC_sd),height=.001,alpha=0.2,color='black')+
  geom_errorbar(aes(ymax=NGM_D+NGM_sd_D,ymin=NGM_D-NGM_sd_D),width=0.001,alpha=0.2,color='black')+
  #geom_text(aes(label=ifelse(MIC>5,Gene,'')),hjust=-0.1, vjust=-0.1,size=3,color='black')+
  geom_abline(intercept=mndli,slope=mndls,alpha=0.5,color='red')+
  geom_abline(intercept=mndui,slope=mndus,alpha=0.5,color='red')+
  annotate("text", 100, mndls*100+mndli+0.005, label = "5%",color='red')+
  annotate("text", 100, mndus*100+mndui+0.005, label = "95%",color='red')+
  ggtitle('Strain growth in NGM 24hr OD - 100uM 5FU')+xlab('MIC [5FU], uM')+
  ylab('OD')+xlim(0,100)+ylim(0,0.25)+
  annotate('text',x = 75, y = 0.25, label = lm_eqn(fitD), parse = TRUE)
mbcD
#dev.copy2pdf(device=cairo_pdf,file="Figures/MIC-NGMTreatment_bac_growth_Starving-color.pdf",width=9,height=9)
#dev.copy2pdf(device=cairo_pdf,file="Figures/MIC-NGMTreatment_bac_growth.pdf",width=9,height=9)
dev.copy2pdf(device=cairo_pdf,file="Figures/MIC-NGMTreatment_bac_growth_NoLabels.pdf",width=9,height=9)

fitC<-lm(NGM_C ~ MIC,qbacmicq)
fitNCqr<-rq(NGM_C ~ MIC,data=qbacmicq,tau=c(0.05,0.95))
mncli<-coefficients(fitNCqr)[1,][[1]]
mncui<-coefficients(fitNCqr)[1,][[2]]
mncls<-coefficients(fitNCqr)[2,][[1]]
mncus<-coefficients(fitNCqr)[2,][[2]]

mbcC<-ggplot(qbacmicq,aes(x=MIC,y=NGM_C))+geom_point(size=1)+stat_smooth(aes(group = 1),method = "lm")+
  xlim(0,100)+ylim(0,0.3)+
  geom_errorbarh(aes(xmax=MIC+MIC_sd,xmin=MIC-MIC_sd),height=.001,alpha=0.2,color='black')+
  geom_errorbar(aes(ymax=NGM_C+NGM_sd_C,ymin=NGM_C-NGM_sd_C),width=0.001,alpha=0.2,color='black')+
  #geom_text(aes(label=Gene),hjust=-0.1, vjust=-0.1,size=3,color='black')+
  geom_abline(intercept=mncli,slope=mncls,alpha=0.5,color='red')+
  geom_abline(intercept=mncui,slope=mncus,alpha=0.5,color='red')+
  annotate("text", 100, mncls*100+mncli+0.005, label = "5%",color='red')+
  annotate("text", 100, mncus*100+mncui+0.005, label = "95%",color='red')+
  ggtitle('Strain growth in NGM 24hr OD - Control')+xlab('MIC [5FU], uM')+ylab('OD')+
  annotate('text',x = 75, y = 0.3, label = lm_eqn(fitC), parse = TRUE)
mbcC
#dev.copy2pdf(device=cairo_pdf,file="Figures/MIC-NGMControl_bac_growth_Starving-color.pdf",width=9,height=9)
#dev.copy2pdf(device=cairo_pdf,file="Figures/MIC-NGMControl_bac_growth.pdf",width=9,height=9)
dev.copy2pdf(device=cairo_pdf,file="Figures/MIC-NGMControl_bac_growth_NoLabels.pdf",width=9,height=9)

fitLB <- lm(LB_22hr ~ MIC, data=qbacmicq)
fitLBqr<-rq(LB_22hr ~ MIC,data=qbacmicq,tau=c(0.05,0.95))
mlbli<-coefficients(fitLBqr)[1,][[1]]
mlbui<-coefficients(fitLBqr)[1,][[2]]
mlbls<-coefficients(fitLBqr)[2,][[1]]
mlbus<-coefficients(fitLBqr)[2,][[2]]

mbcLB<-ggplot(qbacmicq,aes(x=MIC,y=LB_22hr))+geom_point(size=1)+
  stat_smooth(aes(group = 1),method = "lm")+xlim(0,100)+
  geom_errorbarh(aes(xmax=MIC+MIC_sd,xmin=MIC-MIC_sd),height=.001,alpha=0.2,color='black')+
  #geom_text(aes(label=ifelse(MIC>5,Gene,'')),hjust=-0.1, vjust=-0.1,size=3,color='black')+
  geom_abline(intercept=mlbli,slope=mlbls,alpha=0.5,color='red')+
  geom_abline(intercept=mlbui,slope=mlbus,alpha=0.5,color='red')+
  annotate("text", 100, mlbls*100+mlbli+0.02, label = "5%",color='red')+
  annotate("text", 100, mlbus*100+mlbui+0.02, label = "95%",color='red')
mbcLB<-mbcLB+ggtitle('Strain growth in LB 22hr OD - Control')+xlab('MIC [5FU], uM')+ylab('OD')
mbcLB+annotate('text',x = 75, y = 1.1, label = lm_eqn(fitLB), parse = TRUE)
#dev.copy2pdf(device=cairo_pdf,file="Figures/MIC-LB22hr_bac_growth_Starving-color.pdf",width=9,height=9)
#mbcLB<-mbcLB+geom_hline(yintercept=blqq05,color='red',alpha=0.5,linetype='longdash')+geom_hline(yintercept=blqq95,color='red',alpha=0.5,linetype='longdash')
#dev.copy2pdf(device=cairo_pdf,file="Figures/MIC-LB22hr_bac_growth.pdf",width=9,height=9)
dev.copy2pdf(device=cairo_pdf,file="Figures/MIC-LB22hr_bac_growth_NoLabels.pdf",width=9,height=9)

mbcMOPS<-ggplot(qbacmic,aes(x=MIC,y=MOPS_24hr,color=Starving))+geom_point(size=2)+stat_smooth(aes(group = 1),method = "lm")+xlim(0,100)
mbcMOPS<-mbcMOPS+geom_errorbarh(aes(xmax=MIC+MIC_sd,xmin=MIC-MIC_sd),height=.001,alpha=0.2,color='black')#+geom_errorbar(aes(ymax=NGM_C+NGM_sd_C,ymin=NGM_C-NGM_sd_C),width=0.001,alpha=0.2)
mbcMOPS<-mbcMOPS+geom_text(aes(label=Gene),hjust=-0.1, vjust=-0.1,size=3,color='black')
mbcMOPS+ggtitle('Strain growth in MOPS 48hr OD - Control')+xlab('MIC [5FU], uM')+ylab('OD')



alldist<-ggplot(qbacmic,aes(x=MIC,y=reorder(Gene,MIC),color=Starving))+geom_point(size=1) + geom_errorbarh(aes(xmax=MIC+MIC_sd,xmin=MIC-MIC_sd))
alldist<-alldist+theme(axis.text.y = element_text(vjust = 0,size=5))+scale_x_continuous(breaks=seq(0,100,by=10))
alldist+ylab('Gene knockout')+xlab('MIC [5FU], uM')+ggtitle('Protective properties of gene knockouts for C. elegans in 5FU exposure')
#dev.copy2pdf(device=cairo_pdf,file="All_MIC_variation_SD.pdf",width=8,height=40)



#limits <- aes(xmax = MIC_avg + MIC_sd, xmin=MIC_avg - MIC_sd)
#dodge <- position_dodge(width=0)
scr2dist<-ggplot(subset(qbacmicq,MIC>2.5),aes(x=MIC,y=reorder(Gene,MIC)))+geom_point(color='red',size=1) + geom_errorbarh(aes(xmax=MIC+MIC_sd,xmin=MIC-MIC_sd))
scr2dist<-scr2dist+theme(axis.text.y = element_text(vjust = 0,size=5))+scale_x_continuous(breaks=seq(0,100,by=10))
scr2dist+ylab('Gene knockout')+xlab('MIC [5FU], uM')+ggtitle('Protective properties of gene knockouts for C. elegans in 5FU exposure')
dev.copy2pdf(device=cairo_pdf,file="Figures/Secondary-screen_MIC_variation_SD.pdf",width=8,height=20)




######
q90<-quantile(qbacmic$MIC,0.9,na.rm=TRUE)[[1]]
q95<-quantile(qbacmic$MIC,0.95,na.rm=TRUE)[[1]]
q99<-quantile(qbacmic$MIC,0.99,na.rm=TRUE)[[1]]

dist<-ggplot(qbacmic,aes(x=MIC))+stat_ecdf()+ggtitle('Cumulative distribution of MIC values')
dist<-dist+geom_hline(yintercept=0.90,color='green',alpha=0.5,linetype='longdash')+geom_vline(xintercept=q90,color='green',alpha=0.5,linetype='longdash')
dist<-dist+geom_hline(yintercept=0.95,color='blue',alpha=0.5,linetype='longdash')+geom_vline(xintercept=q95,color='blue',alpha=0.5,linetype='longdash')
dist<-dist+geom_hline(yintercept=0.99,color='red',alpha=0.5,linetype='longdash')+geom_vline(xintercept=q99,color='red',alpha=0.5,linetype='longdash')
dist<-dist+annotate("text", 180, 0.91, label = "90%",color='green')
dist<-dist+annotate("text", 180, 0.96, label = "95%",color='blue')
dist<-dist+annotate("text", 180, 1, label = "99%",color='red')
dist<-dist+scale_y_continuous(limits=c(0,1), labels = scales::percent,breaks=seq(0,1,by=0.1))
#dist<-dist+scale_x_continuous(breaks=seq(0,100,by=5))
dist<-dist+scale_x_log10(breaks=c(0,2.5,5,10,20,30,40,50,100))
dist+ylab('')+xlab('log10(MIC [5FU], uM)')+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.copy2pdf(device=cairo_pdf,file="Figures/Cumulative_distribution_of_MIC_adjusted-x-scale.pdf",width=9,height=9)


disth<-ggplot(qbacmic,aes(x=MIC))+geom_histogram(binwidth=5)
disth<-disth+ggtitle('Histogram of MIC values')
#disth<-disth+scale_x_log10(breaks=seq(0,100,by=10))
disth

#Write selected knockouts

write.csv(bacmic,'Data/MICs_and_bacterial_growth-All.csv')
write.csv(qbacmic,'Data/MICs_and_bacterial_growth-Unique.csv')
write.csv(qbacmicq,'Data/MICs_and_bacterial_growth-Unique_clean.csv')

bacmicsen1<-subset(qbacmicq,MIC>1)
bacmicsen25<-subset(qbacmicq,MIC>2.5)
bacmicsen5<-subset(qbacmicq,MIC>5)
bacmicsen10<-subset(qbacmicq,MIC>10)

PT_clean<-subset(qbacmicq,MIC>2.5 & !is.na(Gene))[,c('Gene','MIC')]
PT_all<-subset(qbacmicq, !is.na(Gene) & !is.na(MIC) )[,c('Gene','MIC')]
write.table(PT_clean,file='Data/Gene_MIC_above2.5_ForPT.tab',sep='\t',row.names = FALSE,quote=FALSE)
write.table(PT_all,file='Data_Gene_MIC_all_ForPT.tab',sep='\t',row.names = FALSE,quote=FALSE)

write.csv(bacmicsen1,'Data/MICS_over1.csv')
write.csv(bacmicsen25,'Data/MICS_over25.csv')
write.csv(bacmicsen5,'Data/MICS_over5.csv')
write.csv(bacmicsen10,'Data/MICS_over10.csv')
