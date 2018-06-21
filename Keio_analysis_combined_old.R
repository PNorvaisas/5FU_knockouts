library('ggplot2')
library('gplots')
library('plyr')
library('reshape2')
library(tidyr)
library(quantreg)
library(ellipse)


#Vennerable installation: install.packages("Vennerable", repos="http://R-Forge.R-project.org")

#library(quantreg)

elipsoid=function(df,xvar,yvar,scale=1,groups=''){
  df<-subset(df,!is.na(df[,xvar]) &!is.na(df[,yvar]))
  df_ell <- data.frame()
  if (groups!=''){
    for(g in levels(df[,groups])){
      df_ell <- rbind(df_ell, cbind(as.data.frame(with(df[df$groups==g,],ellipse(cor(xvar, yvar),scale=c(sd(xvar),sd(yvar)),centre=c(mean(x),mean(y))))),
                                    group=g))
    }
  }else {
      df_ell <- as.data.frame( ellipse( cor(df[,c(xvar,yvar)],use='complete.obs'),scale=c( sd(df[,xvar],na.rm=TRUE)*scale ,sd(df[,yvar],na.rm=TRUE)*scale ),centre=c( mean(df[,xvar],na.rm=TRUE),mean(df[,yvar],na.rm=TRUE) )))
  }
  return(df_ell)
}


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



#Supplemented screen 3
scr3fix<-read.table('3rd_screen_location_fix.csv',sep=',',quote = '"',header = TRUE,stringsAsFactors = FALSE)
scr3Scores<-read.table('3rd_screen_Scores.csv',sep=',',quote = '"',header = TRUE,stringsAsFactors = FALSE)
scr3Scores<-rename(scr3Scores, c("Keio.Plate.no."="Plate", "Position"="Well", "Well"="NWell"))

#Fix mistakes with gene-location missmatches
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
scr3Mics[!is.na(scr3Mics$Right_gene),c('Gene')]<-scr3Mics[!is.na(scr3Mics$Right_gene),c('Right_gene')]
scr3Mics$Right_gene<-NULL
#Fix mistakes with gene-location missmatches
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

allscores<-subset(allscores,! is.na(Gene) & Gene!='' & Score!='')
allmics<-subset(allmics,! is.na(Gene) & Gene!='' & MIC!='')

allscores[allscores$Gene=='WT cont','Gene']<-'WT'
allscores[allscores$Gene=='upp cont','Gene']<-'upp'

allmics[allmics$Gene=='WT cont','Gene']<-'WT'
allmics[allmics$Gene=='upp cont','Gene']<-'upp'

allscores[allscores$Gene=='WT','Plate']<-''
allscores[allscores$Gene=='WT','Well']<-''
allmics[allmics$Gene=='WT','Plate']<-''
allmics[allmics$Gene=='WT','Well']<-''




#Real duplicates!
rdupl<-as.factor(unique(keioinfo[which(duplicated(keioinfo$Gene)),]$Gene))


allmics<-merge(allmics,subset(keioinfo,!Gene %in% rdupl)[,c('Gene','Plate','Well')],by='Gene',all.x=TRUE)

allmics[!is.na(allmics$Plate.y) &!is.na(allmics$Well.y),c('Plate.x','Well.x')]<-allmics[!is.na(allmics$Plate.y) &!is.na(allmics$Well.y),c('Plate.y','Well.y')]
allmics<-rename(allmics,c('Plate.x'='Plate','Well.x'='Well'))
allmics$Plate.y<-NULL
allmics$Well.y<-NULL
#Real duplicates in mics
rduplm<-c('dcuC','yedN','yhcE')


scores<-subset(allscores,!(is.na(as.numeric(allscores$Score))))
scores$Score<-as.numeric(scores$Score)

micss<-subset(allmics,!(is.na(as.numeric(allmics$MIC))))
micss$MIC<-as.numeric(micss$MIC)

#Calculate averages for Scores and MICs over replicates
scores_avg<-ddply(scores, .(Gene,Plate,Well,Measure), summarise, Score_avg=mean(Score,na.rm = TRUE),Score_sd=sd(Score,na.rm = TRUE))

#Averaging by gene names, until duplicates are sorted out
mics_avg<-ddply(micss, .(Gene,Plate,Well), summarise, MIC_avg=mean(MIC,na.rm = TRUE),MIC_sd=sd(MIC,na.rm = TRUE))

#duplm should be dcuC and yhcE
duplm<-as.factor(unique(mics_avg[which(duplicated(mics_avg$Gene)),]$Gene))

alls<-dcast(scores_avg,Gene+Plate+Well ~Measure,mean,value.var = c('Score_avg'))


#Evaluate MICS
alls$MIC<-evalmic2(alls)$MIC


#Merging by gene names, until duplicates are sorted out   ,'Plate','Well'
allfull<-merge(alls,mics_avg,by=c('Gene','Plate','Well'),all.x=TRUE,all.y=TRUE)


allfull$MIC<-ifelse(!is.na(allfull$MIC_avg),allfull$MIC_avg,allfull$MIC)
allfull$MIC_avg<-NULL

#Manual fixes
allfull[allfull$Gene=='WT','MIC']<-1
allfull[allfull$Gene=='WT','MIC_sd']<-0

allfull[allfull$Gene=='upp','MIC']<-15
allfull[allfull$Gene=='upp','MIC_sd']<-0


#Find duplicates
dupl<-as.factor(unique(allfull[which(duplicated(allfull$Gene)),]$Gene))




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
allinf<-allinfr[,! colnames(allinfr) %in% c('Gene.y','Row','Column','Comment')]
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


#Bacterial growth
#NGM media start OD=0.057392708
bac<-read.table('Bacteria.csv',sep=',',quote = '"',header = TRUE,stringsAsFactors=FALSE)

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

#Data fully merged!!



#Get only unique instances of knockouts
#How to deal with duplicated entries
#ubacmic<-bacmic[!rev(duplicated(rev(bacmic$Gene))),]
#ubacmic<-ubacmic[,colnames(ubacmic)[c(1:3,6:20)]]
#subset(ubacmic,!is.na(Gene))



#Output folder:
odir<-'Figures_v2'
ddir<-'Data_v2'

#qbacmicq - no outliers
qbacmicq<-subset(bacmic,! Gene %in% c('glnA','aceE','atpB','atpG','atpE','atpF','lpd'))


bacmic<-
qbacmic<-

#
bcq05<-quantile(bacmic$NGM_C,0.05,na.rm=TRUE)[[1]]
bcq95<-quantile(bacmic$NGM_C,0.95,na.rm=TRUE)[[1]]
bdq05<-quantile(bacmic$NGM_D,0.05,na.rm=TRUE)[[1]]
bdq95<-quantile(bacmic$NGM_D,0.95,na.rm=TRUE)[[1]]
blq05<-quantile(bacmic$LB_22hr,0.05,na.rm=TRUE)[[1]]
blq95<-quantile(bacmic$LB_22hr,0.95,na.rm=TRUE)[[1]]

fitbac<-lm(NGM_D ~ NGM_C,data=bacmic)
#confint(fitbac,'(Intercept)',level=0.95)[[2]]
#coefficients(fitbac)[[2]]

fitqr<-rq(NGM_D ~ NGM_C,data=bacmic,tau=c(0.05,0.95))
bgli<-coefficients(fitqr)[1,][[1]]
bgui<-coefficients(fitqr)[1,][[2]]
bgls<-coefficients(fitqr)[2,][[1]]
bgus<-coefficients(fitqr)[2,][[2]]


bacres<-subset(bacmic,NGM_D>NGM_C*bgus+bgui)
bacsens<-subset(bacmic,NGM_D<NGM_C*bgls+bgli)


theme_set(theme_light())

baccor<-ggplot(bacmic,aes(x=NGM_C,y=NGM_D,color=MIC))+
  geom_point(size=1)+ylim(0, .25)+
  ylab(expression(paste('Knockout strain growth OD - 100',mu,'M 5FU')))+
  xlab('Knockout strain growth OD - Control')+
  ggtitle(expression(paste('Growth of knockout strains in control and 100',mu,'M 5FU treatment')))+
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
  labs(color=expression(paste('MIC [5FU], ',mu,'M')))+
  annotate('text',x = 0.125, y = 0.25, label = lm_eqn(fitbac), parse = TRUE)
baccor
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Control-Treatment_NGM_growth.pdf",sep = ''),width=9,height=9)



bacmed<-melt(bacmic[,colnames(bacmic) %in% c('Gene','NGM_C','NGM_D','LB_22hr','MOPS_24hr','MOPS_48hr')],
             id=c('Gene'),variable.name = 'Media',value.name='OD')
bacmed$Media<-factor(bacmed$Media,levels = c('NGM_C','NGM_D','LB_22hr','MOPS_24hr','MOPS_48hr'),
                     labels=c('NGM - 24h','NGM + 100uM 5FU','LB - 22hr','MOPS - 24hr','MOPS - 48hr'))
bacmed<-subset(bacmed,!Media %in% c('MOPS - 48hr')) #, 'MOPS - 48hr'

bachist<-ggplot(bacmed,aes(x=OD,fill=Media))+
  geom_histogram(aes(y=0.01*..density..),position='identity',alpha=0.5,binwidth = 0.01)+
  labs(fill='Media')+xlab('OD')+ylab('')+
  scale_y_continuous(limits=c(0,0.10), labels = scales::percent)+
  ggtitle('Distribution of strain growth')
bachist
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Bac_growth_disribution.pdf",sep=''),width=9,height=9)


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

df_el<-elipsoid(subset(qbacmicq,!is.na(MIC) & ! is.na(NGM_D)),'MIC','NGM_D')

mbcD<-ggplot(qbacmicq,aes(x=MIC,y=NGM_D))+geom_point(size=1)+
  stat_smooth(aes(group = 1),method = "lm")+
  geom_errorbarh(aes(xmax=MIC+MIC_sd,xmin=MIC-MIC_sd),height=.0005,alpha=0.2,color='black')+
  geom_errorbar(aes(ymax=NGM_D+NGM_sd_D,ymin=NGM_D-NGM_sd_D),width=0.0005,alpha=0.2,color='black')+
  geom_text(aes(label=ifelse((NGM_D>MIC*mndus+mndui | NGM_D < MIC*mndls+mndli) | MIC>50 ,Gene,'')),
            hjust=-0.1, vjust=-0.1,size=2)+
  geom_abline(intercept=mndli,slope=mndls,alpha=0.5,color='red')+
  geom_abline(intercept=mndui,slope=mndus,alpha=0.5,color='red')+
  annotate("text", 100, mndls*100+mndli+0.005, label = "5%",color='red')+
  annotate("text", 100, mndus*100+mndui+0.005, label = "95%",color='red')+
  ggtitle(expression(paste('Strain growth in NGM 24hr OD - 100',mu,'M 5FU')))+xlab(expression(paste('MIC [5FU], ',mu,'M')))+
  ylab('OD')+xlim(0,100)+ylim(0,0.25)+
  annotate('text',x = 50, y = 0.25, label = lm_eqn(fitD), parse = TRUE)#+
  #geom_path(data=df_el, aes(x=MIC, y=NGM_D), size=1, linetype=1,color='grey',alpha=0.5)
  #stat_density2d()
mbcD
#dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/MIC-NGMTreatment_bac_growth_Starving-color.pdf',sep=''),width=9,height=9)
#dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/MIC-NGMTreatment_bac_growth.pdf",sep=''),width=9,height=9)
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/MIC-NGMTreatment_bac_growth_NoLabels.pdf",sep=''),width=5,height=5)

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
  geom_text(aes(label=ifelse((NGM_C>MIC*mncus+mncui | NGM_C < MIC*mncls+mncli) | MIC>45,Gene,'')),
            hjust=-0.1, vjust=-0.1,size=2)+
  geom_abline(intercept=mncli,slope=mncls,alpha=0.5,color='red')+
  geom_abline(intercept=mncui,slope=mncus,alpha=0.5,color='red')+
  annotate("text", 100, mncls*100+mncli+0.005, label = "5%",color='red')+
  annotate("text", 100, mncus*100+mncui+0.005, label = "95%",color='red')+
  ggtitle('Strain growth in NGM 24hr OD - Control')+xlab(expression(paste('MIC [5FU], ',mu,'M')))+ylab('OD')+
  annotate('text',x = 50, y = 0.3, label = lm_eqn(fitC), parse = TRUE)
mbcC
#dev.copy2pdf(device=cairo_pdf,file="Figures/MIC-NGMControl_bac_growth_Starving-color.pdf",width=9,height=9)
#dev.copy2pdf(device=cairo_pdf,file="Figures/MIC-NGMControl_bac_growth.pdf",width=9,height=9)
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/MIC-NGMControl_bac_growth_NoLabels.pdf",sep=''),width=6,height=6)

fitLB <- lm(LB_22hr ~ MIC, data=qbacmicq)
fitLBqr<-rq(LB_22hr ~ MIC,data=qbacmicq,tau=c(0.05,0.95))
mlbli<-coefficients(fitLBqr)[1,][[1]]
mlbui<-coefficients(fitLBqr)[1,][[2]]
mlbls<-coefficients(fitLBqr)[2,][[1]]
mlbus<-coefficients(fitLBqr)[2,][[2]]

mbcLB<-ggplot(qbacmicq,aes(x=MIC,y=LB_22hr))+geom_point(size=1)+
  stat_smooth(aes(group = 1),method = "lm")+xlim(0,100)+
  geom_errorbarh(aes(xmax=MIC+MIC_sd,xmin=MIC-MIC_sd),height=.001,alpha=0.2,color='black')+
  geom_text(aes(label=ifelse(((LB_22hr>MIC*mlbus+mlbui | LB_22hr < MIC*mlbls+mlbli) & MIC >1) | MIC>40,Gene,'')),
            hjust=-0.1, vjust=-0.1,size=3)+
  geom_abline(intercept=mlbli,slope=mlbls,alpha=0.5,color='red')+
  geom_abline(intercept=mlbui,slope=mlbus,alpha=0.5,color='red')+
  annotate("text", 100, mlbls*100+mlbli+0.02, label = "5%",color='red')+
  annotate("text", 100, mlbus*100+mlbui+0.02, label = "95%",color='red')+
  ggtitle('Strain growth in LB 22hr OD - Control')+xlab(expression(paste('MIC [5FU], ',mu,'M')))+ylab('OD')+
  annotate('text',x = 50, y = 1.1, label = lm_eqn(fitLB), parse = TRUE)
mbcLB
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/MIC-LB22hr_bac_growth_NoLabels.pdf",sep=''),width=9,height=9)



#All MICs
alldist<-ggplot(qbacmicq,aes(x=MIC,y=reorder(Gene,MIC,max)))+
  geom_point(color='red',size=1) + geom_errorbarh(aes(xmax=MIC+MIC_sd,xmin=MIC-MIC_sd))+
  theme(axis.text.y = element_text(vjust = 0,size=4))+
  scale_x_continuous(breaks=seq(0,100,by=10))+
  ylab('Gene knockout')+
  xlab(expression(paste('MIC [5FU], ',mu,'M')))+
  ggtitle('Protective properties of gene knockouts for C. elegans in 5FU exposure')
alldist
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/MIC_variation_SD_all.pdf",sep=''),width=8,height=150)

#MICs over 1
scr2dist<-ggplot(subset(qbacmicq,MIC>1),aes(x=MIC,y=reorder(Gene,MIC,max)))+
  geom_point(color='red',size=1) + geom_errorbarh(aes(xmax=MIC+MIC_sd,xmin=MIC-MIC_sd))+
  theme(axis.text.y = element_text(vjust = 0,size=4))+
  scale_x_continuous(breaks=seq(0,100,by=10))+
  ylab('Gene knockout')+
  xlab(expression(paste('MIC [5FU], ',mu,'M')))+
  ggtitle('Protective properties of gene knockouts for C. elegans in 5FU exposure')
scr2dist
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/MIC_variation_SD_MIC-over-1.pdf",sep=''),width=8,height=30)



######
q90<-quantile(qbacmicq$MIC,0.9,na.rm=TRUE)[[1]]
q95<-quantile(qbacmicq$MIC,0.95,na.rm=TRUE)[[1]]
q99<-quantile(qbacmicq$MIC,0.99,na.rm=TRUE)[[1]]

dist<-ggplot(qbacmicq,aes(x=MIC))+stat_ecdf()+ggtitle('Cumulative distribution of MIC values')+
  geom_hline(yintercept=0.90,color='green',alpha=0.5,linetype='longdash')+
  geom_vline(xintercept=q90,color='green',alpha=0.5,linetype='longdash')+
  geom_hline(yintercept=0.95,color='blue',alpha=0.5,linetype='longdash')+
  geom_vline(xintercept=q95,color='blue',alpha=0.5,linetype='longdash')+
  geom_hline(yintercept=0.99,color='red',alpha=0.5,linetype='longdash')+
  geom_vline(xintercept=q99,color='red',alpha=0.5,linetype='longdash')+
  annotate("text", 1, 0.92, label = "90%",color='green')+
  annotate("text", 1, 0.97, label = "95%",color='blue')+
  annotate("text", 1, 1, label = "99%",color='red')+
  scale_y_continuous(limits=c(0,1), labels = scales::percent,breaks=seq(0,1,by=0.1))+
  scale_x_log10(breaks=c(0,2.5,5,10,20,30,40,50,75,100))+
  ylab('')+xlab(expression(paste('MIC [5FU], ',mu,'M (log10 scaled)')))+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dist
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Cumulative_distribution_of_MIC_log10-x-scale.pdf",sep=''),width=9,height=9)
#dist<-dist+scale_x_continuous(breaks=seq(0,100,by=5))





write.csv(bacmic,paste(ddir,'/MICs_and_bacterial_growth-All.csv',sep=''))
write.csv(qbacmicq,paste(ddir,'/MICs_and_bacterial_growth-Clean.csv',sep=''))

# 
# #PT lists
# PT_clean<-subset(qbacmicq,MIC>2.5 & !is.na(Gene))[,c('Gene','MIC')]
# PT_all<-subset(qbacmicq, !is.na(Gene) & !is.na(MIC) )[,c('Gene','MIC')]
# write.table(PT_clean,file='Data/Gene_MIC_above2.5_ForPT.tab',sep='\t',row.names = FALSE,quote=FALSE)
# write.table(PT_all,file='Data_Gene_MIC_all_ForPT.tab',sep='\t',row.names = FALSE,quote=FALSE)


