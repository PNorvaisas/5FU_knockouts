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

trim <- function (x) gsub("^\\s+|\\s+$", "", x)


elipsoid=function(df,xvar,yvar,scale=1,groups=''){
  df<-subset(df,!is.na(df[,xvar]) &!is.na(df[,yvar]))
  df_ell <- data.frame()
  #This part not working
  if (groups!=''){
    for(g in levels(df[,groups])){
      df_ell <- rbind(df_ell,
                      cbind(as.data.frame( with(df[df$groups==g,],ellipse(cor(xvar, yvar),scale=c(sd(xvar),sd(yvar)),centre=c(mean(x),mean(y))))),group=g))
    }
  }else {
    cval<-cor(df[,c(xvar,yvar)],use='complete.obs')
    #cval[c(2,3)]<-cval[c(2,3)]*sd(df[,yvar],na.rm=TRUE)/sd(df[,xvar],na.rm=TRUE)
    df_ell <- as.data.frame( ellipse( cval,
                                      scale=c( sd(df[,xvar],na.rm=TRUE)*scale ,sd(df[,yvar],na.rm=TRUE)*scale ),
                                      centre=c( mean(df[,xvar],na.rm=TRUE),mean(df[,yvar],na.rm=TRUE) )))
  }
  return(df_ell)
}


lm_eqn = function(m) {
  fres<-summary(m)
  l <- list(a = format(abs(coef(m)[2]), digits = 2),
            b = format(coef(m)[1], digits = 2),
            r2 = format(summary(m)$r.squared, digits = 2),
            p2 = format(pf(fres$fstatistic[1], fres$fstatistic[2], fres$fstatistic[3],lower.tail = FALSE)[[1]], digits = 2,scientific=TRUE));
  
  if (coef(m)[2] >= 0)  {
    cof <- substitute(italic(y) == b + a %.% italic(x),l)
    full <- substitute(italic(y) == b + a %.% italic(x)*","~~italic(r)^2~"="~r2*","~~italic(p)~"="~p2,l)
  } else {
    cof <- substitute(italic(y) == b - a %.% italic(x),l) 
    full <- substitute(italic(y) == b - a %.% italic(x)*","~~italic(r)^2~"="~r2*","~~italic(p)~"="~p2,l)
  }
  
  stat<-substitute(italic(r)^2~"="~r2*","~~italic(p)~"="~p2,l)
  return(list('Coef'=as.character(as.expression(cof)),
              'Stat'=as.character(as.expression(stat)),
              'Full'=as.character(as.expression(full)),
              'Atop'=as.character(as.expression(paste(cof,'\n',stat,sep='') ))))                 
}

lm_rp = function(m) {
  fres<-summary(m)
  l <- list(a = format(abs(coef(m)[2]), digits = 2),
            b = format(coef(m)[1], digits = 2),
            r2 = format(summary(m)$r.squared, digits = 2),
            p2 = format(pf(fres$fstatistic[1], fres$fstatistic[2], fres$fstatistic[3],lower.tail = FALSE)[[1]], digits = 2,scientific=TRUE));
  eq <- substitute(italic(r)^2~"="~r2*","~~italic(p)~"="~p2,l)
  as.character(as.expression(eq));                 
}

ddir<-'Data_final'
odir<-'Figures_final/Keio_full_screen'

theme_set(theme_light())

setwd("~/Projects/2015-Metformin/Worms")

keioinfo<-read.table('../Keio_library/Keio_library_fully_annotated.csv',sep=',',quote = '"',header = TRUE,stringsAsFactors=FALSE)
keioinfo$X<-NULL
#keioinfo<-subset(keioinfo,!Plate %in% c('91','93','95'))
keio<-read.table('Keio_growth.csv',sep=',',quote = '"',header = TRUE,stringsAsFactors=FALSE)

bacmic<-read.table(paste(ddir,'/MICs_and_bacterial_growth-Complete.csv',sep=''),sep=',',quote = '"',header = TRUE,stringsAsFactors=FALSE)


rep1<-read.table('Keio_NGM_full_screen/Rep1/NGM_Keio_Rep1.csv',
           sep=',',quote = '"',header = TRUE,stringsAsFactors=TRUE)

keio$LB_22hr<-as.numeric(keio$LB_22hr)
keio$MOPS_48hr<-as.numeric(keio$MOPS_48hr)
keio$MOPS_24hr<-as.numeric(keio$MOPS_24hr)


fitMOPS<-lm(MOPS_48hr~LB_22hr,keio)
resMOPS<-summary(fitMOPS)




ggplot(keio,aes(x=LB_22hr,y=MOPS_28hr))+
  geom_abline(intercept=resMOPS$coefficients[[1]],slope=resMOPS$coefficients[[2]],color='red')+
  geom_point()+
  scale_x_continuous(breaks=seq(0,1.2,by=0.1),limits=c(0,1.2))+
  scale_y_continuous(breaks=seq(0,1.2,by=0.1),limits=c(0,1.2))+
  xlab('Growth OD in LB 22hr (Baba2006)')+
  ylab('Growth OD in MOPS 48hr (Baba2006)')+
  annotate('text',x = 0.6, y = 1.2, label = lm_eqn(fitMOPS)$Full, parse = TRUE,color ='red',size=5)
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Growth_MOPS_vs_LB_Baba2006.pdf",sep=''),
             width=5,height=4)

fitMOPS24<-lm(MOPS_24hr~LB_22hr,keio)
resMOPS24<-summary(fitMOPS24)

ggplot(keio,aes(x=LB_22hr,y=MOPS_24hr))+
  geom_abline(intercept=resMOPS24$coefficients[[1]],slope=resMOPS24$coefficients[[2]],color='red')+
  geom_point()+
  scale_x_continuous(breaks=seq(0,1.2,by=0.1),limits=c(0,1.2))+
  scale_y_continuous(breaks=seq(0,1.2,by=0.1),limits=c(0,1.2))+
  xlab('Growth OD in LB 22hr (Baba2006)')+
  ylab('Growth OD in MOPS 24hr (Baba2006)')+
  annotate('text',x = 0.6, y = 1.2, label = lm_eqn(fitMOPS24)$Full, parse = TRUE,color ='red',size=5)
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Growth_MOPS24_vs_LB_Baba2006.pdf",sep=''),
             width=5,height=4)




WT_mean<-subset(bacmic,Gene=='WT')$C_OD_Mean
WT_SD<-subset(bacmic,Gene=='WT')$C_OD_SD


growth<-merge(keioinfo,rep1,by=c('Plate','Well'),all.y=TRUE)

blanks<-subset(growth,is.na(Gene))
realblanks<-subset(blanks,OD<0.1)

ggplot(realblanks,aes(x=OD))+geom_histogram()

background<-mean(realblanks$OD)

growth$growth<-growth$OD-background

write.csv(growth,paste(ddir,'/NGM_OD_new_screen_raw.csv',sep=''),row.names = FALSE)


WT<-subset(growth,Plate %in% c('WT1','WT2'))

WT_mean_new<-mean(WT$growth)
WT_SD_new<-sd(WT$growth)

ggplot(WT,aes(x=growth))+geom_histogram()+
  geom_vline(xintercept = WT_mean,color='red')+
  geom_vline(xintercept = WT_mean+WT_SD,color='red',alpha=0.2)+
  geom_vline(xintercept = WT_mean-WT_SD,color='red',alpha=0.2)


table(growth$Gene)

growthc<-subset(growth,!Plate %in% c('WT1','WT2') & !Gene %in% c('WT','present'))





#Merge with Keio library data
ods<-merge(growthc,keio,by.x=c('JW_id','ECK','Gene'),by.y=c('JW.id','ECK.number','Gene'),all.x=TRUE)
ods[is.na(ods$LB_22hr),c('LB_22hr','MOPS_24hr','MOPS_48hr')]<-keio[match(subset(ods,is.na(LB_22hr))$Gene,keio$Gene),c('LB_22hr','MOPS_24hr','MOPS_48hr')]

ods$LB_22hr<-as.numeric(ods$LB_22hr)

fitLB<-lm(growth~LB_22hr,ods)
resLB<-summary(fitLB)
resLB
resLB$r.squared



ggplot(ods,aes(x=LB_22hr,y=growth))+
  geom_abline(intercept=resLB$coefficients[[1]],slope=resLB$coefficients[[2]],color='red')+
  geom_point()+
  scale_x_continuous(breaks=seq(0,1.2,by=0.1),limits=c(0,1.2))+
  scale_y_continuous(breaks=seq(0,1.2,by=0.1),limits=c(0,1.2))+
  xlab('Growth OD in LB')+
  ylab('Growth OD in NGM - new')+
  annotate('text',x = 0.6, y = 1.2, label = lm_eqn(fitLB)$Full, parse = TRUE,color ='red',size=5)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Growth_NGM-new_vs_LB_forcomparison.pdf",sep=''),
             width=5,height=4)
  

all<-merge(ods,subset(bacmic,Gene!='WT'),by=c('Plate','Well'),all.x=TRUE,all.y=TRUE)


selpl<-subset(all,Plate %in% plates)

fitNGM<-lm(growth~C_OD_Mean,selpl)
resNGM<-summary(fitNGM)
resNGM
resNGM$r.squared

plates<-c('15','17','27','29','31',
          '39','41','45','49','51',
          '55','59','71','73','77',
          '83','87','91','93'
)


errcolor<-'grey90'
pntcolor<-'grey20'
ggplot(selpl,aes(x=C_OD_Mean,y=growth))+
  geom_hline(yintercept = WT_mean_new,color='blue',alpha=0.3)+
  geom_vline(xintercept = WT_mean,color='blue',alpha=0.3)+
  annotate('text',x = 0, y = WT_mean_new, label ='WT',color ='blue',size=5)+
  annotate('text',y = 0, x = WT_mean, label ='WT', color ='blue',size=5)+
  geom_abline(intercept=resNGM$coefficients[[1]],slope=resNGM$coefficients[[2]],color='red')+
  #geom_abline(intercept=0,slope=1,color='grey')+
  geom_errorbarh(aes(xmax=C_OD_Mean+C_OD_SD,xmin=C_OD_Mean-C_OD_SD),color=errcolor)+
  geom_point(color=pntcolor)+
  scale_x_continuous(breaks=seq(0,0.6,by=0.1),limits=c(0,0.6))+
  scale_y_continuous(breaks=seq(0,0.6,by=0.1),limits=c(0,0.6))+
  xlab('Growth OD in NGM - old screen')+
  ylab('Growth OD in NGM - new screen')+
  annotate('text',x = 0.3, y = 0.6, label = lm_eqn(fitNGM)$Full, parse = TRUE,color ='red',size=5)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Growth_NGM-new_vs_NGM-old_selected_plates.pdf",sep=''),
             width=5,height=4)



fitMIC<-lm(MIC~growth,all)
resMIC<-summary(fitMIC)
resMIC
resMIC$r.squared



ggplot(all,aes(x=growth,y=MIC))+
  geom_abline(intercept=resMIC$coefficients[[1]],slope=resMIC$coefficients[[2]],color='red')+
  geom_errorbar(aes(ymax=MIC+MIC_SD,ymin=MIC-MIC_SD),color=errcolor)+
  geom_point()+
  scale_x_continuous(breaks=seq(0,0.6,by=0.1),limits=c(0,0.6))+
  scale_y_continuous(breaks=seq(0,120,by=20),limits=c(0,120))+
  ylab('MIC')+
  xlab('Growth OD in NGM - new screen')+
  annotate('text',x = 0.3, y = 120, label = lm_eqn(fitMIC)$Full, parse = TRUE,color ='red',size=5)


# dev.copy2pdf(device=cairo_pdf,
#              file=paste(odir,"/Growth_MIC_vs_NGM-new.pdf",sep=''),
#              width=5,height=4)


selcol<-c('Plate','Well','Gene.y','MIC','MIC_SD','growth','LB_22hr.x','C_OD_Mean','C_OD_SD','T_OD_Mean','T_OD_SD','GT_Interaction')
allc<-all[,selcol]

allc<-rename(allc,c('Gene.y'='Gene','growth'='NGM_OD_new','LB_22hr.x'='LB_22hr'))

write.csv(allc,paste(ddir,'/NGM_OD_new_screen.csv',sep=''),row.names = FALSE)

