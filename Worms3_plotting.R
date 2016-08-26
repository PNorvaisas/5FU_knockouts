library('ggplot2')
library('gplots')
library('plyr')
library('reshape2')
library(tidyr)
library(quantreg)
library(ellipse)
library(sp)

theme_set(theme_light())

theme_update(panel.background = element_rect(colour = "black"),
      axis.text = element_text(colour = "black"))


setwd("~/Projects/2015-Metformin/Worms")

#Vennerable installation: install.packages("Vennerable", repos="http://R-Forge.R-project.org")

#library(quantreg)

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



#Output folder:
odir<-'Figures_final'
ddir<-'Data_final'

#bacmic - no outliers
bacmic<-read.table(paste(ddir,'/MICs_and_bacterial_growth-Complete.csv',sep=''),sep=',',header=TRUE,stringsAsFactors = FALSE)

#Simple calc

length(subset(bacmic,MIC>5 & Gene!='WT')$Gene)/length(subset(bacmic,Gene!='WT')$Gene)

575/3813


#

fitbac<-lm(T_OD_Mean ~ C_OD_Mean,data=bacmic)

erralpha<-1
errcolor<-'grey80'

showgenes<-c('upp','dcuC','yjjG')
brks<-c(5,10,25,50,100)
micname<-expression(paste('C. elegans\nMIC [5FU], ',mu,'M',sep=''))
mrklist<-c('upp','WT','yjjG')
WTdata<-subset(bacmic,Gene=='WT')
gradcols<-c('black','yellow','red')
#size=0.5
baccor<-ggplot(bacmic,aes(x=C_OD_Mean,y=T_OD_Mean))+
  geom_abline(intercept=fitbac$coefficients[[1]],slope=fitbac$coefficients[[2]],alpha=0.5,color='red')+
  geom_abline(intercept=0,slope=1,alpha=0.1,aes(color='grey'),linetype='longdash',size=0.5)+
  geom_errorbarh(aes(xmax=C_OD_Mean+C_OD_SD,xmin=C_OD_Mean-C_OD_SD),height=0,alpha=erralpha,color=errcolor)+
  geom_errorbar(aes(ymax=T_OD_Mean+T_OD_SD,ymin=T_OD_Mean-T_OD_SD),width=0,alpha=erralpha,color=errcolor)+
  geom_point(aes(size=MIC,colour=MIC),alpha=0.9)+#,alpha=1
  scale_x_continuous(breaks=seq(0,.35,by=.05),limits=c(0,0.35))+
  scale_y_continuous(breaks=seq(0,.25,by=.05),limits=c(0,.25))+
  geom_text(aes(label=ifelse(Gene %in% mrklist,as.character(Gene),'')),
            hjust=-0.1, vjust=-0.1,size=5,colour = "red")+
  scale_colour_gradientn(colours = gradcols,trans = "log",
                         breaks=brks,limits=c(5,100),guide='legend',name=micname)+
  scale_size(range=c(1,6),breaks=brks,name=micname)+
  ylab(expression(paste('E. coli growth OD600nm - 50 ',mu,'M 5FU')))+
  xlab('E. coli growth OD600nm - Control')+
  ggtitle(expression(paste('E. coli growth in NGM 24hr - Control vs 50 ',mu,'M 5FU treatment')))+
  guides(color=guide_legend(), size = guide_legend())
baccor
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Control-Treatment_NGM_growth.pdf",sep = ''),
             width=6.5,height=5)

 

bacmed<-melt(bacmic[,colnames(bacmic) %in% c('Gene','C_OD_Mean','T_OD_Mean','LB_22hr','MOPS_24hr','MOPS_48hr')],
             id=c('Gene'),variable.name = 'Media',value.name='OD')
bacmed$Media<-factor(bacmed$Media,levels = c('C_OD_Mean','T_OD_Mean','LB_22hr','MOPS_24hr','MOPS_48hr'),
                     labels=c('NGM - 24h','NGM + 50uM 5FU','LB - 22hr','MOPS - 24hr','MOPS - 48hr'))
bacmed<-subset(bacmed,!Media %in% c('MOPS - 48hr','MOPS - 24hr')) #, 'MOPS - 48hr'

bachist<-ggplot(bacmed,aes(x=OD,fill=Media))+
  geom_histogram(aes(y=0.01*..density..),position='identity',alpha=0.5,binwidth = 0.1)+
  labs(fill='Media')+xlab('OD')+ylab('')+
  scale_y_continuous(limits=c(0,0.075), labels = scales::percent)+
  scale_x_continuous(breaks=seq(0,1.1,by=.1))+
  ggtitle('Distribution of strain growth')
bachist
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Bac_growth_disribution.pdf",sep=''),width=9,height=9)



#Group Significantly sensitive and resistant bacterial hits
bacmic$BacSensitivity<-ifelse(bacmic$GT_Interaction<0 & bacmic$GT_p.value<0.05,'Synergistic',NA)
bacmic$BacSensitivity<-ifelse(bacmic$GT_Interaction>0 & bacmic$GT_p.value<0.05,'Antagonistic',bacmic$BacSensitivity)
bacmic$BacSensitivity<-factor(bacmic$BacSensitivity,levels=c('Synergistic','Antagonistic'),labels=c('Synergistic','Antagonistic'))



txtsize<-4
txtalpha<-1
txtcolor<-'grey40'


pntsize<-1
pntalpha<-1
pntcolor='grey40'

erralpha<-1
errcolor<-'grey90'

mrkgenes<-c('upp','yjjG','WT')
mrkcolor<-'red'
mrkalpha<-1
mrksize<-4

eqsize<-3


xtitle<-expression(paste('Worm MIC [5FU], ',mu,'M'))
ytitle<-'Bacterial growth OD600nm'
#Control
fitC<-lm(C_OD_Mean ~ MIC,bacmic)
fitC2<-lm(MIC ~ C_OD_Mean,bacmic)


df_MIC_C<-elipsoid(subset(bacmic,!is.na(MIC) & ! is.na(C_OD_Mean)),'MIC','C_OD_Mean')
df_MIC_Cs<-elipsoid(subset(bacmic,!is.na(MIC) & ! is.na(C_OD_Mean)),'MIC','C_OD_Mean',scale=1.2)
df_C_MIC<-elipsoid(subset(bacmic,!is.na(MIC) & ! is.na(C_OD_Mean)),'C_OD_Mean','MIC')
df_C_MICs<-elipsoid(subset(bacmic,!is.na(MIC) & ! is.na(C_OD_Mean)),'C_OD_Mean','MIC',scale=1.2)

sensrescol<-c("blue","green")

#,color=ifelse(Gene %in% mrkgenes,'red',txtcolor)
# mbcC<-ggplot(bacmic,aes(x=MIC,y=C_OD_Mean,color=BacSensitivity))+
#   geom_errorbarh(aes(xmax=MIC+MIC_SD,xmin=MIC-MIC_SD),height=0,alpha=erralpha,color=errcolor)+
#   geom_errorbar(aes(ymax=C_OD_Mean+C_OD_SD,ymin=C_OD_Mean-C_OD_SD),width=0,alpha=erralpha,color=errcolor)+
#   geom_path(data=df_MIC_C, aes(x=MIC, y=C_OD_Mean), size=1, linetype=1,color='red',alpha=0.5)+
#   geom_abline(intercept=fitC$coefficients[[1]],slope=fitC$coefficients[[2]],alpha=0.5,color='red')+
#   geom_point(size=pntsize,alpha=pntalpha)+
#   geom_text(aes(label=ifelse(!point.in.polygon(MIC, C_OD_Mean, df_MIC_Cs$MIC, df_MIC_Cs$C_OD_Mean, mode.checked=FALSE) &
#                                !Gene %in% mrkgenes,as.character(Gene),'')),
#             hjust=-0.1, vjust=-0.75,size=txtsize,alpha=txtalpha,color=txtcolor)+
#   geom_text(aes(label=ifelse(Gene %in% mrkgenes,as.character(Gene),'')),
#             hjust=-0.1, vjust=-0.75,size=mrksize,alpha=mrkalpha,color=mrkcolor)+
#   scale_y_continuous(breaks=seq(0,0.35,by=0.05),limits=c(0,0.35))+
#   scale_x_continuous(breaks=seq(0,125,by=25),limits=c(-12.5,112.5))+
#   scale_color_manual(values=sensrescol,na.value=pntcolor)+
#   labs(color='Knockout\ninteraction with 5-FU')+
#   theme(axis.title=element_blank())
# mbcC+guides(colour=FALSE)
# dev.copy2pdf(device=cairo_pdf,
#              file=paste(odir,"/MIC-NGMControl_bac_growth_SD-elipse_clean.pdf",sep=''),
#              width=5,height=4)
# 
# mbcC+ggtitle('Bacterial growth in NGM 24hr - Control')+
#   theme(axis.title=element_text())+
#   xlab(xtitle)+
#   ylab(ytitle)+
#   annotate('text',x = 60, y = 0.34, label = lm_eqn(fitC)$Full, parse = TRUE,color ='red',size=eqsize)
# dev.copy2pdf(device=cairo_pdf,
#              file=paste(odir,"/MIC-NGMControl_bac_growth_SD-elipse.pdf",sep=''),
#              width=5,height=4)
# dev.copy2pdf(device=cairo_pdf,
#              file=paste(odir,"/MIC-NGMControl_bac_growth_SD-elipse_large.pdf",sep=''),
#              width=9,height=6)

mbcC2<-ggplot(bacmic,aes(y=MIC,x=C_OD_Mean,color=BacSensitivity))+
  geom_errorbar(aes(ymax=MIC+MIC_SD,ymin=MIC-MIC_SD),width=0,alpha=erralpha,color=errcolor)+
  geom_errorbarh(aes(xmax=C_OD_Mean+C_OD_SD,xmin=C_OD_Mean-C_OD_SD),height=0,alpha=erralpha,color=errcolor)+
  geom_path(data=df_C_MIC, aes(x=C_OD_Mean, y=MIC), size=1, linetype=1,color='red',alpha=0.5)+
  geom_abline(intercept=fitC2$coefficients[[1]],slope=fitC2$coefficients[[2]],alpha=0.5,color='red')+
  geom_point(size=pntsize,alpha=pntalpha)+
#   geom_text(aes(label=ifelse(!point.in.polygon(MIC, C_OD_Mean, df_MIC_Cs$MIC, df_MIC_Cs$C_OD_Mean, mode.checked=FALSE) &
#                                !Gene %in% mrkgenes,as.character(Gene),'')),
#             hjust=-0.1, vjust=-0.75,size=txtsize,alpha=txtalpha,color=txtcolor)+
  geom_text(aes(label=ifelse(Gene %in% mrkgenes,as.character(Gene),'')),
            hjust=-0.1, vjust=-0.75,size=mrksize,alpha=mrkalpha,color=mrkcolor)+
  scale_x_continuous(breaks=seq(0,0.35,by=0.05),limits=c(0,0.35))+
  scale_y_continuous(breaks=seq(0,125,by=25),limits=c(-12.5,112.5))+
  scale_color_manual(values=sensrescol,na.value=pntcolor)+
  labs(color='Knockout\ninteraction with 5-FU')+
  theme(axis.title=element_blank())
mbcC2+guides(colour=FALSE)
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/NGMControl-MIC_bac_growth_SD-elipse_clean.pdf",sep=''),
             width=5,height=4)

mbcC2+ggtitle('Bacterial growth in NGM 24hr - Control')+
  theme(axis.title=element_text())+
  ylab(xtitle)+
  xlab(ytitle)+
  annotate('text',x = 0.05, y = 100, label = lm_eqn(fitC2)$Full, parse = TRUE,color ='red',size=eqsize)
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/NGMControl-MIC_bac_growth_SD-elipse.pdf",sep=''),
             width=5,height=4)
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/NGMControl-MIC_bac_growth_SD-elipse_large.pdf",sep=''),
             width=9,height=6)



#Treatment
fitD<-lm(T_OD_Mean ~ MIC,bacmic)
fitD2<-lm(MIC ~ T_OD_Mean,bacmic)


df_MIC_D<-elipsoid(subset(bacmic,!is.na(MIC) & ! is.na(T_OD_Mean)),'MIC','T_OD_Mean')
df_MIC_Ds<-elipsoid(subset(bacmic,!is.na(MIC) & ! is.na(T_OD_Mean)),'MIC','T_OD_Mean',scale=1.2)

df_D_MIC<-elipsoid(subset(bacmic,!is.na(MIC) & ! is.na(T_OD_Mean)),'T_OD_Mean','MIC')
df_D_MICs<-elipsoid(subset(bacmic,!is.na(MIC) & ! is.na(T_OD_Mean)),'T_OD_Mean','MIC',scale=1.2)


# mbcD<-ggplot(bacmic,aes(x=MIC,y=T_OD_Mean,color=BacSensitivity))+
#   geom_errorbarh(aes(xmax=MIC+MIC_SD,xmin=MIC-MIC_SD),height=0,alpha=erralpha,color=errcolor)+
#   geom_errorbar(aes(ymax=T_OD_Mean+T_OD_SD,ymin=T_OD_Mean-T_OD_SD),width=0,alpha=erralpha,color=errcolor)+
#   geom_path(data=df_MIC_D, aes(x=MIC, y=T_OD_Mean), size=1, linetype=1,color='red',alpha=0.5)+
#   geom_abline(intercept=fitD$coefficients[[1]],slope=fitD$coefficients[[2]],alpha=0.5,color='red')+
#   geom_point(size=pntsize,alpha=pntalpha)+
# 
#   geom_text(aes(label=ifelse(!point.in.polygon(MIC, T_OD_Mean, df_MIC_Ds$MIC, df_MIC_Ds$T_OD_Mean, mode.checked=FALSE) &
#                                !Gene %in% mrkgenes,as.character(Gene),'')),
#             hjust=-0.1, vjust=-0.75,size=txtsize,alpha=txtalpha,color=txtcolor)+
#   geom_text(aes(label=ifelse(Gene %in% mrkgenes,as.character(Gene),'')),
#             hjust=-0.1, vjust=-0.75,size=mrksize,alpha=mrkalpha,color=mrkcolor)+
#   labs(color='Knockout\ninteraction with 5-FU')+
#   scale_x_continuous(breaks=seq(0,100,by=25),limits=c(-12.5,112.5))+
#   scale_y_continuous(breaks=seq(0,0.35,by=0.05),limits=c(0,0.35))+
#   scale_color_manual(values=sensrescol,na.value=pntcolor)+
#   theme(axis.title=element_blank())
# mbcD+guides(colour=FALSE)
# dev.copy2pdf(device=cairo_pdf,
#              file=paste(odir,"/MIC-NGMTreatment_bac_growth_SD-elipse_clean.pdf",sep=''),
#              width=5,height=4)
# 
# mbcD+ggtitle(expression(paste('Bacterial growth in NGM 24hr - 50',mu,'M 5FU')))+
#   theme(axis.title=element_text())+
#   xlab(xtitle)+
#   ylab(ytitle)+
#   annotate('text',x = 60, y = 0.34, label = lm_eqn(fitD)$Full, parse = TRUE,color ='red',size=eqsize)
# dev.copy2pdf(device=cairo_pdf,
#              file=paste(odir,"/MIC-NGMTreatment_bac_growth_SD-elipse.pdf",sep=''),
#              width=5,height=4)
# dev.copy2pdf(device=cairo_pdf,
#              file=paste(odir,"/MIC-NGMTreatment_bac_growth_SD-elipse_large.pdf",sep=''),
#              width=9,height=6)
# 

#MIC vs growth

mbcD2<-ggplot(bacmic,aes(y=MIC,x=T_OD_Mean,color=BacSensitivity))+
  geom_errorbar(aes(ymax=MIC+MIC_SD,ymin=MIC-MIC_SD),width=0,alpha=erralpha,color=errcolor)+
  geom_errorbarh(aes(xmax=T_OD_Mean+T_OD_SD,xmin=T_OD_Mean-T_OD_SD),height=0,alpha=erralpha,color=errcolor)+
  geom_path(data=df_D_MIC, aes(x=T_OD_Mean,y=MIC), size=1, linetype=1,color='red',alpha=0.5)+
  geom_abline(intercept=fitD2$coefficients[[1]],slope=fitD2$coefficients[[2]],alpha=0.5,color='red')+
  geom_point(size=pntsize,alpha=pntalpha)+
#   geom_text(aes(label=ifelse(!point.in.polygon(MIC, T_OD_Mean, df_MIC_Ds$MIC, df_MIC_Ds$T_OD_Mean, mode.checked=FALSE) &
#                                !Gene %in% mrkgenes,as.character(Gene),'')),
#             hjust=-0.1, vjust=-0.75,size=txtsize,alpha=txtalpha,color=txtcolor)+
  geom_text(aes(label=ifelse(Gene %in% mrkgenes,as.character(Gene),'')),
            hjust=-0.1, vjust=-0.75,size=mrksize,alpha=mrkalpha,color=mrkcolor)+
  labs(color='Knockout\ninteraction with 5-FU')+
  scale_y_continuous(breaks=seq(0,100,by=25),limits=c(-12.5,112.5))+
  scale_x_continuous(breaks=seq(0,0.35,by=0.05),limits=c(0,0.35))+
  scale_color_manual(values=sensrescol,na.value=pntcolor)+
  theme(axis.title=element_blank())
mbcD2+guides(colour=FALSE)
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/NGMTreatment-MIC_bac_growth_SD-elipse_clean.pdf",sep=''),
             width=5,height=4)

mbcD2+ggtitle(expression(paste('Bacterial growth in NGM 24hr - 50',mu,'M 5FU')))+
  theme(axis.title=element_text())+
  xlab(ytitle)+
  ylab(xtitle)+
  annotate('text',x = 0.05, y = 100, label = lm_eqn(fitD2)$Full, parse = TRUE,color ='red',size=eqsize)
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/NGMTreatment-MIC_bac_growth_SD-elipse.pdf",sep=''),
             width=5,height=4)
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/NGMTreatment-MIC_bac_growth_SD-elipse_large.pdf",sep=''),
             width=9,height=6)


#LB


#Try adding Nassos data

nass<-read.table('Nassos_data/5FU-derivatives.txt',sep='\t',header=TRUE,stringsAsFactors = FALSE)
colnames(nass)<-c('Index','5FU_10uM','5FU_37C','Fluorocytosynes','Fluorodeoxyuridines','Fluorouridines','Bromodeoxyuridine')
nass$GeneID<-transform(nass, Index=colsplit(nass$Index, " - ", c('GeneID','Type')) )$Index$GeneID
nass$Type<-transform(nass, Index=colsplit(nass$Index, " - ", c('GeneID','Type')) )$Index$Type
nass<-transform(nass, GeneID=colsplit(nass$GeneID, "-", c('ECK','Gene')))
nass$Index<-NULL
nass$ECK<-nass$GeneID$ECK
#nass$Gene<-nass$GeneID$Gene
nass$GeneID<-NULL

#nass<-transform(nass, ECKall=colsplit(nass$ECK, "/", c('ECK1','ECK2','ECK3','ECK4')))
# nass$ECK1<-transform(nass, ECK=colsplit(nass$ECK, "/", c('ECK1','ECK2','ECK3','ECK4')))$ECK$ECK1
# nass$ECK2<-transform(nass, ECK=colsplit(nass$ECK, "/", c('ECK1','ECK2','ECK3','ECK4')))$ECK$ECK2
# nass$ECK3<-transform(nass, ECK=colsplit(nass$ECK, "/", c('ECK1','ECK2','ECK3','ECK4')))$ECK$ECK3
# nass$ECK4<-transform(nass, ECK=colsplit(nass$ECK, "/", c('ECK1','ECK2','ECK3','ECK4')))$ECK$ECK3
# nass$ECK<-NULL
# nass$ECK<-nass$ECK1
# nass$ECK1<-NULL
# nassm<-melt(nass,measure.vars = c('ECK','ECK2','ECK3','ECK4'),value.name = 'ECK',variable.name = 'ECKs')
# 


nassECK<-nass$ECK
fcECK<-bacmic$ECK
fco5ECK<-subset(bacmic,MIC>5)$ECK

intersect(nass$ECK,bacmic$ECK)
length(intersect(nass$ECK,fco5ECK))


bacmicLBt<-subset(bacmic,MIC>5 | Gene=='WT')

bacmicLB<-merge(bacmicLBt,nass,by='ECK',all.x=TRUE)

length(subset(nass,ECK %in% bacmicLBt$ECK)$Gene)

GT_mean<-mean(bacmicLB$GT_Interaction,na.rm = TRUE)
GT_SD<-sd(bacmicLB$GT_Interaction,na.rm = TRUE)

bacmicLB$GT_Interaction_norm<-scale(bacmicLB$GT_Interaction,center=TRUE,scale=TRUE)
mean(bacmicLB$GT_Interaction_norm,na.rm = TRUE)
sd(bacmicLB$GT_Interaction_norm,na.rm = TRUE)

ggplot(bacmicLB,aes(x=GT_Interaction_norm,y=X5FU_37C))+
  geom_hline(yintercept = 0,color='red',linetype='longdash',alpha=0.5)+
  geom_vline(xintercept = 0,color='red',linetype='longdash',alpha=0.5)+
  geom_hline(yintercept = 2,color='gray80',linetype='longdash')+
  geom_vline(xintercept = 2,color='gray80',linetype='longdash')+
  geom_hline(yintercept = -2,color='gray80',linetype='longdash')+
  geom_vline(xintercept = -2,color='gray80',linetype='longdash')+
  geom_point(size=0.5)+#,alpha=1
  geom_text(aes(label=ifelse(Gene %in% c('yjjG','upp')|
                               abs(GT_Interaction_norm) >2.5|
                               abs(X5FU_37C) >2.5
                             ,as.character(Gene),'')),
            hjust=-0.1, vjust=-0.1,size=3,colour = "red")+
  scale_y_continuous(breaks=seq(-12,4,by=2))+
  scale_x_continuous(breaks=seq(-8,4,by=2),limits=c(-8,3))+  
  ylab('5-FU and KO interaction in LB')+
  xlab('5-FU and KO interaction in NGM')

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Media_LB-NGM_KO-5FU_Interaction.pdf",sep=''),
             width=5,height=4)



# abs(GT_Interaction_norm) >2 |
# abs(X5FU_37C) >2


bacmicLBm<-melt(bacmicLB,measure.vars = c('X5FU_10uM','X5FU_37C','Fluorocytosynes','Fluorodeoxyuridines','Fluorouridines','Bromodeoxyuridine'),
                variable.name = 'Nassos',value.name = 'Interaction')


ggplot(bacmicLBm,aes(x=GT_Interaction_norm,y=Interaction))+
  geom_hline(yintercept = 0,color='gray',linetype='longdash')+
  geom_vline(xintercept = 0,color='gray',linetype='longdash')+
  geom_hline(yintercept = 2,color='gray80',linetype='longdash')+
  geom_vline(xintercept = 2,color='gray80',linetype='longdash')+
  geom_hline(yintercept = -2,color='gray80',linetype='longdash')+
  geom_vline(xintercept = -2,color='gray80',linetype='longdash')+
  geom_point(aes(size=MIC,colour=MIC),alpha=0.9)+#,alpha=1
  geom_text(aes(label=ifelse(Gene %in% mrklist |
                               abs(GT_Interaction_norm) >2 |
                               abs(Interaction) >2
                             ,as.character(Gene),'')),
            hjust=-0.1, vjust=-0.1,size=5,colour = "red")+
  scale_size(range=c(1,6),breaks=brks,name=micname)+
  scale_colour_gradientn(colours = gradcols,trans = "log",
                         breaks=brks,limits=c(5,100),guide='legend',name=micname)+
  scale_y_continuous(breaks=seq(-8,4,by=2),limits=c(-8,4))+
  scale_x_continuous(breaks=seq(-8,4,by=2),limits=c(-8,4))+  
  ylab('KO Drug Interaction in LB')+
  xlab('KO 5-FU Interaction in NGM')+
  labs(color='Knockout\ninteraction with 5-FU')+
  guides(color=guide_legend(), size = guide_legend())+
  geom_smooth(se = TRUE, method = "lm")+
  facet_wrap(~Nassos,nrow=3,ncol=2)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Media_Knockout_Drug_Interaction.pdf",sep=''),
             width=12,height=16)


ggplot(subset(bacmicLBm,Interaction==''),aes(x=GT_Interaction_norm,y=Interaction))+
  geom_hline(yintercept = 0,color='gray',linetype='longdash')+
  geom_vline(xintercept = 0,color='gray',linetype='longdash')+
  geom_hline(yintercept = 2,color='gray80',linetype='longdash')+
  geom_vline(xintercept = 2,color='gray80',linetype='longdash')+
  geom_hline(yintercept = -2,color='gray80',linetype='longdash')+
  geom_vline(xintercept = -2,color='gray80',linetype='longdash')+
  geom_point(aes(size=MIC,colour=MIC),alpha=0.9)+#,alpha=1
  geom_text(aes(label=ifelse(Gene %in% mrklist |
                               abs(GT_Interaction_norm) >2 |
                               abs(Interaction) >2
                             ,as.character(Gene),'')),
            hjust=-0.1, vjust=-0.1,size=5,colour = "red")+
  scale_size(range=c(1,6),breaks=brks,name=micname)+
  scale_colour_gradientn(colours = gradcols,trans = "log",
                         breaks=brks,limits=c(5,100),guide='legend',name=micname)+
  scale_y_continuous(breaks=seq(-8,4,by=2),limits=c(-8,4))+
  scale_x_continuous(breaks=seq(-8,4,by=2),limits=c(-8,4))+  
  ylab('KO Drug Interaction in LB')+
  xlab('KO 5-FU Interaction in NGM')+
  labs(color='Knockout\ninteraction with 5-FU')+
  guides(color=guide_legend(), size = guide_legend())+
  geom_smooth(se = TRUE, method = "lm")+
  facet_wrap(~Nassos,nrow=3,ncol=2)



fitLB <- lm(LB_22hr ~ MIC, data=bacmicLB)

df_MIC_LB<-elipsoid(subset(bacmicLB,!is.na(LB_22hr)),'MIC','LB_22hr') #!is.na(MIC) & 
df_MIC_LBs<-elipsoid(subset(bacmicLB,!is.na(LB_22hr)),'MIC','LB_22hr',scale=1.8) #!is.na(MIC) & 

# mbcLB<-ggplot(bacmicLB,aes(x=MIC,y=LB_22hr,color=BacSensitivity))+
#   geom_errorbarh(aes(xmax=MIC+MIC_SD,xmin=MIC-MIC_SD),height=0,alpha=erralpha,color=errcolor)+
#   geom_path(data=df_MIC_LB, aes(x=MIC, y=LB_22hr), size=1, linetype=1,color='red',alpha=0.5)+
#   geom_abline(intercept=fitLB$coefficients[[1]],slope=fitLB$coefficients[[2]],alpha=0.5,color='red')+
#   geom_point(size=pntsize,alpha=pntalpha)+
#   geom_text(aes(label=ifelse(!point.in.polygon(MIC, LB_22hr, df_MIC_LBs$MIC, df_MIC_LBs$LB_22hr, mode.checked=FALSE) &
#                                MIC>5 & !Gene %in% mrkgenes,as.character(Gene),'')),
#             hjust=-0.1, vjust=-0.75,size=txtsize,alpha=txtalpha,color=txtcolor)+
#   geom_text(aes(label=ifelse(Gene %in% mrkgenes,as.character(Gene),'')),
#             hjust=-0.1, vjust=-0.75,size=mrksize,alpha=mrkalpha,color=mrkcolor)+
#   scale_x_continuous(breaks=seq(0,100,by=25),limits=c(-12.5,112.5))+
#   scale_y_continuous(breaks=seq(0,1,by=0.2),limits=c(0,1))+
#   scale_color_manual(values=sensrescol,na.value=pntcolor)+
#   labs(color='Knockout\ninteraction with 5-FU')+
#   theme(axis.title=element_blank())
# mbcLB+guides(colour=FALSE)
# #,trans='log'
# 
# dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/MIC-LB22hr_bac_growth_SD-elipse_clean.pdf",sep=''),
#              width=5,height=4)
# 
# 
# mbcLB+ggtitle(expression(paste('Bacterial growth in LB 22hr')))+
#   theme(axis.title=element_text())+
#   xlab(xtitle)+
#   ylab(ytitle)+
#   annotate('text',x = 60, y = 0.95, label = lm_eqn(fitLB)$Full, parse = TRUE,color ='red',size=eqsize)
# 
# dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/MIC-LB22hr_bac_growth_SD-elipse.pdf",sep=''),
#              width=5,height=4)
# dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/MIC-LB22hr_bac_growth_SD-elipse_large.pdf",sep=''),
#              width=9,height=6)



#LB vs NGM_C
fitLBC <- lm(LB_22hr ~ C_OD_Mean, data=bacmicLB)
fitLBCqr<-rq(LB_22hr ~ C_OD_Mean,data=bacmicLB,tau=c(0.05,0.95))
mlbcli<-coefficients(fitLBCqr)[1,][[1]]
mlbcui<-coefficients(fitLBCqr)[1,][[2]]
mlbcls<-coefficients(fitLBCqr)[2,][[1]]
mlbcus<-coefficients(fitLBCqr)[2,][[2]]

df_C_LB<-elipsoid(subset(bacmicLB,!is.na(LB_22hr)),'C_OD_Mean','LB_22hr') #!is.na(MIC) & 
df_C_LBs<-elipsoid(subset(bacmicLB,!is.na(LB_22hr)),'C_OD_Mean','LB_22hr',scale=1.8) #!is.na(MIC) & 

mbcLBC<-ggplot(bacmicLB,aes(x=C_OD_Mean,y=LB_22hr))+
  geom_errorbarh(aes(xmax=C_OD_Mean+C_OD_SD,xmin=C_OD_Mean-C_OD_SD),height=0,alpha=erralpha,color=errcolor)+
  geom_abline(intercept=fitLBC$coefficients[[1]],slope=fitLBC$coefficients[[2]],alpha=0.5,color='red')+
  geom_point(aes(size=MIC,colour=MIC),alpha=0.9)+#,alpha=1
  scale_x_continuous(breaks=seq(0,.35,by=.05),limits=c(0,0.35))+
  scale_y_continuous(breaks=seq(0,1,by=.2),limits=c(0,1))+
  geom_text(aes(label=ifelse(Gene %in% mrklist | LB_22hr<0.4,as.character(Gene),'')),
            hjust=-0.1, vjust=-0.1,size=5,colour = "red")+
  scale_size(range=c(1,6),breaks=brks,name=micname)+
  scale_colour_gradientn(colours = gradcols,trans = "log",
                         breaks=brks,limits=c(5,100),guide='legend',name=micname)+

  ylab('E. coli growth LB OD600nm')+
  xlab('E. coli growth NGM OD600nm')+
  labs(color='Knockout\ninteraction with 5-FU')+
  guides(color=guide_legend(), size = guide_legend())
mbcLBC
+guides(colour=FALSE)
#,trans='log'
#  theme(axis.title=element_blank())+
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/NGM_C-LB22hr_MIC_bac_growth.pdf",sep=''),
             width=5,height=4)


mbcLBC+ggtitle(expression(paste('Bacterial growth in LB 22hr')))+
  theme(axis.title=element_text())+
  xlab('Growth in NGM Control')+
  ylab('Growth in LB')+
  annotate('text',x = 0.2, y = 0.95, label = lm_eqn(fitLBC)$Full, parse = TRUE,color ='red',size=eqsize)

dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/NGM_C-LB22hr_MIC_bac_growth.pdf",sep=''),
             width=5,height=4)
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/NGM_C-LB22hr_MIC_bac_growth_large.pdf",sep=''),
             width=9,height=6)



micname<-expression(paste('C. elegans\nMIC [5FU], ',mu,'M',sep=''))

#Volcano plot C/T
ggplot(bacmic,aes(x=GT_Interaction,y=-log10(GT_p.value)))+
  geom_errorbarh(aes(xmax=GT_Interaction+GT_SD,
                     xmin=GT_Interaction-GT_SD),height=.001,alpha=0.2)+
  geom_point(aes(size=MIC,colour=MIC),alpha=0.7)+
  scale_size(range=c(1,5),breaks = brks,name=micname)+
  scale_colour_gradientn(colours = gradcols,trans = "log",
                         breaks=brks,limits=c(5,100),guide="legend",name=micname)+
  guides(color=guide_legend(), size = guide_legend())+
  geom_text(aes(label=ifelse(GT_p.value<0.01,as.character(Gene),''),color=MIC),
            hjust=-0.1, vjust=-0.1,size=3)+
  ggtitle('Knockout ')+
  xlab('KO 5-FU Interaction')+ylab('-log10(p-val)')+
  ylim(0,5)

dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Bacteria_Keio_Volcano_Control|Treatment_WTnormDiff.pdf",sep = ''),
             width=9,height=6)




#WT relative

wtfit<-lm(T_WTDiff ~ C_WTDiff,data=bacmic)

#brks<-c(0,5,25,50,100)
showKO<-c('upp','yjjG')#,'gcvA','aceE','gltA'

logODcor<-ggplot(bacmic,aes(x=C_WTDiff,y=T_WTDiff))+
  geom_errorbarh(aes(xmax=C_WTDiff+C_WTDiff_SD,xmin=C_WTDiff-C_WTDiff_SD),
                 height=0,alpha=erralpha,color=errcolor)+
  geom_errorbar(aes(ymax=T_WTDiff+T_WTDiff_SD,ymin=T_WTDiff-T_WTDiff_SD),
                width=0,alpha=erralpha,color=errcolor)+
  geom_abline(data=sumfit,aes(intercept=WTDiff_b,slope=WTDiff_a),alpha=0.5,color='red')+
  geom_abline(intercept=0,slope=1,alpha=0.2,aes(color='grey'),linetype='longdash')+
  geom_point(aes(size=MIC,colour=MIC),alpha=0.99)+
  scale_colour_gradientn(colours = gradcols,trans = "log",
                         breaks=brks,limits=c(5,100),guide='legend',name=micname)+
  scale_size(range=c(1,6),breaks=brks,name=micname)+
  ylim(-2.5,0.75)+
  annotate('text',x = -1.5, y = 0.25, label = lm_eqn(wtfit)$Full, parse = TRUE,color='red')+
  geom_text(aes(label=ifelse(Gene %in% showKO,
                             as.character(Gene),'')),
            hjust=-0.1, vjust=-0.1,size=5,colour = "red")+#
  ggtitle('Treatment/Control comparison of bacterial growth logFC\n(Knockout/WT)')+
  xlab('E. coli growth logFC (KO/WT) - Control')+ ylab('E. coli growth logFC (KO/WT) - 5-FU Treatment')+
  guides(color=guide_legend(), size = guide_legend())
logODcor
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Bacteria_Keio_WT_relative_Scatter.pdf",sep = ''),
             width=6,height=5)
# |
#   C_WTDiff< -0.9 |
#   T_WTDiff< -0.9 |
#   C_WTDiff> 0.1 |
#   T_WTDiff> 0.2 |
#   MIC >90


#CTWTDiff_norm_pval<0.05 & abs(CTWTDiff_norm_Mean)>0.75 & !Gene=='WT'




#All MICs
alldist<-ggplot(bacmic,aes(x=MIC,y=reorder(Gene,MIC,max)))+
  geom_point(color='red',size=1) + geom_errorbarh(aes(xmax=MIC+MIC_SD,xmin=MIC-MIC_SD))+
  theme(axis.text.y = element_text(vjust = 0,size=4))+
  scale_x_continuous(breaks=seq(0,100,by=10))+
  ylab('Gene knockout')+
  xlab(expression(paste('MIC [5FU], ',mu,'M')))+
  ggtitle('Protective properties of gene knockouts for C. elegans in 5FU exposure')
alldist
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/MIC_variation_SD_all.pdf",sep=''),width=8,height=150)

#MICs over 5
scr2dist<-ggplot(subset(bacmic,MIC>5),aes(x=MIC,y=reorder(Gene,MIC,max)))+
  geom_errorbarh(aes(xmax=MIC+MIC_SD,xmin=MIC-MIC_SD),height=0,alpha=erralpha,color=errcolor)+
  geom_point(color='red',size=0.1) +
  #scale_x_continuous(breaks=seq(0,100,by=50))+
  scale_x_log10(breaks=c(0,1,10,100))+
  #ylab('Keio Gene Knockouts (MIC>5)')+
  ylab('')+xlab(xtitle)
scr2dist+
  theme(axis.text.y=element_blank(),
        panel.grid.major=element_blank(),
        axis.ticks.y=element_blank())+
  annotation_logticks(base = 10,color='grey50',sides='b')
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/MIC_variation_SD_MIC-over-5_crop.pdf",sep=''),
             width=2,height=6)

#scale_x_log10(breaks=c(0,1,10,100))+

scr2dist+theme(axis.text.y = element_text(vjust = 0,size=4))+
  ggtitle('Protective properties of gene knockouts for C. elegans in 5FU exposure')+
  scale_x_continuous(breaks=seq(0,100,by=25))
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/MIC_variation_SD_MIC-over-5.pdf",sep=''),width=8,height=30)



######
q90<-quantile(bacmic$MIC,0.9,na.rm=TRUE)[[1]]
q95<-quantile(bacmic$MIC,0.95,na.rm=TRUE)[[1]]
q99<-quantile(bacmic$MIC,0.99,na.rm=TRUE)[[1]]

dist<-ggplot(bacmic,aes(x=MIC))+stat_ecdf()+
  geom_hline(yintercept=0.90,color='green',alpha=0.5,linetype='longdash')+
  #geom_vline(xintercept=q90,color='green',alpha=0.5,linetype='longdash')+
  geom_hline(yintercept=0.95,color='blue',alpha=0.5,linetype='longdash')+
  #geom_vline(xintercept=q95,color='blue',alpha=0.5,linetype='longdash')+
  geom_hline(yintercept=0.99,color='red',alpha=0.5,linetype='longdash')+
  #geom_vline(xintercept=q99,color='red',alpha=0.5,linetype='longdash')+
  annotate("text", 2, 0.92, label = "90%",color='green')+
  annotate("text", 2, 0.97, label = "95%",color='blue')+
  annotate("text", 2, 1, label = "99%",color='red')+
  scale_y_continuous(limits=c(0,1), labels = scales::percent,breaks=seq(0,1,by=0.1))+
  scale_x_log10(breaks=c(0,1,10,100))+
  ylab('')+xlab(xtitle)+
  annotation_logticks(base = 10,color='grey50',sides='b')
#  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dist+
  annotate("text", q90+1, 0.9-0.02, label = sprintf("%1.0f", q90),color='green')+
  annotate("text", q95+2, 0.95-0.02, label = sprintf("%1.0f", q95),color='blue')+
  annotate("text", q99+3, 0.99-0.02, label = sprintf("%1.0f", q99),color='red')
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Cumulative_distribution_of_MIC_log10-x-scale_crop.pdf",sep=''),
             width=2,height=6)
dist+ggtitle('Cumulative distribution of MIC values')+
  #scale_x_log10(breaks=c(0,1,5,10,25,50,75,100))+
  annotate("text", q90+1, 0.01, label = sprintf("%1.0f", q90),color='green')+
  annotate("text", q95+2, 0.01, label = sprintf("%1.0f", q95),color='blue')+
  annotate("text", q99+3, 0.01, label = sprintf("%1.0f", q99),color='red')
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Cumulative_distribution_of_MIC_log10-x-scale.pdf",sep=''),
             width=9,height=9)

