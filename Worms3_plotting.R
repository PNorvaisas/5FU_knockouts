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




#Output folder:
odir<-'Figures_v3'
ddir<-'Data_v3'

#qbacmicq - no outliers
bacmic<-read.table(paste(ddir,'/MICs_and_bacterial_growth-All.csv',sep=''),sep=',',header=TRUE,stringsAsFactors = FALSE)
qbacmicq<-read.table(paste(ddir,'/MICs_and_bacterial_growth-Clean.csv',sep=''),sep=',',header=TRUE,stringsAsFactors = FALSE)


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

baccor<-ggplot(bacmic,aes(x=NGM_C,y=NGM_D))+
  geom_point(aes(color=ifelse(Gene %in% bac3a$Gene,'3rd','2nd')),size=1)+ylim(0, .25)+
  ylab(expression(paste('Knockout strain growth OD - 100',mu,'M 5FU')))+
  xlab('Knockout strain growth OD - Control')+
  ggtitle(expression(paste('Growth of knockout strains in control and 100',mu,'M 5FU treatment')))+
  stat_smooth(aes(group = 1),method = "lm")+
  geom_abline(intercept=0,slope=1,alpha=0.5,aes(color='grey'),linetype='longdash')+
  geom_text(aes(label=ifelse(NGM_D>NGM_C*bgus+bgui | NGM_D < NGM_C*bgls+bgli | NGM_C<0.03 ,as.character(Gene),'')),
            hjust=-0.1, vjust=-0.1,size=3)+
  geom_errorbarh(aes(xmax=NGM_C+NGM_sd_C,xmin=NGM_C-NGM_sd_C),height=.001,alpha=0.2)+
  geom_errorbar(aes(ymax=NGM_D+NGM_sd_D,ymin=NGM_D-NGM_sd_D),width=0.001,alpha=0.2)+
  geom_abline(intercept=bgli,slope=bgls,alpha=0.5,color='red')+
  geom_abline(intercept=bgui,slope=bgus,alpha=0.5,color='red')+
  annotate("text", 0.25,0.25*bgls+bgli+0.005, label = "5%",color='red')+
  annotate("text", 0.25,0.25*bgus+bgui+0.005, label = "95%",color='red')+
  scale_x_continuous(breaks=seq(0,.3,by=.05))+
  labs(color='Screen')+
  #labs(color=expression(paste('MIC [5FU], ',mu,'M')))+
  annotate('text',x = 0.125, y = 0.25, label = lm_eqn(fitbac), parse = TRUE)
baccor
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Control-Treatment_NGM_growth_Screen-comparison.pdf",sep = ''),width=9,height=9)




bac3cor<-ggplot(subset(bacmic,Gene %in% bac3m$Gene),aes(x=NGM_C,y=NGM_D,color=MIC))+
  geom_point(size=1)+ylim(0, .25)+
  ylab(expression(paste('Knockout strain growth OD - 100',mu,'M 5FU')))+
  xlab('Knockout strain growth OD - Control')+
  ggtitle(expression(paste('Growth of knockout strains in control and 100',mu,'M 5FU treatment')))+
  #stat_smooth(aes(group = 1),method = "lm")+
  geom_abline(intercept=0,slope=1,alpha=0.5,aes(color='grey'),linetype='longdash')+
  #geom_text(aes(label=ifelse(NGM_D>NGM_C*bgus+bgui | NGM_D < NGM_C*bgls+bgli | NGM_C<0.03 ,as.character(Gene),'')),
  #          hjust=-0.1, vjust=-0.1,size=3)+
  geom_text(aes(label=Gene))+
  geom_errorbarh(aes(xmax=NGM_C+NGM_sd_C,xmin=NGM_C-NGM_sd_C),height=.001,alpha=0.2)+
  geom_errorbar(aes(ymax=NGM_D+NGM_sd_D,ymin=NGM_D-NGM_sd_D),width=0.001,alpha=0.2)+
  #geom_abline(intercept=bgli,slope=bgls,alpha=0.5,color='red')+
  #geom_abline(intercept=bgui,slope=bgus,alpha=0.5,color='red')+
  #annotate("text", 0.25,0.25*bgls+bgli+0.005, label = "5%",color='red')+
  #annotate("text", 0.25,0.25*bgus+bgui+0.005, label = "95%",color='red')+
  scale_x_continuous(breaks=seq(0,.3,by=.05))+
  #labs(color='Screen')+
  labs(color=expression(paste('MIC [5FU], ',mu,'M')))+
  annotate('text',x = 0.125, y = 0.25, label = lm_eqn(fitbac), parse = TRUE)
bac3cor
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Control-Treatment_NGM_growth_3rd_Screen.pdf",sep = ''),width=9,height=9)



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

df_MIC_D<-elipsoid(subset(qbacmicq,!is.na(MIC) & ! is.na(NGM_D)),'MIC','NGM_D')

mbcD<-ggplot(qbacmicq,aes(x=MIC,y=NGM_D))+geom_point(size=1)+
  geom_errorbarh(aes(xmax=MIC+MIC_sd,xmin=MIC-MIC_sd),height=.0001,alpha=0.2,color='black')+
  geom_errorbar(aes(ymax=NGM_D+NGM_sd_D,ymin=NGM_D-NGM_sd_D),width=0.0001,alpha=0.2,color='black')+
  geom_text(aes(label=ifelse((NGM_D>MIC*mndus+mndui | NGM_D < MIC*mndls+mndli) | MIC>50 ,Gene,'')),
            hjust=-0.1, vjust=-0.75,size=3)+
  ggtitle(expression(paste('Strain growth in NGM 24hr OD - 100',mu,'M 5FU')))+
  xlab(expression(paste('MIC [5FU], ',mu,'M')))+
  ylab('OD')+ylim(0,0.25)+
  scale_x_continuous(breaks=seq(0,100,by=25))+
  scale_y_continuous(breaks=seq(0,0.25,by=0.05))+
  annotate('text',x = 50, y = 0.225, label = lm_eqn(fitD), parse = TRUE)
  #++xlim(0,115)
  #geom_path(data=df_el, aes(x=MIC, y=NGM_D), size=1, linetype=1,color='grey',alpha=0.5)
  #stat_density2d()
mbcD

mbcD+geom_path(data=df_MIC_D, aes(x=MIC, y=NGM_D), size=1, linetype=1,color='red',alpha=0.2)
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/MIC-NGMTreatment_bac_growth_NoLabels_SD-elipse.pdf",sep=''),width=5,height=5)

mbcD+stat_smooth(aes(group = 1),method = "lm")+
  geom_abline(intercept=mndli,slope=mndls,alpha=0.5,color='red')+
  geom_abline(intercept=mndui,slope=mndus,alpha=0.5,color='red')+
  annotate("text", 100, mndls*100+mndli+0.005, label = "5%",color='red')+
  annotate("text", 100, mndus*100+mndui+0.005, label = "95%",color='red')
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/MIC-NGMTreatment_bac_growth_NoLabels.pdf",sep=''),width=5,height=5)


fitC<-lm(NGM_C ~ MIC,qbacmicq)
fitNCqr<-rq(NGM_C ~ MIC,data=qbacmicq,tau=c(0.05,0.95))
mncli<-coefficients(fitNCqr)[1,][[1]]
mncui<-coefficients(fitNCqr)[1,][[2]]
mncls<-coefficients(fitNCqr)[2,][[1]]
mncus<-coefficients(fitNCqr)[2,][[2]]

df_MIC_C<-elipsoid(subset(qbacmicq,!is.na(MIC) & ! is.na(NGM_C)),'MIC','NGM_C')

mbcC<-ggplot(qbacmicq,aes(x=MIC,y=NGM_C))+geom_point(size=1)+
  ylim(0,0.3)+
  geom_errorbarh(aes(xmax=MIC+MIC_sd,xmin=MIC-MIC_sd),height=.0001,alpha=0.2,color='black')+
  geom_errorbar(aes(ymax=NGM_C+NGM_sd_C,ymin=NGM_C-NGM_sd_C),width=0.001,alpha=0.2,color='black')+
  geom_text(aes(label=ifelse((NGM_C>MIC*mncus+mncui | NGM_C < MIC*mncls+mncli) | MIC>45,Gene,'')),
            hjust=-0.1, vjust=-0.75,size=3)+

  scale_x_continuous(breaks=seq(0,100,by=25))+
  ggtitle('Strain growth in NGM 24hr OD - Control')+xlab(expression(paste('MIC [5FU], ',mu,'M')))+ylab('OD')+
  annotate('text',x = 50, y = 0.29, label = lm_eqn(fitC), parse = TRUE)
mbcC

mbcC+geom_path(data=df_MIC_C, aes(x=MIC, y=NGM_C), size=1, linetype=1,color='red',alpha=0.2)
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/MIC-NGMControl_bac_growth_NoLabels_SD-elipse.pdf",sep=''),width=5,height=5)

mbcC+geom_abline(intercept=mncli,slope=mncls,alpha=0.5,color='red')+
  stat_smooth(aes(group = 1),method = "lm")+
  geom_abline(intercept=mncui,slope=mncus,alpha=0.5,color='red')+
  annotate("text", 100, mncls*100+mncli+0.005, label = "5%",color='red')+
  annotate("text", 100, mncus*100+mncui+0.005, label = "95%",color='red')
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/MIC-NGMControl_bac_growth_NoLabels.pdf",sep=''),width=5,height=5)


fitLB <- lm(LB_22hr ~ MIC, data=qbacmicq)
fitLBqr<-rq(LB_22hr ~ MIC,data=qbacmicq,tau=c(0.05,0.95))
mlbli<-coefficients(fitLBqr)[1,][[1]]
mlbui<-coefficients(fitLBqr)[1,][[2]]
mlbls<-coefficients(fitLBqr)[2,][[1]]
mlbus<-coefficients(fitLBqr)[2,][[2]]

mbcLB<-ggplot(qbacmicq,aes(x=MIC,y=LB_22hr))+geom_point(size=1)+
  stat_smooth(aes(group = 1),method = "lm")+xlim(0,100)+
  geom_errorbarh(aes(xmax=MIC+MIC_sd,xmin=MIC-MIC_sd),height=.0001,alpha=0.2,color='black')+
  geom_text(aes(label=ifelse(((LB_22hr>MIC*mlbus+mlbui | LB_22hr < MIC*mlbls+mlbli) & MIC >1) | MIC>40,Gene,'')),
            hjust=0.05, vjust=-0.75,size=3)+
  geom_abline(intercept=mlbli,slope=mlbls,alpha=0.5,color='red')+
  geom_abline(intercept=mlbui,slope=mlbus,alpha=0.5,color='red')+
  annotate("text", 100, mlbls*100+mlbli+0.02, label = "5%",color='red')+
  annotate("text", 100, mlbus*100+mlbui+0.02, label = "95%",color='red')+
  ggtitle('Strain growth in LB 22hr OD - Control')+
  xlab(expression(paste('MIC [5FU], ',mu,'M')))+ylab('OD')+
  annotate('text',x = 50, y = 1.1, label = lm_eqn(fitLB), parse = TRUE)
mbcLB
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/MIC-LB22hr_bac_growth_NoLabels.pdf",sep=''),width=5,height=5)



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
  geom_point(color='red',size=0.1) +
  #geom_line(color='red',size=0.1)+
  geom_errorbarh(aes(xmax=MIC+MIC_sd,xmin=MIC-MIC_sd),height=.0001,alpha=0.2)+
  scale_x_continuous(breaks=seq(0,100,by=50))+
  ylab('Keio Gene Knockouts (MIC>1)')+
  xlab(expression(paste('MIC [5FU], ',mu,'M')))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
scr2dist+theme(axis.text.y=element_blank(),panel.grid.major=element_blank(),axis.ticks.y=element_blank())
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/MIC_variation_SD_MIC-over-1_crop.pdf",sep=''),width=2,height=6)

scr2dist+theme(axis.text.y = element_text(vjust = 0,size=4))+
  ggtitle('Protective properties of gene knockouts for C. elegans in 5FU exposure')
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/MIC_variation_SD_MIC-over-1.pdf",sep=''),width=8,height=30)



######
q90<-quantile(qbacmicq$MIC,0.9,na.rm=TRUE)[[1]]
q95<-quantile(qbacmicq$MIC,0.95,na.rm=TRUE)[[1]]
q99<-quantile(qbacmicq$MIC,0.99,na.rm=TRUE)[[1]]

dist<-ggplot(qbacmicq,aes(x=MIC))+stat_ecdf()+
  geom_hline(yintercept=0.90,color='green',alpha=0.5,linetype='longdash')+
  geom_vline(xintercept=q90,color='green',alpha=0.5,linetype='longdash')+
  geom_hline(yintercept=0.95,color='blue',alpha=0.5,linetype='longdash')+
  geom_vline(xintercept=q95,color='blue',alpha=0.5,linetype='longdash')+
  geom_hline(yintercept=0.99,color='red',alpha=0.5,linetype='longdash')+
  geom_vline(xintercept=q99,color='red',alpha=0.5,linetype='longdash')+
  annotate("text", 2, 0.92, label = "90%",color='green')+
  annotate("text", 2, 0.97, label = "95%",color='blue')+
  annotate("text", 2, 1, label = "99%",color='red')+
  scale_y_continuous(limits=c(0,1), labels = scales::percent,breaks=seq(0,1,by=0.1))+
  scale_x_log10(breaks=c(0,1,10,100))+
  ylab('')+xlab(expression(paste('MIC [5FU], ',mu,'M')))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dist
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Cumulative_distribution_of_MIC_log10-x-scale_crop.pdf",sep=''),
             width=2,height=6)
dist+ggtitle('Cumulative distribution of MIC values')
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Cumulative_distribution_of_MIC_log10-x-scale.pdf",sep=''),
             width=6,height=6)

