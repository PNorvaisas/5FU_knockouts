library('ggplot2')
library('gplots')
library('plyr')
library('reshape2')
library(tidyr)
library(quantreg)
library(ellipse)
library(sp)
theme_set(theme_light())


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

#cor(df[,c(xvar,yvar)],use='complete.obs')
# lm(df[,c(xvar,yvar)])$coefficients[[2]]

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



multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}




#Output folder:
odir<-'Figures_final'
ddir<-'Data_final'

#bacmic - no outliers
bacmic<-read.table(paste(ddir,'/MICs_and_bacterial_growth-Complete.csv',sep=''),sep=',',header=TRUE,stringsAsFactors = FALSE)


#
bcq05<-quantile(bacmic$OD_C_Mean,0.05,na.rm=TRUE)[[1]]
bcq95<-quantile(bacmic$OD_C_Mean,0.95,na.rm=TRUE)[[1]]
bdq05<-quantile(bacmic$OD_T_Mean,0.05,na.rm=TRUE)[[1]]
bdq95<-quantile(bacmic$OD_T_Mean,0.95,na.rm=TRUE)[[1]]
blq05<-quantile(bacmic$LB_22hr,0.05,na.rm=TRUE)[[1]]
blq95<-quantile(bacmic$LB_22hr,0.95,na.rm=TRUE)[[1]]

fitbac<-lm(OD_T_Mean ~ OD_C_Mean,data=bacmic)
#confint(fitbac,'(Intercept)',level=0.95)[[2]]
#coefficients(fitbac)[[2]]

fitqr<-rq(OD_T_Mean ~ OD_C_Mean,data=bacmic,tau=c(0.05,0.95))
bgli<-coefficients(fitqr)[1,][[1]]
bgui<-coefficients(fitqr)[1,][[2]]
bgls<-coefficients(fitqr)[2,][[1]]
bgus<-coefficients(fitqr)[2,][[2]]


bacres<-subset(bacmic,OD_T_Mean>OD_C_Mean*bgus+bgui)
bacsens<-subset(bacmic,OD_T_Mean<OD_C_Mean*bgls+bgli)



showgenes<-c('upp','dcuC','yjjG')
brks<-c(0,5,10,25,50,100)
micname<-expression(paste('Worm MIC [5FU], ',mu,'M',sep=''))
mrklist<-c('upp','WT','yjjG')
WTdata<-subset(bacmic,Gene=='WT')
#size=0.5
baccor<-ggplot(bacmic,aes(x=OD_C_Mean,y=OD_T_Mean))+
  geom_abline(intercept=0,slope=1,alpha=0.1,aes(color='grey'),linetype='longdash',size=0.5)+
  geom_errorbarh(aes(xmax=OD_C_Mean+OD_C_SD,xmin=OD_C_Mean-OD_C_SD),height=.001,alpha=0.2)+
  geom_errorbar(aes(ymax=OD_T_Mean+OD_T_SD,ymin=OD_T_Mean-OD_T_SD),width=0.001,alpha=0.2)+
  geom_point(aes(size=MIC,colour=MIC),alpha=0.7)+
  geom_abline(intercept=bgli,slope=bgls,alpha=0.5,color='red')+
  geom_abline(intercept=bgui,slope=bgus,alpha=0.5,color='red')+
  annotate("text", 0.35,0.35*bgls+bgli+0.005, label = "5%",color='red',size=5)+
  annotate("text", 0.35,0.35*bgus+bgui+0.005, label = "95%",color='red',size=5)+
  scale_x_continuous(breaks=seq(0,.35,by=.05),limits=c(0,0.35))+
  scale_y_continuous(breaks=seq(0,.25,by=.05),limits=c(0,.25))+
  annotate('text',x = 0.125, y = 0.24, label = lm_eqn(fitbac)$Full, parse = TRUE)+
  stat_smooth(aes(group = 1),method = "lm")+
  geom_text(aes(label=ifelse(Gene %in% mrklist,as.character(Gene),'')),
            hjust=-0.1, vjust=-0.1,size=5,colour = "red")+
  scale_colour_gradientn(colours = c('black','orange','red','red'),
                         breaks=brks,limits=c(5,100),guide='legend',name=micname)+
  scale_size(range=c(1,6),breaks=brks,name=micname)+
  ylab(expression(paste('Bacterial growth OD600nm - 50',mu,'M 5FU')))+
  xlab('Bacterial growth OD600nm - Control')+
  ggtitle(expression(paste('Bacterial growth in NGM 24hr - Control vs 50',mu,'M 5FU treatment')))+
  guides(color=guide_legend(), size = guide_legend())
baccor
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Control-Treatment_NGM_growth.pdf",sep = ''),
             width=9,height=6)

#OD_T_Mean>OD_C_Mean*bgus+bgui+0.02 |
#OD_T_Mean < OD_C_Mean*bgls+bgli-0.04|
  

bacmed<-melt(bacmic[,colnames(bacmic) %in% c('Gene','OD_C_Mean','OD_T_Mean','LB_22hr','MOPS_24hr','MOPS_48hr')],
             id=c('Gene'),variable.name = 'Media',value.name='OD')
bacmed$Media<-factor(bacmed$Media,levels = c('OD_C_Mean','OD_T_Mean','LB_22hr','MOPS_24hr','MOPS_48hr'),
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


#Quantile ranges for control and treatment

# mbcDl<-mbcD+stat_smooth(aes(group = 1),method = "lm")+
#   geom_abline(intercept=mndli,slope=mndls,alpha=0.5,color='red')+
#   geom_abline(intercept=mndui,slope=mndus,alpha=0.5,color='red')+
#   annotate("text", 100, mndls*100+mndli+0.005, label = "5%",color='red')+
#   annotate("text", 100, mndus*100+mndui+0.005, label = "95%",color='red')
# mbcDl
# dev.copy2pdf(device=cairo_pdf,
#              file=paste(odir,"/MIC-NGMTreatment_bac_growth_NoLabels.pdf",sep=''),width=5,height=5)
# 
# 
# mbcDl+geom_path(data=df_MIC_D, aes(x=MIC, y=OD_T_Mean), size=1, linetype=1,color='red',alpha=0.2)
# dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/MIC-NGMTreatment_bac_growth_Linear_SD-elipse.pdf",sep=''),width=5,height=5)
# 


bdqq05<-quantile(bacmic$OD_T_Mean,0.05,na.rm=TRUE)[[1]]
bdqq95<-quantile(bacmic$OD_T_Mean,0.95,na.rm=TRUE)[[1]]
bcqq05<-quantile(bacmic$OD_C_Mean,0.05,na.rm=TRUE)[[1]]
bcqq95<-quantile(bacmic$OD_C_Mean,0.95,na.rm=TRUE)[[1]]
blqq05<-quantile(bacmic$LB_22hr,0.05,na.rm=TRUE)[[1]]
blqq95<-quantile(bacmic$LB_22hr,0.95,na.rm=TRUE)[[1]]


pntsize<-0.5
pntalpha<-0.5
txtsize<-2
txtalpha<-0.5
txtcolor<-'black'
erralpha<-0.2
errcolor<-'black'

mrkgenes<-c('upp','yjjG','WT')
mrkcolor<-'red'
mrkalpha<-1
mrksize<-3

eqsize<-3

#Control
fitC<-lm(OD_C_Mean ~ MIC,bacmic)
fitNCqr<-rq(OD_C_Mean ~ MIC,data=bacmic,tau=c(0.05,0.95))
mncli<-coefficients(fitNCqr)[1,][[1]]
mncui<-coefficients(fitNCqr)[1,][[2]]
mncls<-coefficients(fitNCqr)[2,][[1]]
mncus<-coefficients(fitNCqr)[2,][[2]]

df_MIC_C<-elipsoid(subset(bacmic,!is.na(MIC) & ! is.na(OD_C_Mean)),'MIC','OD_C_Mean')
df_MIC_Cs<-elipsoid(subset(bacmic,!is.na(MIC) & ! is.na(OD_C_Mean)),'MIC','OD_C_Mean',scale=1.2)

#,color=ifelse(Gene %in% mrkgenes,'red',txtcolor)
mbcC<-ggplot(bacmic,aes(x=MIC,y=OD_C_Mean))+
  geom_errorbarh(aes(xmax=MIC+MIC_SD,xmin=MIC-MIC_SD),height=.0001,alpha=erralpha,color=errcolor)+
  geom_errorbar(aes(ymax=OD_C_Mean+OD_C_SD,ymin=OD_C_Mean-OD_C_SD),width=0.001,alpha=erralpha,color=errcolor)+
  geom_path(data=df_MIC_C, aes(x=MIC, y=OD_C_Mean), size=1, linetype=1,color='red',alpha=0.5)+
  geom_abline(intercept=fitC$coefficients[[1]],slope=fitC$coefficients[[2]],alpha=0.5,color='red')+
  geom_point(size=pntsize,alpha=pntalpha)+
  geom_text(aes(label=ifelse(!point.in.polygon(MIC, OD_C_Mean, df_MIC_Cs$MIC, df_MIC_Cs$OD_C_Mean, mode.checked=FALSE) &
                               !Gene %in% mrkgenes,as.character(Gene),'')),
            hjust=-0.1, vjust=-0.75,size=txtsize,alpha=txtalpha,color=txtcolor)+
  geom_text(aes(label=ifelse(Gene %in% mrkgenes,as.character(Gene),'')),
            hjust=-0.1, vjust=-0.75,size=mrksize,alpha=mrkalpha,color=mrkcolor)+
  scale_y_continuous(breaks=seq(0,0.35,by=0.05),limits=c(0,0.35))+
  scale_x_continuous(breaks=seq(0,125,by=25),limits=c(-12.5,112.5))+
  ggtitle('Bacterial growth in NGM 24hr - Control')+
  xlab(expression(paste('Worm MIC [5FU], ',mu,'M')))+
  ylab('Bacterial growth OD600nm')+
  annotate('text',x = 60, y = 0.34, label = lm_eqn(fitC)$Full, parse = TRUE,color ='red',size=eqsize)
mbcC
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/MIC-NGMControl_bac_growth_SD-elipse.pdf",sep=''),
             width=5,height=5)
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/MIC-NGMControl_bac_growth_SD-elipse_large.pdf",sep=''),
             width=9,height=9)


#Treatment
fitD<-lm(OD_T_Mean ~ MIC,bacmic)
fitNDqr<-rq(OD_T_Mean ~ MIC,data=bacmic,tau=c(0.05,0.95))
mndli<-coefficients(fitNDqr)[1,][[1]]
mndui<-coefficients(fitNDqr)[1,][[2]]
mndls<-coefficients(fitNDqr)[2,][[1]]
mndus<-coefficients(fitNDqr)[2,][[2]]


df_MIC_D<-elipsoid(subset(bacmic,!is.na(MIC) & ! is.na(OD_T_Mean)),'MIC','OD_T_Mean')
df_MIC_Ds<-elipsoid(subset(bacmic,!is.na(MIC) & ! is.na(OD_T_Mean)),'MIC','OD_T_Mean',scale=1.2)


mbcD<-ggplot(bacmic,aes(x=MIC,y=OD_T_Mean))+
  geom_errorbarh(aes(xmax=MIC+MIC_SD,xmin=MIC-MIC_SD),height=.0001,alpha=erralpha,color=errcolor)+
  geom_errorbar(aes(ymax=OD_T_Mean+OD_T_SD,ymin=OD_T_Mean-OD_T_SD),width=0.0001,alpha=erralpha,color=errcolor)+
  geom_path(data=df_MIC_D, aes(x=MIC, y=OD_T_Mean), size=1, linetype=1,color='red',alpha=0.5)+
  geom_abline(intercept=fitD$coefficients[[1]],slope=fitD$coefficients[[2]],alpha=0.5,color='red')+
  geom_point(size=pntsize,alpha=pntalpha)+

  geom_text(aes(label=ifelse(!point.in.polygon(MIC, OD_T_Mean, df_MIC_Ds$MIC, df_MIC_Ds$OD_T_Mean, mode.checked=FALSE) &
                               !Gene %in% mrkgenes,as.character(Gene),'')),
            hjust=-0.1, vjust=-0.75,size=txtsize,alpha=txtalpha,color=txtcolor)+
  geom_text(aes(label=ifelse(Gene %in% mrkgenes,as.character(Gene),'')),
            hjust=-0.1, vjust=-0.75,size=mrksize,alpha=mrkalpha,color=mrkcolor)+
  ggtitle(expression(paste('Bacterial growth in NGM 24hr - 50',mu,'M 5FU')))+
  xlab(expression(paste('Worm MIC [5FU], ',mu,'M')))+
  ylab('Bacterial growth OD600nm')+
  scale_x_continuous(breaks=seq(0,100,by=25),limits=c(-12.5,112.5))+
  scale_y_continuous(breaks=seq(0,0.35,by=0.05),limits=c(0,0.35))+
  annotate('text',x = 60, y = 0.34, label = lm_eqn(fitD)$Full, parse = TRUE,color ='red',size=eqsize)
mbcD
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/MIC-NGMTreatment_bac_growth_SD-elipse.pdf",sep=''),
             width=5,height=5)
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/MIC-NGMTreatment_bac_growth_SD-elipse_large.pdf",sep=''),
             width=9,height=9)

multiplot(mbcC, mbcD,  cols=2)
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/MIC-NGM-All_bac_growth_SD-elipse_large.pdf",sep=''),
             width=9,height=6)



#LB
fitLB <- lm(LB_22hr ~ MIC, data=bacmic)
fitLBqr<-rq(LB_22hr ~ MIC,data=bacmic,tau=c(0.05,0.95))
mlbli<-coefficients(fitLBqr)[1,][[1]]
mlbui<-coefficients(fitLBqr)[1,][[2]]
mlbls<-coefficients(fitLBqr)[2,][[1]]
mlbus<-coefficients(fitLBqr)[2,][[2]]

df_MIC_LB<-elipsoid(subset(bacmic,!is.na(MIC) & ! is.na(LB_22hr)),'MIC','LB_22hr')
df_MIC_LBs<-elipsoid(subset(bacmic,!is.na(MIC) & ! is.na(LB_22hr)),'MIC','LB_22hr',scale=1.8)

mbcLB<-ggplot(bacmic,aes(x=MIC,y=LB_22hr))+
  geom_errorbarh(aes(xmax=MIC+MIC_SD,xmin=MIC-MIC_SD),height=.00001,alpha=erralpha,color=errcolor)+
  geom_path(data=df_MIC_LB, aes(x=MIC, y=LB_22hr), size=1, linetype=1,color='red',alpha=0.5)+
  geom_abline(intercept=fitLB$coefficients[[1]],slope=fitLB$coefficients[[2]],alpha=0.5,color='red')+
  geom_point(size=pntsize,alpha=pntalpha)+
  geom_text(aes(label=ifelse(!point.in.polygon(MIC, LB_22hr, df_MIC_LBs$MIC, df_MIC_LBs$LB_22hr, mode.checked=FALSE) &
                               MIC>5 & !Gene %in% mrkgenes,as.character(Gene),'')),
            hjust=-0.1, vjust=-0.75,size=txtsize,alpha=txtalpha,color=txtcolor)+
  geom_text(aes(label=ifelse(Gene %in% mrkgenes,as.character(Gene),'')),
            hjust=-0.1, vjust=-0.75,size=mrksize,alpha=mrkalpha,color=mrkcolor)+
  ggtitle(expression(paste('Bacterial growth in LB 22hr')))+
  xlab(expression(paste('Worm MIC [5FU], ',mu,'M')))+
  ylab('Bacterial growth OD600nm')+
  scale_x_continuous(breaks=seq(0,100,by=25),limits=c(-12.5,112.5))+
  scale_y_continuous(breaks=seq(0,1,by=0.2),limits=c(0,1))+
  annotate('text',x = 60, y = 0.95, label = lm_eqn(fitLB)$Full, parse = TRUE,color ='red',size=eqsize)
mbcLB
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/MIC-LB22hr_bac_growth_SD-elipse.pdf",sep=''),
             width=5,height=5)
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/MIC-LB22hr_bac_growth_SD-elipse_large.pdf",sep=''),
             width=9,height=9)



micname<-expression(paste('Worm MIC [5FU], ',mu,'M',sep=''))

#Volcano plot C/T
ggplot(bacmic,aes(x=CTODDiff_norm_Mean,y=-log10(CTODDiff_norm_pval)))+
  geom_errorbarh(aes(xmax=CTODDiff_norm_Mean+CTODDiff_norm_SD,
                     xmin=CTODDiff_norm_Mean-CTODDiff_norm_SD),height=.001,alpha=0.2)+
  geom_point(aes(size=MIC,colour=MIC),alpha=0.7)+
  scale_size(range=c(1,5),breaks = brks,name=micname)+
  scale_colour_gradientn(colours = c('black','orange','red'),
                         breaks=brks,limits=c(5,100),guide="legend",name=micname)+
  guides(color=guide_legend(), size = guide_legend())+
  geom_text(aes(label=ifelse(CTODDiff_norm_pval<0.01,as.character(Gene),''),color=MIC),
            hjust=-0.1, vjust=-0.1,size=3)+
  ggtitle('Treatment - Control bacteria growth OD difference')+
  xlab('OD difference')+ylab('-log10(p-val)')+
  ylim(0,5)

dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Bacteria_Keio_Volcano_Control|Treatment_ODnormDiff.pdf",sep = ''),
             width=9,height=6)


#Volcano plot C/T WT relative
ggplot(bacmic,aes(x=CTWTDiff_norm_Mean,y=-log10(CTWTDiff_norm_pval)))+
  geom_errorbarh(aes(xmax=CTWTDiff_norm_Mean+CTWTDiff_norm_SD,
                     xmin=CTWTDiff_norm_Mean-CTWTDiff_norm_SD),height=.001,alpha=0.2)+
  geom_point(aes(size=MIC,colour=MIC),alpha=0.7)+
  scale_size(range=c(1,5),breaks = brks,name=micname)+
  scale_colour_gradientn(colours = c('black','orange','red'),
                         breaks=brks,limits=c(5,100),guide="legend",name=micname)+
  guides(color=guide_legend(), size = guide_legend())+
  geom_text(aes(label=ifelse(CTWTDiff_norm_pval<0.01,as.character(Gene),''),color=MIC),
            hjust=-0.1, vjust=-0.1,size=3)+
  ggtitle('Treatment/Control bacteria growth logFC\n(WT relative)')+
  xlab('logFC')+ylab('-log10(p-val)')+
  ylim(0,5)

dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Bacteria_Keio_Volcano_Control|Treatment_WTnormDiff.pdf",sep = ''),
             width=9,height=6)


#WT relative

repfit<-read.table(paste(ddir,'/Bacterial_growth_linear_fits.csv',sep=''),sep=',',header=TRUE,stringsAsFactors = FALSE)
sumfit<-data.frame(WTDiff_a=lmsum(lm(WTDiff_T_Mean ~ WTDiff_C_Mean,data=bacmic))$a,
                    WTDiff_b=lmsum(lm(WTDiff_T_Mean ~ WTDiff_C_Mean,data=bacmic))$b,
                    WTDiff_r2=lmsum(lm(WTDiff_T_Mean ~ WTDiff_C_Mean,data=bacmic))$r2,
                    WTDiff_pval=lmsum(lm(WTDiff_T_Mean ~WTDiff_C_Mean,data=bacmic))$p2)





brks<-c(0,5,10,25,50,100)



ggplot(bacmic,aes(x=WTDiff_C_Mean,y=WTDiff_T_Mean))+
  geom_errorbarh(aes(xmax=WTDiff_C_Mean+WTDiff_C_SD,xmin=WTDiff_C_Mean-WTDiff_C_SD),
                 height=.001,alpha=0.2)+
  geom_errorbar(aes(ymax=WTDiff_T_Mean+WTDiff_T_SD,ymin=WTDiff_T_Mean-WTDiff_T_SD),
                width=0.001,alpha=0.2)+
  geom_abline(data=sumfit,aes(intercept=WTDiff_b,slope=WTDiff_a),alpha=0.5,color='red')+
  geom_point(aes(size=MIC,colour=MIC),alpha=0.6)+
  scale_colour_gradientn(colours = c('black','orange','red','red'),
                         breaks=brks,limits=c(5,100),guide='legend',name=micname)+
  scale_size(range=c(1,6),breaks=brks,name=micname)+
  geom_abline(intercept=0,slope=1,alpha=0.2,aes(color='grey'),linetype='longdash')+
  
  geom_text(aes(label=ifelse(Gene %in% c('upp','yjjG'),
                             as.character(Gene),'')),hjust=-0.1, vjust=-0.1,size=4,color=mrkcolor)+
#  xlim(-2.5,0.5)+ylim(-2.5,0.5)+
  ggtitle('Treatment/Control comparison of bacterial growth logFC\n(Knockout/WT)')+
  xlab('Bacteria growth logFC - Control')+ ylab('Bacteria growth logFC - 5-FU Treatment')+
  guides(color=guide_legend(), size = guide_legend())
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Bacteria_Keio_WT_relative_Scatter.pdf",sep = ''),
             width=9,height=9)



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
  geom_point(color='red',size=0.1) +
  #geom_line(color='red',size=0.1)+
  geom_errorbarh(aes(xmax=MIC+MIC_SD,xmin=MIC-MIC_SD),height=.0001,alpha=0.2)+
  scale_x_continuous(breaks=seq(0,100,by=50))+
  ylab('Keio Gene Knockouts (MIC>5)')+
  xlab(expression(paste('MIC [5FU], ',mu,'M')))
 # theme(axis.text.x = element_text(angle = 90, hjust = 1))
scr2dist+
  theme(axis.text.y=element_blank(),
        panel.grid.major=element_blank(),
        axis.ticks.y=element_blank())
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/MIC_variation_SD_MIC-over-5_crop.pdf",sep=''),width=2,height=6)

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
  ylab('')+xlab(expression(paste('MIC [5FU], ',mu,'M')))
#  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dist
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Cumulative_distribution_of_MIC_log10-x-scale_crop.pdf",sep=''),
             width=2,height=6)
dist+ggtitle('Cumulative distribution of MIC values')+
  scale_x_log10(breaks=c(0,1,5,10,25,50,75,100))+
  annotate("text", q90+1, 0.05, label = sprintf("%1.0f", q90),color='green')+
  annotate("text", q95+2, 0.05, label = sprintf("%1.0f", q95),color='blue')+
  annotate("text", q99+3, 0.05, label = sprintf("%1.0f", q99),color='red')
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Cumulative_distribution_of_MIC_log10-x-scale.pdf",sep=''),
             width=9,height=9)

