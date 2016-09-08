library('ggplot2')
library('gplots')
library('plyr')
library('reshape2')
library(tidyr)
library(quantreg)
library(ellipse)
library(sp)
library(car)

library('gtools')


getinfo<-function(cof) {
  df<-data.frame(cof)
  df$Comparisons<-rownames(df)
  df<-rename(df,c('Std..Error'='SE','Pr...t..'='p.value'))
  df$Stars<-stars.pval(df$p.value)
  dfm<-melt(df,id.vars=c('Comparisons'),variable.name='Stat',value.name = 'Value')
  dfs<-dcast(dfm,.~Comparisons+Stat,value.var = 'Value')
  return(dfs)
}

getstats<-function(xvar,yvar) {
  fitbq<-lm(yvar~xvar)
  resbq<-summary(fitbq)
  coefs<-getinfo(resbq$coefficients)
  return(coefs)
}



GFP<-read.csv('Leo_GFP/GFP_normalised.csv',sep=',',header=TRUE)

GFPm<-melt(GFP,id.vars = c('Replicate','Strain','PL'),variable.name = 'Plasmid',value.name = 'GFPnorm')
GFPm$logGFP<-log(GFPm$GFPnorm,2)
GFPm$PL<-as.factor(GFPm$PL)

ggplot(subset(GFPm,!Plasmid=='pUA66'),aes(x=Plasmid,y=logGFP,fill=PL))+
  geom_boxplot()+
  facet_grid(.~Strain)


GFPStat<-ddply(GFPm,.(Strain,Plasmid,PL),
                  summarise,GFP_mean=mean(logGFP),GFP_SD=sd(logGFP))


GFPmS<-ddply(GFPm,.(Strain,Plasmid),
                 summarise,PL_p.value=getstats(PL,logGFP)[,'xvar1_p.value'],
                 PL_Stars=getstats(PL,logGFP)[,'xvar1_Stars'])

GFPmS$PL<-'1'




GFPComplete<-merge(GFPStat,GFPmS,by=c('Strain','PL','Plasmid'),all.x=TRUE,all.y=TRUE)

GFPC<-subset(GFPComplete,Plasmid!='pUA66')



ggplot(GFPC,aes(y=GFP_mean,x=Plasmid,fill=PL))+
  geom_errorbar(aes(ymin=GFP_mean-GFP_SD,ymax=GFP_mean+GFP_SD),position="dodge", width=.5)+
  geom_bar(stat="identity", position="dodge", width=.5)+
  geom_text(aes(label=PL_Stars),position='dodge',size=5,color='black')+
  labs(fill='Pyradoxal, uM')+
  xlab('Plasmid')+
  ylab('GFP logFC')+
  facet_grid(.~'Strain'+Strain)

dev.copy2pdf(device=cairo_pdf,
             file='Leo_GFP/GFP_Stats.pdf',
             width=9,height=6)

write.csv(GFPC,paste('Leo_GFP/GFP_Stats.csv',sep=''),row.names = FALSE)
