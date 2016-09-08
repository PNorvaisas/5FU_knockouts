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

library(ellipse)

theme_set(theme_light())

setwd("~/Projects/2015-Metformin/Worms")

#Output folder:
odir<-'Figures_final/Biolog'
ddir<-'Data_final'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")



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

lmsum<-function(m){fres<-summary(m)
l <- list(b = as.double(coef(m)[1], digits = 2),
          a = as.double(coef(m)[2], digits = 2),
          r2 = as.double(summary(m)$r.squared, digits = 3),
          p2 = as.double(pf(fres$fstatistic[1], fres$fstatistic[2], fres$fstatistic[3],lower.tail = FALSE)[[1]], digits = 5))
return(l)
}

#Plotting


allwb<-read.table(paste(ddir,'/Biolog_Combined_Summary_Statistics.csv',sep=''),sep=',',header=TRUE,stringsAsFactors = FALSE)



clean<-subset(allwb,Well!='A1' & !UniqueName %in% c('Positive Control-PM5','2-Hydroxy Benzoic Acid','L-Leucine') )



gradcolours<-c('black','yellow','orange','red')

exclude<-c('2-Hydroxy Benzoic Acid','L-Leucine')

fitNGMDiff<-lm(NGMDiff_T_Mean ~ NGMDiff_C_Mean,data=subset(allwb,!Name %in% exclude))
NGMa<-fitNGMDiff$coefficients[[2]]
NGMb<-fitNGMDiff$coefficients[[1]]

fitODI<-lm(ODInt_T_Mean ~ ODInt_C_Mean,data=allwb)
ODIa<-fitODI$coefficients[[2]]
ODIb<-fitODI$coefficients[[1]]

df_NGMT_NGMC<-elipsoid(subset(allwb,!is.na(NGMDiff_C_Mean) & ! is.na(NGMDiff_T_Mean)),
                       'NGMDiff_C_Mean','NGMDiff_T_Mean',scale=0.7)



df_ODIT_ODIC<-elipsoid(subset(allwb,!is.na(ODInt_C_Mean) & ! is.na(ODInt_T_Mean)),
                       'ODInt_C_Mean','ODInt_T_Mean',scale=0.7)

wmedname<-'C. elegans\nphenotype'
#mrknames<-c('L-Arabinose')
mrknames<-c('Adenosine','Adenosine-PM5','2-Deoxy Adenosine',
            'Thymidine','Inosine',
            '2-Deoxycytidine-PM5','Cytidine-PM5',
            'Uridine','Uracil-PM5',
            'Uridine-PM5','2-Deoxyuridine-PM5',
            'L-Arabinose')
txtalpha<-0.7
txtsize<-3
mrksize<-4
mrkcolor<-'red'
mrkalpha<-1
pntalpha<-0.8


erralpha<-1
errcolor<-'grey90'


BiologNGM<-ggplot(clean,aes(x=NGMDiff_C_Mean,y=NGMDiff_T_Mean))+
  geom_abline(aes(intercept=NGMb,slope=NGMa),alpha=0.5,color='red')+
  geom_errorbarh(aes(xmax=NGMDiff_C_Mean+NGMDiff_C_SD,xmin=NGMDiff_C_Mean-NGMDiff_C_SD),height=0,alpha=erralpha,color=errcolor)+
  geom_errorbar(aes(ymax=NGMDiff_T_Mean+NGMDiff_T_SD,ymin=NGMDiff_T_Mean-NGMDiff_T_SD),width=0,alpha=erralpha,color=errcolor)+
  geom_point(aes(size=W_Median,colour=W_Median),alpha=pntalpha)+
  scale_size(range=c(1,5),name=wmedname)+
  scale_colour_gradientn(colours = gradcolours,
                         breaks=c(0,1,2,3,4),limits=c(0,4),guide="legend",name=wmedname)+
  geom_abline(intercept=0,slope=1,alpha=0.2,color='grey',linetype='longdash')+
  ylim(-2.6,2.1)+xlim(-1.1,1.9)+
  ggtitle('Treatment/Control comparison of bacterial growth logFC\n(NGM+Metabolite/NGM)')+
  xlab('Bacteria growth logFC - Control')+ ylab('Bacteria growth logFC - 5-FU Treatment')+
  #annotate('text',x = -0.25, y =2.1, label = eqn, parse = TRUE,color ='red',size=4)+
  guides(color=guide_legend(), size = guide_legend())
BiologNGM
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Bacteria_Scatter_Control|Treatment_NGMsubs_clean.pdf",sep = ''),
             width=7,height=4)
morenames<-c('L-Lyxose','Orotic acid-PM5','alpha-Methyl-D-Glucoside','L-Sorbose')
BiologNGM+
  geom_text(aes(label=ifelse(UniqueName %in% mrknames|
                             UniqueName %in% morenames,as.character(UniqueName),'')),
            hjust=-0.1, vjust=-0.75,size=mrksize,alpha=mrkalpha,color=mrkcolor)
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Bacteria_Scatter_Control|Treatment_NGMsubs.pdf",sep = ''),
             width=7,height=4)

BiologODI<-ggplot(allwb,aes(x=ODInt_C_Mean,y=ODInt_T_Mean))+
  geom_abline(aes(intercept=ODIb,slope=ODIa),alpha=0.5,color='red')+
  geom_errorbarh(aes(xmax=ODInt_C_Mean+ODInt_C_SD,xmin=ODInt_C_Mean-ODInt_C_SD),height=0,alpha=erralpha,color=errcolor)+
  geom_errorbar(aes(ymax=ODInt_T_Mean+ODInt_T_SD,ymin=ODInt_T_Mean-ODInt_T_SD),width=0,alpha=erralpha,color=errcolor)+
  geom_point(aes(size=W_Median,colour=W_Median),alpha=pntalpha)+
  scale_size(range=c(1,5),name=wmedname)+
  scale_colour_gradientn(colours = gradcolours,
                         breaks=c(0,1,2,3,4),limits=c(0,4),guide="legend",name=wmedname)+
  geom_abline(intercept=0,slope=1,alpha=0.2,color='grey',linetype='longdash')+
  scale_x_continuous(breaks=c(0,4,8,12,16),limits = c(0,16))+
  ylim(0,8)+
  #xlim(0,16)+
  ggtitle('Treatment/Control comparison of bacterial growth OD intergral')+
  xlab('Bacteria growth OD integral - Control')+ ylab('Bacteria growth OD integral - 5-FU Treatment')+
  #annotate('text',x = -0.25, y =2.1, label = eqn, parse = TRUE,color ='red',size=4)+
  guides(color=guide_legend(), size = guide_legend())
BiologODI
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Bacteria_Scatter_ControlvsTreatment_Raw_clean.pdf",sep = ''),
             width=7,height=4)
BiologODI+
#   geom_text(aes(label=ifelse(!point.in.polygon(ODInt_T_Mean, ODInt_C_Mean, df_ODIT_ODIC$ODInt_T_Mean, df_ODIT_ODIC$ODInt_C_Mean, mode.checked=FALSE) &
#                                          !UniqueName %in% mrknames &
#                                          (W_Median>3|MT_p.value<0.05),as.character(UniqueName),''),colour=W_Median),
#                       hjust=-0.1, vjust=-0.75,size=txtsize,alpha=txtalpha)+
  geom_text(aes(label=ifelse(UniqueName %in% c(mrknames,c('Positive Control-PM5',
                                                          'Negative Control',
                                                          'Negative Control-PM5')),as.character(UniqueName),'')),
            hjust=-0.1, vjust=-0.75,size=mrksize,alpha=mrkalpha,color=mrkcolor)
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Bacteria_Scatter_ControlvsTreatment_Raw.pdf",sep = ''),
             width=7,height=4)

#Relationship between W_median and MT_Interaction
ggplot(allwb,aes(x=as.factor(W_Median),y=MT_Interaction))+geom_boxplot()


#Relationship between W_median and MT_Interaction
ggplot(allwb,aes(x=as.factor(W_Median),y=NGMDiff_C_Mean))+geom_boxplot()


ggplot(allwb,aes(y=W_Median,x=NGMDiff_C_Mean))+geom_point()

ggplot(allwb,aes(x=as.factor(W_Median),y=NGMDiff_T_Mean))+geom_boxplot()


#


#Volcano plot Control
volC<-ggplot(clean,aes(x=NGMDiff_C_Mean,y=-log10(NGMDiff_C_p.value)))+
  geom_hline(yintercept = -log(0.05,10),color='red')+
  geom_errorbarh(aes(xmax=NGMDiff_C_Mean+NGMDiff_C_SD,xmin=NGMDiff_C_Mean-NGMDiff_C_SD),height=.001,alpha=0.2)+
  geom_point(aes(size=W_Median,colour=W_Median),alpha=pntalpha)+
  scale_size(range=c(1,5),name=wmedname)+
  scale_colour_gradientn(colours = gradcolours,
                         breaks=c(0,1,2,3,4),limits=c(0,4),guide="legend",name=wmedname)+
  guides(color=guide_legend(), size = guide_legend())+
  geom_text(aes(label=ifelse(NGMDiff_C_p.value<0.01 | W_Median>2,as.character(UniqueName),''),colour=W_Median),hjust=-0.1, vjust=-0.1,size=3)+
  ggtitle('Control: Metabolite+NGM/NGM bacteria growth logFC')+
  xlab('logFC')+ylab('-log10(p-value)')+
  xlim(-3,3)
volC
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Bacteria_Volcano_Control.pdf",sep = ''),
             width=9,height=6)


#Volcano plot Treatment
volT<-ggplot(clean,aes(x=NGMDiff_T_Mean,y=-log10(NGMDiff_T_p.value),color=W_Median))+
  geom_hline(yintercept = -log(0.05,10),color='red')+
  geom_errorbarh(aes(xmax=NGMDiff_T_Mean+NGMDiff_T_SD,xmin=NGMDiff_T_Mean-NGMDiff_T_SD),height=.001,alpha=0.2)+
  geom_point(aes(size=W_Median,colour=W_Median),alpha=pntalpha)+
  scale_size(range=c(1,5),name=wmedname)+
  scale_colour_gradientn(colours = gradcolours,
                         breaks=c(0,1,2,3,4),limits=c(0,4),guide="legend",name=wmedname)+
  guides(color=guide_legend(), size = guide_legend())+
  geom_text(aes(label=ifelse(NGMDiff_T_p.value<0.01 | W_Median>2,as.character(UniqueName),''),colour=W_Median),hjust=-0.1, vjust=-0.1,size=3)+
  ggtitle('Treatment: Metabolite+NGM/NGM bacteria growth logFC')+
  xlab('logFC')+ylab('-log10(p-value)')+
  xlim(-3,3)
volT
#,size=-log10(NGMDiff_T_pval)
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Bacteria_Volcano_Treatment.pdf",sep = ''),
             width=9,height=6)

multiplot(volC,volT,cols=2)
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Bacteria_Volcano_ControlandTreatment.pdf",sep = ''),
             width=12,height=6)


#Volcano plot Metabolite Treatment interaction
ggplot(clean,aes(x=MT_Interaction,y=-log10(MT_p.value)))+
  geom_hline(yintercept = -log(0.05,10),color='red')+
  geom_errorbarh(aes(xmax=MT_Interaction+MT_SD,xmin=MT_Interaction-MT_SD),height=.001,alpha=0.2)+
  geom_point(aes(size=W_Median,colour=W_Median),alpha=pntalpha)+
  scale_size(range=c(1,5),name=wmedname)+
  scale_colour_gradientn(colours = gradcolours,
                         breaks=c(0,1,2,3,4),limits=c(0,4),guide="legend",name=wmedname)+
  guides(color=guide_legend(), size = guide_legend())+
  geom_text(aes(label=ifelse(MT_p.value<0.01| W_Median>2,as.character(UniqueName),''),color=W_Median),hjust=-0.1, vjust=-0.1,size=3)+
  ggtitle('Treatment/Control bacteria growth logFC')+
  xlab('logFC')+ylab('-log10(p-val)')+
  #ylim(0,15)+
  xlim(-3,3)
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Bacteria_Volcano_Control|Treatment_MT_Interaction_p.pdf",sep = ''),
             width=9,height=6)

#Volcano plot Metabolite Treatment interaction FDR
ggplot(clean,aes(x=MT_Interaction,y=-log10(MT_FDR)))+
  geom_hline(yintercept = -log(0.05,10),color='red')+
  geom_errorbarh(aes(xmax=MT_Interaction+MT_SD,xmin=MT_Interaction-MT_SD),height=.001,alpha=0.2)+
  geom_point(aes(size=W_Median,colour=W_Median),alpha=pntalpha)+
  scale_size(range=c(1,5),name=wmedname)+
  scale_colour_gradientn(colours = gradcolours,
                         breaks=c(0,1,2,3,4),limits=c(0,4),guide="legend",name=wmedname)+
  guides(color=guide_legend(), size = guide_legend())+
  geom_text(aes(label=ifelse(MT_FDR<0.01| W_Median>2,as.character(UniqueName),''),color=W_Median),hjust=-0.1, vjust=-0.1,size=3)+
  ggtitle('Treatment/Control bacteria growth logFC')+
  xlab('logFC')+ylab('-log10(p-val)')+
  #ylim(0,15)+
  xlim(-3,3)
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Bacteria_Volcano_Control|Treatment_MT_Interaction_FDR.pdf",sep = ''),
             width=9,height=6)







#Distribution of differences
dist<-ggplot(subset(allwb,Well!='A1' & !Name %in% c('Positive Control','Negative Control')),
             aes(x=MT_Interaction))+
  geom_histogram(binwidth = 0.2)+
  ggtitle('Distribution of Treatment over Control logFC distribution\n(NGM bacground substracted)')+
  xlab('Treatment/Control logFC')+
  facet_grid(PlateGroup~.)
dist

dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Bacteria_Distribution_logFC.pdf",sep = ''),
             width=12,height=6)







#Heatmap bacteria

bgg <- colorRampPalette(c("blue", "gray", "green"))(n = 32)

reorderfun_mean = function(d,w) { reorder(d, w, agglo.FUN = mean) }
reorderfun_median = function(d,w) { reorder(d, w, agglo.FUN = median) }



#Heatmaps worms
allwbh<-subset(allwb, Well!='A1' &
                 !Name %in% c('Positive Control','Negative Control','2-Hydroxy Benzoic Acid') &
                 (W_Median>1 | MT_FDR<0.05))


heat<-allwbh[order(allwbh$W_Mean,decreasing=TRUE),]
#heat$Namep<-paste(heat$UniqueName,ifelse(stars.pval(heat$MT_p.value)!=' ',paste(' (',stars.pval(heat$MT_p.value),')',sep=''),''),sep='')
#heat$NameFDR<-paste(heat$UniqueName,ifelse(stars.pval(heat$MT_FDR)!=' ',paste(' (',stars.pval(heat$MT_FDR),')',sep=''),''),sep='')

Names<-heat$UniqueName
comp<-c('W1','W2','W3','W4')
st<-heat[,comp]
data<-data.matrix(st)

# data        = data.matrix(allwb[,c('W1','W2','W3','W4')])
# distance    = dist(data,method="minkowski")
# cluster     = hclust(distance, method="ward.D")#, method="ward.D"
# dendrogram  = as.dendrogram(cluster)
# Rowv        = rowMeans(data, na.rm = T)
# RowMedians  = apply(data,1,median)
# dendrogram  = reorder(dendrogram, Rowv)
# 

bor <- colorRampPalette(c("black",'yellow', "orange", "red"))(n = 5)

hmap<-heatmap.2(data,key=uk,Colv=FALSE,Rowv=FALSE,trace='none',labRow=Names,
                labCol=comp,col=bor,
                xlab='Replicate',
                dendrogram="none",scale="none",na.color="grey",
                cexRow=0.6,cexCol=1,margin=c(8,16),
                lwid=c(0.2,0.8),reorderfun=reorderfun_mean,symkey=FALSE)

#Bacteria heatmap

dat<-heat
#dat$Namepnorm<-paste(dat$UniqueName,ifelse(stars.pval(dat$MT_p.value)!=' ',paste(' (',stars.pval(dat$MT_p.value),')',sep=''),''),sep='')

#dat$Namepnorm<-paste(dat$UniqueName,ifelse(stars.pval(dat$MT_FDR)!=' ',paste(' (',stars.pval(dat$MT_FDR),')',sep=''),''),sep='')

dat$pnorm<-stars.pval(dat$MT_FDR)

Names<-dat$UniqueName
Names<-dat$pnorm

comp<-c('MT_Interaction','MT_Interaction')
dats<-dat[,comp]

#,Rowv=FALSE
heatmap.2(as.matrix(dats),key=TRUE,Colv=FALSE,trace='none',labRow=Names,
          labCol=comp,col=bgg,
          xlab='Replicate',Rowv=FALSE,
          dendrogram="none",scale="none",na.color="grey",
          cexRow=0.4,cexCol=1,margin=c(8,16),
          lwid=c(0.2,0.8),reorderfun=reorderfun_mean,symkey=TRUE)

dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Bacteria_Heatmap_logFC-mean_Worm-ordered_WMo1_BacSig_names_new.pdf',sep = ''),
             width=4,height=12)
