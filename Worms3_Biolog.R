library('ggplot2')
library('gplots')
library('plyr')
library('reshape2')
library('tidyr')
library('gtools')

library(ellipse)
library(sp)
library(xlsx)

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

trim <- function (x) gsub("^\\s+|\\s+$", "", x)

setwd("~/Projects/2015-Metformin/Worms")

ddir<-'Data_final'
odir<-'Figures_final/Biolog/'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")

#Please add descriptors used in Design file to this list!
#Plase add all non-numeric or descriptory columns in Summary file to this list!
design<-c('File','Plate','Strain','Type','Media','Replicate','Well','Index',
          'Data','Name','EcoCycID','Group','Description','Metformin_mM','Sugar_mM','Descriptor')

data<-read.table('./Biolog/Bacteria_data/Bacteria_collected/Summary.csv',sep=',',quote = "\"",header=TRUE)


cols<-colnames(data)
realanot<-intersect(design,cols)

wellspec<-c('Well','Index','Name','EcoCycID','Group','Description')


#File	Plate	Strain	Type	Media	Replicate	Well	Index	Data	Name	EcoCycID	Group	Description

dm<-melt(data,id=realanot,
         variable.name='Descriptor',value.name='Value')

dm$ReplicateB<-paste('B',as.character(dm$Replicate),sep='')

dm$Type<-ifelse(dm$Type=='Control','C','T')

byrepsided<-dcast(dm,Plate+Strain+Well+Index+Name+EcoCycID+Descriptor~Type+ReplicateB,
             fun.aggregate = NULL,value.var = c('Value'),fill = as.numeric(NA))


#Summary
dm_sum<-ddply(dm,.(Plate,Type,Strain,Well,Index,Name,EcoCycID,Descriptor),
              summarise,Mean=mean(Value,na.rm=TRUE),SD=sd(Value,na.rm=TRUE),N=length(Value))

#Give quick summary
ref<-subset(dm_sum,Well=='A1')
ref<-rename(ref, c("Mean"="A1_Mean", "SD"="A1_SD", "N"="A1_N"))
ref<-ref[,!colnames(ref) %in% c('Well','Index','Name','EcoCycID')]

dm_ref<-merge(dm_sum,ref,by=c('Plate','Type','Strain','Descriptor'),all.x=TRUE)



dmrefs<-melt(dm_ref,measure.vars = c('Mean','SD','N','A1_Mean','A1_SD','A1_N'),
         variable.name='Stats',value.name='Value')

sumsided<-dcast(dmrefs,Plate+Strain+Well+Index+Name+EcoCycID+Descriptor~Type+Stats,
                  fun.aggregate = NULL,value.var = c('Value'),fill = as.numeric(NA))

sumsidedw<-sumsided[, -grep("_N", colnames(sumsided))]

write.csv(byrepsided,paste(ddir,'/Biolog_Bacteria_Summary_replicates_by-treatment.csv',sep=''))
write.csv(sumsidedw,paste(ddir,'/Biolog_Bacteria_Summary_averages_by-treatment.csv',sep=''))


#Worm data
wormd<-read.table('./Biolog/Worm_data/Biolog_Worm_linearised.csv',sep=',',quote = "\"",header=TRUE)
bioinfo<-read.table('../Biolog_results/Biolog_metabolites_Ecocyc.csv',sep=',',quote = "\"",header=TRUE)
bioinfo$Metabolite<-trim(bioinfo$Metabolite)

wormdata<-merge(wormd,bioinfo,by=c('Plate','Well'),all.x=TRUE)
wormdata$Name<-NULL
wormdata<-rename(wormdata,c('Metabolite'='Name'))

wormbyrep<-dcast(wormdata,Plate+Well+Index+Name+EcoCycID~Replicate,
                 fun.aggregate = NULL,value.var = c('Value'),fill = as.numeric(NA))


wormsum<-ddply(wormdata,.(Plate,Well),
               summarise,W_Mean=mean(Value,na.rm=TRUE),W_SD=sd(Value,na.rm=TRUE),
               W_Median=median(Value,na.rm=TRUE))

wormall<-merge(wormbyrep,wormsum,by=c('Plate','Well'))
wormall<-rename(wormall,c('1'='W1','2'='W2','3'='W3','4'='W4'))
write.csv(wormall,paste(ddir,'/Biolog_Worms_Summary.csv',sep=''))




#Calculate stats

ints<-subset(dm,Descriptor %in% c('Int_750nm_log'))#,'Int_750nm','Max_750nm','Max_750nm_log','X24h_750nm','X24h_750nm_log'

intrefs<-subset(ints,Well=='A1')
intrefs<-rename(intrefs, c("Value"="A1"))
intrefs<-intrefs[,union(setdiff(realanot,wellspec),c('Descriptor','A1'))]

in_ref<-merge(ints,intrefs,by=c('File','Plate','Type','Data','Replicate','Descriptor'),
              all.x = TRUE,all.y=TRUE)
in_ref$NGMDiff<-in_ref$Value-in_ref$A1
in_ref$UniqueName<-as.factor(paste(as.character(in_ref$Name),'-PM5',sep=''))


in_ref<-rename(in_ref,c('Value'='logODInt'))

inrefm<-melt(in_ref[,!colnames(in_ref) %in% c('A1')],measure.vars = c('logODInt','NGMDiff'),
             variable.name='Measure',value.name='Value')

insided<-dcast(inrefm,Plate+Well+Replicate+ReplicateB+Name+UniqueName+EcoCycID+Group+Description~Measure+Type,
                fun.aggregate = NULL,value.var = c('Value'),fill = as.numeric(NA))



PM5Names<-unique(subset(insided,Plate=='PM5')$Name)
PM12Names<-unique(subset(insided,Plate %in% c('PM1','PM2A'))$Name)
PM1Names<-unique(subset(insided,Plate %in% c('PM1'))$Name)
PM2Names<-unique(subset(insided,Plate %in% c('PM2A'))$Name)

dubl<-setdiff(intersect(PM5Names,PM12Names),c('Negative Control','Positive Control'))


insided$UniqueName<-ifelse(insided$Plate=='PM5',
                           paste(as.character(insided$Name),'-PM5',sep=''),
                           as.character(insided$Name))
insided$PlateGroup<-ifelse(insided$Plate=='PM5','PM5','PM1&PM2A')

insidedA<-insided
insidedA$PlateGroup<-'All'

exclude<-c('2-Hydroxy Benzoic Acid','L-Leucine')#,'L-Leucine'

insidedd<-merge(insided,insidedA,all.x=TRUE,all.y=TRUE)
repfitsjoined<-ddply(subset(insidedd,!UniqueName %in% exclude &
                              !Name %in% c('Positive Control','Negative Control')), .(PlateGroup), summarise,
               logODInt_a=lmsum(lm(`logODInt_T` ~ `logODInt_C`))$a,
               logODInt_b=lmsum(lm(`logODInt_T` ~ `logODInt_C`))$b,
               logODInt_r2=lmsum(lm(`logODInt_T` ~ `logODInt_C`))$r2,
               logODInt_p=lmsum(lm(`logODInt_T` ~ `logODInt_C`))$p2,
               NGMDiff_a=lmsum(lm(`NGMDiff_T` ~ `NGMDiff_C`))$a,
               NGMDiff_b=lmsum(lm(`NGMDiff_T` ~ `NGMDiff_C`))$b,
               NGMDiff_r2=lmsum(lm(`NGMDiff_T` ~ `NGMDiff_C`))$r2,
               NGMDiff_p=lmsum(lm(`NGMDiff_T` ~ `NGMDiff_C`))$p2)
repfitsjoined$Type<-'Replicates'
#Join raw data and fits for normalisation
insided_f<-merge(insided,repfitsjoined[,c('PlateGroup','NGMDiff_a','NGMDiff_b')],by='PlateGroup',all.x=TRUE,all.y=TRUE)
insided_f$CTDiff<-insided_f$NGMDiff_T-insided_f$NGMDiff_C
#insided_f$NGMDiff_T_norm<-(insided_f$NGMDiff_T-insided_f$NGMDiff_b)/insided_f$NGMDiff_a
insided_f$CTDiff_norm<-insided_f$NGMDiff_T-(insided_f$NGMDiff_C*insided_f$NGMDiff_a+insided_f$NGMDiff_b)

insided_f<-subset(insided_f,PlateGroup!='All')
insided_f<-insided_f[,!colnames(insided_f) %in% c('NGMDiff_a','NGMDiff_b')]

insidedm<-melt(insided_f,id.vars = c('PlateGroup','Plate','Well','Replicate','ReplicateB',
                                     'Name','UniqueName','EcoCycID','Group','Description'),
             variable.name='Measure',value.name='Value')


insidstat<-ddply(insidedm,.(PlateGroup,Plate,Well,Name,UniqueName,EcoCycID,Group,Description,Measure),
                summarise,
                Mean=mean(Value,na.rm=TRUE),SD=sd(Value,na.rm=TRUE),
                pval=ifelse(sd(Value,na.rm=TRUE)==0,NA,t.test(Value,mu=0)$p.value))


insidstatm<-melt(insidstat,measure.vars = c('Mean','SD','pval'),
               variable.name='Stat',value.name='Value')



insidsum<-dcast(insidstatm,PlateGroup+Plate+Well+Name+UniqueName+EcoCycID+Group+Description~Measure+Stat,
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
               NGMDiff_p=lmsum(lm(`NGMDiff_T_Mean` ~ `NGMDiff_C_Mean`))$p2)
sumfitsjoined$Type<-'Summary'
allfits<-merge(repfitsjoined,sumfitsjoined,all.x=TRUE,all.y=TRUE)


insidsumc<-insidsum[, -grep("logODInt_._pval", colnames(insidsum))]


allrep<-merge(wormall[,c('Plate','Well','W1','W2','W3','W4','W_Mean','W_SD','W_Median')],
              subset(byrepsided,Descriptor=='Int_750nm_log')[,c('Plate','Well',
                                                                'C_B1','C_B2','C_B3','C_B4',
                                                                'T_B1','T_B2','T_B3','T_B4')],
              by=c('Plate','Well'))

  
allwbr<-merge(allrep,insidsumc,by=c('Plate','Well'))
colnames(allwbr)
allwb<-allwbr[,c(18,1:2,19:23,3:17,24:39)]
colnames(allwb)



bioexpl<-read.table(paste(ddir,'/Biolog_column_explanation.csv',sep=''),
                    sep=',',quote = '"',header = TRUE,stringsAsFactors=FALSE)


#Officia file with explanations
explrd<-subset(bioexpl,Column %in% colnames(allwb))
explrd<-explrd[match(colnames(allwb),explrd$Column),]


write.csv(allwb,paste(ddir,'/Biolog_Combined_Summary_Statistics.csv',sep=''),row.names = FALSE)
write.xlsx2(explrd, file=paste(ddir,'/Biolog_Combined_Summary_Statistics.xlsx',sep=''), sheetName="Readme",row.names = FALSE,showNA=FALSE)
write.xlsx2(allwb, file=paste(ddir,'/Biolog_Combined_Summary_Statistics.xlsx',sep=''), sheetName="Data", append=TRUE,row.names = FALSE)#showNA=FALSE

write.xlsx2(explrd, file='/Users/Povilas/Projects/B-D-H paper/figures and data/figure 5/final files/table S5.xlsx', sheetName="Biolog_Readme", append=TRUE,row.names = FALSE,showNA=FALSE)
write.xlsx2(allwb, file='/Users/Povilas/Projects/B-D-H paper/figures and data/figure 5/final files/table S5.xlsx', sheetName="Biolog_Data", append=TRUE,row.names = FALSE)#showNA=FALSE








#Plotting


clean<-subset(allwb,Well!='A1' & !UniqueName %in% c('Positive Control-PM5','2-Hydroxy Benzoic Acid','L-Leucine') )

pmfit<-subset(allfits,PlateGroup=='All'& Type=='Replicates')


gradcolours<-c('black','yellow','orange','red')


linfit <- list(a = format(pmfit$NGMDiff_a, digits = 2),
          b = format(pmfit$NGMDiff_b, digits = 2),
          r2 =format(pmfit$NGMDiff_r2, digits = 2),
          p2 =format(pmfit$NGMDiff_p, digits = 2,scientific=TRUE));

eqnexpr <- substitute(italic(y) == b + a %.% italic(x)*","~~italic(r)^2~"="~r2*","~~italic(p)~"="~p2,linfit)
eqn<-as.character(as.expression(eqnexpr))

lsugs<-as.character(subset(allwb,Description=='carbohydrate' & (grepl('L-', Name) | grepl('Glucoside', Name)))$Name)



df_NGMT_NGMC<-elipsoid(subset(allwb,!is.na(NGMDiff_C_Mean) & ! is.na(NGMDiff_T_Mean)),
                       'NGMDiff_C_Mean','NGMDiff_T_Mean',scale=0.7)

wmedname<-'C. elegans\nphenotype'
#mrknames<-c('L-Arabinose')
mrknames<-c('Adenosine','Adenosine-PM5','2-Deoxy Adenosine',
            'Thymidine','Inosine',
            '2`-Deoxycytidine-PM5','Cytidine-PM5',
            'Uridine','Uracil-PM5',
            'Uridine-PM5','2`-Deoxyuridine-PM5',
            'L-Arabinose')
txtalpha<-0.7
txtsize<-3
mrksize<-4
mrkcolor<-'red'
mrkalpha<-1
pntalpha<-0.8


erralpha<-1
errcolor<-'grey90'

theme_set(theme_light())
ggplot(clean,aes(x=NGMDiff_C_Mean,y=NGMDiff_T_Mean))+
  geom_abline(data=pmfit,aes(intercept=NGMDiff_b,slope=NGMDiff_a),alpha=0.5,color='red')+
  geom_errorbarh(aes(xmax=NGMDiff_C_Mean+NGMDiff_C_SD,xmin=NGMDiff_C_Mean-NGMDiff_C_SD),height=0,alpha=erralpha,color=errcolor)+
  geom_errorbar(aes(ymax=NGMDiff_T_Mean+NGMDiff_T_SD,ymin=NGMDiff_T_Mean-NGMDiff_T_SD),width=0,alpha=erralpha,color=errcolor)+
  geom_point(aes(size=W_Median,colour=W_Median),alpha=pntalpha)+
  scale_size(range=c(1,5),name=wmedname)+
  scale_colour_gradientn(colours = gradcolours,
                         breaks=c(0,1,2,3,4),limits=c(0,4),guide="legend",name=wmedname)+
  geom_abline(intercept=0,slope=1,alpha=0.2,color='grey',linetype='longdash')+
#   geom_text(aes(label=ifelse(!point.in.polygon(NGMDiff_T_Mean, NGMDiff_C_Mean, df_NGMT_NGMC$NGMDiff_T_Mean, df_NGMT_NGMC$NGMDiff_C_Mean, mode.checked=FALSE) &
#                                !UniqueName %in% mrknames &
#                                (W_Median>3|CTDiff_norm_pval<0.05),as.character(UniqueName),''),colour=W_Median),
#             hjust=-0.1, vjust=-0.75,size=txtsize,alpha=txtalpha)+
#   geom_text(aes(label=ifelse(UniqueName %in% mrknames,as.character(UniqueName),'')),
#             hjust=-0.1, vjust=-0.75,size=mrksize,alpha=mrkalpha,color=mrkcolor)+
  ylim(-2.6,2.1)+xlim(-1.1,1.9)+
  ggtitle('Treatment/Control comparison of bacterial growth logFC\n(NGM+Metabolite/NGM)')+
  xlab('Bacteria growth logFC - Control')+ ylab('Bacteria growth logFC - 5-FU Treatment')+
  #annotate('text',x = -0.25, y =2.1, label = eqn, parse = TRUE,color ='red',size=4)+
  guides(color=guide_legend(), size = guide_legend())
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Bacteria_Scatter_Control|Treatment_NGMsubs_clean.pdf",sep = ''),
             width=7,height=4)
#,colour=W_Median
#CTDiff_norm_pval<0.01 | W_Median>2 || Name %in% c('Uridine','Adenosine')

#+geom_abline(data=allfits,aes(intercept=logODInt_b,slope=logODInt_a),alpha=0.5,color='red')




#Volcano plot Control
volC<-ggplot(clean,aes(x=NGMDiff_C_Mean,y=-log10(NGMDiff_C_pval)))+
  geom_errorbarh(aes(xmax=NGMDiff_C_Mean+NGMDiff_C_SD,xmin=NGMDiff_C_Mean-NGMDiff_C_SD),height=.001,alpha=0.2)+
  geom_point(aes(size=W_Median,colour=W_Median),alpha=pntalpha)+
  scale_size(range=c(1,5),name=wmedname)+
  scale_colour_gradientn(colours = gradcolours,
                         breaks=c(0,1,2,3,4),limits=c(0,4),guide="legend",name=wmedname)+
  guides(color=guide_legend(), size = guide_legend())+
  geom_text(aes(label=ifelse(NGMDiff_C_pval<0.01 | W_Median>2,as.character(UniqueName),''),colour=W_Median),hjust=-0.1, vjust=-0.1,size=3)+
  ggtitle('Control: Metabolite+NGM/NGM bacteria growth logFC')+
  xlab('logFC')+ylab('-log10(p-value)')+
  xlim(-3,3)
volC
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Bacteria_Volcano_Control.pdf",sep = ''),
             width=9,height=6)


#Volcano plot Treatment
volT<-ggplot(clean,aes(x=NGMDiff_T_Mean,y=-log10(NGMDiff_T_pval),color=W_Median))+
  geom_errorbarh(aes(xmax=NGMDiff_T_Mean+NGMDiff_T_SD,xmin=NGMDiff_T_Mean-NGMDiff_T_SD),height=.001,alpha=0.2)+
  geom_point(aes(size=W_Median,colour=W_Median),alpha=pntalpha)+
  scale_size(range=c(1,5),name=wmedname)+
  scale_colour_gradientn(colours = gradcolours,
                         breaks=c(0,1,2,3,4),limits=c(0,4),guide="legend",name=wmedname)+
  guides(color=guide_legend(), size = guide_legend())+
  geom_text(aes(label=ifelse(NGMDiff_T_pval<0.01 | W_Median>2,as.character(UniqueName),''),colour=W_Median),hjust=-0.1, vjust=-0.1,size=3)+
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


#Volcano plot C/T
ggplot(clean,aes(x=CTDiff_Mean,y=-log10(CTDiff_pval)))+
  geom_errorbarh(aes(xmax=CTDiff_Mean+CTDiff_SD,xmin=CTDiff_Mean-CTDiff_SD),height=.001,alpha=0.2)+
  geom_point(aes(size=W_Median,colour=W_Median),alpha=pntalpha)+
  scale_size(range=c(1,5),name=wmedname)+
  scale_colour_gradientn(colours = gradcolours,
                         breaks=c(0,1,2,3,4),limits=c(0,4),guide="legend",name=wmedname)+
  guides(color=guide_legend(), size = guide_legend())+
  geom_text(aes(label=ifelse(CTDiff_pval<0.01| W_Median>2,as.character(UniqueName),''),color=W_Median),hjust=-0.1, vjust=-0.1,size=3)+
  ggtitle('Treatment/Control bacteria growth logFC')+
  xlab('logFC')+ylab('-log10(p-val)')+
  xlim(-3,3)
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Bacteria_Volcano_Control|Treatment_NGMsubs.pdf",sep = ''),
             width=9,height=6)

#Volcano plot C/T normalised
ggplot(clean,aes(x=CTDiff_norm_Mean,y=-log10(CTDiff_norm_pval)))+
  geom_hline(yintercept = -log10(0.05),color='red',linetype='longdash',alpha=0.2)+
  annotate("text", 3,-log10(0.05)+0.05, label = "p=0.05",color='red',size=4,alpha=0.5)+
  geom_errorbarh(aes(xmax=CTDiff_norm_Mean+CTDiff_norm_SD,xmin=CTDiff_norm_Mean-CTDiff_norm_SD),height=.001,alpha=0.2)+
  geom_point(aes(size=W_Median,colour=W_Median),alpha=pntalpha)+
  scale_size(range=c(1,5),name=wmedname)+
  scale_colour_gradientn(colours = gradcolours,
                         breaks=c(0,1,2,3,4),limits=c(0,4),guide="legend",name=wmedname)+
  guides(color=guide_legend(), size = guide_legend())+
  geom_text(aes(label=ifelse((CTDiff_norm_pval<0.01 | W_Median>2) |
                               (CTDiff_norm_pval<0.05 & W_Median>1),
                             as.character(UniqueName),''),color=W_Median),hjust=-0.1, vjust=-0.1,size=3)+
  ggtitle('Treatment/Trend bacteria growth logFC')+
  xlab('logFC')+ylab('-log10(p-value)')+ xlim(-3.5,3.5)
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Bacteria_Volcano_Control|Treatment_NGMsubsTreatNorm.pdf",sep = ''),
             width=9,height=6)






#Distribution of differences
dist<-ggplot(subset(insidsum,Well!='A1' & !Name %in% c('Positive Control','Negative Control')),
       aes(x=CTDiff_Mean))+
  geom_histogram(binwidth = 0.2)+
  ggtitle('Distribution of Treatment over Control logFC distribution\n(NGM bacground substracted)')+
  xlab('Treatment/Control logFC')+
  facet_grid(PlateGroup~.)
dist

distnorm<-ggplot(subset(insidsum,Well!='A1' & !Name %in% c('Positive Control','Negative Control')),
       aes(x=CTDiff_norm_Mean))+
  geom_histogram(binwidth = 0.2)+
  ggtitle('Distribution of Treatment over Control logFC distribution\n(NGM bacground substracted and treatment effect normalised)')+
  xlab('Treatment/Control logFC')+
  facet_grid(PlateGroup~.)

multiplot(dist,distnorm,cols=2)
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Bacteria_Distribution_logFC.pdf",sep = ''),
             width=12,height=6)







#Heatmap bacteria

bgg <- colorRampPalette(c("blue", "gray", "green"))(n = 32)

reorderfun_mean = function(d,w) { reorder(d, w, agglo.FUN = mean) }
reorderfun_median = function(d,w) { reorder(d, w, agglo.FUN = median) }

for (tp in c('CTDiff','CTDiff_norm')) {
  insided_f_rep<-dcast(insided_f,Plate+Well+Name~Replicate,
                       fun.aggregate = NULL,value.var = c(tp),fill = as.numeric(NA))
  ins_f_repc<-subset(insided_f_rep,Well!='A1' & !Name %in% c('Positive Control','Negative Control'))
  heat<-merge(ins_f_repc,allwb[,c('Plate','Well','Name','UniqueName','W1','W2','W3','W4','CTDiff_pval','CTDiff_norm_pval')],by=c('Plate','Well','Name'))
  heat$Namepnorm<-paste(heat$UniqueName,ifelse(stars.pval(heat$CTDiff_norm_pval)!=' ',paste(' (',stars.pval(heat$CTDiff_norm_pval),')',sep=''),''),sep='')
  heat$Namep<-paste(heat$UniqueName,ifelse(stars.pval(heat$CTDiff_pval)!=' ',paste(' (',stars.pval(heat$CTDiff_pval),')',sep=''),''),sep='')
  if (tp=='CTDiff') {
    Names<-heat$Namep
    tpname<-''
  } else {
    Names<-heat$Namepnorm
    tpname<-'-normalised'
  }
  
  comp<-c('1','2','3','4')
  st<-heat[,comp]
  
  for (uk in c(TRUE,FALSE)) {
    if (uk){
      usekey<-'_Key'
      wd<-12
      he<-10
    } else {
      usekey<-''
      wd<-4
      he<-40
    }
    heatmap.2(as.matrix(st),key=uk,Colv=FALSE,trace='none',labRow=Names,
              labCol=comp,col=bgg,
              xlab='Replicate',
              dendrogram="row",scale="none",na.color="grey",
              cexRow=0.4,cexCol=1,
              margin=c(8,16),
              lwid=c(0.2,0.8),
              reorderfun=reorderfun_mean,symkey=TRUE)
    dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Bacteria_Heatmap_logFC',tpname,usekey,'.pdf',sep = ''),
                 width=wd,height=he)
  }
}



#Heatmaps worms
allwbh<-subset(allwb, Well!='A1' &
                 !Name %in% c('Positive Control','Negative Control','2-Hydroxy Benzoic Acid') &
                 (W_Median>1 | CTDiff_norm_pval<0.05))


heat<-allwbh[order(allwbh$W_Mean,decreasing=TRUE),]
heat$Namep<-paste(heat$UniqueName,ifelse(stars.pval(heat$CTDiff_pval)!=' ',paste(' (',stars.pval(heat$CTDiff_pval),')',sep=''),''),sep='')
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

for (uk in c(TRUE,FALSE)) {
  if (uk){
    usekey<-'_Key'
    wd<-12
    he<-10
  } else {
    usekey<-''
    wd<-4
    he<-12
  }
  #hmap<-
  heatmap.2(data,key=uk,Colv=FALSE,Rowv=FALSE,trace='none',labRow=Names,
            labCol=comp,col=bor,
            xlab='Replicate',
            dendrogram="none",scale="none",na.color="grey",
            cexRow=0.6,cexCol=1,margin=c(8,16),
            lwid=c(0.2,0.8),reorderfun=reorderfun_mean,symkey=FALSE)
  #hmap ,reorderfun=reorderfun_mean
  #reorderfun=reorderfun_median
  dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Worms_Heatmap_Score_WMo1_BacSig',usekey,'.pdf',sep = ''),
               width=wd,height=he)
}

hmap<-heatmap.2(data,key=uk,Colv=FALSE,Rowv=FALSE,trace='none',labRow=Names,
                labCol=comp,col=bor,
                xlab='Replicate',
                dendrogram="none",scale="none",na.color="grey",
                cexRow=0.6,cexCol=1,margin=c(8,16),
                lwid=c(0.2,0.8),reorderfun=reorderfun_mean,symkey=FALSE)

#Bacteria heatmap

dat<-heat
dat$Namepnorm<-paste(dat$UniqueName,ifelse(stars.pval(dat$CTDiff_norm_pval)!=' ',paste(' (',stars.pval(dat$CTDiff_norm_pval),')',sep=''),''),sep='')
dat$pnorm<-stars.pval(dat$CTDiff_norm_pval)

#Names<-dat$UniqueName
Names<-dat$pnorm
comp<-c('CTDiff_norm_Mean','CTDiff_norm_Mean')
dats<-dat[,comp]

#,Rowv=FALSE
heatmap.2(as.matrix(dats),key=TRUE,Colv=FALSE,trace='none',labRow=Names,
          labCol=comp,col=bgg,
          xlab='Replicate',Rowv=FALSE,
          dendrogram="none",scale="none",na.color="grey",
          cexRow=0.4,cexCol=1,margin=c(8,16),
          lwid=c(0.2,0.8),reorderfun=reorderfun_mean,symkey=TRUE)

dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Bacteria_Heatmap_logFC-mean_Worm-ordered_WMo1_BacSig_pval.pdf',sep = ''),
             width=4,height=12)

