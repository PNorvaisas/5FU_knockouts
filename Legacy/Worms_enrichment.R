library('ggplot2')
library('gplots')
library('plyr')
library('reshape2')
library(tidyr)
library('Vennerable')
library('GOplot')
library(quantreg)

library('KEGG.db')
library('KEGGREST')
library('biomaRt')



# source("https://bioconductor.org/biocLite.R")
# biocLite("biomaRt")
# listDatasets(useMart("ENSEMBL_MART_ENSEMBL", host = "www.ensembl.org"))
# mart=useMart("bacterial_mart_5")
# listDatasets(mart)



keggInfo('pathway')
keggLink("pathway", "eco")

keggLink("pathway", "eco")
keggList("pathway")
keggLink("pathway", "eco")

keggeco<-keggList("eco")

keggmec<-keggLink("pathway", "eco")

url <- color.pathway.by.objects("path:eco00260",
                                c("eco:b0002", "eco:c00263"),
                                c("#ff0000", "#00ff00"), c("#ffff00", "yellow"))
if(interactive()){
  browseURL(url)
}

keggGet(c('path:eco00260'),option=c('kgml'))
keggGet(c('path:eco00260'),option=c('kgml'))




GOsplit<-function(x,rtrn="GO"){
  if (grepl('~',x)) {
    sep='~'
  } else if (grepl(':',x)) {
    sep=':'
  } else {
    sep=' '
  }
  splt<-strsplit(x,split=sep)[[1]]
  res<-list()
  if (length(splt)==2){
    res$GO<-as.character(splt[1])
    res$Term<-as.character(splt[2])
  } else {
    res$GO<-NA
    res$Term<-as.character(splt[1])
  }
  if (rtrn=="GO"){
    return(res$GO)
  } else if (rtrn=="Term"){
    return(res$Term)
  }
  
}

read.GO<-function(flnm,filter=TRUE,pval='Benjamini'){
  GO<-read.table(flnm,sep='\t',quote = '"',header = TRUE,stringsAsFactors=FALSE)
  as.data.frame(GO)
  #rownames(GO)<-str(rownames(GO))
  GO$ID<-as.character(lapply(GO$Term,GOsplit,rtrn="GO"))
  GO$Term<-as.character(lapply(GO$Term,GOsplit,rtrn="Term"))
  if (pval=='Benjamini'){
    GO$adj_pval<-GO$Benjamini
  } else if (pval=='Bonferoni'){
    GO$adj_pval<-GO$Bonferoni
  } else if (pval=='FDR'){
    GO$adj_pval<-GO$FDR
  } else if (pval=='PValue'){
    GO$adj_pval<-GO$PValue
  } else {
    GO$adj_pval<-GO$PValue
  }
  
  GO[GO$Category=="KEGG_PATHWAY",]$Category<-"KP"
  GO[GO$Category=="GOTERM_BP_FAT",]$Category<-"BP"
  GO[GO$Category=="GOTERM_CC_FAT",]$Category<-"CC"
  GO[GO$Category=="GOTERM_MF_FAT",]$Category<-"MF"
  if (filter){
    GO<-subset(GO,Category %in% c("BP","MF","CC"))
  }
  
  return(GO)
}

qbacmiq<-read.table('Data/MICs_and_bacterial_growth-Unique_clean.csv',sep=',',header=TRUE)
unicom<-read.table('Data/MICs_Unknown_function.csv',sep=',',header=TRUE,stringsAsFactors = FALSE)






PLP<-read.table('Data/Genes_using_PLP.csv',sep=',',header=FALSE)
all<-subset(qbacmic,!is.na(Gene))$Gene
PLPg<-PLP$V1
MIC1<-subset(qbacmic,!is.na(Gene) & MIC>1)$Gene
MIC25<-subset(qbacmic,!is.na(Gene) & MIC>2.5)$Gene


PLPl<-list('PLP using'=as.character(PLPg),'Whole Keio library'=all,'MIC>1'=MIC1,'MIC>2.5'=MIC25)

plot(Venn(PLPl),type='ellipses', doWeights = FALSE)#,type='ellipses',type='ellipses'
dev.copy2pdf(device=cairo_pdf,file="Figures/PLP_using_enzymes.pdf",width=9,height=9)


PLPc<-PLPl[c('PLP using','MIC>1','MIC>2.5')]
PLP_df<-plyr::ldply(PLPc, cbind)
PLP_df<-rename(PLP_df,c('.id'='List','1'='Gene'))


plps<-subset(qbacmic, Gene %in% unique(PLP_df$Gene))
plps$'PLP using'<-ifelse(plps$Gene %in% PLPc$'PLP using',TRUE,FALSE)
plps$'MIC>1'<-ifelse(plps$Gene %in% PLPc$'MIC>1',TRUE,FALSE)
plps$'MIC>2.5'<-ifelse(plps$Gene %in% PLPc$'MIC>2.5',TRUE,FALSE)
plps<-plps[,c(1,19:21,2:7,9:18,8)]
write.csv(plps,'Data/PLP_use_MIC_overlap.csv')



unw<-list('Whole Keio library'=all,'Unknown function'=subset(unicom,Unknown.function=='TRUE')$Gene,'MIC>2.5'=MIC25)


plot(Venn(unw),doEuler=TRUE, doWeights = FALSE)#,type='ellipses',type='ellipses' ,type='ellipses'  
dev.copy2pdf(device=cairo_pdf,file="Figures/Unknown_enzymes_3sets.pdf",width=9,height=9)





fitqr<-rq(NGM_D ~ NGM_C,data=qbacmicq,tau=c(0.05,0.95))
bgli<-coefficients(fitqr)[1,][[1]]
bgui<-coefficients(fitqr)[1,][[2]]
bgls<-coefficients(fitqr)[2,][[1]]
bgus<-coefficients(fitqr)[2,][[2]]

fitqr2<-rq(NGM_D ~ NGM_C,data=qbacmicq,tau=c(0.10,0.90))
bgli2<-coefficients(fitqr2)[1,][[1]]
bgui2<-coefficients(fitqr2)[1,][[2]]
bgls2<-coefficients(fitqr2)[2,][[1]]
bgus2<-coefficients(fitqr2)[2,][[2]]


q05<-quantile(qbacmic$MIC,0.05,na.rm=TRUE)[[1]]
q10<-quantile(qbacmic$MIC,0.1,na.rm=TRUE)[[1]]
q90<-quantile(qbacmic$MIC,0.9,na.rm=TRUE)[[1]]
q95<-quantile(qbacmic$MIC,0.95,na.rm=TRUE)[[1]]
q99<-quantile(qbacmic$MIC,0.99,na.rm=TRUE)[[1]]



bacres95<-subset(qbacmicq,NGM_D>NGM_C*bgus+bgui)$Gene
bacsens05<-subset(qbacmicq,NGM_D<NGM_C*bgls+bgli)$Gene

bacres90<-subset(qbacmicq,NGM_D>NGM_C*bgus2+bgui2)$Gene
bacsens10<-subset(qbacmicq,NGM_D<NGM_C*bgls2+bgli2)$Gene

wormres95<-subset(qbacmicq,MIC>q95)$Gene
wormsens05<-subset(qbacmicq,MIC<q05)$Gene

wormres90<-subset(qbacmicq,MIC>q90)$Gene
wormres25<-subset(qbacmicq,MIC>2.5)$Gene
wormresall<-subset(qbacmicq,MIC>1)$Gene
wormresall5<-subset(qbacmicq,MIC>5)$Gene
wormsens10<-subset(qbacmicq,MIC<q10)$Gene

resmic5=list('Bacteria sensitive bottom 5%'=bacsens05,
            'Bacteria resistant top 5%'=bacres95,
            'Worms resistant top 5%'=wormres95,
            'Worms sensitive bottom 5%'=wormsens05)

resmic10=list('Bacteria sensitive bottom 10%'=bacsens10,
              'Bacteria resistant top 10%'=bacres90,
             'Worms resistant top 10%'=wormres90) # ,'Worms sensitive bottom 10%'=wormsens10 ,'Worms resistant all'=wormresall5

resmic_df<-plyr::ldply(resmic10, cbind)
resmic_df<-rename(resmic_df,c('.id'='List','1'='Gene'))


rmc<-subset(qbacmic, Gene %in% unique(resmic_df$Gene))
rmc$'Bacteria resistant top 5%'<-ifelse(rmc$Gene %in% subset(resmic_df,List=="Bacteria resistant top 5%")$Gene,TRUE,FALSE)
rmc$'Bacteria sensitive bottom 5%'<-ifelse(rmc$Gene %in% subset(resmic_df,List=="Bacteria sensitive bottom 5%")$Gene,TRUE,FALSE)
rmc$'Worms resistant top 5%'<-ifelse(rmc$Gene %in% subset(resmic_df,List=="Worms resistant top 5%")$Gene,TRUE,FALSE)
rmc$'Worms sensitive bottom 5%'<-ifelse(rmc$Gene %in% subset(resmic_df,List=="Worms sensitive bottome 5%%")$Gene,TRUE,FALSE)
rmc<-rmc[,c(1,19:21,2:7,9:18,8)]
write.csv(rmc,'Data/Venn_Worm-Bacteria_resistance_overlap.csv')

plot(Venn(resmic10), doWeights = FALSE,type='ellipses')#,type='ellipses'  ,type='ellipses' ,type='ellipses'
dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Venn_Worm-Bacteria_resistance_overlap_MICover2.5_10perc.pdf",sep=''),width=9,height=9)



enbct<-qbacmicq[,c('Gene','UniACC','MIC','NGM_D')]
enbc<-rename(enbct, c('MIC'="logFC",'UniACC'='ID'))



Over25EC<-read.GO("Enrichment/Over25_chart_EC-bac.txt",pval='Benjamini',filter=FALSE)

circ<-circle_dat(Over25EC,enbc)

circm<-merge(circ,qbacmicq, by.x='genes',by.y='UniACC',all.x = TRUE)
circm$genes<-NULL
circm<-rename(circm, c('Gene'="genes"))
circm$logp<--log10(circm$adj_pval)
circm$term<-as.factor(circm$term)
circm$category<-as.factor(circm$category)
circm$strain<-ifelse(circm$category=='KP',substr(circm$ID,1,3),'eco')
circm$strain<-as.factor(circm$strain)

circms<-ddply(circm,.(category,strain,ID,term),summarise,
               MIC_avg_sum=mean(MIC,na.rm = TRUE),MIC_sd_sum=sd(MIC,na.rm = TRUE),
               MIC_med_sum=median(MIC,na.rm = TRUE),
               MIC_05_sum=quantile(MIC,probs=c(0.05),na.rm = TRUE),MIC_95_sum=quantile(MIC,probs=c(0.95),na.rm = TRUE),
               MIC_25_sum=quantile(MIC,probs=c(0.25),na.rm = TRUE),MIC_75_sum=quantile(MIC,probs=c(0.75),na.rm = TRUE),
               NGM_D_avg_sum=mean(NGM_D,na.rm = TRUE),NGM_D_sd_sum=sd(NGM_D,na.rm = TRUE),
               NGM_D_med_sum=median(NGM_D,na.rm = TRUE),
               NGM_D_05_sum=quantile(NGM_D,probs=c(0.05),na.rm = TRUE),NGM_D_95_sum=quantile(NGM_D,probs=c(0.95),na.rm = TRUE),
               NGM_D_25_sum=quantile(NGM_D,probs=c(0.25),na.rm = TRUE),NGM_D_75_sum=quantile(NGM_D,probs=c(0.75),na.rm = TRUE),
               NGM_C_avg_sum=mean(NGM_C,na.rm = TRUE),NGM_C_sd_sum=sd(NGM_C,na.rm = TRUE),
               NGM_C_med_sum=median(NGM_C,na.rm = TRUE),
               NGM_C_05_sum=quantile(NGM_C,probs=c(0.05),na.rm = TRUE),NGM_C_95_sum=quantile(NGM_C,probs=c(0.95),na.rm = TRUE),
               NGM_C_25_sum=quantile(NGM_C,probs=c(0.25),na.rm = TRUE),NGM_C_75_sum=quantile(NGM_C,probs=c(0.75),na.rm = TRUE),
               LB_sum=mean(LB_22hr,na.rm = TRUE),LB_sd_sum=sd(LB_22hr,na.rm = TRUE),
               LB_med_sum=median(LB_22hr,na.rm = TRUE),
               LB_05_sum=quantile(LB_22hr,probs=c(0.05),na.rm = TRUE),LB_95_sum=quantile(LB_22hr,probs=c(0.95),na.rm = TRUE),
               LB_25_sum=quantile(LB_22hr,probs=c(0.25),na.rm = TRUE),LB_75_sum=quantile(LB_22hr,probs=c(0.75),na.rm = TRUE),
               gcount=as.numeric(length(unique(genes))),logp_sum=mean(logp) )
enr<-merge(circm,circms,by=c('category','strain','ID','term'))

enrexp<-rename(enr,c('genes'='Gene'))[,!colnames(enr) %in% c('Comment','JW_id','Starving','logFC','zscore')]
write.csv(enrexp,'Data/Over2.5_EC-bac_Enrichment.csv')

Slevel=c(2,3,0)
SAbr=c('alpha=.01','alpha=.001','all')
Stat<-data.frame(slevel=Slevel,sabr=SAbr)


TName=c('BP','KP')
TTitle=c('Biological process', 'KEGG pathway')
TAbr=c('BP','KEGG-PWY')
Types<-data.frame(tname=TName,tabr=TAbr,ttitle=TTitle)

Media=c('LB_22hr','NGM_D','NGM_C')
MTitle=c('LB 22hr growth','NGM + 100uM 5FU growth','NGM growth')
MAbr=c('LB','NGM-5FU','NGM')
MOrder=c('LB_med_sum','NGM_D_med_sum','NGM_C_med_sum')
MQvars=c('LB','NGM_D','NGM_C')
MScale=c(1,0.3,0.3)
Medias<-data.frame(media=Media,mtitle=MTitle,mabr=MAbr,morder=MOrder,mscale=MScale,mqvars=MQvars)

##Growth distributions
for(m in 1:nrow(Medias)) {
  mr<-Medias[m,]
  for(s in 1:nrow(Stat)) {
    sr <- Stat[s,]
    for(t in 1:nrow(Types)) {
      tr <- Types[t,]
      sel<-subset(enr,category==as.character(tr$tname) & strain=='eco' & logp_sum>sr$slevel)
      #& logp_sum>2 & logp_sum>2  & logp_sum>2  & logp_sum>2
      gtitle<-paste('Enriched E. Coli MG1655 ',tr$ttitle,' - ',mr$mtitle,' distribution, ',sr$sabr,sep = '')
      fname<-paste('Figures/Enrichment/OD/',tr$tabr,'_BG_',mr$mabr,'_',sr$sabr,'.pdf',sep = '')
      print(gtitle)
      print(fname)
      termBG<-ggplot(sel,
                     aes(y=sel[,as.character(mr$media)],x=reorder(sel$term,sel[,as.character(mr$morder)]),fill=gcount,color=logp_sum))+
        geom_boxplot(position='identity')+coord_flip()
      termBG<-termBG+labs(fill='Number of genes',color='-log10(adj. p-value)')
      termBG<-termBG+scale_fill_gradient( high="red",low='white')
      termBG<-termBG+scale_color_gradient(limits=c(0,8), high="blue",low='gray')
      termBG<-termBG+ylab('OD')+ylim(0,mr$mscale)
      termBG<-termBG+xlab(tr$ttitle)+ggtitle(gtitle)
      #termBG
      cairo_pdf(fname,width=9,height=9)
      print(termBG)
      dev.off()
      #ggsave(plot=termBG,file=fname,width=9,height=9)
      #dev.copy2pdf(device=cairo_pdf,file=fname,width=12,height=12)
    }
  }
}

##MIC distributions
for(s in 1:nrow(Stat)) {
  sr <- Stat[s,]
  for(t in 1:nrow(Types)) {
    tr <- Types[t,]
    sel<-subset(enr,category==as.character(tr$tname) & strain=='eco' & logp_sum>sr$slevel)
    #& logp_sum>2 & logp_sum>2  & logp_sum>2  & logp_sum>2
    gtitle<-paste('Enriched E. Coli MG1655 ',tr$ttitle,' - MIC distribution, ',sr$sabr,sep = '')
    fname<-paste('Figures/Enrichment/MIC/',tr$tabr,'_MIC_',sr$sabr,'.pdf',sep = '')
    print(gtitle)
    print(fname)
    termMIC<-ggplot(sel,
                    aes(y=MIC,x=reorder(term,MIC_med_sum),fill=gcount,color=logp_sum))+
      geom_boxplot(position='identity')+coord_flip()
    termMIC<-termMIC+labs(fill='Number of genes',color='-log10(adj. p-value)')
    termMIC<-termMIC+scale_fill_gradient( high="red",low='white')
    termMIC<-termMIC+scale_color_gradient(limits=c(0,8), high="blue",low='gray')
    termMIC<-termMIC+ylab('MIC [5FU], uM')+ylim(0,100)
    termMIC<-termMIC+xlab(tr$ttitle)+ggtitle(gtitle)
    termMIC
    cairo_pdf(fname,width=9,height=9)
    print(termMIC)
    dev.off()
  }
}



Media2=c('LB_22hr','NGM_D','NGM_C')
MTitle2=c('LB 22hr growth','NGM + 100uM 5FU growth','NGM growth')
MAbr2=c('LB','NGM-5FU','NGM')
MOrder2=c('LB_med_sum','NGM_D_med_sum','NGM_C_med_sum')
MQvars2=c('LB','NGM_D','NGM_C')
MScale2=c(1,0.2,0.2)
Medias2<-data.frame(media=Media2,mtitle=MTitle2,mabr=MAbr2,morder=MOrder2,mscale=MScale2,mqvars=MQvars2)

##MIC vs Media
for(m in 1:nrow(Medias2)) {
  mr<-Medias2[m,]
  for(s in 1:nrow(Stat)) {
    sr <- Stat[s,]
    for(t in 1:nrow(Types)) {
      tr <- Types[t,]
#       if (mr$mqvars=='LB') {
#         sel<-subset(circms,category==as.character(tr$tname) & strain=='eco' & logp_sum>sr$slevel)
#       } else {
#         sel<-subset(circms,category==as.character(tr$tname) & strain=='eco' & logp_sum>sr$slevel)
#       }
      sel<-subset(enr,category==as.character(tr$tname) & strain=='eco' & logp_sum>sr$slevel)
      #& logp_sum>2 & logp_sum>2  & logp_sum>2  & logp_sum>2
      gtitle<-paste('Enriched E. Coli MG1655 ',tr$ttitle,': MIC',' - ',mr$mtitle,' distribution, ',sr$sabr,sep = '')
      fname<-paste('Figures/Enrichment/MICvsOD/',tr$tabr,'_MIC_vs_BG_',mr$mabr,'_',sr$sabr,'.pdf',sep = '')
      print(gtitle)
      print(fname)
      kpcor<-ggplot(sel,
                    aes(x=MIC_med_sum,y=sel[,paste(mr$mqvars,'_med_sum',sep='')],color=logp_sum))+
        geom_point(aes(size=gcount))
      kpcor<-kpcor+geom_errorbarh(aes(xmax=MIC_75_sum,
                                      xmin=MIC_25_sum),height=mr$mscale*0.0125,alpha=0.5)
      kpcor<-kpcor+geom_errorbar(aes(ymax=sel[,paste(mr$mqvars,'_75_sum',sep='')],
                                     ymin=sel[,paste(mr$mqvars,'_25_sum',sep='')]),width=1,alpha=0.5)
      kpcor<-kpcor+geom_text(aes(label=term),color='black',hjust=-0.1, vjust=-0.5,size=3)+
        scale_color_gradient(limits=c(0,8), high="blue",low='gray')+
        labs(size='Number of genes',color='-log10(adj. p-value)')+
        xlab('MIC [5FU], uM')+ylab('OD')+ylim(0,mr$mscale)+xlim(0,100)+
        ggtitle(gtitle)
      cairo_pdf(fname,width=9,height=9)
      print(kpcor)
      dev.off()
      #ggsave(plot=termBG,file=fname,width=9,height=9)
      #dev.copy2pdf(device=cairo_pdf,file=fname,width=12,height=12)
    }
  }
}

fitbac<-lm(NGM_D ~ NGM_C,data=qbacmicq)
confint(fitbac,'NGM_C',level=0.95)

for(s in 1:nrow(Stat)) {
  sr <- Stat[s,]
  for(t in 1:nrow(Types)) {
    tr <- Types[t,]
        
    sel<-subset(enr,category==as.character(tr$tname) & strain=='eco' & logp_sum>sr$slevel)
    #& logp_sum>2 & logp_sum>2  & logp_sum>2  & logp_sum>2
    gtitle<-paste('Enriched E. Coli MG1655 ',tr$ttitle,': NGM - NGM + 100um 5FU correlation, ',sr$sabr,sep = '')
    fname<-paste('Figures/Enrichment/NGMvsNGM-5FU/',tr$tabr,'_NGM_vs_NGM-5FU_',sr$sabr,'.pdf',sep = '')
    print(gtitle)
    print(fname)
    kpcor<-ggplot(sel,aes(x=NGM_C_med_sum,y=NGM_D_med_sum,color=logp_sum))+
      geom_point(aes(size=gcount))+
      geom_abline(aes(fill='1:1'),intercept=0,slope=1,alpha=0.5,color='grey',linetype='longdash')+
      geom_abline(aes(fill='Trend for all knockouts'),intercept=coef(fitbac)[[1]],slope=coef(fitbac)[[2]],alpha=0.5,color='red')
    kpcor<-kpcor+geom_errorbarh(aes(xmax=NGM_C_75_sum,xmin=NGM_C_25_sum),height=0.15*0.0125,alpha=0.5)
    kpcor<-kpcor+geom_errorbar(aes(ymax=NGM_D_75_sum,ymin=NGM_D_25_sum),width=0.2*0.0125,alpha=0.5)
    kpcor<-kpcor+geom_text(aes(label=term),color='black',hjust=-0.1, vjust=-0.5,size=3)+
      scale_color_gradient(limits=c(0,8), high="blue",low='gray')+
      labs(size='Number of genes',color='-log10(adj. p-value)',fill='Guides')+
      xlab('OD, NGM - Control')+ylab('OD, NGM + 100uM 5FU')+ylim(0,0.15)+xlim(0,0.2)+
      ggtitle(gtitle)
    cairo_pdf(fname,width=9,height=9)
    print(kpcor)
    dev.off()
  }
}



#ggsave(plot=termBG,file=fname,width=9,height=9)
#dev.copy2pdf(device=cairo_pdf,file=fname,width=12,height=12)



circgo<-subset(circm,category %in% c('KP','BP'))

GOBar(subset(circgo, category == 'KP'))
GOBar(circgo, display = 'multiple')
GOBar(subset(circgo,adj_pval>0.01), display = 'multiple')
GOBubble(circgo, title = 'Bubble plot', color = c('deeppink', 'blue', 'chartreuse1'), display = 'multiple', labels = 4)
GOBubble(circgo)
circgo<-subset(circgo,adj_pval<0.01 & count>4)
GOCircle(subset(circgo, category == 'BP'))
#dev.copy2pdf(device=cairo_pdf,file="GO/GO_AUCnorm_All.pdf",width=12,height=8)

BP<-unique(subset(circgo, category == 'BP' & count>4)[,c('ID','term','adj_pval','zscore')])$ID
CC<-unique(subset(circgo, category == 'CC' & count>4)[,c('ID','term','adj_pval','zscore')])$ID
GOCircle(circgo,nsub=c(CC[1],BP[1:9]))

GOCircle(subset(circgo,zscore<0 & category == 'BP'))
#dev.copy2pdf(device=cairo_pdf,file="GO/GO_AUCnorm_All.pdf",width=12,height=8)
GOCircle(subset(circgo,zscore>0 & category == 'BP'))
#dev.copy2pdf(device=cairo_pdf,file="GO/GO_Peak_UP.pdf",width=12,height=8)

#circ<-subset(circ,category=='BP')


MIC_top<-unique(circ[with(circ,order(logFC,decreasing=TRUE)),]$genes)[1:10]
genes_bottom<-unique(circ[with(circ,order(logFC,decreasing=FALSE)),]$genes)[1:10]


BP_top<-as.character(unique(subset(circgo, category == 'BP' & count>4)[,c('ID','term','adj_pval','zscore')]$term))
KP_top<-as.character(unique(subset(circgo, category == 'KP' & count>4)[,c('ID','term','adj_pval','zscore')]$term))


chord <- chord_dat(circgo, process=KP_top[1:10])
chord <- chord_dat(circgo, process=BP_top[1:10])
GOChord(chord, space = 0.02, gene.space = 0.25, gene.size = 5)
#dev.copy2pdf(device=cairo_pdf,file="GO/GO_AUCnorm_All_circle.pdf",width=16,height=16)

GOCluster(circ, KP_top[1:4], clust.by = 'term', term.width = 2)




