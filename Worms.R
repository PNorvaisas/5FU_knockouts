library('ggplot2')
library('gplots')
library('plyr')
library('reshape2')
library(tidyr)

scr2<-read.table('Secondary_screen_PN.csv',sep=',',quote = '"',header = TRUE)
scr2<-subset(scr2, ! Gene %in% c('XXXXXXX','no bacteria','empty'))

scr2$MIC.1_factor<-as.factor(scr2$MIC.1)
scr2$MIC.1<-as.numeric(as.character(scr2$MIC.1))
scr2$MIC.2<-as.numeric(as.character(scr2$MIC.2))
scr2$MIC_avg<-apply(scr2[,c('MIC.1','MIC.2')],1,mean)

corlin<-ggplot(scr2,aes(x=MIC.1,y=MIC.2))+geom_point(alpha=0.2,color='red',size=5)+stat_smooth(aes(group = 1),method = "lm")+xlab('MIC replicate 1')+ylab('MIC replicate 2')
corlin

scr2avg<-melt(scr2,id=c('MIC_avg'),measure.vars = c('MIC.1','MIC.2'),value.name='MIC_values')
scr2avg$Replicates<-factor(scr2avg$variable,
                       levels = c('MIC.1','MIC.2'),
                       labels = c("1","2")) 


scr2avg<-subset(scr2avg,!is.na(MIC_avg))
comparison<-ggplot(scr2avg,aes(x=MIC_avg,y=MIC_values,color=Replicates))+geom_point(alpha=0.2,size=5)+stat_smooth(aes(group = 1),method = "lm")+xlab('MIC average')+ylab('Experimental MIC measurements')+xlim(0,100)+ylim(0,100)
comparison<-comparison+ggtitle('Consistency of replicates in secondary screen (5-100uM 5FU)')
comparison
dev.copy2pdf(device=cairo_pdf,file="Consistency_scr1.pdf",width=11.69,height=8.27)

comphist<-ggplot(scr2avg,aes(x=as.factor(MIC_avg),y=MIC_values))+geom_boxplot()+xlab('MIC average')+ylab('Experimental MIC measurements')
comphist
dev.copy2pdf(device=cairo_pdf,file="Consistency_scr1_hist.pdf",width=11.69,height=8.27)


scr1<-read.table('Primary_screen_PN_clean.csv',sep=',',quote = '"',header = TRUE,stringsAsFactors=FALSE)
scr1<-subset(scr1, ! Gene %in% c('XXXXXXX','no bact','empty','',NA))
scr1[scr1$Gene %in% c('WT?','WT control', 'dodgy "WT"'),'Gene']<-'WT'

allg<-unique(scr1$Gene)
concs<-c(0,1,2.5,5)

mat<-matrix(ncol=5)
for (sdrug in c('X0','X1','X2.5','X5')){
  sconc=as.numeric(gsub('X','',sdrug))
  for (sworm in c(0,3,6,9)){
    for (edrug in c('X0','X1','X2.5','X5')){
      econc=as.numeric(gsub('X','',edrug))
      if (edrug!=sdrug & match(econc,concs)-match(sconc,concs)==1){
        for (eworm in c(0,3,6,9)){
          ln<-length(subset(scr1,scr1[,sdrug]==sworm & scr1[,edrug]==eworm)$Gene)
          mat<-rbind(mat,c(sconc,sworm,econc,eworm,ln))
        }
      }
    }
  }
}
links<-as.data.frame(mat)
colnames(links)<-c('Sdrug','Sworm','Edrug','Eworm','Links')
links<-links[!is.na(links$Links),]

occur<-sort(table(scr1$Gene))
dupl<-occur[occur > 1]



scr1m<-melt(scr1,id=c("Gene","Keio.Plate.no.","Position","Details","Faults"),measure.vars=c("X0","X1","X2.5","X5"),value.name='Worm.development')
scr1m$FFU<-as.numeric(gsub('X','',scr1m$variable))
scr1m$variable<-NULL
scr1m$Gene<-as.factor(scr1m$Gene)

#selection<-subset(scr1m, ! Gene %in% scr2$Gene | Gene=='WT')
#selection<-subset(scr1m,Gene=='WT')
#selection<-subset()
#selection<-subset(selection, Gene %in% scr1$Gene[1:100])
selection<-scr1m

#selection<-subset(scr1m,!Gene %in% rownames(dupl))
linkage<-ggplot(selection[!is.na(selection$Worm.development),],aes(x=FFU,y=Worm.development,color=Gene))+geom_line(alpha=0.5,size=1)
linkage+theme(panel.background = element_rect(fill = 'white', colour = 'white'))
dev.copy2pdf(device=cairo_pdf,file="Screen1_links.pdf",width=11.69,height=8.27)


t2<-scr2[,c('Gene','MIC_avg')]
colnames(t2)<-c('Gene','FFU')
t2$Worm.development<-0
t1<-scr1m[,c('Gene','FFU','Worm.development')]

joined<-merge(t1,t2,all.x=TRUE,all.y=TRUE)


fulllinkage<-ggplot(joined[!is.na(joined$Worm.development),],aes(x=FFU,y=Worm.development,color=Gene))+geom_line(alpha=0.5,size=1)
fulllinkage+theme(panel.background = element_rect(fill = 'white', colour = 'white'))
dev.copy2pdf(device=cairo_pdf,file="Joined_links.pdf",width=11.69,height=8.27)




chosen<-subset(joined,FFU>5)

fin<-subset(joined,Gene %in% subset(joined,FFU>5)$Gene & FFU==5)

comp<-merge(fin,chosen,by=c('Gene'),all.x=TRUE,all.y=TRUE)
comp$FFU.x<-NULL
comp$Worm.development.y<-NULL
comp$Worm.development.x<-factor(comp$Worm.development.x,levels=c('6','9'),labels=c('++','+++'))

sixandnine<-ggplot(comp[!is.na(comp$Worm.development.x),],aes(x=Worm.development.x,y=FFU.y))+geom_boxplot(notch=TRUE)+xlab('Worm development at [5FU]=5uM')+ylab('Final MIC')
sixandnine<-sixandnine+ggtitle('Final MIC distributions based on worm development at [5FU]=5uM')
sixandnine
dev.copy2pdf(device=cairo_pdf,file="Final_MIC_by_5FUat5.pdf",width=11.69,height=8.27)


allg<-unique(scr1$Gene)

zeros<-subset(scr1,scr1$'X0'==9 & scr1$'X1'==0 & (scr1$'X2.5'==0 | scr1$'X2.5'==3)&  scr1$'X5'==0)$Gene #scr1$'X2.5'==0 & 
zeros<-unique(zeros)
twofives<-subset(scr1,scr1$'X0'==9 & (scr1$'X1'==3 | scr1$'X1'==6 | scr1$'X1'==9 )& scr1$'X2.5'==0 &  scr1$'X5'==0)$Gene
twofives<-unique(twofives)
fivers<-subset(scr1,scr1$'X0'==9 & scr1$'X1'>0 & scr1$'X2.5'>0 &  scr1$'X5'==0)$Gene
fivers<-unique(fivers)


mc<-data.frame('Gene'=unique(scr1$Gene))

mics<-merge(mc,scr2[,c('Gene','MIC_avg')],by=c('Gene'),all.x=TRUE,all.y=TRUE)
colnames(mics)<-c('Gene','MIC')
mics[mics$Gene %in% zeros,'MIC' ]<-1
mics[mics$Gene %in% twofives,'MIC' ]<-2.5
mics[mics$Gene %in% fivers,'MIC' ]<-5

dist<-ggplot(mics,aes(x=MIC))+stat_ecdf()+geom_hline(yintercept=0.90,color='green')+geom_hline(yintercept=0.95,color='blue')+geom_hline(yintercept=0.98,color='red')+ scale_x_continuous(breaks=1:100)#+geom_histogram(binwidth=5)
dist<-dist+theme(axis.text.x = element_text(angle = 90, hjust = 1))+ggtitle('Cumulative histogram of MIC values')
dist<-dist+annotate("text", 110, 0.91, label = "90%",color='green')
dist<-dist+annotate("text", 110, 0.96, label = "95%",color='blue')
dist<-dist+annotate("text", 110, 0.99, label = "98%",color='red')
dist+scale_x_continuous(breaks=c(0,10,20,30,40,50,60,70,80,90,100))
dev.copy2pdf(device=cairo_pdf,file="Cumulative_distribution_of_MIC.pdf",width=11.69,height=8.27)

atop5<-mics[mics$MIC>=5,]$Gene
atop5<-unique(atop5[!is.na(atop5)])

atop10<-mics[mics$MIC>=10,]$Gene
atop10<-unique(atop10[!is.na(atop10)])

atop25<-mics[mics$MIC>=25,]$Gene
atop25<-unique(atop25[!is.na(atop25)])

atop30<-mics[mics$MIC>=30,]$Gene
atop30<-unique(atop30[!is.na(atop30)])

atop40<-mics[mics$MIC>=40,]$Gene
atop40<-unique(atop40[!is.na(atop40)])

write.csv(allg,'All.csv')
write.csv(atop10,'Top_above10.csv')
write.csv(atop25,'Top_above25.csv')
write.csv(atop30,'Top_above30.csv')
write.csv(atop40,'Top_above40.csv')

zeroto6<-unique(subset(scr1,scr1$'X5'==6 & scr1$'X2.5'==0)$Gene)
zeroto3<-unique(subset(scr1,scr1$'X5'==3 & scr1$'X2.5'==0)$Gene)

intersect(atop30,zeroto6)
intersect(atop25,zeroto6)

naat5<-unique(subset(scr1,is.na(scr1$'X5'))$Gene)
zerozeroand6<-unique(subset(scr1,scr1$'X5'==6 & scr1$'X1'==0 & scr1$'X2.5'==0))
intersect(scr2$Gene,naat5)


scr169<-unique(subset(scr1,scr1$'X5'==9 | scr1$'X5'==6)$Gene)
scr2all<-unique(scr2$Gene)
length(scr169)
length(scr2all)
length(intersect(scr169,scr2all))

na25s<-unique(subset(scr1,scr1$'X5'==0 & is.na(scr1$'X1') & is.na(scr1$'X2.5') & scr1$'X0'==9)$Gene)

zero125<-unique(subset(scr1,scr1$'X1'==0 & scr1$'X2.5'==0 & is.na(scr1$'X5'))$Gene)
zero15<-unique(subset(scr1,scr1$'X1'==0 & scr1$'X5'==0 & is.na(scr1$'X2.5'))$Gene)
zero255<-unique(subset(scr1,scr1$'X2.5'==0 & scr1$'X5'==0)$Gene)

redoNA<-setdiff(unique(subset(scr1,is.na(scr1$'X5'))$Gene), zero125)
zeroby2<-union(zero125,zero15 )
zeroby2<-unique(union(zeroby2,zero255 ))
redoNS<-setdiff(unique(subset(scr1,scr1$'X5'==0)$Gene),zeroby2)
add3<-unique(subset(scr1,scr1$'X5'==3)$Gene)
extension<-union(redoNS,redoNA)
extension<-union(extension,missed)
extension<-unique(union(extension,add3))
missed<-setdiff(scr169,scr2all)

write.csv(missed,'Missed69.csv')
write.csv(redoNA,'Redo_NA.csv')
write.csv(redoNS,'Redo_NS.csv')
write.csv(add3,'Onestar.csv')
write.csv(extension,'All_extension.csv')

