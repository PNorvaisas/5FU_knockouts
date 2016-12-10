packs<-c('Biobase','base','car','multcomp','ggplot2','gplots',"org.EcK12.eg.db","KEGGREST","GO.db")


write('',file="R_citations.bib",append=FALSE)
for (pac in packs) {
  x<-citation(pac)
  toBibtex(x)
  write(toBibtex(x),file="R_citations.bib",append=TRUE)
}