library(stringr)
library(ggplot2)
library("dplyr")

DMR<-data.table::fread("DMR-DMRcate-Annotation-10abr-singuiones.tsv",sep="\t")

DMR$seqnames<-factor(str_sub(DMR$seqnames,4),levels = paste(1:22,sep=""))

ggplot(DMR,aes(x=start,y=meanbetafc,colour =ifelse(abs(meanbetafc)>0.01,"A","B")))+
  geom_point(size=1, shape=1)+ facet_grid(~seqnames,scale="free",switch="both")+
  theme(
    axis.line.y = element_line(color="black", size = 1),
    strip.text.x = element_text(size = 18),
    axis.text.y = element_text(face = "bold",                                                                                                 size = 18),panel.border = element_blank(), 
    panel.background = element_blank(), 
    panel.grid = element_blank(), 
    panel.spacing.x = unit(0.1,"line"),
    axis.text.x = element_blank(),
    axis.title.x=element_blank(),
    axis.ticks.x=element_blank()) +
  geom_hline(yintercept=0.01, linetype="solid", color = "red", size=1)+ 
  geom_hline(yintercept=-0.01,linetype="solid",color="red",size=1)+
  scale_color_manual(guide=FALSE, values=c("darkgreen", "darkblue"))