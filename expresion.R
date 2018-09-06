# ==============================================================
# Código para análisis de datos Mirna
# microarreglos de expresión Illumina
# Claudia Rangel Escareño
# Genómica Computacional 4 Agosto 2017
# Ultima modificación: 27 Agosto 2017 
# ==============================================================

library(affy) 
library(limma)
library(gtools)
library(gplots)
library(GO.db)
library(annotate)

data = read.table(file="DataNormalized.txt", header=TRUE, sep="\t")
etiquetas = read.table(file="etiquetas.txt", header=TRUE)
sampleInfo = data.table::fread(file="SamplesR.txt", header=TRUE)

sampleInfoexpresion = data.table::fread(file="../../Expresion/ComparacionesMetilacion-Expresion.csv",sep=",",header=T)
sampleInfoexpresion<-merge(sampleInfoexpresion,sampleInfo,by.x=c("ID Expresion"),by.y=c("Expediente"))


# Colores por grupo
colores = sampleInfoexpresion[,]
colores = apply(colores,2,as.character)

# Checando la correcta lectura de los datos
# ----------------------------------------------------------------
dim(data)

# Seleccionamos sólo aquellos genes cuyo p-value<0.5    
# de un total de 34656 tenemos p.values<.5 -> 27453
# ---------------------------------------------------
pval.ind = 55:104
pvalues = data[,pval.ind]
pvalues.avg = apply(pvalues,1,mean)
hist(pvalues.avg, col="steelblue", main="Histogram of p-values<0.5")
sum(pvalues.avg<.5)



# Con base a esos seleccionados ahora tomamos los 
# datos de niveles de expresión que estarán en la
# variable exp.matrix con columnas y renglones etiquetados
# ---------------------------------------------------
indx.rowpval = which(pvalues.avg<.5)
entID = data[indx.rowpval,"ENTREZ_GENE_ID"]
sym = data[indx.rowpval,1]
exp.matrix = data[indx.rowpval,-55:-104]
colnames(exp.matrix) = as.character(sampleInfoexpresion$Etiqueta)
rownames(exp.matrix) = sym

expval = apply(exp.matrix,2,min)
min(expval)
positive.exp.matrix = log2(exp.matrix+25)
exp.matrix = positive.exp.matrix
boxplot(exp.matrix, las=2,col=colores ,cex.axis=.7, main="Datos Normalizados", pch=".")


# Initial data after filtering out the p-values
# ----------------------------------------------
NoTranscripts = 27499
cont.list = read.table(file="groupsFile.txt", sep="\t", header=TRUE)
C = read.table(file="contrastsFile.txt", sep="\t", header=TRUE)


# initialize variables
# ---------------------
nr = dim(cont.list)[1]
nc = dim(cont.list)[2]
lfc.status = c()
B.status = c()
pVal.status = c()
all.contrasts = c()
B = 0
M = 0.3
pv = 2

for (i in 1:(dim(C)[1]) )
{
  # design matrix
  # --------------
  grpH = C[i,1]
  grpK = C[i,2]
  H = cbind(cont.list[grpH,3:nc])
  if (length(which(H==0))>0) {H = H[-which(H==0)]}
  K = cbind(cont.list[grpK,3:nc])
  if (length(which(K==0))>0) {K = K[-which(K==0)]}
  
  rws = (length(H) + length(K))
  design = matrix(rep(0,2*rws), nrow=rws)
  colnames(design) = c(as.character(cont.list$Grupos[grpH]),as.character(cont.list$Grupos[grpK]))
  design[1:(length(H)),1]=1
  design[(length(H)+1):rws,2]=1
  
  # contrasts matrix
  # -----------------
  myContrast = paste(as.character(cont.list$Grupos[grpH]), as.character(cont.list$Grupos[grpK]), sep =" - ")
  cont.matrix = makeContrasts(myContrast, levels=design)
  contraste = paste(as.character(cont.list$Grupos[grpH]), as.character(cont.list$Grupos[grpK]), sep ="_vs_")
  all.contrasts = c(all.contrasts, contraste)
  
  # linear model fit
  # -----------------
  array.indx = stack(cbind(H,K))[,-2,drop=FALSE]
  data.Cont = exp.matrix[,c(array.indx)$values]
  fit = lmFit(data.Cont,design)
  fitC = contrasts.fit(fit, cont.matrix)
  fitCB = eBayes(fitC)
  TT = topTable(fitCB, coef=1, adjust = "fdr", sort.by = "logFC", number=NoTranscripts, genelist=fit$genes)
  lfc.status = cbind(lfc.status,TT$logFC)
  B.status = cbind(B.status,TT$B)
  pVal.status = cbind(pVal.status,-log10(TT$P.Value))
  print(paste(contraste, "... done"))
  
  
  # Selection point
  # --------------------------------
  #selected = TT[(TT$B>B & abs(TT$logFC)>M),]
  indx.DE = (TT$P.Value<0.01 & abs(TT$logFC)>M)
  selected = TT[indx.DE,]
  if (nrow(selected)>1)
  {
    
    # outputs: Data
    # -----------------
    file.name = paste("Stats-",contraste,".csv", sep="")
    write.csv(selected, file=file.name)
    # documento .HTML
    file.name = paste(contraste,"html", sep=".")
    htmlpage(list(entID[indx.DE]), filename = file.name, title=paste("Differentially Expressed Genes:",contraste,sep=" "), othernames = cbind(rownames(selected),selected), table.head = c("Entrez ID", "Gene symbol",colnames(selected)), tabletable.center = TRUE)
    
    # output: Graphics
    # -----------------
    file.name = paste("VolcanoPlot-",contraste,".pdf", sep="")
    pdf(file=file.name,  height=6.5, width=11)
    x0 = min(lfc.status[,i]) -.5
    x1 = max(lfc.status[,i]) +.5
    y0 = min(pVal.status[,i]) -.5
    y1 = max(pVal.status[,i]) +.5
    plot(fitCB$coefficients, -log10(fitCB$p.value), pch=".", col="blue", ylim=c(y0,y1), xlim=c(x0,x1), main=contraste, cex.lab=1.3, ylab="p-value", xlab="log fold-change")
    par(new=T)
    abline(v=-M, col="brown", ylab="", xlab="")
    par(new=T)
    abline(v=M, col="brown", ylab="", xlab="")
    par(new=T)
    abline(h=pv, col="black", ylab="", xlab="")
    par(new=T)
    ind1 = abs(fitCB$coefficients)>M
    ind2 = (-log10(fitCB$p.value)) > pv
    ind3 = (fitCB$coefficients>M & (-log10(fitCB$p.value)) > pv)
    ind4 = (fitCB$coefficients< -M & (-log10(fitCB$p.value))> pv)
    x = as.matrix(fitCB$coef[ind1])
    y = as.matrix(-log10(fitCB$p.value)[ind1])
    plot(x, y, col="magenta",ylim=c(y0,y1), xlim=c(x0,x1),main="", pch = 20 ,xlab="", ylab="",cex.lab=1.2)
    x = as.matrix(fitCB$coef[ind2])
    y = as.matrix(-log10(fitCB$p.value)[ind2])
    par(new=T)
    plot(x, y, col="orange",  ylim=c(y0,y1), xlim=c(x0,x1), main="", pch = 20,xlab="", ylab="",cex.lab=1.2)
    x = as.matrix(fitCB$coef[ind3])
    y = as.matrix(-log10(fitCB$p.value)[ind3])
    par(new=T)
    plot(x, y, col="red",  ylim=c(y0,y1), xlim=c(x0,x1), main="", pch = 20,xlab="", ylab="",cex.lab=1.2)
    x = as.matrix(fitCB$coef[ind4])
    y = as.matrix(-log10(fitCB$p.value)[ind4])
    par(new=T)
    plot(x, y, col="darkgreen", ylim=c(y0,y1), xlim=c(x0,x1), main="", pch = 20,xlab="", ylab="",cex.lab=1.2)
    dev.off()
    
    # -------------------------------------			
    # Hierarchical Clustering and Heatmaps
    # -------------------------------------
    
    # heatmap with the only the microarrays involved 
    # variable data.Cont has ALL expression data rows but only the microarrays in play
    
    datos.clus = data.Cont[rownames(selected),]
    file.name = paste("Heatmap-",contraste,".pdf", sep="")
    pdf(file=file.name, height=14, width=8)
    par(oma = c(3,1,3,4),mar=c(12,5,2,2)+0.1)
    ind.hmap = heatmap.2(as.matrix(datos.clus), col=greenred(400), scale="row", key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=.55, cexCol=.8, main = contraste, labRow = rownames(selected))
    dev.off()
    Idx = ind.hmap$rowInd[length(ind.hmap$rowInd):1] 
    genes.hmap = cbind(selected[Idx,], datos.clus[Idx,])
    archivo = paste("GenesHeatmap",contraste,".csv", sep="_")
    write.csv(genes.hmap, file=archivo)
    
    print(paste(contraste, " plots and data ... done"))
  }
  else
  {print(paste(contraste, " empty contrast"))}
}

# Expression Matrix
# --------------------------------------------
write.csv(exp.matrix, file="MatrizExpresion_exp.csv")


# analysis of summary statistics per contrast
# --------------------------------------------
print("Summary statistcs... done")
colnames(lfc.status) = all.contrasts
colnames(B.status) = all.contrasts
colnames(pVal.status) = all.contrasts
pdf(file="SummaryStatistics.pdf", height=9, width=13)
par(mar= c(14,4,4,2) + 0.1)
boxplot(lfc.status, main="Log Fold-change distribution", las =2, col="wheat", cex = 0.5, cex.axis =.8)
boxplot(B.status, main="B-statistic distribution", las = 2, col="salmon",cex = 0.5, cex.axis =.8)
boxplot(pVal.status, main="P.Value distribution", las = 2, col="steelblue",cex = 0.5, cex.axis =.8)
par(mfrow=c(1,2))
for (k in 1:(dim(C)[1]))
{
  hist(lfc.status[,k], xlab="Log Fold-change values", main = all.contrasts[k], col="wheat")
  hist(B.status[,k], xlab="B-statistic values", main = "", col="salmon")
}
dev.off()
print("Process completed")









