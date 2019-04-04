#methylation analysis
library(ggplot2)
library(limma)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(RColorBrewer)
library(missMethyl)
library(matrixStats)
library(minfiData)
library(Gviz)
library(DMRcate)
library(stringr)
library(wateRmelon)

targets<-read.csv("Samples.csv")
#targets_21mayo<-read.csv("Sample_Sheet-21Mayo.csv")
#targets_21mayo<-targets_21mayo[match(targets$Sample_Name.1,targets_21mayo$Sample_Name),]

rgSet <- read.metharray.exp(targets=targets)

targets$ID <- paste(targets$Condicion,targets$Sample_Name,sep=".")


sampleNames(rgSet) <- targets$ID

detP <- detectionP(rgSet)

pal <- brewer.pal(8,"Dark2")

pdf("pvalues_methylation.pdf")
par(mfrow=c(1,1),las=2)
barplot(colMeans(detP), col=pal[factor(targets$Sample_Group)],
 cex.names=0.8,ylab="Mean detection p-values",ylim=c(0,0.015))

abline(h=0.01,col="red")
legend("topleft", legend=levels(factor(targets$Sample_Group)), fill=pal,
 bg="white")


barplot(colMeans(detP), col=pal[factor(targets$Sample_Group)], 
 cex.names=0.8, ylim = c(0,0.002), ylab="Mean detection p-values")
legend("topleft", legend=levels(factor(targets$Sample_Group)), fill=pal,
 bg="white")

dev.off()

qcReport(rgSet, sampNames=targets$ID, sampGroups=targets$Sample_Group,
 pdf="qcReport.pdf")

mSetSq <- preprocessQuantile(rgSet)

mSetRaw <- preprocessRaw(rgSet)

pdf("density.pdf")
par(mfrow=c(1,1))
densityPlot(rgSet, sampGroups=targets$Sample_Group,main="Raw", legend=FALSE)
legend("top", legend = levels(factor(targets$Sample_Group)),
 text.col=brewer.pal(8,"Dark2"))
densityPlot(getBeta(mSetSq), sampGroups=targets$Sample_Group,
 main="Normalized", legend=FALSE)
legend("top", legend = levels(factor(targets$Sample_Group)),
 text.col=brewer.pal(8,"Dark2"))
dev.off()


pdf("mds.pdf")
par(mfrow=c(1,1))
plotMDS(getM(mSetSq), top=1000, gene.selection="common",pch=10,
 col=pal[factor(targets$Sample_Group)])
legend("top", legend=levels(factor(targets$Sample_Group)), text.col=pal,
 bg="white", cex=0.7)
plotMDS(getM(mSetSq), top=1000, gene.selection="common",pch=10,
 col=pal[factor(targets$Slide)])
legend("top", legend=levels(factor(targets$Slide)), text.col=pal,
 bg="white", cex=0.7)

dev.off()


plotMDS(getM(mSetSq), top=1000, gene.selection="common",pch=10,
 col=pal[targets_new$TS1])
legend("top", legend=levels(targets_new$TS1), text.col=pal,
 bg="white", cex=0.7)




pdf("mds.3.pdf")
par(mfrow=c(1,3))
plotMDS(getM(mSetSq), top=1000, gene.selection="common", pch=10,
 col=pal[factor(targets$Sample_Group)], dim=c(1,3))
legend("top", legend=levels(factor(targets$Sample_Group)), text.col=pal,
 cex=0.7, bg="white")
plotMDS(getM(mSetSq), top=1000, gene.selection="common", pch=10,
 col=pal[factor(targets$Sample_Group)], dim=c(2,3))
legend("topleft", legend=levels(factor(targets$Sample_Group)), text.col=pal,
 cex=0.7, bg="white")
plotMDS(getM(mSetSq), top=1000, gene.selection="common",pch=10,
 col=pal[factor(targets$Sample_Group)], dim=c(3,4))
legend("topright", legend=levels(factor(targets$Sample_Group)), text.col=pal,
 cex=0.7, bg="white")
dev.off()

detP <- detP[match(featureNames(mSetSq),rownames(detP)),]

keep <- rowSums(detP < 0.01) == ncol(mSetSq)
table(keep)

mSetSqFlt <- mSetSq[keep,]

mSetSqFlt <- dropLociWithSnps(mSetSqFlt)

xReactiveProbes <- read.csv(file="~/Downloads/illumina450k_filtering-master/48639-non-specific-probes-Illumina450k.csv")
keep <- !(featureNames(mSetSqFlt) %in% xReactiveProbes$TargetID)
table(keep)

mSetSqFlt <- mSetSqFlt[keep,]

mSetSqFlt
#remove non-specific

mVals <- getM(mSetSqFlt)

bVals <- getBeta(mSetSqFlt)

#Annotation process
ann450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450kSub <- ann450k[match(rownames(bVals),ann450k$Name), c(1:4,12:19,24:ncol(ann450k))]



pdf("B-M_Values.pdf")
par(mfrow=c(1,1))
densityPlot(bVals, sampGroups=targets$Sample_Group, main="Beta values",
 legend=FALSE, xlab="Beta values")
legend("top", legend = levels(factor(targets$Sample_Group)),
 text.col=brewer.pal(8,"Dark2"))
densityPlot(mVals, sampGroups=targets$Sample_Group, main="M-values",
 legend=FALSE, xlab="M values")
legend("topleft", legend = levels(factor(targets$Sample_Group)),
 text.col=brewer.pal(8,"Dark2"))
dev.off()



###################DMPs analysis

cpg.annotate.pvalue<- function(datatype = c("array", "sequencing"), object, annotation = c(array = "IlluminaHumanMethylation450k",
    annotation = "ilmn12.hg19"), analysis.type = c("differential",
    "variability"), design, contrasts = FALSE, cont.matrix = NULL,
    fdr = 0.05, coef, ...){ 
    if (datatype == "array") {
        stopifnot(is.matrix(object))
        analysis.type <- match.arg(analysis.type)
        switch(analysis.type, differential = {
            stopifnot(is.matrix(design))
            if (!contrasts) {
                stopifnot(colnames(design)[1] == "(Intercept)")
            } else {
                stopifnot(!is.null(cont.matrix))
            }
            fit <- lmFit(object, design, ...)
            if (contrasts) {
                stopifnot(coef %in% colnames(cont.matrix))
                fit <- contrasts.fit(fit, cont.matrix)
            }
            fit <- eBayes(fit)
            tt <- topTable(fit, coef = coef, number = nrow(object))
            nsig <- sum(tt$P.Value < fdr)
            if (nsig == 0) {
                message("Your contrast returned no individually significant probes. Set pcutoff manually in dmrcate() to return DMRs, but be warned there is an 
increased risk of Type I errors.")
            }
            if (nsig > 0 & nsig <= 100) {
                message(paste("Your contrast returned", nsig,
                  "individually significant probes; a small but real effect. Consider manually setting the value of pcutoff to return more DMRs, but be warned t
hat doing this increases the risk of Type I errors."))
            }
            if (nsig > 100) {
                message(paste("Your contrast returned", nsig,
                  "individually significant probes. We recommend the default setting of pcutoff in dmrcate()."))
            }
            betafit <- lmFit(ilogit2(object), design, ...)
            if (contrasts) {
                betafit <- contrasts.fit(betafit, cont.matrix)
            }
            betafit <- eBayes(betafit)
            betatt <- topTable(betafit, coef = coef, number = nrow(object))
            m <- match(rownames(tt), rownames(betatt))

             tt$betafc <- betatt$logFC[m]
            m <- match(rownames(tt), rownames(object))
            object <- object[m, ]
            RSobject <- RatioSet(object, annotation = annotation)
            RSanno <- getAnnotation(RSobject)
            stat <- tt$t
            annotated <- data.frame(ID = rownames(object), stat = stat,
                CHR = RSanno$chr, pos = RSanno$pos, betafc = tt$betafc,
                indfdr = tt$P.Value, is.sig = tt$P.Value <
                  fdr)
        }, variability = {
            RSobject <- RatioSet(object, annotation = annotation)
            RSanno <- getAnnotation(RSobject)
            wholevar <- var(object)
            weights <- apply(object, 1, var)
            weights <- weights/mean(weights)
            annotated <- data.frame(ID = rownames(object), stat = weights,
                CHR = RSanno$chr, pos = RSanno$pos, betafc = rep(0,
                  nrow(object)), indfdr = rep(0, nrow(object)),
                is.sig = weights > quantile(weights, 0.95))
        })
        annotated <- annotated[order(annotated$CHR, annotated$pos),
            ]
        class(annotated) <- "annot"
        return(annotated)
    }
    if (datatype == "sequencing") {
        if (!all(c("stat", "chr", "pos", "diff", "fdr") %in%
            colnames(object)))
            stop("Error: object does not contain all required columns, was it created by DSS::DMLtest()? Must contain colNames 'stat', 'chr', 'pos', 'diff' and 'fdr'.")
        annotated <- data.frame(ID = rownames(object), stat = object$stat,
            CHR = object$chr, pos = object$pos, betafc = object$diff,
            indfdr = object$fdr, is.sig = object$fdr < fdr)
        annotated <- annotated[order(annotated$CHR, annotated$pos),
            ]
        class(annotated) <- "annot"
    }
    else {
message("Error: datatype must be one of 'array' or 'sequencing'")
    }
    return(annotated)
} 

contraste <- "Suicida-Control"
casos <- targets$Condicion=="Suicida"
control <-targets$Condicion=="Control"
data_matrix <- bVals

design_matrix<-matrix(0,dim(data_matrix)[2],2)
design_matrix[control,2]<-1
design_matrix[casos,1]<-1

colnames(design_matrix)<-strsplit(contraste,"-")[[1]]

rownames(design_matrix)<-colnames(data_matrix)

fit <- lmFit(data_matrix, design_matrix)

contMatrix <- makeContrasts(contrasts=contraste, levels=design_matrix)
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)
# look at the numbers of DM CpGs at FDR < 0.05
print(summary(decideTests(fit2)))
DMPs <- topTable(fit2, num=Inf, coef=1, genelist=ann450kSub,adjust.method = "fdr")

head(DMPs)

myAnnotation <- cpg.annotate(object = bVals, datatype = "array", what = "Beta",
 analysis.type = "differential", design = design_matrix,
 contrasts = TRUE, cont.matrix = contMatrix,
 coef = "Suicida-Control", arraytype = "450K",fdr=0.1)




DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
head(DMRs$results)



results.ranges <- extractRanges(DMRs, genome = "hg19")
# set up the grouping variables and colours
groups <- pal[1:length(unique(targets_new$Condicion))]
names(groups) <- levels(factor(targets_new$Condicion))
cols <- groups[as.character(factor(targets_new$Condicion))]
samps <- 1:nrow(targets_new)
# draw the plot for the top DMR
DMR.plot(ranges=results.ranges, dmr=1, CpGs=bVals, phen.col=cols, what = "Beta",
 arraytype = "450K", pch=10, toscale=TRUE, plotmedians=TRUE,
 genome="hg19", samps=samps)



## GO term enrichment

sigCpGs <- DMPs$Name[DMPs$P.Value<0.00005]

length(sigCpGs)

all <- DMPs$Name

length(all)

par(mfrow=c(1,1))
gst <- gometh(sig.cpg=sigCpGs, all.cpg=all, plot.bias=TRUE)


sigCpGs <- DMPs$Name[DMPs$P.Value<0.05]
all <- DMPs$Name
par(mfrow=c(1,1))
gst <- gometh(sig.cpg=sigCpGs, all.cpg=all, plot.bias=TRUE)



topGO(gst, number=10)



##Age calculation
edad<-agep(bVals)


df.edad<-data.frame(targets_21mayo,DNAm=edad)

pdf("edad.pdf")
p.edad<-ggplot(aes(x=Edad,y=DNAm,fill=Condicion),data=df.edad)+
    geom_point(aes(colour=Condicion))+theme_classic()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + guides(fill=FALSE)


print(p.edad+geom_smooth(method='lm',se=FALSE,aes(colour=Condicion))+ labs(colour = "Condition",x="Age"))

p.edad<-ggplot(aes(x=Edad,y=DNAm,fill=Inflamacion),data=df.edad)+
    geom_point(aes(colour=Inflamacion))+theme_classic()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + guides(fill=FALSE)


print(p.edad+geom_smooth(method='lm',se=FALSE,aes(colour=Inflamacion))+ labs(colour = "Medical Condition",x="Age"))

p.edad<-ggplot(aes(x=Edad,y=DNAm,fill=Tpatolo),data=df.edad)+
    geom_point(aes(colour=Tpatolo))+theme_classic()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + guides(fill=FALSE)


print(p.edad+geom_smooth(method='lm',se=FALSE,aes(colour=Tpatolo))+ labs(colour = "Pathology",x="Age"))

p.edad<-ggplot(aes(x=Edad,y=DNAm,fill=Tsus),data=df.edad)+
    geom_point(aes(colour=Tsus))+theme_classic()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + guides(fill=FALSE)


print(p.edad+geom_smooth(method='lm',se=FALSE,aes(colour=Tsus))+ labs(colour = "Substance abuse",x="Age"))

p.edad<-ggplot(aes(x=Edad,y=DNAm,fill=Tpsic),data=df.edad)+
    geom_point(aes(colour=Tpsic))+theme_classic()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + guides(fill=FALSE)


print(p.edad+geom_smooth(method='lm',se=FALSE,aes(colour=Tpsic))+ labs(colour = "Psychiatric conditions",x="Age"))


p.edad<-ggplot(aes(x=Edad,y=DNAm,fill=Condicion),data=df.edad)+
    geom_point(aes(colour=Tpatolo))+theme_classic()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + guides(fill=FALSE)


print(p.edad+geom_smooth(method='lm',se=FALSE)+ labs(x="Age"))


dev.off()


df.edad%>% group_by(Tpsic) %>% do(lm(DNAm~Edad,data=.) %>% coef %>% as_data_frame)


##Cell count
cellCounts <- estimateCellCounts(rgSet)

df.cell<-melt(cbind(cellCounts,targets_21mayo[,-4]))

pdf("cell.pdf")
p.cell<-ggplot(aes(x=variable,y=value,fill=Condicion),data=df.cell)+
    geom_point(aes(colour=Condicion))+theme_classic()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + guides(fill=FALSE)


print(p.cell+labs(colour = "Condition",x="Cell type"))

p.cell<-ggplot(aes(x=variable,y=value,fill=Condicion),data=df.cell)+
    geom_boxplot(aes(colour=Condicion))+theme_classic()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + guides(fill=FALSE)


print(p.cell+labs(colour = "Condition",x="Cell type"))


p.cell<-ggplot(aes(x=variable,y=value,fill=Inflamacion),data=df.cell)+
    geom_boxplot(aes(colour=Inflamacion))+theme_classic()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + guides(fill=FALSE)


print(p.cell+labs(colour = "Medical Condition",x="Cell type"))


p.cell<-ggplot(aes(x=variable,y=value,fill=Tpatolo),data=df.cell)+
    geom_boxplot(aes(colour=Tpatolo))+theme_classic()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + guides(fill=FALSE)


print(p.cell+labs(colour = "Pathology",x="Cell type"))


p.cell<-ggplot(aes(x=variable,y=value,fill=Tsus),data=df.cell)+
    geom_boxplot(aes(colour=Tsus))+theme_classic()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + guides(fill=FALSE)


print(p.cell+labs(colour = "Substance abuse",x="Cell type"))


p.cell<-ggplot(aes(x=variable,y=value,fill=Tpsic),data=df.cell)+
    geom_boxplot(aes(colour=Tpsic))+theme_classic()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + guides(fill=FALSE)


print(p.cell+labs(colour = "Psychiatric conditions", x="Cell type"))

dev.off()