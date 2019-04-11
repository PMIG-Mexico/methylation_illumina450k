library("ChAMP")
library("minfi")
library("wateRmelon")

myLoad <- ChAMP::champ.load("./data/methylation/",arraytype="450K")

myLoad$pd$Slide<-as.factor(myLoad$pd$Slide)

champ.QC(beta = myLoad$beta,
         pheno=myLoad$pd$Sample_Group,
         mdsPlot=TRUE,
         densityPlot=TRUE,
         dendrogram=TRUE,
         PDFplot=TRUE,
         Rplot=TRUE,
         Feature.sel="None",
         resultsDir="./CHAMP_QCimages_prenorm_/")

#Normalization
myNorm=champ.norm(beta=myLoad$beta,rgSet=myLoad$rgSet,mset=myLoad$mset,resultsDir="./CHAMP_Normalization/",arraytype="450K", cores=4)

#QC post norm
champ.QC(beta = myNorm,
         pheno=myLoad$pd$Sample_Group,
         mdsPlot=TRUE,
         densityPlot=TRUE,
         dendrogram=TRUE,
         PDFplot=TRUE,
         Rplot=TRUE,
         Feature.sel="None",
         resultsDir="./CHAMP_QC_postnorm/")

#batcheffects 
#plot SVD plot to analyze batch effects
champ.SVD(beta= myNorm, rgSet= myLoad$rgSet, pd= myLoad$pd, resultsDir = "./CHAMP_SVDimages_norm/")

#array
combat_a=champ.runCombat(beta=myNorm, pd= myLoad$pd, batchname=c("Array","Slide"))
champ.SVD(beta= combat_a, rgSet= myLoad$rgSet, pd= myLoad$pd, resultsDir = "./CHAMP_SVDimages_array/")
#slide
combat_a_b1=champ.runCombat(beta=combat_a, pd=myLoad$pd, batchname=c("Slide"))
champ.SVD(beta= combat_a_b1, rgSet= myLoad$rgSet, pd= myLoad$pd, resultsDir = "./CHAMP_SVDimages_array_slide/")

#findingDMPs
myDMP=champ.DMP(beta=combat_a_b1, arraytype="450K")

#Finding DMRs
DMRsimple=champ.DMR(beta=combat_a_b1,pheno=myLoad$pd$Sample_Group)

myDMR=champ.DMR(beta=combat_a_b1,pheno=myLoad$pd$Sample_Group,arraytype="450K",method = "Bumphunter",minProbes=7, adjPvalDmr= 0.1, cores=2, smooth=TRUE, useWeights=FALSE,permutations=NULL,B=1000, nullMethod="bootstrap",meanLassoRadius=375,minDmrSep=1000,minDmrSize=50,adjPvalProbe=0.05,Rplot=T,PDFplot=T,resultsDir="./CHAMP_ProbeLasso/")
DMR=myDMR$BumphunterDMR
DMR$seqnames <- as.factor(substr(DMR$seqnames, 
                                 4, 100))
index <- apply(DMR, 1, function(x) which(probe.features$CHR == 
                                           x[1] & probe.features$MAPINFO >= as.numeric(x[2]) & probe.features$MAPINFO <= 
                                           as.numeric(x[3])))
Anno <- data.frame(DMRindex = unname(unlist(sapply(names(index), 
                                                   function(x) rep(x, length(index[[x]]))))), probe.features[do.call(c, 
                                                                                                                     index), 1:8])


#DMR.plot(ranges, dmr, CpGs, what=c("Beta", "M"),
#         arraytype=c("EPIC"), phen.col,
#         genome = c("hg19"),
#         samps = NULL, ...)
#write.csv(myDMPs, "myDMPS_0.1_BH_second")


write.csv(myDMRS_saas, "myDMRS_saas")
myCombatmyCombat...[,grep("RES",pheno$group)]
pheno=myLoad$pd$Sample_Group[,grep("RES",myLoad$pd$Sample_Group)]
beta=myCombat[,grep("RES",myLoad$pd$Sample_Group)]

champ.EpiMod(beta=myCombat_slide_array_plate_age,pheno=myLoad_NRES$pd$Sample_Group)
saveRDS(myLoad, "myLoad_NRES_March22")

champ.EpiMod(beta=myCombat_slide_array_plate_age_NRES[,grep("RES",myLoad$pd$Sample_Group)],pheno=myLoad$pd[grep("RES",myLoad$pd$Sample_Group), ]$Sample_Group, arraytype="EPIC")



GSEA_CD_goseq=champ.GSEA(beta=myCombat_slide_array_plate_age,DMP=myDMPS_C_D, DMR=myDMRs_C_D_Bump, method="goseq", arraytype="EPIC")
write.csv(GSEA_CD_goseq, "GSEA_CD_goseq")

GSEA_RES_goseq=champ.GSEA(beta=myCombat_slide_array_plate_age_NRES[,grep("RES",myLoad$pd$Sample_Group])],DMP=myDMP, DMR=myDMR, method="goseq", arraytype="EPIC")
write.csv(GSEA_RES_goseq, "GSEA_RES_goseq")devtools::install_github('reptalex/phylofactor')
