##10_importKangArray.R


rm(list=ls())
options(stringsAsFactors = F)

library(biomaRt); library(WGCNA); library(pSI)

setwd("~/Github/BrainSpan//")


##Get Kang et al.,

if(file.exists("./raw_data/Kang_GSE25219/Kang_GSE25219_raw.RData")) {
  load("./raw_data/Kang_GSE25219/Kang_GSE25219_raw.RData")
} else {
  library(GEOquery)
  gsedat = getGEO("GSE25219", GSEMatrix = T)
  datExpr = exprs(gsedat[[1]])
  datMeta = pData(phenoData(gsedat[[1]]))
  datProbes = pData(featureData(gsedat[[1]]))
  save(file="./raw_data/Kang_GSE25219/Kang_GSE25219_raw.RData", datExpr, datMeta,datProbes)
}


for (i in 1:nrow(datProbes)) {
  grout <- regexpr("gene:",as.character(datProbes[i,"mrna_assignment"]),fixed=TRUE)
  datProbes$ensg[i] <- substr(as.character(datProbes[i,"mrna_assignment"]),grout[1]+5,grout[1]+19)
}

keep <- substr(datProbes$ensg,1,4)=="ENSG"
datExpr <- datExpr[keep,]
datProbes <- datProbes[keep,]

collapseDat <- collapseRows(datET = datExpr,rowGroup = datProbes$ensg,rowID = rownames(datExpr))
datExpr <- collapseDat$datETcollapsed
datProbes <- datProbes[match(rownames(datExpr), datProbes$ensg),]
all(colnames(datExpr)==rownames(datMeta)) ## Samples match
all(rownames(datExpr)==datProbes$ensg) ## Samples match




datMeta$Subject = gsub("brain code: ", "", datMeta$characteristics_ch1)
datMeta$Region = gsub("region: ", "", datMeta$characteristics_ch1.1)
datMeta$Sex[datMeta$characteristics_ch1.3 == "Sex: M"] = "M"; datMeta$Sex[datMeta$characteristics_ch1.3 == "Sex: F"] = "F";  datMeta$Sex = factor(datMeta$Sex)
datMeta$Age = (gsub("age: ", "", datMeta$characteristics_ch1.4))
datMeta$Stage = as.numeric(gsub("Stage: ", "", datMeta$characteristics_ch1.5))
datMeta$PMI = as.numeric(gsub("postmortem interval: ", "", datMeta$characteristics_ch1.6))
datMeta$pH = as.numeric(gsub("ph: ", "", datMeta$characteristics_ch1.7))
datMeta$RIN = as.numeric(gsub("rna integrity number: ", "", datMeta$characteristics_ch1.8))


num=as.numeric(substr(datMeta$Age,0,2))
prenatal = grepl("PCW", datMeta$Age)
postnatal = !prenatal
yrs = grepl("Y", datMeta$Age)
mos = grepl("PMonth",datMeta$Age)

datMeta$Age_years[yrs] = num[yrs]
datMeta$Age_years[prenatal] = -(40 - num[prenatal])/40
datMeta$Age_years[mos] = num[mos]/12

num[mos] = num[mos]/12
datMeta$Age_unit = num
datMeta$Unit = "Postnatal (Years)"; datMeta$Unit[prenatal] = "Prenatal (Post-Conception Week)"


br = datMeta$Region
datMeta$Region_broad[br=="AMY"] = "Amygdala"
datMeta$Region_broad[br=="A1C" | br == "ITC"| br == "STC" | br == "TC"] = "Temporal Cortex"
datMeta$Region_broad[br=="CB" | br=="CBC"| br == "URL"] = "Cerebellum"
datMeta$Region_broad[br=="CGE" | br=="LGE" | br =="MGE" | br == "STR" | br=="VF"] = "Striatum"
datMeta$Region_broad[br=="DFC" | br == "MFC" | br == "PC"| br == "VFC" | br=="M1C" | br == "OFC" | br=="FC"] = "Frontal Cortex"
datMeta$Region_broad[br=="DTH" | br == "MD" | br=="DIE"] = "Thalamus"
datMeta$Region_broad[br=="HIP"] = "Hippocampus"
datMeta$Region_broad[br=="IPC"| br == "S1C" | br == "PC"] = "Parietal Cortex"
datMeta$Region_broad[br=="OC"| br == "V1C"] = "Occipital Cortex"
datMeta$Region_broad = as.factor(datMeta$Region_broad)

datMeta$Period = NA
datMeta$Period[datMeta$Stage <=7 & datMeta$Stage >= 3] = "fetal"
datMeta$Period[datMeta$Stage <=12 & datMeta$Stage >= 8] = "postnatal"
datMeta$Period[datMeta$Stage <=15 & datMeta$Stage >= 13] = "adult"
datMeta$Period = factor(datMeta$Period, levels=c("fetal", "postnatal", "adult"))


kang = vector(mode="list", length=3)
kang$datExpr = datExpr
kang$datMeta = datMeta
kang$datProbes = datProbes


regions = c("Cortex", "Amygdala","Hippocampus", "Thalamus", "Striatum", "Cerebellum")
datExpr.regional = as.data.frame(matrix(NA,ncol=length(regions),nrow=nrow(datExpr)))
colnames(datExpr.regional)= regions
rownames(datExpr.regional) = rownames(datExpr)

datExpr.regional = as.data.frame(matrix(NA,ncol=length(regions),nrow=nrow(datExpr)))
colnames(datExpr.regional)= regions
rownames(datExpr.regional) = rownames(datExpr)

datExpr.regionAndPeriod = as.data.frame(matrix(NA,ncol=length(regions),nrow=nrow(datExpr)))

datExpr.qn = normalizeBetweenArrays(datExpr,"quantile")

for(r in regions) {
  datExpr.regional[,r] = apply(datExpr.qn[,which(grepl(r, datMeta$Region_broad) & datMeta$Period!="fetal")],1,mean)
  
}

datExpr.regional.pSI = specificity.index(datExpr.regional)
pSI.count(datExpr.regional.pSI)
kang$datExpr.regional.pSI = datExpr.regional.pSI
save(kang,file="./working_data/Kang_CollapsedRows_16874genes_1340samples.Rdata")


datMeta$PMI[is.na(datMeta$PMI)] = mean(datMeta$PMI,na.rm=T)
datMeta$pH[is.na(datMeta$pH)] = mean(datMeta$pH,na.rm=T)
datMeta$RIN[is.na(datMeta$RIN)] = mean(datMeta$RIN,na.rm=T)

mod = model.matrix(~as.factor(Stage) + RIN + Sex + PMI + pH, data=datMeta)
beta = (solve(t(mod)%*%mod)%*%t(mod)) %*% t(datExpr)

kang$datExpr.regress = datExpr - t(as.matrix(mod[,16:19]) %*% as.matrix(beta[16:19,]))


dE = read.csv("./raw_data/BrainSpan_RNAseq_genes_matrix_csv/expression_matrix.csv")
dM =  read.csv("./raw_data/BrainSpan_RNAseq_genes_matrix_csv/columns_metadata.csv",row.names=1)
dP = read.csv("./raw_data/BrainSpan_RNAseq_genes_matrix_csv/rows_metadata.csv",row.names=1)
rownames(dE) = dP$ensembl_gene_id
