#11_importRNAseqFromWebsite.R

rm(list=ls())
setwd("~/Github/BrainSpan//")


datExpr = read.csv("./raw_data/BrainSpan_RNAseq_genes_matrix_csv/expression_matrix.csv",row.names=1,head=F)
datMeta = read.csv("./raw_data/BrainSpan_RNAseq_genes_matrix_csv/columns_metadata.csv")
datProbes = read.csv("./raw_data/BrainSpan_RNAseq_genes_matrix_csv/rows_metadata.csv")

rownames(datExpr) = datProbes$ensembl_gene_id




num=as.numeric(substr(datMeta$age,0,2))
prenatal = grepl("pcw", datMeta$age); postnatal = !prenatal
yrs = grepl("yrs", datMeta$age); mos = grepl("mos",datMeta$age)

datMeta$Age_years[yrs] = num[yrs]
datMeta$Age_years[prenatal] = -(40 - num[prenatal])/40
datMeta$Age_years[mos] = num[mos]/12

num[mos] = num[mos]/12
datMeta$Age_unit = num
datMeta$Unit = "Postnatal (Years)"; datMeta$Unit[prenatal] = "Prenatal (Post-Conception Week)"


br = datMeta$structure_acronym
datMeta$Region_broad[br=="AMY"] = "Amygdala"
datMeta$Region_broad[br=="A1C" | br == "ITC"| br == "STC" | br == "TCx"] = "Temporal Cortex"
datMeta$Region_broad[br=="CB" | br=="CBC"| br == "URL"] = "Cerebellum"
datMeta$Region_broad[br=="CGE" | br=="LGE" | br =="MGE" | br == "STR" | br=="VF"] = "Striatum"
datMeta$Region_broad[br=="DFC" | br == "MFC" | br == "PC"| br == "VFC" | br=="M1C" | br == "OFC" | br=="FC"] = "Frontal Cortex"
datMeta$Region_broad[br=="DTH" | br == "MD" | br=="DIE"] = "Thalamus"
datMeta$Region_broad[br=="HIP"] = "Hippocampus"
datMeta$Region_broad[br=="IPC"| br == "S1C" | br == "PCx"] = "Parietal Cortex"
datMeta$Region_broad[br=="Ocx"| br == "V1C"] = "Occipital Cortex"
datMeta$Region_broad = as.factor(datMeta$Region_broad)


datExpr = log2(datExpr+1)

#Filter out low expressed genes
to_keep = (apply(datExpr>.5,1,sum) > ncol(datExpr)/2); table(to_keep)
datExpr = datExpr[to_keep,]
datProbes = datProbes[to_keep,]


brainspan = vector(mode="list", length=3)
brainspan$datExpr = datExpr
brainspan$datMeta = datMeta
brainspan$datProbes = datProbes



save(file="./working_data/BrainSpan_RNAseqFromWebsite.Rdata", brainspan)
