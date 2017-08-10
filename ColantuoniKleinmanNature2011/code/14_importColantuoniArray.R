#14_importColantuoniArray.R

rm(list=ls())


datExpr= read.delim("./raw_data/Colantuoni_GSE30272/GSE30272_ExprsMtxCleanedN269_31SVN.txt")
rownames(datExpr) = datExpr$ID; datExpr=datExpr[,-1]



datMeta = read.delim("./raw_data/Colantuoni_GSE30272/GSE30272_ModelMatrix.txt")
datMeta2 = read.delim("./raw_data/Colantuoni_GSE30272/GSE30272_SVN31Matrix.txt")
datMeta3 = data.frame(t(read.delim("./raw_data/Colantuoni_GSE30272/GSE30272_series_matrix.txt", nrow=56,skip=55)))
datMeta3 = datMeta3[-1,]
datMeta3$Batch = factor(gsub("array batch: ", "", datMeta3$X9))
datMeta3$Age = as.numeric(gsub("age: ", "", datMeta3$X10))
datMeta3$Sex = as.factor(gsub("Sex: ", "", datMeta3$X11))
datMeta3$PMI = as.numeric(gsub("postmortem interval [(]pmi[)]:", "", datMeta3$X13))
datMeta3$RIN = as.numeric(gsub("rna integrity number [(]rin[)]: ", "", datMeta3$X15))
datMeta = datMeta3[match(colnames(datExpr), rownames(datMeta3)),]



datProbes = read.delim("./raw_data/Colantuoni_GSE30272/GPL4611_noParents.an.txt/GPL4611_noParents.an.txt",skip=7)
datProbes= datProbes[match(rownames(datExpr), datProbes$ProbeName),]

bm = useMart("ENSEMBL_MART_ENSEMBL", "hsapiens_gene_ensembl", host="grch37.ensembl.org")
A = listAttributes(bm); F = listFilters(bm)
datProbes2 = getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
                   filters = "external_gene_name", values = datProbes$GeneSymbols, mart=bm)
datProbes = cbind(datProbes, datProbes2[match(datProbes$GeneSymbols, datProbes2$external_gene_name),])


collapseDat <- collapseRows(datET = datExpr,rowGroup = datProbes$ensembl_gene_id,rowID = rownames(datExpr))
datExpr <- collapseDat$datETcollapsed
datProbes <- datProbes[match(rownames(datExpr), datProbes$ensembl_gene_id),]
all(rownames(datExpr)==datProbes$ensg) ## Samples match
all(colnames(datExpr)==rownames(datMeta))


CP = read.delim("~/Github/CrossDisorder_GWAS_annotation/data/FBD_combinedgenes.txt",head=FALSE)
dE = scale(datExpr)

idx = na.omit(match(CP$V1, rownames(datExpr)))
t.test(apply(datExpr[idx,],2,mean), datMeta$Age>0)
plot(apply(dE[idx,],2,mean)~datMeta$Age)
t.test(apply(dE[idx,],2,mean) ~ (datMeta$Age>0))

a = datMeta$Age
a[a<0] = a[a<0] * 40
plot(apply(dE[idx,a<0],2,mean) ~ a[a<0])
