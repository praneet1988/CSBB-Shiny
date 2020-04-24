options(warn=-1)
##############################  RUVSEQ  Normalization ##############################
NumberofControl <- lengthcontrol
NumberofTreatment <- lengthtreatment
Counts <- 10
NumberofSamples <- 1
##############################################


countData <- data_temp
filter <- apply(countData, 1, function(x) length(x[x>Counts])>=NumberofSamples)
filtered <- countData[filter,]
genes <- rownames(filtered)
aa <- length(genes)
Controlrep <- rep("Ctl",NumberofControl)
Treatmentrep <- rep("Trt",NumberofTreatment)
x <- as.factor(c(Controlrep,Treatmentrep))
set <- newSeqExpressionSet(as.matrix(filtered), phenoData = data.frame(x, row.names=colnames(filtered)))


##upper quartlie normalization
set <- betweenLaneNormalization(set, which="upper")

##Differential expression analysis using Empirical RUVg
design <- model.matrix(~x, data=pData(set))
y <- DGEList(counts=counts(set), group=x)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, method="deviance", robust=TRUE, subset=NULL)
y <- estimateGLMTagwiseDisp(y)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
topTags(lrt)
result <- topTags(lrt, n=dim(y)[1]+1, adjust.method="BH", sort.by="logFC")
DEresult <- as.data.frame(result)
return(DEresult)