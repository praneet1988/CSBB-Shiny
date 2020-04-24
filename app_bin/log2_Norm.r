data <- as.matrix(read.table(inFile$datapath, header = TRUE, sep="\t", row.names=1, check.names=F))
UData <- log2(data+1)
UData <- data.frame(UData)