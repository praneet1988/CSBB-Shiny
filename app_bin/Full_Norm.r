data <- as.matrix(read.table(inFile$datapath, header = TRUE, sep="\t", row.names=1, check.names=F))
UData <- betweenLaneNormalization(data, which="full", round=FALSE)
UData <- data.frame(UData)