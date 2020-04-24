data <- as.matrix(read.table(inFile$datapath, header = TRUE, sep="\t", row.names=1, check.names=F))
a <- c()
for(f in 1:nrow(data))
{
 data_rescale <- rescale(data[f,], to = c(0, 1), from = range(data[f,], na.rm = TRUE, finite = TRUE))
 a <- rbind(a,data_rescale)
}
UData <- cbind(rownames(data),a)
UData <- as.matrix(UData, rownames=FALSE)