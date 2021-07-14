data <- as.matrix(DataUse())
data <- data[apply(data[,-1], 1, function(x) !all(x==0)),]
data.t <- t(data)
data.t <- data.t[ , which(apply(data.t, 2, var) != 0)]
pca <- prcomp(data.t, center=T, scale. = T)
pc1 <- round(pca$sdev[1]^2/sum(pca$sdev^2)*100,2)
pc2 <- round(pca$sdev[2]^2/sum(pca$sdev^2)*100,2)
PC1_use <- paste0("PC1", "(", pc1, "%)")
PC2_use <- paste0("PC2", "(", pc2, "%)")
Samples_temp <- rownames(data.t)
Samples <- factor(Samples_temp)
scores <- data.frame(Samples_temp, pca$x[,1:3])
MIN_X <- min(scores$PC1)
Max_X <- max(scores$PC1)
header <- "Principal Component Analysis"
plot <- qplot(x=PC1, y=PC2, data=scores, colour=Samples, xlim=c(MIN_X-75,Max_X+75)) + xlab(PC1_use) + ylab(PC2_use) + geom_point(shape=1) + geom_text(aes(label=Samples_temp), hjust=0, vjust=0) + scale_size_area() + theme(axis.text = element_text(size = 14),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"),legend.key = element_rect(fill = "white"),legend.background = element_rect(fill = "white"),panel.grid.major = element_line(),panel.grid.minor = element_blank(),panel.background = element_rect(fill = "white")) + ggtitle(header)
