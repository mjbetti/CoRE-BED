col1 <- c("ABC", "EpiMap")
col2 <- c(25, 25)
col3 <- c(11, 7)

data_a <- data.frame(col1, col2, col3)
names(data_a) <- c("method", "num_variants", "num_overlap")

library(RColorBrewer)
colpal = brewer.pal(2,"Set1")
mycol  = colpal[c(2)] 

pdf("Figure7b.pdf", height=2,width=7)
layout(matrix(c(1,1,1,2,2), nrow = 3, ncol = 1, byrow = TRUE))
par(mar=c(5,15,4,2))
barplot(data_a[,2],horiz=TRUE,names=data_a$method,border=0,col="gray90",xlab="Number of SNPs",las=2,xaxt="n",xlim=c(0,40),cex.names=0.9)
barplot(data_a[,3],horiz=TRUE,names=""        ,border=0,col=mycol,add=T)
legend("bottomright",c("Regulatory overlap","No regulatory overlap"), fill=c(colpal[2],"gray90"), border=0, bty="n", cex=0.9)
#mtext(side=3,"a",adj=0,cex=1.25,font=2,las=1)
dev.off()