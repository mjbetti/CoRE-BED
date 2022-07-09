data_a = read.table("figure.data.txt",h=T,sep="\t")
data_a = data_a[order(data_a[,3]),]

library(RColorBrewer)
colpal = brewer.pal(2,"Set1")
mycol  = colpal[c(1)] 

pdf("figure.pdf", height=2,width=7)
layout(matrix(c(1,1,1,2,2), nrow = 3, ncol = 1, byrow = TRUE))
par(mar=c(5,15,4,2))
barplot(data_a[,2],horiz=TRUE,names=data_a$core_bed,border=0,col="gray90",xlab="Number of SNPs",las=2,xaxt="n",xlim=c(0,25),cex.names=0.9)
barplot(data_a[,3],horiz=TRUE,names=""        ,border=0,col=mycol,add=T)
legend("bottomright",c("Within a CoRE-BED element","Not within a CoRE-BED element"), fill=c(colpal[1],"gray90"), border=0, bty="n", cex=0.9)
#mtext(side=3,"a",adj=0,cex=1.25,font=2,las=1)
dev.off()