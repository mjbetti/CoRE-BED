########## Figure 2 ##########

data   = read.table("figure2.data.txt",h=T,sep="\t")
legend = read.table("figure2.legend.txt",h=T,sep="\t",fill=T)

library(RColorBrewer)
colpal = brewer.pal(9,"Set1")
mycol  = c(subset(legend$COL,legend$TOPLOT==1),1)
mypch  = c(subset(legend$PCH,legend$TOPLOT==1),19)

pdf("figure2.pdf",height=4,width=7)
layout(matrix(c(1,1,1,2,2), nrow = 1, ncol = 5, byrow = TRUE))
par(mar=c(5,6,4,0))
plot(data$S2G_recall    , data$S2G_precision,pch=mypch,xlim=c(0,0.4),ylim=c(0,1),xlab="Recall",ylab="Precision",cex=1.25,col=colpal[mycol],main="",bty="l",cex.axis=1.25,cex.lab=1.25); 
text(data$S2G_recall[14], data$S2G_precision[14], "  cS2G"    , cex=1.25, col=colpal[mycol[14]], adj=c(0,0.5))
text(data$S2G_recall[ 1], data$S2G_precision[ 1], "  Exon"    , cex=1.00, col=colpal[mycol[ 1]], adj=c(0,0.5))
text(data$S2G_recall[ 2], data$S2G_precision[ 2], "  Promoter", cex=1.00, col=colpal[mycol[ 2]], adj=c(0,1.0))
text(data$S2G_recall[ 6], data$S2G_precision[ 6], "  GTEx"    , cex=1.00, col=colpal[mycol[ 6]], adj=c(0,0.5))
text(data$S2G_recall[ 7], data$S2G_precision[ 7], "  eQTLGen" , cex=1.00, col=colpal[mycol[ 7]], adj=c(0,0.0))
text(data$S2G_recall[ 9], data$S2G_precision[ 9], "  EpiMap"  , cex=1.00, col=colpal[mycol[ 9]], adj=c(0,0.5))
text(data$S2G_recall[10], data$S2G_precision[10], "  ABC"     , cex=1.00, col=colpal[mycol[10]], adj=c(0,0.5))
text(data$S2G_recall[13], data$S2G_precision[13], "  Cicero"  , cex=1.00, col=colpal[mycol[13]], adj=c(0,1.0))
text(data$S2G_recall[ 5], data$S2G_precision[ 5], "\n(Closest TSS)", cex=1.00, col=colpal[mycol[ 5]], adj=c(0.5,0.6))
par(mar=c(1,1,4,0))
plot(0,0,col=0,xlab="",ylab="",xaxt="n",yaxt="n",bty="n")
legend("topleft",legend=as.character(legend$LEGEND),pch=legend$PCH,col=colpal[legend$COL],text.font=legend$BOLD.LEGEND,bty="n",cex=0.8)
dev.off()


########## Figure 3 ##########

data_a = read.table("figure3a.data.txt",h=T,sep="\t")
data_a = data_a[c(1,2,5,6,7,9,10,13,14),]
data_a = data_a[order(data_a[,3]),]
data_b = read.table("figure3b.data.txt",h=F,sep="\t")

data_b_names = data_b[-which(duplicated(data_b[,1])),1]
mywidth      = table(data_b[,1])[data_b_names]
mywidth      = mywidth / max(mywidth) 

library(RColorBrewer)
colpal = brewer.pal(9,"Set1")
mycol  = colpal[c(4,5,2,8,5,8,2,9,1)] 

pdf("figure3.pdf",width=8,height=4)
layout(matrix(c(1,1,1,2,2), nrow = 1, ncol = 5, byrow = TRUE))
par(mar=c(5,18,4,1))
barplot(data_a[,2],horiz=TRUE,names=data_a$S2G,border=0,col="gray90",xlab="Number of SNP-gene-disease triplets",las=2,xaxt="n",xlim=c(0,10000),cex.names=0.9)
barplot(data_a[,3],horiz=TRUE,names=""        ,border=0,col=mycol,add=T)
legend("bottomright",c("Estimated number of correct triplets","Estimated number of incorrect triplets"), fill=c(colpal[1],"gray90"), border=0, bty="n", cex=0.9)
mtext(side=3,"a",adj=0,cex=1.25,font=2,las=1)
#
par(mar=c(5,1,4,5))
boxplot(
subset(data_b[,2],data_b[,1]==data_b_names[1]),
subset(data_b[,2],data_b[,1]==data_b_names[2]),
subset(data_b[,2],data_b[,1]==data_b_names[3]),
subset(data_b[,2],data_b[,1]==data_b_names[4]),
subset(data_b[,2],data_b[,1]==data_b_names[5]),
subset(data_b[,2],data_b[,1]==data_b_names[6]),
subset(data_b[,2],data_b[,1]==data_b_names[7]),
subset(data_b[,2],data_b[,1]==data_b_names[8]),
subset(data_b[,2],data_b[,1]==data_b_names[9]),
border=mycol,lwd=1.0,xlab="Triplet confidence score",width=mywidth,yaxt="n",horizontal=TRUE,ylim=c(0,1),frame=F,outline=F)
mtext(side=3,"b",adj=0,cex=1.25,font=2,las=1)
dev.off()


########## Figure 4 ##########

library(RColorBrewer)

data_a = read.table("figure4.pval.txt",h=T,sep="\t")
data_b = read.table("figure4.pip.txt",h=T,sep="\t")
data_c = read.table("figure4.gene.txt",h=T,sep="\t")

myplot <- function(TRAIT,CHR,POS,GENE,dist,OUTNAME){
    mydata_a = data_a[which(data_a$TRAIT==TRAIT & data_a$CHR==CHR),]
    mydata_b = data_b[which(data_b$TRAIT==TRAIT & data_b$CHR==CHR & data_b$POS>(POS-dist) & (data_b$POS<POS+dist)),]

    mycol = rep(1,nrow(mydata_a))
    rg    = which(mydata_a$POS %in% unique(mydata_b$POS))
    mycol[rg] = 10
    mycol = brewer.pal(12,"Paired")[mycol]
    
    pdf(OUTNAME,width=7,height=5)
    par(mfrow=c(2,1)) 
    ##############################
    par(mar=c(0,5,4,1))
    plot(mydata_a$POS/1000000,-log10(mydata_a$P),col=mycol,pch=16,xlab="",ylab=expression(-log[10]*"(p-value)"),xlim=(POS+dist*c(-1,1))/1000000,bty="l",xaxt="n",cex.axis=1.5,cex.lab=1.5,ylim=c(0,1.05*max(-log10(mydata_a$P),na.rm=T)))
    if (length(rg)>1){eps=750} else {eps=0}
    for (i in 1:length(rg)){
        lines(rep((eps*(i*2-3)+mydata_a$POS[rg[i]])/1000000,2),c(-9999,-log10(mydata_a$P[rg[i]])),col=brewer.pal(12,"Paired")[4],lwd=2)
    }
    for (i in 1:length(rg)){
        points(mydata_a$POS[rg[i]]/1000000,-log10(mydata_a$P[rg[i]]),col=mycol[rg[i]],pch=15,cex=1.25)
        text  (mydata_a$POS[rg[i]]/1000000,-log10(mydata_a$P[rg[i]]),paste0(" ",unique(mydata_b$SNP)[i]," "),col=mycol[rg[i]],cex=1.25,adj=c(2-i,0.5))
    }
    axis(3,at=(POS-dist)/1000000,TRAIT,hadj=0,tick=F,font=2,cex.axis=1.5)
    ##############################
    par(mar=c(5,5,0,1))
    genes = data_c[which(data_c[,1]==CHR & data_c$END>=(POS-dist) & data_c$START<=(POS+dist)),-1]
    genes[1,1] = max(genes[1,1],POS-dist)
    genes[nrow(genes),2] = min(genes[nrow(genes),2],POS+dist)
    mycol2 = rep("grey",nrow(genes)); mycol2[which(genes[,3]==GENE)]=brewer.pal(12,"Paired")[4]
    plot(0,0,col=0,xlab=paste0("Position on chromosome ",CHR," (Mb)"),ylab="",xlim=(POS+dist*c(-1,1))/1000000,ylim=c(0.25,2),bty="l",yaxt="n",cex.axis=1.5,cex.lab=1.5)   #,xaxt="n"
    for (i in 1:length(rg)){
        lines(rep((eps*(i*2-3)+mydata_a$POS[rg[i]])/1000000,2),c(-9999,-log10(mydata_a$P[rg[i]])),col=brewer.pal(12,"Paired")[4],lwd=2)
    }
    myy=rep(c(1.5,1,0.5),nrow(genes))
    for (i in 1:nrow(genes)){
        polygon(genes[i,c(1,2,2,1)]/1000000,myy[i]+0.1*c(-1,-1,1,1),col=mycol2[i],border=0)
        text(sum(genes[i,c(1,2)])/2000000,myy[i]-0.1,genes[i,3],adj=c(0.5,1),col=1,cex=1.1)
    }
    myannot = subset(mydata_b$annot,mydata_b$gene==GENE)
    for (i in 1:length(rg)){
        text((eps*(i*2-3)+mydata_a$POS[rg[i]])/1000000,2,paste0(" ",myannot[i]," "),col=mycol[rg[i]],adj=c(2-i,1),cex=1.25)
    }
    dev.off()
}

myplot("Type 2 diabetes",  11,  2858636, "CDKN1C", 100000,"figure4a.pdf")
myplot(         "Asthma",  3, 187907615,   "BCL6", 600000,"figure4b.pdf")
myplot(         "Eczema",  2, 242698640,  "PDCD1", 150000,"figure4c.pdf")
myplot(            "HDL", 13, 113927208,  "LAMP1", 100000,"figure4d.pdf")


########## Figure 5 ##########

library(RColorBrewer)
colpal = brewer.pal(12,"Set3")[-c(2,9)]

data_a = read.table("figure5a.data.txt",h=T,fill=T)
y0     = data_a[,1]
y0.se  = data_a[,2]; y0.se[8] = 0
y1     = data_a[,3]
y1.se  = data_a[,4]; y1.se[8] = 0
y2     = data_a[,5]
nbgenes      = data_a[,6]# c(100,200,500,1000,2000,5000,10000,19995)
nbgenes_axis = data_a[,7]# c(100,"",500,"1K","2K","5K","10K","20K")

data_bc = read.table("figure5bc.data.txt",h=T)
traits  = data_bc$trait
Me      = data_bc$Me
Ge      = data_bc$Ge
Ge_c    = data_bc$Ge_c
Ge_lf   = data_bc$Ge_lf
logscaletick = rep(1:10,6)*c(rep(10**0,10),rep(10**1,10),rep(10**2,10),rep(10**3,10),rep(10**4,10),rep(10**5,10))
indep_traits = c(45,12,1,11,13,18,19,21,9,41,28,23,48,42,44,49)
indep_names  = c("Hair Color", "MPV", "Cholesterol", "MCH", "MC", "BMD", "Balding", "Height", "HLSRC", "FEV1/FVC", "Eczema", "DBP", "Age Menarche", "FVC", "Chronotype", "#Children")
tokeep  = c(1,20,21,27,28,43,39,45); 
mynames = c("Cholesterol","BMI","Height","AID","Eczema","Neuroticism","T2D","Hair color")
adj0 = c(0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0)
adj1 = c(0.5,0.5,0.5,0.5,1.0,0.0,0.5,0.5)
adj2 = c(0.5,0.0,0.0,0.5,0.5,0.5,0.5,0.5)

pdf("figure5.pdf",height=8,width=8)
layout(matrix(c(1,1,2,3), nrow = 2, ncol = 2, byrow = TRUE))
#
par(mar=c(5,5,3,2))
plot (c(0,nbgenes),c(0,y1),type="l",xlab="\nNo. of top-ranked genes",ylab=expression("Proportion of "*italic("h")["gene"]**2),xlim=c(1,20000),ylim=c(0,1.15),xaxt="n",col="red",cex.lab=1.25,bty="l")
polygon(c(c(0,nbgenes),rev(c(0,nbgenes))),c(0,y1-1.96*y1.se,rev(y1+1.96*y1.se),0),col="lightgrey",border=0)
polygon(c(c(0,nbgenes),rev(c(0,nbgenes))),c(0,y0-1.96*y0.se,rev(y0+1.96*y0.se),0),col="grey",border=0)
#abline(h=1,col=2,lty=2)
lines(c(0,nbgenes),c(0,y0),type="l",col="darkblue",lwd=2)
lines(c(0,nbgenes),c(0,y1),type="l",col="red",lwd=2)
lines(c(0,nbgenes),c(0,y2),type="l",col="red",lty=2,lwd=2)
axis(1,at=nbgenes,nbgenes_axis,las=2)
legend("bottomright",c("cS2G - Validation (S-LDSC; 122K)","cS2G - Training (PolyFun; 337K)","Closest TSS - Validation (S-LDSC; 122K)"),col=c("red","red","darkblue"),lty=c(1,2,1),lwd=2,bty="n")
mtext(side=3,"a\n",adj=0,cex=1.5,font=2)
#
par(mar=c(6,5,3,2))
plot(Me,Ge,pch=16,col="lightgrey",log="xy",xlab="",ylab="",xlim=c(min(Me),100000),ylim=c(min(c(Ge_c,Ge_lf)),5000),xaxt="n",yaxt="n",bty="l")
mtext(side=1, text=bquote(atop("Effective number of causal SNPs ("*italic("M")["e"]*")")), line = 4)
mtext(side=2, text=bquote(atop("Effective number of causal genes ("*italic("G")["e"]*")")), line = 1.5)
axis(1,at=logscaletick,rep("",60),tck=-0.02,lwd.ticks=0.5); axis(2,at=logscaletick,rep("",60),tck=-0.02,lwd.ticks=0.5)
axis(1,at=c(10**2,10**3,10**4,10**5,10**6,10**7),c(expression(10**2),expression(10**3),expression(10**4),expression(10**5),expression(10**6),expression(10**7)))
axis(2,at=c(1,10,100,1000),c(1,10,100,1000))
points(median(Me[indep_traits]),median(Ge[indep_traits]),pch=15,col="red",cex=1.5)
text(min(Me),5000,paste0("Correlation = ",round(cor(log(Me[indep_traits]),log(Ge[indep_traits])),2)),adj=c(0,1))
for (i in 1:length(tokeep)){
    points(Me[traits[tokeep[i]]],Ge[traits[tokeep[i]]],pch=16,col=colpal[i])
    text  (Me[traits[tokeep[i]]],Ge[traits[tokeep[i]]],paste0("  ",mynames[i],"  "),adj=c(adj0[i],adj1[i]),col=colpal[i],font=2)
}
mtext(side=3,"b\n",adj=0,cex=1.5,font=2)
#
plot(Ge_c,Ge_lf,pch=16,xlim=c(min(c(Ge_c,Ge_lf)),5000),ylim=c(min(c(Ge_c,Ge_lf)),5000),log="xy",col="lightgrey",xaxt="n",yaxt="n",xlab="",ylab="",bty="l")
mtext(side=1, text=bquote(atop("Effective number of causal genes","for common SNPs ("*italic("G")["e,common"]*")")), line = 4)
mtext(side=2, text=bquote(atop("Effective number of causal genes","for low-freq. SNPs ("*italic("G")["e,low-freq"]*")")), line = 2)
axis(1,at=logscaletick,rep("",60),tck=-0.02,lwd.ticks=0.5); axis(2,at=logscaletick,rep("",60),tck=-0.02,lwd.ticks=0.5)
axis(1,at=c(1,10,100,1000),c(1,10,100,1000))
axis(2,at=c(1,10,100,1000),c(1,10,100,1000))
abline(0,1,col="grey")
points(median(Ge_c[indep_traits]),median(Ge_lf[indep_traits]),pch=15,col=2,cex=1.5)
text(min(c(Ge_c,Ge_lf)),5000,paste0("Correlation = ",round(cor(log(Ge_c[indep_traits]),log(Ge_lf[indep_traits])),2)),adj=c(0,1))
for (i in 1:length(tokeep)){
    points(Ge_c[traits[tokeep[i]]],Ge_lf[traits[tokeep[i]]],pch=16,col=colpal[i])
    text  (Ge_c[traits[tokeep[i]]],Ge_lf[traits[tokeep[i]]],paste0("  ",mynames[i],"  "),adj=c(adj0[i],adj2[i]),col=colpal[i],font=2)
}
mtext(side=3,"c\n",adj=0,cex=1.5,font=2)
dev.off()



