library(ggplot2)

class1_df <- data.frame(matrix(nrow = 0, ncol = 3))
names(class1_df) <- c("gwas", "h2", "se")
class2_df <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(class2_df) <- c("gwas", "h2", "se")
class3_df <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(class3_df) <- c("gwas", "h2", "se")
class4_df <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(class4_df) <- c("gwas", "h2", "se")
class5_df <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(class5_df) <- c("gwas", "h2", "se")
class6_df <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(class6_df) <- c("gwas", "h2", "se")
class7_df <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(class7_df) <- c("gwas", "h2", "se")
class8_df <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(class8_df) <- c("gwas", "h2", "se")
any_df <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(any_df) <- c("gwas", "h2", "se")

for (result in list.files("/home/bettimj/gamazon_rotation/mod_core-bed/uk_biobank/h2/anno_1000G", pattern = "*.results")) {
	print(result)
	file <- read.table(paste0("/home/bettimj/gamazon_rotation/mod_core-bed/uk_biobank/h2/anno_1000G/", result), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
	df <- as.data.frame(file)
	df$Prop._h2[df$Prop._h2 < 0] <- 0
	class1_vec <- c(result, df[1,3], df[1,4])
	class1_df <- rbind(class1_df, class1_vec)
	class2_vec <- c(result, df[2,3], df[2,4])
	class2_df <- rbind(class2_df, class2_vec)
	class3_vec <- c(result, df[3,3], df[3,4])
	class3_df <- rbind(class3_df, class3_vec)
	class4_vec <- c(result, df[4,3], df[4,4])
	class4_df <- rbind(class4_df, class4_vec)
	class5_vec <- c(result, df[5,3], df[5,4])
	class5_df <- rbind(class5_df, class5_vec)
	class6_vec <- c(result, df[6,3], df[6,4])
	class6_df <- rbind(class6_df, class6_vec)
	class7_vec <- c(result, df[7,3], df[7,4])
	class7_df <- rbind(class7_df, class7_vec)
	class8_vec <- c(result, df[8,3], df[6,4])
	class8_df <- rbind(class8_df, class8_vec)
	any_vec <- c(result, df[9,3], df[9,4])
	any_df <- rbind(any_df, any_vec)
}

#Class 1
names(class1_df) <- c("gwas", "h2", "se")
class1_df$h2 <- as.numeric(class1_df$h2)
class1_df$se <- as.numeric(class1_df$se)
class1_df$gwas <- strsplit(class1_df$gwas, "_h")
class1_df$gwas <- unlist(lapply(class1_df$gwas, `[`, 1))
dodge <- position_dodge(width = 0.9)
limits <- aes(ymax = class1_df$h2 + class1_df$se,
              ymin = class1_df$h2 - class1_df$se)

png(file = "class1_69_traits_h2.png",
height=4000, width=5000, res = 300)
#margins below in inches: mai=c(xbot,xlef,xtop,xrig), mgp=axis title,axis labels,axis line
par(pty="m",mai=c(1,1.5,.5,.5),mgp=c(3,1,0))
ggplot(data = class1_df, aes(x = gwas, y = h2, fill = "darkorchid")) + 
geom_bar(stat = "identity", position = dodge) +
  geom_errorbar(limits, position = dodge, width = 0.25) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.ticks.x=element_blank(),
        axis.title.x=element_blank(), legend.position = "none")
dev.off()

#Class 2
names(class2_df) <- c("gwas", "h2", "se")
class2_df$h2 <- as.numeric(class2_df$h2)
class2_df$se <- as.numeric(class2_df$se)
class2_df$gwas <- strsplit(class2_df$gwas, "_h")
class2_df$gwas <- unlist(lapply(class2_df$gwas, `[`, 1))
dodge <- position_dodge(width = 0.9)
limits <- aes(ymax = class2_df$h2 + class2_df$se,
              ymin = class2_df$h2 - class2_df$se)

png(file = "class2_69_traits_h2.png",
height=4000, width=5000, res = 300)
#margins below in inches: mai=c(xbot,xlef,xtop,xrig), mgp=axis title,axis labels,axis line
par(pty="m",mai=c(1,1.5,.5,.5),mgp=c(3,1,0))
ggplot(data = class2_df, aes(x = gwas, y = h2, fill = "darkorchid")) + 
geom_bar(stat = "identity", position = dodge) +
  geom_errorbar(limits, position = dodge, width = 0.25) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.ticks.x=element_blank(),
        axis.title.x=element_blank(), legend.position = "none")
dev.off()

#Class 3
names(class3_df) <- c("gwas", "h2", "se")
class3_df$h2 <- as.numeric(class3_df$h2)
class3_df$se <- as.numeric(class3_df$se)
class3_df$gwas <- strsplit(class3_df$gwas, "_h")
class3_df$gwas <- unlist(lapply(class3_df$gwas, `[`, 1))
dodge <- position_dodge(width = 0.9)
limits <- aes(ymax = class3_df$h2 + class3_df$se,
              ymin = class3_df$h2 - class3_df$se)

png(file = "class3_69_traits_h2.png",
height=4000, width=5000, res = 300)
#margins below in inches: mai=c(xbot,xlef,xtop,xrig), mgp=axis title,axis labels,axis line
par(pty="m",mai=c(1,1.5,.5,.5),mgp=c(3,1,0))
ggplot(data = class3_df, aes(x = gwas, y = h2, fill = "darkorchid")) + 
geom_bar(stat = "identity", position = dodge) +
  geom_errorbar(limits, position = dodge, width = 0.25) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.ticks.x=element_blank(),
        axis.title.x=element_blank(), legend.position = "none")
dev.off()

#Class 4
names(class4_df) <- c("gwas", "h2", "se")
class4_df$h2 <- as.numeric(class4_df$h2)
class4_df$se <- as.numeric(class4_df$se)
class4_df$gwas <- strsplit(class4_df$gwas, "_h")
class4_df$gwas <- unlist(lapply(class4_df$gwas, `[`, 1))
dodge <- position_dodge(width = 0.9)
limits <- aes(ymax = class4_df$h2 + class4_df$se,
              ymin = class4_df$h2 - class4_df$se)

png(file = "class4_69_traits_h2.png",
height=4000, width=5000, res = 300)
#margins below in inches: mai=c(xbot,xlef,xtop,xrig), mgp=axis title,axis labels,axis line
par(pty="m",mai=c(1,1.5,.5,.5),mgp=c(3,1,0))
ggplot(data = class4_df, aes(x = gwas, y = h2, fill = "darkorchid")) + 
geom_bar(stat = "identity", position = dodge) +
  geom_errorbar(limits, position = dodge, width = 0.25) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.ticks.x=element_blank(),
        axis.title.x=element_blank(), legend.position = "none")
dev.off()

#Class 5
names(class5_df) <- c("gwas", "h2", "se")
class5_df$h2 <- as.numeric(class5_df$h2)
class5_df$se <- as.numeric(class5_df$se)
class5_df$gwas <- strsplit(class5_df$gwas, "_h")
class5_df$gwas <- unlist(lapply(class5_df$gwas, `[`, 1))
dodge <- position_dodge(width = 0.9)
limits <- aes(ymax = class5_df$h2 + class5_df$se,
              ymin = class5_df$h2 - class5_df$se)

png(file = "class5_69_traits_h2.png",
height=4000, width=5000, res = 300)
#margins below in inches: mai=c(xbot,xlef,xtop,xrig), mgp=axis title,axis labels,axis line
par(pty="m",mai=c(1,1.5,.5,.5),mgp=c(3,1,0))
ggplot(data = class5_df, aes(x = gwas, y = h2, fill = "darkorchid")) + 
geom_bar(stat = "identity", position = dodge) +
  geom_errorbar(limits, position = dodge, width = 0.25) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.ticks.x=element_blank(),
        axis.title.x=element_blank(), legend.position = "none")
dev.off()

#Class 6
names(class6_df) <- c("gwas", "h2", "se")
class6_df$h2 <- as.numeric(class6_df$h2)
class6_df$se <- as.numeric(class6_df$se)
class6_df$gwas <- strsplit(class6_df$gwas, "_h")
class6_df$gwas <- unlist(lapply(class6_df$gwas, `[`, 1))
dodge <- position_dodge(width = 0.9)
limits <- aes(ymax = class6_df$h2 + class6_df$se,
              ymin = class6_df$h2 - class6_df$se)

png(file = "class6_69_traits_h2.png",
height=4000, width=5000, res = 300)
#margins below in inches: mai=c(xbot,xlef,xtop,xrig), mgp=axis title,axis labels,axis line
par(pty="m",mai=c(1,1.5,.5,.5),mgp=c(3,1,0))
ggplot(data = class6_df, aes(x = gwas, y = h2, fill = "darkorchid")) + 
geom_bar(stat = "identity", position = dodge) +
  geom_errorbar(limits, position = dodge, width = 0.25) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.ticks.x=element_blank(),
        axis.title.x=element_blank(), legend.position = "none")
dev.off()

#Class 7
names(class7_df) <- c("gwas", "h2", "se")
class7_df$h2 <- as.numeric(class7_df$h2)
class7_df$se <- as.numeric(class7_df$se)
class7_df$gwas <- strsplit(class7_df$gwas, "_h")
class7_df$gwas <- unlist(lapply(class7_df$gwas, `[`, 1))
dodge <- position_dodge(width = 0.9)
limits <- aes(ymax = class7_df$h2 + class7_df$se,
              ymin = class7_df$h2 - class7_df$se)

png(file = "class7_69_traits_h2.png",
height=4000, width=5000, res = 300)
#margins below in inches: mai=c(xbot,xlef,xtop,xrig), mgp=axis title,axis labels,axis line
par(pty="m",mai=c(1,1.5,.5,.5),mgp=c(3,1,0))
ggplot(data = class7_df, aes(x = gwas, y = h2, fill = "darkorchid")) + 
geom_bar(stat = "identity", position = dodge) +
  geom_errorbar(limits, position = dodge, width = 0.25) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.ticks.x=element_blank(),
        axis.title.x=element_blank(), legend.position = "none")
dev.off()

#Class 8
names(class8_df) <- c("gwas", "h2", "se")
class8_df$h2 <- as.numeric(class8_df$h2)
class8_df$se <- as.numeric(class8_df$se)
class8_df$gwas <- strsplit(class8_df$gwas, "_h")
class8_df$gwas <- unlist(lapply(class8_df$gwas, `[`, 1))
dodge <- position_dodge(width = 0.9)
limits <- aes(ymax = class8_df$h2 + class8_df$se,
              ymin = class8_df$h2 - class8_df$se)

png(file = "class8_69_traits_h2.png",
height=4000, width=5000, res = 300)
#margins below in inches: mai=c(xbot,xlef,xtop,xrig), mgp=axis title,axis labels,axis line
par(pty="m",mai=c(1,1.5,.5,.5),mgp=c(3,1,0))
ggplot(data = class8_df, aes(x = gwas, y = h2, fill = "darkorchid")) + 
geom_bar(stat = "identity", position = dodge) +
  geom_errorbar(limits, position = dodge, width = 0.25) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.ticks.x=element_blank(),
        axis.title.x=element_blank(), legend.position = "none")
dev.off()

#Any
names(any_df) <- c("gwas", "h2", "se")
any_df$h2 <- as.numeric(any_df$h2)
any_df$se <- as.numeric(any_df$se)
any_df$gwas <- strsplit(any_df$gwas, "_h")
any_df$gwas <- unlist(lapply(any_df$gwas, `[`, 1))
dodge <- position_dodge(width = 0.9)
limits <- aes(ymax = any_df$h2 + any_df$se,
              ymin = any_df$h2 - any_df$se)

png(file = "any_69_traits_h2.png",
height=4000, width=5000, res = 300)
#margins below in inches: mai=c(xbot,xlef,xtop,xrig), mgp=axis title,axis labels,axis line
par(pty="m",mai=c(1,1.5,.5,.5),mgp=c(3,1,0))
ggplot(data = any_df, aes(x = gwas, y = h2, fill = "darkorchid")) + 
geom_bar(stat = "identity", position = dodge) +
  geom_errorbar(limits, position = dodge, width = 0.25) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.ticks.x=element_blank(),
        axis.title.x=element_blank(), legend.position = "none")
dev.off()

print(mean(class1_df$h2))
print(mean(class2_df$h2))
print(mean(class3_df$h2))
print(mean(class4_df$h2))
print(mean(class5_df$h2))
print(mean(class6_df$h2))
print(mean(class7_df$h2))
print(mean(class8_df$h2))
print(mean(any_df$h2))

#[1] 0.4875247
#[1] 0.02504062
#[1] 0.4794796
#[1] 0.04241531
#[1] 0.3038595
#[1] 0.04055155
#[1] 0.4564202
#[1] 0.2697938
#[1] 1
