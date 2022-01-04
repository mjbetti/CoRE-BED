##Enhancers only
#Construct a data table containing the number of tissues used and regulatory elements for CoRE-BED Impute, CoRE-BED, ABC, EpiMap ChromHMM, and GTEx v8 eQTLs
num_tissues <- c(827, 28, 131, 833, 49)
num_enhancers <- c(571226, 413840, 267370, 510904, 1851432)
enhancer_names <- c("CoRE-BED Impute enhancers", "CoRE-BED enhancers", "ABC enhancers", "EpiMap ChromHMM enhancers", "GTEx v8 fine-mapped eQTLs")

enhancer_df <- data.frame(num_tissues, num_enhancers, enhancer_names)

pdf("enhancers_tissues_vs_elements_scatter.pdf")
colors <- c("blue",
            "red",
            "gold",
            "green",
            "orange")
reordered_groups <- factor(enhancer_names, levels = c("CoRE-BED Impute enhancers",
                                             "CoRE-BED enhancers",
                                             "ABC enhancers",
                                             "EpiMap ChromHMM enhancers",
                                             "GTEx v8 fine-mapped eQTLs"))
plot(num_tissues, num_enhancers,
     pch = 19,
     xlab = "Number of cell/tissue types",
     ylab = "Number of cis-regulatory elements",
     col = colors)

# Legend
legend("topright",
       legend = enhancer_names,
       pch = 19,
       col = colors)
dev.off()

##All regulatory elements
#Construct a data table containing the number of tissues used and regulatory elements for CoRE-BED Impute, CoRE-BED, ABC, EpiMap ChromHMM, and GTEx v8 eQTLs
num_tissues <- c(827, 28, 131, 833, 49)
num_regs <- c(698982, 529303, 267370, 510904, 1851432)
enhancer_names <- c("CoRE-BED Impute regulatory elements", "CoRE-BED regulatory elements", "ABC enhancers", "EpiMap ChromHMM enhancers", "GTEx v8 fine-mapped eQTLs")

enhancer_df <- data.frame(num_tissues, num_regs, enhancer_names)

pdf("all_regs_tissues_vs_elements_scatter.pdf")
colors <- c("blue",
            "red",
            "gold",
            "green",
            "orange")
reordered_groups <- factor(enhancer_names, levels = c("CoRE-BED Impute regulatory elements",
                                             "CoRE-BED regulatory elements",
                                             "ABC enhancers",
                                             "EpiMap ChromHMM enhancers",
                                             "GTEx v8 fine-mapped eQTLs"))
plot(num_tissues, num_regs,
     pch = 19,
     xlab = "Number of cell/tissue types",
     ylab = "Number of cis-regulatory elements",
     col = colors)

# Legend
legend("topright",
       legend = enhancer_names,
       pch = 19,
       col = colors)
dev.off()

##Combined
##All regulatory elements
#Construct a data table containing the number of tissues used and regulatory elements for CoRE-BED Impute, CoRE-BED, ABC, EpiMap ChromHMM, and GTEx v8 eQTLs
num_tissues <- c(827, 827, 28, 28, 131, 833, 49)
num_regs <- c(571226, 698982, 413840, 529303, 267370, 510904, 1851432)
enhancer_names <- c("CoRE-BED Impute enhancers", "CoRE-BED Impute all regulatory elements", "CoRE-BED enhancers", "CoRE-BED all regulatory elements", "ABC enhancers", "EpiMap ChromHMM enhancers", "GTEx v8 fine-mapped eQTLs")

enhancer_df <- data.frame(num_tissues, num_regs, enhancer_names)

pdf("combined_all_regs_tissues_vs_elements_scatter.pdf")
colors <- c("blue",
            "red",
            "gold",
            "green",
            "pink",
            "gray",
            "orange")
reordered_groups <- factor(enhancer_names, levels = c("CoRE-BED Impute enhancers",
"CoRE-BED Impute all regulatory elements",
"CoRE-BED enhancers",
                                             "CoRE-BED all regulatory elements",
                                             "ABC enhancers",
                                             "EpiMap ChromHMM enhancers",
                                             "GTEx v8 fine-mapped eQTLs"))
plot(num_tissues, num_regs,
     pch = 19,
     xlab = "Number of cell/tissue types",
     ylab = "Number of cis-regulatory elements",
     col = colors)

# Legend
legend("topright",
       legend = enhancer_names,
       pch = 19,
       col = colors)
dev.off()