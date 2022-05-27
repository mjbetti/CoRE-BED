##Combined
##All regulatory elements
#Construct a data table containing the number of tissues used and regulatory elements for CoRE-BED-test, CoRE-BED-train, ABC, EpiMap ChromHMM, and GTEx v8 eQTLs
num_tissues <- c(827, 827, 28, 28, 131, 827, 49)
num_regs <- c(571226, 698982, 413840, 529303, 267370, 510904, 1851432)
enhancer_names <- c("CoRE-BED-training enhancers", "CoRE-BED-training all regulatory elements", "CoRE-BED-test enhancers", "CoRE-BED-test all regulatory elements", "ABC enhancers", "EpiMap ChromHMM enhancers", "GTEx v8 fine-mapped eQTLs")

enhancer_df <- data.frame(num_tissues, num_regs, enhancer_names)

pdf("combined_all_regs_tissues_vs_elements_scatter_train_test.pdf")
colors <- c("blue",
            "red",
            "gold",
            "green",
            "pink",
            "gray",
            "orange")
reordered_groups <- factor(enhancer_names, levels = c("CoRE-BED-training enhancers",
"CoRE-BED-training all regulatory elements",
"CoRE-BED-test enhancers",
                                             "CoRE-BED-test all regulatory elements",
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