#PREPARATION FOR DIFFERENTIAL EXPRESSION ANALYSIS
library(tximport)
samples <- c("GBMEC5785", "GBMEC5441", "GBMEC5433", "GBMEC5377", "EC-1", "EC-2", "EC-3")
files <- file.path("../results/salmon", samples, "quant.sf")
names(files) <- samples
txi <- tximport(files, type = "salmon", txOut = TRUE)

condition <- factor(c("Tumor", "Tumor", "Tumor", "Tumor", "Control", "Control", "Control"))
coldata <- data.frame(row.names = samples, condition)

#DIFFERENCIAL EXPRESSION ANALYSIS
library(DESeq2)
dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "Tumor", "Control"))

write.csv(as.data.frame(res), file = "../results/deseq2_results.csv")
saveRDS(res, file = "../results/DESeq2/res.rds")

#HEATMAP
vsd <- vst(dds, blind = FALSE) #variance stabilizing transformation to make the data suitable for the heatmap or clustering
res_ordered <- res[order(res$padj), ] #ordered by p-value adjusted
res_ordered <- res_ordered[!is.na(res_ordered$padj), ] #omit any NA to avoid problems
top20_genes <- rownames(res_ordered)[1:20] #generate the top 20 genes to visualize
mat <- assay(vsd)[top20_genes, ] 

heatmap(mat, 
        scale = "row", 
        col = colorRampPalette(c("blue", "white", "red"))(100), 
        margins = c(10, 10), 
        main = "Top 20 genes más significativos")

#VOLCANO PLOT
library(volcanoPlot)
library(ggplot2)
library(ggrepel)

res_df <- as.data.frame(res)
res_df <- na.omit(res_df)
saveRDS(res_df, file = "../results/DESeq2/res_df.rds")

res_df$comp_grp <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) >= 1,
                          ifelse(res_df$log2FoldChange > 0, "Upregulated", "Downregulated"),
                          "Not Significant") #Adding labels 

res_df$gene_name <- sapply(strsplit(rownames(res_df), "\\|"), `[`, 6) #specifying gene names so it doesn't take transcripts names

top_genes <- res_df[res_df$comp_grp != "Not Significant", ] #adjusting threshold
top_genes <- top_genes[order(top_genes$padj), ][1:10, ] #generating top 10 genes to visualize
saveRDS(top_genes, file = "../results/DESeq2/top_genes.rds")

top_genes$gene_name <- sapply(strsplit(rownames(top_genes), "\\|"), `[`, 6) #specifying gene names so it doesn't take transcripts names

res_df$threshold <- as.factor(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1) #setting it as factor


ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = comp_grp)) +
  geom_point(alpha = 0.4, size = 1.75) +
  geom_text_repel(data = top_genes, aes(label = gene_name), size = 3) +
  theme_minimal() +
  scale_color_manual(values = c("Upregulated" = "red", 
                                "Downregulated" = "blue", 
                                "Not Significant" = "gray")) +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
  theme(legend.position = "none")
ggsave("../results/DESeq2/volcano_plot.png")