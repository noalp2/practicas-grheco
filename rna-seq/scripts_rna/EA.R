#PREPARATION FOR ENRICHMENT ANALYSIS
res_df <- readRDS("../results/DESeq2/res_df.rds")
res_df$Gene.symbol <- sapply(strsplit(rownames(res_df), "\\|"), `[`, 6)

pathfindR_input <- res_df[, c("Gene.symbol", "log2FoldChange", "padj")]
colnames(pathfindR_input) <- c("Gene.symbol", "logFC", "adj.P.Val")
pathfindR_input <- pathfindR_input[pathfindR_input$adj.P.Val < 0.05, ]

#ENRICHMENT ANALYSIS + GRAPHICAL SUMMARY
library(pathfindR)
output_res_df <- run_pathfindR(pathfindR_input, 
                        gene_sets = "KEGG", 
                        p_val_threshold = 0.05, 
                        output_dir = "../results/pathfindR", 
                        convert2alias = TRUE, 
                        n_processes = 1)

write.csv(pathfindR_input, file = "../results/pathfindR/pathfindR_input.csv", row.names = FALSE)
saveRDS(pathfindR_input, file = "../results/pathfindR/pathfindR_input.rds")
write.csv(output_res_df, file = "../results/pathfindR/pathfindR_results.csv", row.names = FALSE)
saveRDS(output_res_df, file = "../results/pathfindR/pathfindR_results.rds")

#DENDROGRAM 
top_terms <- output_res_df[order(output_res_df$lowest_p), ][1:20, ]
top_terms$ID <- top_terms$Term_Description
output_res_df_clustered <- cluster_enriched_terms(top_terms, plot_dend = TRUE)

write.csv(output_res_df_clustered, file = "../results/pathfindR/pathfindR_results_clustered.csv", row.names = FALSE)
saveRDS(output_res_df_clustered, file = "../results/pathfindR/pathfindR_results_clustered.rds")
