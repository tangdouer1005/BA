library(ggcorrplot)
library(tidyr)
#Import data
DEG_expr <- read.csv("DEG_expr.csv",row.names = 1)
#CIBERSORT results
cibersort_result <- read.csv("cibersort_result.csv",row.names = 1)
data <- t(DEG_expr)
cibersort_result <- cibersort_result[,1:22]
if(max(data) < 50) {data <- 2^data}
data <- t(DEG_expr)
identical(rownames(data),rownames(cibersort_result))
#Spearman correlation analysis
calculate_correlation <- function(Gene_expr, cibersort_result) {
  cor_matrix <- data.frame("Gene" = character(),"im_cell" = character(),"Cor" = numeric(),"p-value" = numeric(), stringsAsFactors = FALSE)
  
  for (i in 1:ncol(cibersort_result)) {
    result <- cor.test(Gene_expr, cibersort_result[, i], method = "spearman")
    new_row <- data.frame("Gene" = "Gene", "im_cell" = colnames(cibersort_result)[i], "Cor" = result$estimate, "p-value" = result$p.value)
    cor_matrix <- rbind(cor_matrix, new_row)
  }
  return(cor_matrix)
}
#Choose gene
data_to_calculate <- data[,"COL1A1"]
data_to_calculate <- data.frame(data_to_calculate)
results <- data.frame("Gene" = character(),"im_cell" = character(),"Cor" = numeric(),"p-value" = numeric(), stringsAsFactors = FALSE)
for (i in 1:ncol(data_to_calculate)) {
  print(i)
  gene_expr <- data_to_calculate[, i]
  corr_result <- calculate_correlation(gene_expr, cibersort_result)
  results <- rbind(results, corr_result)
}
#The gene names were modified manually
results$Gene <- c(rep("COL1A1", 22))
results <- results[,c(1,2,3)]
colnames(results)
#The long data was converted to a wide format using the spread function
results_wide <- spread(results, key = im_cell, value = Cor)
rownames(results_wide) <- results_wide$Gene
results_wide$Gene <- NULL
#Plot the correlation matrix
ggcorrplot(t(results_wide), 
           show.legend = T, 
           colors = c("#2166AC", "white", "#B2182B"), 
           digits = 2, 
           lab = T)