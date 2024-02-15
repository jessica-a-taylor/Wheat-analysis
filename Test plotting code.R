# Load required packages.
source("Wheat-analysis/Functions/loadLibraries.R")
loadLibraries()
rm(loadLibraries)

# Specify paths for files needed for the analysis.
path_to_significant_DEGs <- "Wheat-analysis/Data/Significant DEGs/"
path_to_all_DEGs <- "Wheat-analysis/Data/All DEGs/"
path_to_gene_list <- "Wheat-analysis/Data/Wheat protein coding genes.csv"

path_to_data_output <- "Wheat-analysis/Data/Enrichment results/"
path_to_graph_output <- "Wheat-analysis/Graphs/"

source("Wheat-analysis/Functions/getChIP_seq_data.R")
source("Wheat-analysis/Functions/getGeneCoordinates.R")
source("Wheat-analysis/Functions/Gene width&range.R")
source("Wheat-analysis/Functions/AxisGroup column.R")

# Import ChIP-seq data using 'getChIP_seq_data.R function.
nextflowOutput <- getChIP_seq_data()
rm(getChIP_seq_data)

# Plot results.
# Import 'proportions' and 'frequencies' data for all gene sets into one dataframe.
allFrequency_results <- data.frame()

for (mod in unique(nextflowOutput$Mod.TF)) {
  for (genotype in c("Asory", "Gabo_1RS", "Campesino", "Rialto", "Savannah")) {
    
    allProportions_results <- data.frame()
    averageProportions_results <- data.frame()
    
    for (set in which(grepl(genotype, list.files(paste(path_to_data_output, "All proportions", sep = "")))==TRUE)) {
      
      allProportions <- read.xlsx(paste(path_to_data_output, "All proportions/", list.files(paste(path_to_data_output, "All proportions", sep = ""))[set], sep = ""), sheet = mod)
      allProportions_results <- rbind(allProportions_results, allProportions)
    }
    for (set in which(grepl(genotype, list.files(paste(path_to_data_output, "Average proportions", sep = "")))==TRUE)) {
      
      averageProportions <- read.xlsx(paste(path_to_data_output, "Average proportions/", list.files(paste(path_to_data_output, "Average proportions", sep = ""))[set], sep = ""), sheet = mod)
      averageProportions_results <- rbind(averageProportions_results, averageProportions)
    }
    for (set in which(grepl(genotype, list.files(paste(path_to_data_output, "All frequencies", sep = "")))==TRUE)) {
      
      allFrequency <- read.xlsx(paste(path_to_data_output, "All frequencies/", list.files(paste(path_to_data_output, "All frequencies", sep = ""))[set], sep = ""), sheet = mod)
      allFrequency_results <- rbind(allFrequency_results, allFrequency)
    }
    
    # Specify the pair-wise comparisons to be used for the t-tests.
    my_comparisons <- list(c(unique(allProportions_results$geneSet)[which(grepl("control", unique(allProportions_results$geneSet))==TRUE)], 
                             unique(allProportions_results$geneSet)[which(grepl("upreg", unique(allProportions_results$geneSet))==TRUE)]), 
                           c(unique(allProportions_results$geneSet)[which(grepl("control", unique(allProportions_results$geneSet))==TRUE)], 
                             unique(allProportions_results$geneSet)[which(grepl("downreg", unique(allProportions_results$geneSet))==TRUE)]),
                           c(unique(allProportions_results$geneSet)[which(grepl("upreg", unique(allProportions_results$geneSet))==TRUE)], 
                             unique(allProportions_results$geneSet)[which(grepl("downreg", unique(allProportions_results$geneSet))==TRUE)]))
    
    # Perform t-tests and specify the y-axis coordinates for plotting the results.
    stat.test <- allProportions_results %>% group_by(axisGroup) %>% 
      t_test(overlap ~ geneSet, comparisons = my_comparisons) %>% 
      mutate(y.position = rep(c(1.08, 0.96, 1.01), times = 10))
    
    # Plot the average proportion of overlap for each gene set.
    plot <- ggbarplot(averageProportions_results, x = "geneSet", y="mean.overlap", ylab = "Average Enrichment",
                      color = "black", fill = "geneSet", lab.size = 12,
                      palette = c("bisque2", "darksalmon", "coral3"), 
                      title = mod) + theme_bw() +
      # Add the results of the t-tests.
      stat_pvalue_manual(
        stat.test, 
        label = "p.adj.signif", size = 4,
        tip.length = 0.01, hide.ns = FALSE) +
      
      # Set the limits of the y-axis.
      coord_cartesian(ylim= c(0,1.1), clip = "off") +
      
      # Set font sizes.
      font("title", size = 14) +
      font("ylab", size = 12) +
      font("legend.title", size = 12) +
      font("legend.text", size = 10) +
      font("caption", size = 12) +
      font("axis.text", size = 14)
    
    # Use the 'facet' function to plot separate graphs for each genomic region.
    plot <- facet(plot, facet.by = "axisGroup", nrow = 1, panel.labs.font = list(size = 10),
                  panel.labs = list(axisGroup = c("Intergenic","Promotor \n(1kb)","Promotor \n(500bp)", "20%",            
                                                  "40%","60%","80%","100%","Downstream \n(200bp)","Intergenic")))
    
    # Adjust the appearance of the plot e.g. legend position.
    plot <- ggpar(plot, font.xtickslab = FALSE, ticks = FALSE, legend = "bottom", xlab = FALSE, legend.title = "",
                  font.ytickslab = 8)
    
    # Export each plot as a PDF.
    pdf(paste(path_to_graph_output, "Enrichment per region/", genotype, " ", mod, ".pdf", sep = ""), width = 12, height = 5)
    print(plot)
    dev.off()
  }
}

# Plot the frequency of overlap across all regions, except the intergenic space
# for each modification.
for (mod in unique(nextflowOutput$Mod.TF)) {
  df <- allFrequency_results[allFrequency_results$Mod.TF==mod,]
  
  plot <- ggbarplot(df, x = "geneSet", y="frequency", ylab = "% of modified genes",
                    color = "black", fill = "geneSet", lab.size = 12,
                    palette = c("bisque2", "darksalmon", "coral3"), 
                    title = mod) + theme_bw() +
    coord_cartesian(ylim= c(0,100), clip = "off") +
    
    font("title", size = 14) +
    font("ylab", size = 12) +
    font("legend.title", size = 12) +
    font("legend.text", size = 10) +
    font("caption", size = 12) +
    font("axis.text", size = 14) 
  
  plot <- facet(plot, facet.by = "genotype", nrow = 1)
  plot <- ggpar(plot, font.xtickslab = FALSE, ticks = FALSE, legend = "bottom", xlab = FALSE, legend.title = "",
                font.ytickslab = 8)
  
  pdf(paste(path_to_graph_output, "Percent modified genes/", mod, ".pdf", sep = ""), width = 6, height = 5)
  print(plot)
  dev.off()
}