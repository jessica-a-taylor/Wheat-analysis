# Import list of all protein coding genes.
genomeData <- as.data.frame(read.csv(path_to_gene_list))
genomeData <- genomeData[,-1]

# Perform enrichment analysis on each gene set.
for (geneSet in names(sampleGenes)) {
  print(geneSet)
  
  # For each gene set...
  df <- sampleGenes[[geneSet]]
  
  # Create excel workbook into which the output data will be saved.
  allProportions_wb <- createWorkbook()
  averageProportions_wb <- createWorkbook()
  allFrequency_wb <- createWorkbook()
  
  # Use the 'getGeneCoordinates' function to get the coordinates for each 
  # genomic region for the current gene set.
  geneRegions <- getGeneCoordinates(df, genomeData)
  
  # Use bedtools intersect to find the overlap between ChIP-seq peaks and each gene region.
  for (mod in unique(nextflowOutput$Mod.TF)) {
    
    # Add a worksheet to each workbook in which the output data will be saved.
    addWorksheet(allProportions_wb, sheetName = mod)
    addWorksheet(averageProportions_wb, sheetName = mod)
    addWorksheet(allFrequency_wb, sheetName = mod)
    
    # Create empty dataframes for temporarily storing the output.
    allProportions_data <- data.frame()
    averageProportions_data <- data.frame()
    frequency_data <- data.frame()
    
    # Filter for the current modification/TF.
    peaksPerModification <- nextflowOutputBed[which(nextflowOutputBed$Mod.TF==mod),]
    
    for (region in names(geneRegions)) {
      # Create a bed file with the coordinates for the current gemonic region.
      queryBed <- GRanges(geneRegions[[region]][,c("Gene","seqnames","start","end","width")])
      
      # Use bedtools intersect function to find the overlap between the genomic
      # region and ChIP-seq peaks for the current modification/TF.
      overlap <- bt.intersect(peaksPerModification, queryBed, wo = TRUE) 
      
      # Merge multiple peaks overlapping the same region.
      # Determine the proportion of overlap, and store in 'mergedOverlap'
      mergedOverlap <- data.frame()
      
      for (gene in unique(geneRegions[[region]]$Gene)) {
        row <- which(geneRegions[[region]]$Gene == gene)
        
        # If the gene overlaps with one or more ChIP-seq peaks in the current 
        # region, find the sum of the peaks.
        if (gene %in% overlap$V12) {
          peakOverlaps <- overlap[overlap$V12==gene,]
          
          mergedOverlap <- rbind(mergedOverlap, data.frame(Gene = gene,
                                                           seqnames = geneRegions[[region]][row, "seqnames"],
                                                           start = geneRegions[[region]][row, "start"],
                                                           end = geneRegions[[region]][row, "end"],
                                                           width = geneRegions[[region]][row, "width"],
                                                           overlap = sum(peakOverlaps$V13)/peakOverlaps$V10[1]))
          
        # If the gene does not overlap with a ChIP-seq peak in the current 
        # region, overlap equals zero.
        } else mergedOverlap <- rbind(mergedOverlap, data.frame(Gene = gene,
                                                                seqnames = geneRegions[[region]][row, "seqnames"],
                                                                start = geneRegions[[region]][row, "start"],
                                                                end = geneRegions[[region]][row, "end"],
                                                                width = geneRegions[[region]][row, "width"],
                                                                overlap = 0))
        
      }
      # Maximum proportion of overlap equals 1.
      mergedOverlap[which(mergedOverlap$overlap > 1), "overlap"] <- 1
      
      # Add the overlaps data for all genes in all regions to 'allProportions_data'.
      allProportions_data <- rbind(allProportions_data, data.frame(Gene = mergedOverlap$Gene,
                                                                   seqnames = mergedOverlap$seqnames,
                                                                   start = mergedOverlap$start,
                                                                   end = mergedOverlap$end, 
                                                                   width = mergedOverlap$width,
                                                                   region = rep(region, times = nrow(mergedOverlap)),
                                                                   Mod.TF = rep(mod, times = nrow(mergedOverlap)),
                                                                   overlap = mergedOverlap$overlap,
                                                                   geneSet = rep(geneSet, times = nrow(mergedOverlap)),
                                                                   axisGroup = rep(geneRegionAxisLocations(region), times = nrow(mergedOverlap))))
      
      # Add the average overlap across genes in each region to 'averageProportions_data'.
      averageProportions_data <- rbind(averageProportions_data, data.frame(region = region,
                                                                           Mod.TF = mod,
                                                                           mean.overlap = mean(mergedOverlap$overlap),
                                                                           sd.overlap = sd(mergedOverlap$overlap),
                                                                           frequency = signif((nrow(mergedOverlap[which(mergedOverlap$overlap>0.3),])/nrow(mergedOverlap))*100, digits = 3),
                                                                           geneSet = geneSet,
                                                                           axisGroup = geneRegionAxisLocations(region)))
    }
    # Add the 'proportions' datasets to the corresponding excel worksheet.
    writeData(allProportions_wb, sheet = mod, allProportions_data)
    writeData(averageProportions_wb, sheet = mod, averageProportions_data)
    
    # Determine the frequency of overlap across all regions, except
    # the intergenic space.
    frequency_data <- rbind(frequency_data, 
                            data.frame(Mod.TF = mod,
                                       frequency = length(unique(allProportions_data[which(allProportions_data$overlap > 0 &
                                                                                             grepl("Intergenic", allProportions_data$region)==FALSE), 
                                                                                     "Gene"]))/length(unique(allProportions_data$Gene))*100,
                                       geneSet = str_match(geneSet, "^.* ([a-z]+)$")[,2],
                                       genotype = str_match(geneSet, "^(.*) [a-z]+$")[,2]))
    
    # Add the 'frequencies' data to the corresponding excel worksheet.
    writeData(allFrequency_wb, sheet = mod, frequency_data)
  }
  # Save each worksheet to the corresponding workbook and export as an .xlsx file.
  #saveWorkbook(allProportions_wb, paste(path_to_data_output, "All proportions/", geneSet, ".xlsx", sep = ""), overwrite = TRUE) 
  #saveWorkbook(averageProportions_wb, paste(path_to_data_output, "Average proportions/", geneSet, ".xlsx", sep = ""), overwrite = TRUE)
  #saveWorkbook(allFrequency_wb, paste(path_to_data_output, "All frequencies/", geneSet, ".xlsx", sep = ""), overwrite = TRUE)
}

rm(df, nextflowOutputBed, queryBed, allProportions_wb, averageProportions_wb, allProportions_data, allFrequency_wb, allFrequency_data, averageProportions_data, overlap, region, geneSet)

# Plot results.
# Import 'proportions' and 'frequencies' data fo all gene set into one dataframe.
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
    
  plot <- facet(plot, facet.by = "Genotype", nrow = 1)
  plot <- ggpar(plot, font.xtickslab = FALSE, ticks = FALSE, legend = "bottom", xlab = FALSE, legend.title = "",
                font.ytickslab = 8)
  
  pdf(paste(path_to_graph_output, "Percent modified genes/", mod, ".pdf", sep = ""), width = 6, height = 5)
  print(plot)
  dev.off()
}