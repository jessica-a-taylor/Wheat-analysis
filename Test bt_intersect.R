# Test bedtools intersect.

library(rtracklayer)
library(bedtoolsr)
library(bedr)

# Specify paths for files needed for the analysis.
path_to_significant_DEGs <- "Wheat-analysis/Data/Significant DEGs/"
path_to_all_DEGs <- "Wheat-analysis/Data/All DEGs/"
path_to_gene_list <- "Wheat-analysis/Data/Wheat protein coding genes.csv"

path_to_data_output <- "Wheat-analysis/Data/Enrichment results/"

# If not existing, create folders for the different data outputs of the enrichment analysis.
for (folder in c("All frequencies", "All proportions", "Average proportions")) {
  if (file.exists(paste(path_to_data_output, folder, sep = ""))==FALSE) {
    dir.create(file.path(paste(path_to_data_output, folder, sep = "")))
  } else next
}

path_to_graph_output <- "Wheat-analysis/Graphs/"

# If not existing, create folders for the graphical outputs of the enrichment analysis.
for (folder in c("Enrichment per region", "Percent modified genes")) {
  if (file.exists(paste(path_to_graph_output, folder, sep = ""))==FALSE) {
    dir.create(file.path(paste(path_to_graph_output, folder, sep = "")))
  } else next
}

rm(folder)
# Load required functions.
source("Wheat-analysis/Functions/getGeneCoordinates.R")
source("Wheat-analysis/Functions/Gene width&range.R")
source("Wheat-analysis/Functions/AxisGroup column.R")

# Import ChIP-seq data using 'getChIP_seq_data.R function.
source("Wheat-analysis/Functions/getChIP_seq_data.R")
nextflowOutput <- getChIP_seq_data()
write.csv(nextflowOutput, "Wheat-analysis/Data/Nextflow output summary.csv")

# Convert 'nextflowOutput dataset into bedfile for the bt.intersect function.
nextflowOutput <- nextflowOutput[order(nextflowOutput$seqnames, nextflowOutput$start, nextflowOutput$end),]
nextflowOutputBed <- GRanges(nextflowOutput[,c(1:3,8)])

# Import sample gene sets from 'Significant DEGs' folder. Store in a list.
sampleGenes <- list()

for (file in list.files(path_to_significant_DEGs)) {
  # Use 'str_match' function to name each dataset stored in the list e.g. "Asory upreg".
  sampleGenes[[str_match(file, "^([a-zA-Z]+.*).csv$")[,2]]] <- as.data.frame(read.csv(paste(path_to_significant_DEGs, file, sep = "")))
  
  colnames(sampleGenes[[str_match(file, "^([a-zA-Z]+.*).csv$")[,2]]]) <- c("Gene", "seqnames")
}

# Import control gene sets from 'All DEGs' folder. Add to the list.
for (file in list.files(path_to_all_DEGs)) {
  df <- as.data.frame(read.csv(paste(path_to_all_DEGs, file, sep = "")))
  
  # Remove genes not assigned to a chromosome.
  df <- data.frame(Gene = df[-which(str_match(df$Gene, "^TraesCS(U)02G[0-9]+$")[,2]=="U"),])
  
  # Extract the chromosome from the gene name and create a 'seqnames' column.
  df$seqnames <- str_match(df$Gene, "^TraesCS([0-9]+[a-zA-Z]+)02G[0-9]+$")[,2]
  
  # Sample 10000 random control genes from the unfiltered list of DEGs.
  # Use 'str_match' function to name each dataset stored in the list e.g. "Asory control".
  sampleGenes[[paste(str_match(file, "^([a-zA-Z]+.*) vs [a-zA-Z]+.*.csv$")[,2], "control", sep = " ")]] <- df[sample(1:nrow(df), 10000),]
}

rm (file, df)

# Import list of all protein coding genes.
genomeData <- as.data.frame(read.csv(path_to_gene_list))
genomeData <- genomeData[,-1]

# Perform enrichment analysis on each gene set.
geneSet <- "Asory downreg"
mod <- "H3K4me3"
region <- "Gene20"

df <- sampleGenes[[geneSet]]
  
  # Create excel workbook into which the output data will be saved.
  allProportions_wb <- createWorkbook()
  averageProportions_wb <- createWorkbook()
  allFrequency_wb <- createWorkbook()
  
  # Use the 'getGeneCoordinates' function to get the coordinates for each 
  # genomic region for the current gene set.
  geneRegions <- getGeneCoordinates(df, genomeData)
  
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
    
      # Create a bed file with the coordinates for the current gemonic region.
    queryBed <- geneRegions[[region]][,c("Gene","seqnames","start","end","width")]
    #queryBed <- queryBed[order(queryBed$seqnames, queryBed$start, queryBed$Gene),]
    
    queryBed <- bedr.sort.region(queryBed, check.chr = FALSE)

      queryBed <- GRanges(queryBed)

      # Use bedtools intersect function to find the overlap between the genomic
      # region and ChIP-seq peaks for the current modification/TF.
      overlap <- bt.intersect(peaksPerModification, queryBed, wo = TRUE, sorted = TRUE) 
      
      # Merge multiple peaks overlapping the same region.
      # Determine the proportion of overlap, and store in 'mergedOverlap'
      #mergedOverlap <- data.frame()
      #
      #for (gene in unique(geneRegions[[region]]$Gene)) {
      #  row <- which(geneRegions[[region]]$Gene == gene)
      #  
      #  # If the gene overlaps with one or more ChIP-seq peaks in the current 
      #  # region, find the sum of the peaks.
      #  if (gene %in% overlap$V12) {
      #    peakOverlaps <- overlap[overlap$V12==gene,]
      #    
      #    mergedOverlap <- rbind(mergedOverlap, data.frame(Gene = gene,
      #                                                     seqnames = geneRegions[[region]][row, "seqnames"],
      #                                                     start = geneRegions[[region]][row, "start"],
      #                                                     end = geneRegions[[region]][row, "end"],
      #                                                     width = geneRegions[[region]][row, "width"],
      #                                                     overlap = sum(peakOverlaps$V13)/peakOverlaps$V10[1]))
      #    
      #    # If the gene does not overlap with a ChIP-seq peak in the current 
      #    # region, overlap equals zero.
      #  } else mergedOverlap <- rbind(mergedOverlap, data.frame(Gene = gene,
      #                                                          seqnames = geneRegions[[region]][row, "seqnames"],
      #                                                          start = geneRegions[[region]][row, "start"],
      #                                                          end = geneRegions[[region]][row, "end"],
      #                                                          width = geneRegions[[region]][row, "width"],
      #                                                          overlap = 0))
      #  
      #}
      ## Maximum proportion of overlap equals 1.
      #mergedOverlap[which(mergedOverlap$overlap > 1), "overlap"] <- 1
      #
      ## Add the overlaps data for all genes in all regions to 'allProportions_data'.
      #allProportions_data <- rbind(allProportions_data, data.frame(Gene = mergedOverlap$Gene,
      #                                                             seqnames = mergedOverlap$seqnames,
      #                                                             start = mergedOverlap$start,
      #                                                             end = mergedOverlap$end, 
      #                                                             width = mergedOverlap$width,
      #                                                             region = rep(region, times = nrow(mergedOverlap)),
      #                                                             Mod.TF = rep(mod, times = nrow(mergedOverlap)),
      #                                                             overlap = mergedOverlap$overlap,
      #                                                             geneSet = rep(geneSet, times = nrow(mergedOverlap)),
      #                                                             axisGroup = rep(geneRegionAxisLocations(region), times = nrow(mergedOverlap))))
      #
      ## Add the average overlap across genes in each region to 'averageProportions_data'.
      #averageProportions_data <- rbind(averageProportions_data, data.frame(region = region,
      #                                                                     Mod.TF = mod,
      #                                                                     mean.overlap = mean(mergedOverlap$overlap),
      #                                                                     sd.overlap = sd(mergedOverlap$overlap),
      #                                                                     frequency = signif((nrow(mergedOverlap[which(mergedOverlap$overlap>0.3),])/nrow(mergedOverlap))*100, digits = 3),
      #                                                                     geneSet = geneSet,
      #                                                                     axisGroup = geneRegionAxisLocations(region)))
    #
    # Add the 'proportions' datasets to the corresponding excel worksheet.
    #writeData(allProportions_wb, sheet = mod, allProportions_data)
    #writeData(averageProportions_wb, sheet = mod, averageProportions_data)
    #
    ## Determine the frequency of overlap across all regions, except
    ## the intergenic space.
    #frequency_data <- rbind(frequency_data, 
    #                        data.frame(Mod.TF = mod,
    #                                   frequency = length(unique(allProportions_data[which(allProportions_data$overlap > 0 &
    #                                                                                         grepl("Intergenic", allProportions_data$region)==FALSE), 
    #                                                                                 "Gene"]))/length(unique(allProportions_data$Gene))*100,
    #                                   geneSet = str_match(geneSet, "^.* ([a-z]+)$")[,2],
    #                                   genotype = str_match(geneSet, "^(.*) [a-z]+$")[,2]))
    #
    ## Add the 'frequencies' data to the corresponding excel worksheet.
    #writeData(allFrequency_wb, sheet = mod, frequency_data)
  
  # Save each worksheet to the corresponding workbook and export as an .xlsx file.
  #saveWorkbook(allProportions_wb, paste(path_to_data_output, "All proportions/", geneSet, ".xlsx", sep = ""), overwrite = TRUE) 
  #saveWorkbook(averageProportions_wb, paste(path_to_data_output, "Average proportions/", geneSet, ".xlsx", sep = ""), overwrite = TRUE)
  #saveWorkbook(allFrequency_wb, paste(path_to_data_output, "All frequencies/", geneSet, ".xlsx", sep = ""), overwrite = TRUE)
