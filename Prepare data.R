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
write.csv(nextflowOutput, "Wheat-introgression-analysis/Data/Nextflow output summary.csv")

# Convert 'nextflowOutput dataset into bedfile for the bt.intersect function.
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
