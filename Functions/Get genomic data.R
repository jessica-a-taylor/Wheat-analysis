library(biomaRt)

# Import genome.
genomeData <- as.data.frame(read.csv("Wheat data/IWGSC_v1.1_HC_20170706_cds.csv"))

# Filter for protein-coding genes.
geneTypes <- useEnsemblGenomes(biomart = "plants_mart", dataset = "taestivum_eg_gene")
geneTypes <- as.data.frame(getBM(attributes = listAttributes(geneTypes)[c(1,6:9,20),"name"], mart = geneTypes))
geneTypes <- geneTypes[geneTypes$gene_biotype=="protein_coding",]

genomeData <- genomeData[which(genomeData$Gene %in% geneTypes$ensembl_gene_id),]

# Remove duplicate genes.
genomeData <- genomeData[which(genomeData$GeneID == str_match(genomeData$GeneID, "^([a-zA-Z0-9]*.1)$")[,2]),]

# Add 'width' column.
source("Functions/Gene width&range.R")
genomeData$width <- getWidth(genomeData)

# Change the '+/-' strand to 'positive/negative'.
genomeData[which(genomeData$strand=="+"), "strand"] <- "positive"
genomeData[which(genomeData$strand=="-"), "strand"] <- "negative"

# Add 'expressionLevel' column to each dataset.
expressionData <- as.data.frame(read.csv("Wheat data/Wheat expression data.csv"))

genomeData <- genomeData[which(genomeData$Gene %in% expressionData$transcript),]
genomeData$TPM <- expressionData$Average

write.csv(genomeData, file = "Wheat data/Wheat protein coding genes.csv")

