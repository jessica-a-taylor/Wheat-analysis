# Import expression data.
expressionData <- as.data.frame(read.csv("Wheat data/PRJNA289545 wheat expression.csv"))

expressionData$transcript <- paste(str_extract(expressionData$transcript, "TraesCS[0-9][A-Z]0"), 2,
                                   str_extract(expressionData$transcript, "G[0-9]+"), sep = "")

# Filter out 'low confidence' (LC) genes and those not assigned to a chromosome.
expressionData <- expressionData[-which(grepl("U", expressionData$transcript) | grepl("LC", expressionData$transcript)),]

for (gene in unique(str_match(expressionData$transcript, "^(.*).[0-9]+")[,2])) {
  if (nrow(expressionData[str_match(expressionData$transcript, "^(.*).[0-9]+")[,2]==gene,]) > 1) {

    # For genes with multiple splice variants, use the variant with the highest expression.
    # If both variants' expression have equal expression, just use one.
    # Remove the other variants.
    if (mean(expressionData[str_match(expressionData$transcript, "^(.*).[0-9]+")[,2]==gene,"Average"]) == expressionData[str_match(expressionData$transcript, "^(.*).[0-9]+")[,2]==gene,"Average"][1]) {
      expressionData <- expressionData[-which(str_match(expressionData$transcript, "^(.*).[0-9]+")[,2] == gene &
                                                 expressionData$transcript != paste(gene, ".1", sep = "")),]
    } 
    else if (sum(expressionData[str_match(expressionData$transcript, "^(.*).[0-9]+")[,2]==gene,"Average"]) > 0) {
      expressionData <- expressionData[-which(str_match(expressionData$transcript, "^(.*).[0-9]+")[,2]==gene &
                                                expressionData$Average != max(expressionData[str_match(expressionData$transcript, "^(.*).[0-9]+")[,2]==gene,
                                                                                             "Average"])),]
    }
  } else next
}

expressionData$transcript <- str_match(expressionData$transcript, "^(.*).1$")[,2]

expressionData <- expressionData[which(expressionData$transcript %in% genomeData$Gene),]

write.csv(expressionData, "Wheat data/Wheat expression data.csv")