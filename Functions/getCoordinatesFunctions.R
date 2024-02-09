# Function of determining the coordinates of the intergenic regions between the current and adjacent gene, 
# accounting for their orientation on the DNA strands.

intergenicCoordinatesFunction <- function(df, genomeData, region) {
  
  # Variations of the functions depending on gene orientations.
  getIntergenicCoordinates <- list(UpstreamIntergenic = 
                                     list(current_positive =
                                            list(getAdjacentGene = function(currentGene) {return(currentGene - 1)},
                                                 
                                                 adjacent_positive =
                                                   list(getCoordinates = function(currentGene, adjacentGene, geneCoordinates, gene, genomeData) {
                                                     if ((genomeData[currentGene, "start"] - 1001)-(genomeData[adjacentGene, "end"] + 201) > 0) {
                                                       geneCoordinates$start <- genomeData[adjacentGene, "end"] + 201
                                                       geneCoordinates$end <- genomeData[currentGene, "start"] - 1001
                                                       geneCoordinates$width <- geneCoordinates$end - geneCoordinates$start
                                                       geneCoordinates$ranges <- paste(geneCoordinates$start, "-", geneCoordinates$end, sep = "")
                                                     } else geneCoordinates <- geneCoordinates[-1,]
                                                     return(geneCoordinates)
                                                     }),
                                                 
                                                 adjacent_negative =
                                                   list(getCoordinates = function(currentGene, adjacentGene, geneCoordinates, gene, genomeData) {
                                                     if ((genomeData[currentGene, "start"] - 1001) - (genomeData[adjacentGene, "end"] + 1001) > 0) {
                                                       geneCoordinates$start <- genomeData[adjacentGene, "end"] + 1001
                                                       geneCoordinates$end <- genomeData[currentGene, "start"] - 1001
                                                       geneCoordinates$width <- geneCoordinates$end - geneCoordinates$start
                                                       geneCoordinates$ranges <- paste(geneCoordinates$start, "-", geneCoordinates$end, sep = "")
                                                     } else geneCoordinates <- geneCoordinates[-1,]
                                                     return(geneCoordinates)
                                                     })),
                                          
                                          current_negative =
                                            list(getAdjacentGene = function(currentGene) {return(currentGene + 1)},
                                                 
                                                 adjacent_positive =
                                                   list(getCoordinates = function(currentGene, adjacentGene, geneCoordinates, gene, genomeData) {
                                                     if ((genomeData[adjacentGene, "start"] - 1001) - (genomeData[currentGene, "end"] + 1001) > 0) {
                                                       geneCoordinates$start <- genomeData[currentGene, "end"] + 1001
                                                       geneCoordinates$end <- genomeData[adjacentGene, "start"] - 1001
                                                       geneCoordinates$width <- geneCoordinates$end - geneCoordinates$start
                                                       geneCoordinates$ranges <- paste(geneCoordinates$start, "-", geneCoordinates$end, sep = "")
                                                     } else geneCoordinates <- geneCoordinates[-1,]
                                                     return(geneCoordinates)
                                                     }),
                                                 
                                                 adjacent_negative =
                                                   list(getCoordinates = function(currentGene, adjacentGene, geneCoordinates, gene, genomeData) {
                                                     if ((genomeData[adjacentGene, "start"] - 201) - (genomeData[currentGene, "end"] + 1001) > 0) {
                                                       geneCoordinates$start <- genomeData[currentGene, "end"] + 1001
                                                       geneCoordinates$end <- genomeData[adjacentGene, "start"] - 201
                                                       geneCoordinates$width <- geneCoordinates$end - geneCoordinates$start
                                                       geneCoordinates$ranges <- paste(geneCoordinates$start, "-", geneCoordinates$end, sep = "")
                                                     } else geneCoordinates <- geneCoordinates[-1,]
                                                     return(geneCoordinates)
                                                     }))),
                                   
                                   DownstreamIntergenic = 
                                     list(current_positive =
                                            list(getAdjacentGene = function(currentGene) {return(currentGene + 1)},
                                                 
                                                 adjacent_positive =
                                                   list(getCoordinates = function(currentGene, adjacentGene, geneCoordinates, gene, genomeData) {
                                                     if ((genomeData[adjacentGene, "start"] - 1001) - (genomeData[currentGene, "end"] + 201) > 0) {
                                                       geneCoordinates$start <- genomeData[currentGene, "end"] + 201
                                                       geneCoordinates$end <- genomeData[adjacentGene, "start"] - 1001
                                                       geneCoordinates$width <- geneCoordinates$end - geneCoordinates$start
                                                       geneCoordinates$ranges <- paste(geneCoordinates$start, "-", geneCoordinates$end, sep = "")
                                                     } else geneCoordinates <- geneCoordinates[-1,]
                                                     return(geneCoordinates)
                                                     }),
                                                 
                                                 adjacent_negative =
                                                   list(getCoordinates = function(currentGene, adjacentGene, geneCoordinates, gene, genomeData) {
                                                     if ((genomeData[adjacentGene, "start"] - 201) - (genomeData[currentGene, "end"] + 201) > 0) {
                                                       geneCoordinates$start <- genomeData[currentGene, "end"] + 201
                                                       geneCoordinates$end <- genomeData[adjacentGene, "start"] - 201
                                                       geneCoordinates$width <- geneCoordinates$end - geneCoordinates$start
                                                       geneCoordinates$ranges <- paste(geneCoordinates$start, "-", geneCoordinates$end, sep = "")
                                                     } else geneCoordinates <- geneCoordinates[-1,]
                                                     return(geneCoordinates)
                                                     })),
                                          
                                          current_negative =
                                            list(getAdjacentGene = function(currentGene) {return(currentGene - 1)},
                                                 
                                                 adjacent_positive =
                                                   list(getCoordinates = function(currentGene, adjacentGene, geneCoordinates, gene, genomeData) {
                                                     if ((genomeData[currentGene, "start"] - 201) - (genomeData[adjacentGene, "end"] + 201) > 0) {
                                                       geneCoordinates$start <- genomeData[adjacentGene, "end"] + 201
                                                       geneCoordinates$end <- genomeData[currentGene, "start"] - 201
                                                       geneCoordinates$width <- geneCoordinates$end - geneCoordinates$start
                                                       geneCoordinates$ranges <- paste(geneCoordinates$start, "-", geneCoordinates$end, sep = "")
                                                     } else geneCoordinates <- geneCoordinates[-1,]
                                                     return(geneCoordinates)
                                                     }),
                                                 
                                                 adjacent_negative =
                                                   list(getCoordinates = function(currentGene, adjacentGene, geneCoordinates, gene, genomeData) {
                                                     if ((genomeData[currentGene, "start"] - 201) - (genomeData[adjacentGene, "end"] + 1001) > 0) {
                                                       geneCoordinates$start <- genomeData[adjacentGene, "end"] + 1001
                                                       geneCoordinates$end <- genomeData[currentGene, "start"] - 201
                                                       geneCoordinates$width <- geneCoordinates$end - geneCoordinates$start
                                                       geneCoordinates$ranges <- paste(geneCoordinates$start, "-", geneCoordinates$end, sep = "")
                                                     } else geneCoordinates <- geneCoordinates[-1,]
                                                     return(geneCoordinates)
                                                     })))) 
  
  # Determine coordinates of intergenic regions.
  if (nrow(df) >= 1) {
    
    # Create dataframe for storing the information for each gene.
    newGeneCoordinates <- data.frame()
    
    df <- genomeData[which(genomeData$Gene %in% df$Gene),]
  
    for (gene in df$Gene) {
      
      geneCoordinates <- df[df$Gene==gene,]
      
      currentGene <- which(genomeData$Gene==gene)
      currentStrand <- paste("current_", genomeData[currentGene,"strand"], sep = "")
      
      adjacentGene <- getIntergenicCoordinates[[region]][[currentStrand]]$getAdjacentGene(currentGene)
      adjacentStrand <- paste("adjacent_", genomeData[adjacentGene,"strand"], sep = "")
      
      newGeneCoordinates <- rbind(newGeneCoordinates,
                                  getIntergenicCoordinates[[region]][[currentStrand]][[adjacentStrand]]$getCoordinates(currentGene, adjacentGene, geneCoordinates, gene, genomeData))
    }
  }
  return(newGeneCoordinates)
}

# Function for determining the coordinates of the promotor regions 500 bp and 1000 bp upstream of the TSS.
promotorCoordinatesFunction <- function(df, region) {
  
  # Variations of the functions depending on gene orientations.
  getPromotorCoordinates <- list(Promotor500 =
                                   list(current_positive = 
                                          list(getCoordinates = function(data) {
                                            data$start <- data$start-500
                                            data$end <- data$start+500
                                            data$width <- data$end - data$start
                                            data$ranges <- paste(data$start, "-", data$end, sep = "")
                                            return(data)}),
                                     
                                     current_negative = 
                                          list(getCoordinates = function(data) {
                                            data$start <- data$end
                                            data$end <- data$end+500
                                            data$width <- data$end - data$start
                                            data$ranges <- paste(data$start, "-", data$end, sep = "")
                                            return(data)})),
                                 
                                 Promotor1000 =
                                   list(current_positive = 
                                          list(getCoordinates = function(data) {
                                            data$start <- data$start-1000
                                            data$end <- data$start+1000
                                            data$width <- data$end - data$start
                                            data$ranges <- paste(data$start, "-", data$end, sep = "")
                                            return(data)}),
                                     
                                        current_negative = 
                                          list(getCoordinates = function(data) {
                                            data$start <- data$end
                                            data$end <- data$end+1000
                                            data$width <- data$end - data$start
                                            data$ranges <- paste(data$start, "-", data$end, sep = "")
                                            return(data)})))
  
  
  # Determine coordinates of promotor regions.
  if (nrow(df) >= 1) {
    
    # Create dataframe for storing the information for each gene.
    newGeneCoordinates <- data.frame()
    
    df <- genomeData[which(genomeData$Gene %in% df$Gene),]
    
    for (row in 1:nrow(df)) {
      geneCoordinates <- df[row,]
      
      currentStrand <- paste("current_", df[row,"strand"], sep = "")
      
      newGeneCoordinates <- rbind(newGeneCoordinates, getPromotorCoordinates[[region]][[currentStrand]]$getCoordinates(geneCoordinates))
    }
  }
  return(newGeneCoordinates)
}

# Function for determining the coordinates of the 200 bp downstream regions.
downstreamCoordinatesFunction <- function(df) {
 
  # Variations of the functions depending on gene orientations.
  getDownstreamCoordinates <- list(current_positive = 
                                   list(getCoordinates = function(data) {
                                     data$start <- data$end
                                     data$end <- data$end + 200
                                     data$width <- data$end - data$start
                                     data$ranges <- paste(data$start, "-", data$end, sep = "")
                                     return(data)}),
                                 
                                 current_negative = 
                                   list(getCoordinates = function(data) {
                                     data$start <- data$start - 200
                                     data$end <- data$start + 200
                                     data$width <- data$end - data$start
                                     data$ranges <- paste(data$start, "-", data$end, sep = "")
                                     return(data)}))
  
  
  # Determine coordinates of downstream regions.
  if (nrow(df) >= 1) {
    
    # Create dataframe for storing the information for each gene.
    newGeneCoordinates <- data.frame()
    
    df <- genomeData[which(genomeData$Gene %in% df$Gene),]
    
    for (row in 1:nrow(df)) {
      geneCoordinates <- df[row,]
      
      currentStrand <- paste("current_", df[row,"strand"], sep = "")
      
      newGeneCoordinates <- rbind(newGeneCoordinates, getDownstreamCoordinates[[currentStrand]]$getCoordinates(geneCoordinates))
    }
  }
  return(newGeneCoordinates)
}


# Function for determining the coordinates of the gene body in 20% intervals.
genebodyCoordinatesFunction <- function(df, geneRegions) {
  
  if (nrow(df) >= 1) {
    
    # Create dataframe for storing the information for each gene.
    for (region in c("Gene20", "Gene40", "Gene60", "Gene80", "Gene100")) {
      geneRegions[[region]] <- data.frame()
    }
    
    df <- genomeData[which(genomeData$Gene %in% df$Gene),]
    
    for (row in 1:nrow(df)) {
      for (region in c("Gene20", "Gene40", "Gene60", "Gene80", "Gene100")) {
        geneRegions[[region]] <- rbind(geneRegions[[region]], df[row,])
      }
      
      if (df[row, "strand"]=="positive") {
        
        geneRegions[["Gene20"]][row,c("start","end","width","ranges")] <- c(as.numeric(df[row,"start"]),
                                                     as.numeric(df[row,"start"]) + (as.numeric(df[row,"width"])*0.2),
                                                     (as.numeric(df[row,"start"]) + (as.numeric(df[row,"width"])*0.2) - as.numeric(df[row,"start"])),
                                                     paste(df[row,"start"], "-", df[row,"start"] + df[row,"width"]*0.2, sep = ""))
        
        geneRegions[["Gene40"]][row,c("start","end","width","ranges")] <- c(as.numeric(geneRegions[["Gene20"]][row,c("start","end","width","ranges")]$end) + 1, 
                                                     as.numeric(geneRegions[["Gene20"]][row,c("start","end","width","ranges")]$end) + 1 + (as.numeric(df[row,"width"])*0.2), 
                                                     (as.numeric(geneRegions[["Gene20"]][row,c("start","end","width","ranges")]$end) + 1 + (as.numeric(df[row,"width"])*0.2) - (as.numeric(geneRegions[["Gene20"]][row,c("start","end","width","ranges")]$end) + 1)),
                                                     paste(as.numeric(geneRegions[["Gene20"]][row,c("start","end","width","ranges")]$end) + 1, "-", as.numeric(geneRegions[["Gene20"]][row,c("start","end","width","ranges")]$end) + 1 + df[row,"width"]*0.2, sep = ""))
        
        geneRegions[["Gene60"]][row,c("start","end","width","ranges")] <- c(as.numeric(geneRegions[["Gene40"]][row,c("start","end","width","ranges")]$end) + 1,
                                                     as.numeric(geneRegions[["Gene40"]][row,c("start","end","width","ranges")]$end) + 1 + (as.numeric(df[row,"width"])*0.2),
                                                     (as.numeric(geneRegions[["Gene40"]][row,c("start","end","width","ranges")]$end) + 1 + (as.numeric(df[row,"width"])*0.2) - (as.numeric(geneRegions[["Gene40"]][row,c("start","end","width","ranges")]$end) + 1)),
                                                     paste(as.numeric(geneRegions[["Gene40"]][row,c("start","end","width","ranges")]$end) + 1, "-", as.numeric(geneRegions[["Gene40"]][row,c("start","end","width","ranges")]$end) + 1 + df[row,"width"]*0.2, sep = ""))
        
        geneRegions[["Gene80"]][row,c("start","end","width","ranges")] <- c(as.numeric(geneRegions[["Gene60"]][row,c("start","end","width","ranges")]$end) + 1,
                                                     as.numeric(geneRegions[["Gene60"]][row,c("start","end","width","ranges")]$end) + 1 + (as.numeric(df[row,"width"])*0.2),
                                                     (as.numeric(geneRegions[["Gene60"]][row,c("start","end","width","ranges")]$end) + 1 + (as.numeric(df[row,"width"])*0.2) - (as.numeric(geneRegions[["Gene60"]][row,c("start","end","width","ranges")]$end) + 1)),
                                                     paste(as.numeric(geneRegions[["Gene60"]][row,c("start","end","width","ranges")]$end) + 1, "-", as.numeric(geneRegions[["Gene60"]][row,c("start","end","width","ranges")]$end) + 1 + df[row,"width"]*0.2, sep = ""))
        
        geneRegions[["Gene100"]][row,c("start","end","width","ranges")] <- c(as.numeric(geneRegions[["Gene80"]][row,c("start","end","width","ranges")]$end) + 1, 
                                                      as.numeric(df[row,"end"]), 
                                                      as.numeric(df[row,"end"]) - (as.numeric(geneRegions[["Gene80"]][row,c("start","end","width","ranges")]$end) + 1),
                                                      paste(as.numeric(geneRegions[["Gene80"]][row,c("start","end","width","ranges")]$end) + 1, "-", as.numeric(df[row,"end"]), sep = ""))
      }
      else if (df[row, "strand"]=="negative"){
        geneRegions[["Gene20"]][row,c("start","end","width","ranges")] <- c(as.numeric(df[row,"end"]) - (as.numeric(df[row,"width"])*0.2),
                                                     as.numeric(df[row,"end"]),
                                                     as.numeric(df[row,"end"] - (df[row,"end"] - df[row,"width"]*0.2)),
                                                     paste(df[row,"end"] - df[row,"width"]*0.2, "-", df[row,"end"], sep = ""))
        
        geneRegions[["Gene40"]][row,c("start","end","width","ranges")] <- c(as.numeric(geneRegions[["Gene20"]][row,c("start","end","width","ranges")]$start) - (as.numeric(df[row,"width"])*0.2),
                                                     as.numeric(geneRegions[["Gene20"]][row,c("start","end","width","ranges")]$start) - 1,
                                                     (as.numeric(geneRegions[["Gene20"]][row,c("start","end","width","ranges")]$start) - 1) - (as.numeric(geneRegions[["Gene20"]][row,c("start","end","width","ranges")]$start) - (as.numeric(df[row,"width"])*0.2)),
                                                     paste(as.numeric(geneRegions[["Gene20"]][row,c("start","end","width","ranges")]$start) - (as.numeric(df[row,"width"])*0.2), "-", as.numeric(geneRegions[["Gene20"]][row,c("start","end","width","ranges")]$start) - 1, sep = ""))
        
        geneRegions[["Gene60"]][row,c("start","end","width","ranges")] <- c(as.numeric(geneRegions[["Gene40"]][row,c("start","end","width","ranges")]$start) - (as.numeric(df[row,"width"])*0.2),
                                                     as.numeric(geneRegions[["Gene40"]][row,c("start","end","width","ranges")]$start) - 1, 
                                                     (as.numeric(geneRegions[["Gene40"]][row,c("start","end","width","ranges")]$start) - 1) - (as.numeric(geneRegions[["Gene40"]][row,c("start","end","width","ranges")]$start) - (as.numeric(df[row,"width"])*0.2)),
                                                     paste(as.numeric(geneRegions[["Gene40"]][row,c("start","end","width","ranges")]$start) - (as.numeric(df[row,"width"])*0.2), "-", as.numeric(geneRegions[["Gene40"]][row,c("start","end","width","ranges")]$start) - 1, sep = ""))
        
        geneRegions[["Gene80"]][row,c("start","end","width","ranges")] <- c(as.numeric(geneRegions[["Gene60"]][row,c("start","end","width","ranges")]$start) - (as.numeric(df[row,"width"])*0.2),
                                                     as.numeric(geneRegions[["Gene60"]][row,c("start","end","width","ranges")]$start) - 1,
                                                     (as.numeric(geneRegions[["Gene60"]][row,c("start","end","width","ranges")]$start) - 1) - (as.numeric(geneRegions[["Gene60"]][row,c("start","end","width","ranges")]$start) - (as.numeric(df[row,"width"])*0.2)),
                                                     paste(as.numeric(geneRegions[["Gene60"]][row,c("start","end","width","ranges")]$start) - (as.numeric(df[row,"width"])*0.2), "-", as.numeric(geneRegions[["Gene60"]][row,c("start","end","width","ranges")]$start) - 1, sep = ""))
        
        geneRegions[["Gene100"]][row,c("start","end","width","ranges")] <- c(as.numeric(df[row,"start"]),
                                                      as.numeric(geneRegions[["Gene80"]][row,c("start","end","width","ranges")]$start) - 1,
                                                      (as.numeric(geneRegions[["Gene80"]][row,c("start","end","width","ranges")]$start) - 1) - as.numeric(df[row,"start"]),
                                                      paste(as.numeric(df[row,"start"]), "-", as.numeric(geneRegions[["Gene80"]][row,c("start","end","width","ranges")]$start) - 1, sep = ""))
      }
    }
  }
  return(geneRegions)
}