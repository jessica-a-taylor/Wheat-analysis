# Get the coordinates of each gene region.
getGeneCoordinates <- function(df, genomeData) {
  source("Wheat-analysis/Functions/getCoordinatesFunctions.R")
  
  geneRegions <- list()
  
  # Determine the coordinates of the intergenic regions.
  for (region in c("UpstreamIntergenic", "DownstreamIntergenic")) {
    geneRegions[[region]] <- intergenicCoordinatesFunction(df, genomeData, region)
  }
  
  # Determine the coordinates of the promotor regions.
  for (region in c("Promotor500", "Promotor1000")) {
    geneRegions[[region]] <- promotorCoordinatesFunction(df, region)
  }
  
  # Determine the coordinates of the 200 bp downstream region.
  geneRegions[["Downstream"]] <- downstreamCoordinatesFunction(df)
  
  # Determine the coordinates of the gene body in 20% intervals.
  geneRegions <- genebodyCoordinatesFunction(df, geneRegions)
}
