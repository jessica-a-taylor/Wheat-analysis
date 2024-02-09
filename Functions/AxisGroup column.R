# Function to add a column to the dataframe with the numbers on the x axis that will correspond with each gene region.
geneRegionAxisLocations <- function(region) {
  grouping <- c(seq(from = -60, to = -20, by = 20), seq(from = 20, to = 100, by = 20), seq(from = 140, to = 160, by = 20))
  axisGroup <- c("UpstreamIntergenic", "Promotor1000", "Promotor500", "Gene20", "Gene40", "Gene60",
                 "Gene80", "Gene100", "Downstream", "DownstreamIntergenic")
  
  axisGroupColumn <- grouping[which(axisGroup == region)]

  return(axisGroupColumn) 
}
