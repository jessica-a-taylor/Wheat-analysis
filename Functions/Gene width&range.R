# Function to get gene width and range.
getWidth <- function(dataToUse) {
  return(dataToUse$end - dataToUse$start)
}

getRange <- function(dataToUse) {
  return(paste(dataToUse$start,"-",dataToUse$end, sep = ""))
}
