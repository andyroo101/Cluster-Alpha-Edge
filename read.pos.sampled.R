# Copyright 2019 A. London, B. Jenkins
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
#
# Read a random sample of a pos file on opening
#   0 < sampleFraction < 1 = what fraction of data to read
#   Chunks of data are read at time, then sampled and stored in memory
#   and then returned
read.pos.sampled <- function(posFileName, sampleFraction){
  require(dplyr)
  
  # Check inputs
  if(!file.exists(posFileName)){
    stop("File does not exist:", posFileName)
  }
  
  if (sampleFraction<0 | sampleFraction>1){
    stop("Invalid sample fraction, should be >0 and <1")
  }
  
  # size of floats = 4 bytes
  sizeOfFloat <- 4
  # Number of fields (x,y,z,m, etc.)
  numOfFields <- 4
  # names of the fields
  fieldNames <- c("x","y","z","m")
  # Read in chunks at a time (in number of floats)
  chunkSize <- 32000*numOfFields
  # get file size
  fileInfo<-file.info(posFileName)
  fileSize<-fileInfo$size/sizeOfFloat
  
  # Sample number
  totalNumberToReturn <- sampleFraction*(fileSize/numOfFields)
  
  if (totalNumberToReturn < 1) {
    # return an empty data frame with a warning
    warning("Sampling resulted in no ions being returned")
    return(data.frame(x=single(), y=single(), z=single(), m=single()))
  }
  
  sampleChunkN <- ceiling(chunkSize*sampleFraction)
  
  to.read = file(posFileName, "rb")
  currentIndex <- 1
  i <- 1
  listData <- vector("list", ceiling(fileSize/chunkSize))
  while(currentIndex < fileSize) {
    # Read data into one big vector
    temp <- readBin(to.read,double(), size=sizeOfFloat, n=chunkSize, endian="big")
    # Sample that vector
    if ((currentIndex+chunkSize)<fileSize) {
      # Not at the end of the file
      # Sample ceil(chunkSize*sampleFraction) from the loaded data (sampleChunkN)
      temp <- temp[(sort(rep(sample(chunkSize/numOfFields,sampleChunkN), numOfFields))-1)*numOfFields+1 + rep(0:(numOfFields-1))]
    } else {
      # At the end of the file
      # temp length will be less than chunkSize
      sampleChunkN <- ceiling(length(temp)*sampleFraction/numOfFields)
      temp <- temp[(sort(rep(sample(length(temp)/numOfFields,sampleChunkN), numOfFields))-1)*numOfFields+1 + rep(0:(numOfFields-1))]
    }
    listData[[i]] <- as.data.frame(t(matrix(temp, nrow=numOfFields, dimnames=list(fieldNames))))
    
    i <- i + 1
    currentIndex <- currentIndex + chunkSize
  }
  close(to.read)
  
  posFile<-bind_rows(listData)
  
  return(posFile)
}

# Same as above, but returns a matrix, not a data.frame
read.pos.sampled.matrix <- function(posFileName, sampleFraction){
  require(dplyr)
  
  # Check inputs
  if(!file.exists(posFileName)){
    stop("File does not exist:", posFileName)
  }
  
  if (sampleFraction<0 | sampleFraction>1){
    stop("Invalid sample fraction, should be >0 and <1")
  }
  
  # size of floats = 4 bytes
  sizeOfFloat <- 4
  # Number of fields (x,y,z,m, etc.)
  numOfFields <- 4
  # names of the fields
  fieldNames <- c("x","y","z","m")
  # Read in chunks at a time (in number of floats)
  chunkSize <- 32000*numOfFields
  # get file size
  fileInfo<-file.info(posFileName)
  fileSize<-fileInfo$size/sizeOfFloat
  
  # Sample number
  totalNumberToReturn <- sampleFraction*(fileSize/numOfFields)
  
  if (totalNumberToReturn < 1) {
    # return an empty data frame with a warning
    warning("Sampling resulted in no ions being returned")
    return(data.frame(x=single(), y=single(), z=single(), m=single()))
  }
  
  sampleChunkN <- ceiling(chunkSize*sampleFraction)
  
  to.read = file(posFileName, "rb")
  currentIndex <- 1
  i <- 1
  listData <- vector("list", ceiling(fileSize/chunkSize))
  while(currentIndex < fileSize) {
    # Read data into one big vector
    temp <- readBin(to.read,double(), size=sizeOfFloat, n=chunkSize, endian="big")
    # Sample that vector
    if ((currentIndex+chunkSize)<fileSize) {
      # Not at the end of the file
      # Sample ceil(chunkSize*sampleFraction) from the loaded data (sampleChunkN)
      temp <- temp[(sort(rep(sample(chunkSize/numOfFields,sampleChunkN), numOfFields))-1)*numOfFields+1 + rep(0:(numOfFields-1))]
    } else {
      # At the end of the file
      # temp length will be less than chunkSize
      sampleChunkN <- ceiling(length(temp)*sampleFraction/numOfFields)
      temp <- temp[(sort(rep(sample(length(temp)/numOfFields,sampleChunkN), numOfFields))-1)*numOfFields+1 + rep(0:(numOfFields-1))]
    }
    listData[[i]] <- t(matrix(temp, nrow=numOfFields, dimnames=list(fieldNames)))
    
    i <- i + 1
    currentIndex <- currentIndex + chunkSize
  }
  close(to.read)
  
  posFile <- do.call(rbind,listData)
  
  return(posFile)
}