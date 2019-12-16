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
# Find clusters on the edge of an alpha hull
#   Input: pos file and cluster stats file path
#   Cluster stats may be exported from IVAS (csv)
#   or from posgen*
#     *where the bounding box of the cluster is exported which requires
#		a patched version of posgen
#
#   AtomicDensity = Atomic density of material (atoms/nm3)
#   DetectionEfficiency = Detection efficiency of leap used
#   SamplingFraction = How much one would like to sample the POS file by
#   Returns a list of which cluster IDs are deemed to be on the edge of the data
findEdgeClusters <- function (posFileName, 
                              clusterStatsFile, 
                              AtomicDensity = 86, 
                              DetectionEfficiency = 0.37,
                              SamplingFraction = 0.005,
                              AlphaValue = 0) {
  
  library("tidyverse")
  library("geometry")
  library("alphashape3d")
  
  # get time for timing execution
  StartTime <- Sys.time()
  
  #cluster file path
  # Check for IVAS or posgen format
  cfHeader <- read.csv(clusterStatsFile, nrows=1)
  if (ncol(cfHeader)==2){
    ClusterFormat = "IVAS"
  }else{
    ClusterFormat = "posgen"
  }  
  
  if (ClusterFormat == "IVAS") {
    ClusterImport <- read.csv(clusterStatsFile, skip = 10)
    
    #### Cluster Import ####
    ClusterImport <-
      ClusterImport %>% select(
        cid = X,
        Solute.Ions,
        Ranged.Ions,
        Total.Ions,
        Center_x..nm..Ranged,
        Center_y..nm..Ranged,
        Center_z..nm..Ranged,
        Extent_x..nm..Ranged,
        Extent_y..nm..Ranged,
        Extent_z..nm..Ranged
      ) %>%
      mutate(
        x_max = Center_x..nm..Ranged + Extent_x..nm..Ranged,
        x_min = Center_x..nm..Ranged - Extent_x..nm..Ranged,
        y_max = Center_y..nm..Ranged + Extent_y..nm..Ranged,
        y_min = Center_y..nm..Ranged - Extent_y..nm..Ranged,
        z_max = Center_z..nm..Ranged + Extent_z..nm..Ranged,
        z_min = Center_z..nm..Ranged - Extent_z..nm..Ranged,
        xpos = Center_x..nm..Ranged,
        ypos = Center_y..nm..Ranged,
        zpos = Center_z..nm..Ranged
      ) %>%
      filter(grepl("Cluster", X))
    
  } else {
    ClusterImport <- read_delim(clusterStatsFile,
                                "\t", escape_double = FALSE, trim_ws = TRUE)
    
    ClusterImport <- ClusterImport %>%
      mutate(xpos=X, ypos=Y, zpos=Z, x_max = Xmax, x_min = Xmin, y_max = Ymax, y_min = Ymin, z_max = Zmax, z_min = Zmin)
    
    ClusterImport$cid <- paste("Cluster ",1:nrow(ClusterImport))
  }
  
  ClusterImportedData <- ClusterImport
  
  #### Parameters For Calculating Alpha Value####
  if (AlphaValue == 0) {
    AlphaValue <- 1 / (AtomicDensity * DetectionEfficiency * SamplingFraction)
  }
  
  
  #### Load filtered pos file to improve speed ####
  FilterPosFile <- read.pos.sampled(posFileName, SamplingFraction)
  
  print(paste0(
    "Alpha Value: ",
    AlphaValue,
    " Sampling fraction: ",
    SamplingFraction
  ))
  
  
  #### Find extreme co-ords for each cluster (modelling as cuboids) ####
  # FIXME: probably should be getting cluster id from "X"
  ClusterImport <- ClusterImport %>%
    mutate(
      ClusterID = cid,
      Point1 = paste0(x_max, ",", y_max, ",", z_max),
      Point2 = paste0(x_max, ",", y_max, ",", z_min),
      Point3 = paste0(x_max, ",", y_min, ",", z_max),
      Point4 = paste0(x_max, ",", y_min, ",", z_min),
      Point5 = paste0(x_min, ",", y_max, ",", z_max),
      Point6 = paste0(x_min, ",", y_max, ",", z_min),
      Point7 = paste0(x_min, ",", y_min, ",", z_max),
      Point8 = paste0(x_min, ",", y_min, ",", z_min)
    ) %>%
    select(ClusterID,
           Point1,
           Point2,
           Point3,
           Point4,
           Point5,
           Point6,
           Point7,
           Point8) %>%
    gather(Point, value,-ClusterID) %>%
    arrange(ClusterID) %>%
    separate(value, c("x", "y", "z"), sep = ",") %>%
    mutate(x = as.numeric(x),
           y = as.numeric(y),
           z = as.numeric(z))
  
  ClusterLocationMatrix <-
    as.matrix(ClusterImport %>% select(x, y, z))
  
  
  #### Using alpha-shape 3d ####
  M <- as.matrix(FilterPosFile[, 1:3])
  
  # free some memory
  rm(FilterPosFile)
  gc()
  
  ahullstart <- Sys.time()
  AlphaHullShape <- ashape3d(M, alpha =  AlphaValue)
  #plot(AlphaHullShape, indexAlpha = 1)
  print(Sys.time()-ahullstart)
  #### Ensuring there is only one connected volume ####
  #comp <- components_ashape3d(AlphaHullShape, indexAlpha = "all")
  #NumberVolumes <- as.data.frame(table(comp[[1]]))
  NumberVolumes <- data.frame(x=1) # assume one volume is created = much faster
  
  #### Determing which clusters are edge clusters and returning values ####
  if (nrow(NumberVolumes) != 1) {
    stop(paste0("Error. Multiple volumes created"))
  } else{
    print(paste0("Alpha value has created one volume"))
    
    findTime <- Sys.time()
    # determining which clusters are edge clusters
    ClustersInHull <-
      inashape3d(AlphaHullShape, indexAlpha = 1, ClusterLocationMatrix)
    print(Sys.time()-findTime)
    
    # getting names of clusters that are edge
    ClusterLocation2 <- cbind(ClusterImport, ClustersInHull) %>%
      dplyr::filter(ClustersInHull == "FALSE")
    
    TotalEdgeClusters <-
      as.data.frame(unique(ClusterLocation2$ClusterID))
    
    TotalEdgeClusters <-
      TotalEdgeClusters %>% transmute(ID = gsub("Cluster ", "", unique(ClusterLocation2$ClusterID)))
    
    TotalEdgeClusters <- sort(as.numeric(TotalEdgeClusters$ID))
  }
  
  
  print(paste0(
    "Toal number of clusters detected: ",
    length(TotalEdgeClusters)
  ))
  
  # Pritn total exe time
  EndTime <- Sys.time()
  TotalTime <- EndTime - StartTime
  print(TotalTime)
  
  return(list(TotalEdgeClusters=TotalEdgeClusters, ClusterImportedData=ClusterImportedData))
}