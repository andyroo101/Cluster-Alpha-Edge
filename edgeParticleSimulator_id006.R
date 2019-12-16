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


# Edge particle simulation
StartTime <- Sys.time()
# Generate a set of non-colliding particles (possibly with a core-shell structure)
setwd(dir = "../posgen2")
source('../AlphaEdge/particlePlace.R')
source('../AlphaEdge/edgeParticleSimulator_xmlGen_all.R')
source('../AlphaEdge/findEdgeClusters.R')
source('../AlphaEdge/read.pos.sampled.R')

AtomicDensity <- 20
boxSize <- 40
runs <- 1000
expId <- "id006"

# Generate pos file using posgen xml script
#   Read template, write particles

summaryList <- vector("list", length = runs)

for (fileId in 1:runs) {
  print(fileId)
  particleCompositionCore  <- t(matrix(c(1, 1), nrow = 2))
  particleCompositionShell <- t(matrix(c(2, 1), nrow = 2))
  
  particles <- edgeParticleSim_coreShell(fileName = sprintf("%03d_gen.xml",fileId),
                                      posFileName = sprintf("%03d_full.pos",fileId),
                                      ionDensity = AtomicDensity, 
                                      numParticles = 40,
                                      boxSize = boxSize, 
                                      boxShape = "sphere",
                                      rCore = 1.587,
                                      rShell = 2,
                                      compShell = particleCompositionShell,
                                      compCore = particleCompositionCore,
                                      matrixComp = data.frame(mass=c(1,2,3),count=c(.1,.1,99.8)),
                                      clusterstatsFileName   = sprintf("%03d_stats.txt",fileId), 
                                      unclusterstatsFileName = sprintf("%03d_stats_matrix.txt",fileId), 
                                      sizedistFileName       = sprintf("%03d_sizeDist.txt",fileId), 
                                      clusteredposFileName   = sprintf("%03d_cluster.pos",fileId), 
                                      unclusteredposFileName = sprintf("%03d_matrix.pos",fileId), 
                                      clusteridposFileName   = sprintf("%03d_cluster.index.pos",fileId))
  #   Run posgen
  # Perform cluster selection as part of posgen
  # This generates the files listed in the input XML file
  posgenOut <- system(sprintf("posgen.exe %03d_gen.xml",fileId), wait = TRUE, intern=TRUE)
  
  # Find edge clusters
  edgeClustersAlpha <- findEdgeClusters(posFileName = sprintf("%03d_full.pos",fileId), 
                                        clusterStatsFile = sprintf("%03d_stats.txt",fileId),
                                        AtomicDensity = AtomicDensity, 
                                        DetectionEfficiency = 1, 
                                        SamplingFraction = 0.01,
                                        AlphaValue = 12)
  
  
  # Analyse results
  # Make data frame of: testNo, clusterNo, simSize, detected?, alphaEdge?, actualEdge, composition, measuredRg, soluteCount, matrixCount,
  # Need to analyse composion
  # Which clusters were identified, and which ones were edge clusters (make truth table)
  
  # match generated to detected clusters (by nearest centre distance)
  df1 <- select(particles,"x","y","z")
  df2 <- select(edgeClustersAlpha$ClusterImportedData, x=xpos, y=ypos, z=zpos)
  distances <- as.matrix(dist(bind_rows(df1, df2)))
  
  row.start <- nrow(df1)+1
  row.end <- nrow(df1) + nrow(df2)
  col.start <- 1
  col.end <- nrow(df1)
  # This \/ returns a list as long as the detected clusters, with which number it matches in the original list
  distanceIndex<-apply(distances[row.start:row.end, col.start:col.end], 1, which.min)
  # This is the distance between those points (should be within original cluster radius)
  #d<-apply(distances[row.start:row.end, col.start:col.end], 1, min)
  
  # Construct data frame of results
  summaryResults <- particles
  colnames(summaryResults) <- paste("sim", colnames(summaryResults), sep = "_")
  
  temp <- bind_cols(summaryResults[distanceIndex,], edgeClustersAlpha$ClusterImportedData)
  summaryResults <- bind_rows(temp, summaryResults[-distanceIndex,])
  
  # Define if they were edge or not based on the simulated positions
  # If any of (x,y,z) +/- size exceeds the upper or lower bounding box limit then they are edge clusters
  summaryResults$sim_edge <-  ((summaryResults$sim_R+summaryResults$sim_r)> boxSize/2)
  summaryResults$detected <- !is.na(summaryResults$Y)
  summaryResults$edge <- FALSE
  summaryResults$edge[edgeClustersAlpha$TotalEdgeClusters] <- TRUE
  
  summaryResults$ID <- fileId
  summaryList[[fileId]] <- summaryResults
  
  # clean up
  # copy generated file to separate folder "results/"
  if (!dir.exists(paste0("results/results_", expId))){
    # create results sub-directory
    dir.create(paste0("results/results_", expId))
  }
  
  file.copy(from = c(
    sprintf("%03d_stats.txt",fileId), 
    sprintf("%03d_sizeDist.txt",fileId), 
    sprintf("%03d_gen.xml",fileId)
  ),
  to = paste0("results/results_", expId),
  overwrite = TRUE)
  
  # delete files to free up space
  if (TRUE) {
    file.remove(
      sprintf("%03d_full.pos",fileId),
      sprintf("%03d_stats_matrix.txt",fileId), 
      sprintf("%03d_cluster.pos",fileId), 
      sprintf("%03d_matrix.pos",fileId), 
      sprintf("%03d_cluster.index.pos",fileId),
      sprintf("%03d_stats.txt",fileId), 
      sprintf("%03d_sizeDist.txt",fileId), 
      sprintf("%03d_gen.xml",fileId))
  }
  
}

EndTime <- Sys.time()
TotalTime <- EndTime - StartTime
print(TotalTime)

# save the data as an RData format
assign(expId, bind_rows(summaryList))
save("id006", file="id006_data")