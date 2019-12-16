# Cluster Alpha Edge
Edge detection of clusters for point-cloud data using an alpha-hull. Scripts written in R are used for the simulation of test data (see `edgeParticleSimulator.R` ) and finding the edge clusters from the associated cluster pos file results (see `findEdgeClusters.R`). The "full" pos file, cluster data or statistics (including the data required to get the cluster bounding box) as well as some settings (atomic density, detection efficiency and sampling fraction).

# Required R libraries
```
  library("tidyverse")
  library("geometry")
  library("alphashape3d")
```

# Installation
Copy or extract the contents of the repository into a local folder called “AlphaEdge”. The experiment scripts, detailed below, assume the [posgen exe](https://apttools.sourceforge.net) is available in one folder above the R scripts directory:
```
folder/
      ├── ClusterAlphaEdge
      ├── posgen
```
The `findEdgeClustersConvex.R` script can use an unmodified version of posgen(https://apttools.sourceforge.net), however the `findEdgeClusters.R` script requires the bounding box of each cluster in the cluster statistics file produced by the cluster search step. A source code patch is available to add this functionality to posgen (`posgenBoundingBoxPatch.txt`). This will require the program to be (compiled from scratch)[http://apttools.sourceforge.net/compile-posgen.html] after the patch has been applied.
Then you should be able to run any of the `edgeParticleSimulator` scripts which generate a simulated dataset, cluster search and then find edge clusters. Or you can use the `findEdgeClustersConvex` function on your own data to identify the edge clusters. See the `findEdgeClustersConvex` file for usage.

 
# Experiments
Different simulation experiments are identified by an id, a description of each can be found in the script ending in the id number, also below for reference:
## id002
`edgeParticleSimulator_coreShell_id002.R`

1000 simulations.

20 core/shell particles in a 20 nm3 box, 1.58 core radius, 2 nm shell radius (4 nm total diameter) composition as follows:

| Part | H    | C       | Fe    |
|------|-----|--------|-----|
| Core  | 99.9| 0.01  | 0  |
| Shell | 0.01  | 99.9 | 0  |
| Matrix  | 0.01  | 0.01  | 99.8  |

## id004
`edgeParticleSimulator_id004.R`

## id007
`edgeParticleSimulator_id007.R`

## id008
`edgeParticleSimulator_id008.R`


# Authors: 
Andrew London
Ben Jenkins


# License 
Copyright 2019 A. London, B. Jenkins
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
