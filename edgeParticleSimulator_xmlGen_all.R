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


# Generate an XML script for posgen that will create a 
# number of random or core-shell particles, and then cropped by a cuboid or sphere
# edgeParticleSim_sphere creates clusters with no structure
# edgeParticleSim_coreshell creates clusters with a core-shell structure
edgeParticleSim_sphere <- function(fileName, 
                                   posFileName, 
                                   ionDensity, 
                                   boxSize, 
                                   boxShape = "cube",
                                   numParticles, 
                                   particleValue, 
                                   particleComp,
                                   matrixComp = data.frame(mass=c(1,2,3),count=c(1,1,98)),
                                   clusterstatsFileName = "cluster-stats.txt", 
                                   unclusterstatsFileName = "uncluster-stats.txt",
                                   sizedistFileName = "sizeDist.txt", 
                                   clusteredposFileName = "cluster.pos", 
                                   unclusteredposFileName = "matrix.pos", 
                                   clusteridposFileName = "cluster.indexed.pos",
                                   dmaxSettings = list(dclassify = "0.0",
                                                       knn = "1",
                                                       dmax = "0.5", 
                                                       dbulk = "0.2",
                                                       derode = "0.2")){
  # ** Range file and core ions are hard coded at the moment **
  # example inputs
  # fileName <- "sim-edge-particles2.xml"
  # posFileName <- "test2.pos"
  # ionDensity <- 20
  # boxSize <- 20 (centered on the origin)
  # numParticles <- 20 number of particles to generate
  # particleValue the value passed to the size generator in placeParticles
  # particleComp <- t(matrix(c(1, 0.5, 2, 0.5), nrow = 2))
  #   50% mass value of 1 and 50% mass value of 2
  
  particles <- particlePlace(boxSize=boxSize,
                             boxShape=boxShape, 
                             numParticles = numParticles, 
                             particleValue = particleValue)
  
  fileParts <- vector(length = nrow(particles)+2, mode= "list")
  
  fileParts[[1]] <- 
    c('<!DOCTYPE posscript SYSTEM "posscript.dtd">',
      '<posscript>',
      '	<version value="0.0.1"/>',
      '',
      '	<!-- store the matrix data -->',
      ' <spatrand>',
      sprintf(' \t<bound x="%f" y="%f" z="%f"/>',boxSize,boxSize,boxSize),
      sprintf(' \t<density value="%2.3f"/>', ionDensity),
      '',
      sprintf('		<atom index="(000)" mass="%f"/>\n		<count value="%f"/>\n', matrixComp$mass, matrixComp$count),
      '',
      ' </spatrand>',
      '',
      sprintf(' <geotransform> <translate><x value="%2.5f"/><y value="%2.5f"/><z value="%2.5f"/></translate></geotransform>',
              -boxSize/2,
              -boxSize/2,
              -boxSize/2),
      '')
  
  if (boxShape == "sphere") {
    # clip the cube into a sphere
    temp <- sprintf('	<clip invert="false">  <sphere radius="%2.5f"> <x value="%2.5f"/> <y value="%2.5f"/> <z value="%2.5f"/> </sphere>  </clip>',
                    boxSize/2,
                    0,
                    0,
                    0)
    fileParts[[1]] <- c(fileParts[[1]], 
                        '	<!-- clip matrix data -->',
                        '',
                        temp, 
                        '')
  }
  
  # clip out the matrix from the particle regions
  for (i in 1:nrow(particles)) {
    temp <- sprintf('	<clip invert="true">  <sphere radius="%2.5f"> <x value="%2.5f"/> <y value="%2.5f"/> <z value="%2.5f"/> </sphere>  </clip>',
                    particles$r[i],
                    particles$x[i],
                    particles$y[i],
                    particles$z[i])
    fileParts[[1]] <- c(fileParts[[1]], temp)
  }
  # add the store command at the end, to store the matrix data
  fileParts[[1]] <- c(fileParts[[1]], c('	<store name="running_dataset" forceram="true" clear="true"/>',
                                        '',
                                        ''))
  
  # place particles here
  for (i in 1:nrow(particles)) {
    fileParts[[i+1]] <- c(
      sprintf('\t<!-- START PPT  %d -->',i),
      particleGenSphere(cx=particles$x[i],cy = particles$y[i], cz = particles$z[i], r = particles$r[i], density = ionDensity, composition = particleComp),
      '\t        <recall name="running_dataset"/>',
      '\t        <store name="running_dataset" clear="true" forceram="true"/>',
      sprintf('\t<!-- END PPT  %d -->',i),
      '')
  }
  
  
  tempStr <- c('	<!-- restore the running dataset -->',
               '	<recall name="running_dataset"/>',
               '')
  
  if (boxShape == "sphere") {
    tempStr <- c(tempStr, '\t<!-- clip to the bounding box (sphere) -->',
                 '\t<clip>',
                 sprintf('\t\t<sphere radius="%f">', boxSize/2),
                 sprintf('\t\t\t<%s value="%f"/>', c("x","y","z"), rep(0,3)),
                 '\t\t</sphere>',
                 '\t</clip>')
  } else {
    tempStr <- c(tempStr, '\t<!-- clip to the bounding box -->',
                 '\t<clip>',
                 '\t\t<box>',
                 sprintf('\t\t\t<%s value="%f"/>', c("x","y","z"), rep(-boxSize/2,3)),
                 sprintf('\t\t\t<%s value="%f"/>', c("x","y","z"), rep(boxSize/2,3)),
                 '\t\t</box>',
                 '\t</clip>')
  }
  
  tempStr <- c(tempStr, '',
               sprintf('\t<possave file="%s"/>', posFileName),
               '',
               '<!-- CLUSTER ANALYSIS -->',
               '	<cluster>',
               '		<algorithm value="maxsep">',
               sprintf('			<dclassify value="%s" knn="%s"/> <!-- Coring (pre-cluster) dist; zero to disable this step-->', dmaxSettings$dclassify, dmaxSettings$knn),
               sprintf('			<dmax value="%s"/> <!-- Max sep dist -->', dmaxSettings$dmax),
               sprintf('			<dbulk value="%s"/> <!-- AKA envelope-->', dmaxSettings$dbulk),
               sprintf('			<derode value="%s"/> <!-- erosion distance-->', dmaxSettings$derode),
               '		</algorithm>',
               '		<range file="simple.RRNG"/>',
               '		<!--select ions for cluster core (aka solute)-->',
               '		<core>',
               '			<typelist>',
               '				<atomtype symbol="H"/>',
               '				<atomtype symbol="C"/>',
               '			</typelist>',
               '		</core>',
               '		<!--select ions for bulk (aka matrix)-->',
               '		<bulk>',
               '			<typelist>',
               '				<atomtype symbol="Fe"/> ',
               '			</typelist>',
               '		</bulk>',
               '',
               '		<sizeclip nmin="20"/> ',
               '		<unranged foroutput="true" forstats="true"/>',
               '',
               sprintf('		<clusterstats core="true" bulk="true" percluster="true" file="%s"/>', clusterstatsFileName),
               sprintf('		<!--unclusterstats file="%s"/-->',unclusterstatsFileName),
               '',
               sprintf('		<sizedist file="%s"/>', sizedistFileName),
               sprintf('		<!--clustered-pos file="%s" retain="true"/-->', clusteredposFileName),
               sprintf('		<!--unclustered-pos file="%s" retain="true"/-->', unclusteredposFileName),
               '',
               sprintf('		<clusterid file="%s" offset="1"/>', clusteridposFileName),
               '		',
               '	</cluster>',
               '</posscript>')
  
  fileParts[[length(fileParts)]] <- tempStr
  
  
  # write XML file
  writeLines(unlist(fileParts),fileName)
  
  # return the generated particle positions and sizes
  return(particles)
}


edgeParticleSim_coreShell <- function(fileName, 
                                      posFileName, 
                                      ionDensity, 
                                      boxSize, 
                                      boxShape = "cube",
                                      numParticles, 
                                      rCore, 
                                      rShell, 
                                      compShell,
                                      compCore,
                                      matrixComp = data.frame(mass=c(1,2,3),count=c(1,1,98)),
                                      clusterstatsFileName = "cluster-stats.txt", 
                                      unclusterstatsFileName = "uncluster-stats.txt",
                                      sizedistFileName = "sizeDist.txt", 
                                      clusteredposFileName = "cluster.pos", 
                                      unclusteredposFileName = "matrix.pos", 
                                      clusteridposFileName = "cluster.indexed.pos",
                                      dmaxSettings = list(dclassify = "0.0",
                                                          knn = "1",
                                                          dmax = "0.5", 
                                                          dbulk = "0.2",
                                                          derode = "0.2")){
  # ** Range file and core ions are hard coded at the moment **
  # example inputs
  # fileName <- "sim-edge-particles2.xml"
  # posFileName <- "test2.pos"
  # ionDensity <- 20
  # boxSize <- 20 (centered on the origin)
  # numParticles <- 20 number of particles to generate
  # particleValue the value passed to the size generator in placeParticles
  # particleComp <- t(matrix(c(1, 0.5, 2, 0.5), nrow = 2))
  #   50% mass value of 1 and 50% mass value of 2
  
  particles <- particlePlace(boxSize=boxSize, boxShape = boxShape, numParticles = numParticles, particleValue = rShell)
  
  fileParts <- vector(length = nrow(particles)+2, mode= "list")
  
  fileParts[[1]] <- 
    c('<!DOCTYPE posscript SYSTEM "posscript.dtd">',
      '<posscript>',
      '	<version value="0.0.1"/>',
      '',
      '	<!-- store the matrix data -->',
      ' <spatrand>',
      sprintf(' \t<bound x="%f" y="%f" z="%f"/>',boxSize,boxSize,boxSize),
      sprintf(' \t<density value="%2.3f"/>', ionDensity),
      '',
      sprintf('		<atom index="(000)" mass="%f"/>\n		<count value="%f"/>\n', matrixComp$mass, matrixComp$count),
      '',
      ' </spatrand>',
      '',
      sprintf(' <geotransform> <translate><x value="%2.5f"/><y value="%2.5f"/><z value="%2.5f"/></translate></geotransform>',
              -boxSize/2,
              -boxSize/2,
              -boxSize/2),
      '')
  
  if (boxShape == "sphere") {
    # clip the cube into a sphere
    temp <- sprintf('	<clip invert="false">  <sphere radius="%2.5f"> <x value="%2.5f"/> <y value="%2.5f"/> <z value="%2.5f"/> </sphere>  </clip>',
                    boxSize/2,
                    0,
                    0,
                    0)
    fileParts[[1]] <- c(fileParts[[1]], 
                        '	<!-- clip matrix data -->',
                        '',
                        temp, 
                        '')
  }
  
  # clip out the matrix from the particle regions
  for (i in 1:nrow(particles)) {
    temp <- sprintf('	<clip invert="true">  <sphere radius="%2.5f"> <x value="%2.5f"/> <y value="%2.5f"/> <z value="%2.5f"/> </sphere>  </clip>',
                    particles$r[i],
                    particles$x[i],
                    particles$y[i],
                    particles$z[i])
    fileParts[[1]] <- c(fileParts[[1]], temp)
  }
  # add the store command at the end, to store the matrix data
  fileParts[[1]] <- c(fileParts[[1]], c('	<store name="running_dataset" forceram="true" clear="true"/>',
                                        '',
                                        ''))
  
  # place particles here
  for (i in 1:nrow(particles)) {
    fileParts[[i+1]] <- c(
      sprintf('\t<!-- START PPT  %d -->',i),
      particleGenCoreShell(cx=particles$x[i],
                           cy = particles$y[i], 
                           cz = particles$z[i], 
                           rCore = rCore,
                           rShell = particles$r[i], 
                           density = ionDensity, 
                           compositionCore = compCore,
                           compositionShell = compShell),
      '\t        <recall name="running_dataset"/>',
      '\t        <store name="running_dataset" clear="true" forceram="true"/>',
      sprintf('\t<!-- END PPT  %d -->',i),
      '')
  }
  
  fileParts[[length(fileParts)]] <- c('	<!-- restore the running dataset -->',
                                      '	<recall name="running_dataset"/>',
                                      '')
  
  
  if (boxShape == "sphere") {
    tempStr <- c('\t<!-- clip to the bounding box (sphere) -->',
                 '\t<clip>',
                 sprintf('\t\t<sphere radius="%f">', boxSize/2),
                 sprintf('\t\t\t<%s value="%f"/>', c("x","y","z"), rep(0,3)),
                 '\t\t</sphere>',
                 '\t</clip>')
  } else {
    tempStr <- c('\t<!-- clip to the bounding box -->',
                 '\t<clip>',
                 '\t\t<box>',
                 sprintf('\t\t\t<%s value="%f"/>', c("x","y","z"), rep(-boxSize/2,3)),
                 sprintf('\t\t\t<%s value="%f"/>', c("x","y","z"), rep(boxSize/2,3)),
                 '\t\t</box>',
                 '\t</clip>')
  }
  
  fileParts <- c(fileParts, tempStr, sprintf('\t<possave file="%s"/>', posFileName),
                 '',
                 '<!-- CLUSTER ANALYSIS -->',
                 '	<cluster>',
                 '		<algorithm value="maxsep">',
                 sprintf('			<dclassify value="%s" knn="%s"/> <!-- Coring (pre-cluster) dist; zero to disable this step-->', dmaxSettings$dclassify, dmaxSettings$knn),
                 sprintf('			<dmax value="%s"/> <!-- Max sep dist -->', dmaxSettings$dmax),
                 sprintf('			<dbulk value="%s"/> <!-- AKA envelope-->', dmaxSettings$dbulk),
                 sprintf('			<derode value="%s"/> <!-- erosion distance-->', dmaxSettings$derode),
                 '		</algorithm>',
                 '		<range file="simple.RRNG"/>',
                 '		<!--select ions for cluster core (aka solute)-->',
                 '		<core>',
                 '			<typelist>',
                 '				<atomtype symbol="H"/>',
                 '				<atomtype symbol="C"/>',
                 '			</typelist>',
                 '		</core>',
                 '		<!--select ions for bulk (aka matrix)-->',
                 '		<bulk>',
                 '			<typelist>',
                 '				<atomtype symbol="Fe"/> ',
                 '			</typelist>',
                 '		</bulk>',
                 '',
                 '		<sizeclip nmin="20"/> ',
                 '		<unranged foroutput="true" forstats="true"/>',
                 '',
                 sprintf('		<clusterstats core="true" bulk="true" percluster="true" file="%s"/>', clusterstatsFileName),
                 sprintf('		<!--unclusterstats file="%s"/-->',unclusterstatsFileName),
                 '',
                 sprintf('		<sizedist file="%s"/>', sizedistFileName),
                 sprintf('		<!--clustered-pos file="%s" retain="true"/-->', clusteredposFileName),
                 sprintf('		<!--unclustered-pos file="%s" retain="true"/-->', unclusteredposFileName),
                 '',
                 sprintf('		<clusterid file="%s" offset="1"/>', clusteridposFileName),
                 '		',
                 '	</cluster>',
                 '</posscript>')
  
  
  # write XML file
  writeLines(unlist(fileParts),fileName)
  
  # return the generated particle positions and sizes
  return(particles)
}

## function to place a sphere of specific size
particleGenSphere <- function(cx, cy, cz, r, density, composition){
  # cx, cy, cz are centre coordinates
  # r = radius
  # density = ions/volume
  # composition is the mass and count value as a matrix, e.g.
  # t(matrix(c(1, 0.5, 2, 0.5), nrow = 2)) would be half mass 1
  # and half mass 2, posgen normalises the count values
  
  # make the atom indicies
  compList <- sprintf('		<atom index="(000)" mass="%f"/>\n		<count value="%f"/>',composition[,1], composition[,2])

  # combine the spatrand, translate, clip and translate together
  return(c('	<spatrand>',              # generate random box
           sprintf('    <bound x="%f" y="%f" z="%f"/>',r*2,r*2,r*2),
           sprintf('		<density value="%2.3f"/>', density),
           '',
           compList,
           '	</spatrand> ',
           '	<geotransform>',
           '		<translate>',           # centre it on the origin
           sprintf('			<x value="%f"/>',-r),
           sprintf('			<y value="%f"/>',-r),
           sprintf('			<z value="%f"/>',-r),
           '		</translate>',
           '	</geotransform>',
           '	<clip>',                  # clip out a sphere
           sprintf('		<sphere radius="%f">',r),
           '			<x value="0"/>',
           '			<y value="0"/>',
           '			<z value="0"/>',
           '		</sphere>',
           '	</clip>',
           '	<geotransform>',
           '		<translate>',           # translate it to the required point
           sprintf('			<x value="%f"/>',cx),
           sprintf('			<y value="%f"/>',cy),
           sprintf('			<z value="%f"/>',cz),
           '		</translate>',
           '	</geotransform>'))
}

## function to place a core-shell sphere of specific size (shell and core same density)
particleGenCoreShell <- function(cx, cy, cz, rCore, rShell, density, compositionCore, compositionShell) {
  # cx, cy, cz are centre coordinates
  # rCore = radius of the core
  # rShell = outer radius of the shell
  #   therefore, shell thickness = rShell-rCore
  # density = ions/volume
  # composition is the mass and count value as a matrix, e.g.
  #   t(matrix(c(1, 0.5, 2, 0.5), nrow = 2)) would be half mass 1
  #   and half mass 2, posgen normalises the count values
  #   The composition of the core and shell are specified separately
  
  # make the atom indicies of the core
  compListCore <- sprintf('		<atom index="(000)" mass="%f"/>\n		<count value="%f"/>',
                          compositionCore[,1], compositionCore[,2])
  
  # make the atom indicies of the shell
  compListShell <- sprintf('		<atom index="(000)" mass="%f"/>\n		<count value="%f"/>',
                           compositionShell[,1], compositionShell[,2])
  
  # combine the spatrand, translate, clip and translate together
  return(c('	<spatrand>',              # generate random box
           sprintf('    <bound x="%f" y="%f" z="%f"/>',rCore*2,rCore*2,rCore*2),
           sprintf('		<density value="%2.3f"/>', density),
           '',
           compListCore,
           '	</spatrand> ',
           '	<geotransform>',
           '		<translate>',           # centre it on the origin
           sprintf('			<x value="%f"/>',-rCore),
           sprintf('			<y value="%f"/>',-rCore),
           sprintf('			<z value="%f"/>',-rCore),
           '		</translate>',
           '	</geotransform>',
           '	<clip>',                  # clip out a sphere
           sprintf('		<sphere radius="%f">',rCore),
           '			<x value="0"/>',
           '			<y value="0"/>',
           '			<z value="0"/>',
           '		</sphere>',
           '	</clip>',
           '<store name="core" clear="true" forceram="true"/>',
           '',
           '<!-- Generate shell atoms -->',
           '	<spatrand>',              # generate random box
           sprintf('    <bound x="%f" y="%f" z="%f"/>',rShell*2,rShell*2,rShell*2),
           sprintf('		<density value="%2.3f"/>', density),
           '',
           compListShell,
           '	</spatrand> ',
           '	<geotransform>',
           '		<translate>',           # centre it on the origin
           sprintf('			<x value="%f"/>',-rShell),
           sprintf('			<y value="%f"/>',-rShell),
           sprintf('			<z value="%f"/>',-rShell),
           '		</translate>',
           '	</geotransform>',
           '	<clip invert="true">',   # clip out core
           sprintf('		<sphere radius="%f">',rCore),
           '			<x value="0"/>',
           '			<y value="0"/>',
           '			<z value="0"/>',
           '		</sphere>',
           '	</clip>',
           '	<clip>',                  # clip outer
           sprintf('		<sphere radius="%f">',rShell),
           '			<x value="0"/>',
           '			<y value="0"/>',
           '			<z value="0"/>',
           '		</sphere>',
           '	</clip>',
           '',
           '<!-- merge it with its core -->',
           '<recall name="core"/>',
           '',
           '<!-- Translate particle -->',
           '	<geotransform>',
           '		<translate>',           # translate it to the required point
           sprintf('			<x value="%f"/>',cx),
           sprintf('			<y value="%f"/>',cy),
           sprintf('			<z value="%f"/>',cz),
           '		</translate>',
           '	</geotransform>'))
  
}