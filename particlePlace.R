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
# place particles of different sizes without overlapping
# inputs:
# box size
# number of particles required
# particle value (radius)
# separationMargin extra space around each cluster
# type = cube or sphere
# returns data.frame of positions and radius
particlePlace <- function (boxSize = 20, numParticles = 20, particleValue = 2, separationMargin = 1, boxShape="cube") {
  
  require(dplyr)
  
  # function which defines what the particle size distribution will be
  # it takes the number of sizes required n, and returns an array of sizes
  sizeDist <- function(n=20, v=2){
    # mono-disperse with size v
    return(vector(mode = "numeric", length = n)+v)
  }
  
  particleRadiusDist <- sizeDist(numParticles, particleValue)
  particleMaxRadius <- max(particleRadiusDist)
  
  if (boxShape=="cube") {
    # make a check for the volume of the box and the particles requested (to avoid very slow computation)
    # box volume
    boxVolume <- boxSize^3
    particleVolume <- 4*pi*sum(particleRadiusDist^3)/3
    if (particleVolume/boxVolume > 0.4) {
      warning(sprintf("Particles occupy %2.0f%% of the box volume",100*particleVolume/boxVolume))
      if (particleVolume/boxVolume > 1) {
        stop("Particle volume greater than 100%")
      }
    }
    
    
    particles <- data.frame(x = runif(numParticles, min = -boxSize/2 - particleMaxRadius, max = boxSize/2 + particleMaxRadius),
                            y = runif(numParticles, min = -boxSize/2 - particleMaxRadius, max = boxSize/2 + particleMaxRadius),
                            z = runif(numParticles, min = -boxSize/2 - particleMaxRadius, max = boxSize/2 + particleMaxRadius),
                            r = particleRadiusDist)
    # try and place each particle, check for collisions
    for (i in 1:numParticles){
      okToPlace <- FALSE
      while (!okToPlace) {
        # make new point
        p <- particles[i,1:3]
        # check is placement is ok
        d <- sqrt(rowSums(sweep(as.matrix(particles[,1:3]),2,as.numeric(p))^2)) < (particles$r + particles$r[i] + separationMargin^2)
        # exclude self matching
        d[i] <- FALSE
        if (sum(d)==0){
          okToPlace <- TRUE
        } else {
          # make a new entry for this row
          particles[i,1:3] <- runif(3, min = -boxSize/2 - particleMaxRadius, 
                                    max =  boxSize/2 + particleMaxRadius)
        }
      }
    }
    return(particles)
    
  } else if (boxShape == "sphere") {
    # make a check for the volume of the box and the particles requested (to avoid very slow computation)
    # box volume
    boxVolume <- (4*pi/3)*(boxSize/2)^3
    particleVolume <- 4*pi*sum(particleRadiusDist^3)/3
    if (particleVolume/boxVolume > 0.4) {
      warning(sprintf("Particles occupy %2.0f%% of the box volume",100*particleVolume/boxVolume))
      if (particleVolume/boxVolume > 1) {
        stop("Particle volume greater than 100%")
      }
    }
    
    particles <- data.frame(phi = runif(numParticles, min = 0, max = 2*pi),
                            costheta = runif(numParticles, min = -1, max = 1),
                            u = runif(numParticles, min = 0, max = 1),
                            r = particleRadiusDist)
    
    particles <- particles %>%
      mutate(theta = acos(costheta),
             R = (boxSize/2 + particleMaxRadius) * (u^(1/3)),
             x = R * sin( theta) * cos( phi ),
             y = R * sin( theta) * sin( phi ),
             z = R * cos( theta )) %>%
      select(x,y,z,r,R)

    # try and place each particle, check for collisions
    for (i in 1:numParticles){
      okToPlace <- FALSE
      while (!okToPlace) {
        # make new point
        p <- particles[i,1:3]
        # check is placement is ok
        d <- sqrt(rowSums(sweep(as.matrix(particles[,1:3]),2,as.numeric(p))^2)) < (particles$r + particles$r[i] + separationMargin^2)
        # exclude self matching
        d[i] <- FALSE
        if (sum(d)==0){
          okToPlace <- TRUE
        } else {
          # make a new entry for this row
          temp <- data.frame(phi = runif(1, min = 0, max = 2*pi),
                             costheta = runif(1, min = -1, max = 1),
                             u = runif(1, min = 0, max = 1),
                             r = sizeDist(1, particleValue))
          
          particles[i,1:5] <- temp %>%
            mutate(theta = acos(costheta),
                   R = (boxSize/2 + particleMaxRadius) * (u^(1/3)),
                   x = R * sin( theta) * cos( phi ),
                   y = R * sin( theta) * sin( phi ),
                   z = R * cos( theta )) %>%
            select(x,y,z,r,R)
          
        }
      }
    }
    return(particles)
  }
}