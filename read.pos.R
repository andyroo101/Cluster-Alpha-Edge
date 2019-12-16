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
# Read a binary 'pos' of (x,y,z,m) data for atom probe in to a data frame
# Supply path to file as posFileName
# This function is not very memory efficient and cannot be used on large files
read.pos <- function(posFileName){
  # size of floats = 4 bytes
  sizeOfFloat = 4
  # get file size
  if (!file.exists(posFileName)) {
    stop(paste0("Pos file:'", posFileName,"' does not exist"))
  }
  fileInfo<-file.info(posFileName)
  fileSize<-fileInfo$size/sizeOfFloat
  to.read = file(posFileName, "rb")
  posFile<-readBin(to.read,double(),size=4,n=fileSize,endian="big")
  close(to.read)
  posFile<-t(matrix(posFile,nrow=4,dimnames=list(c("x","y","z","m"))))
  posFile<-as.data.frame(posFile)
  return(posFile)
}