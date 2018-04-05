#########################################################################
#     rMSI - R package for MSI data handling and visualization
#     Copyright (C) 2014 Pere Rafols Soler
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.
############################################################################

#' export_imzML
#' 
#' Exports and rMSI object to an imzML formated file.
#'
#' @param img an rMSI object.
#' @param save_path the output path where the .imzML and .ibd files will be stored.
#'
#' @export
#'
export_imzML <- function(img, save_path)
{
  dir.create(path.expand(save_path), showWarnings = F, recursive = T)
  baseFileName <- unlist(strsplit( basename(img$name), "\\."))[1]
  
  #1- Create the offset matrix
  offMat <- list()
  offMat$UUID <- uuid_timebased()
  offMat$continuous_mode <- TRUE
  offMat$compression_mz <- FALSE
  offMat$compression_int <- FALSE
  offMat$mz_dataType <- "double"
  offMat$int_dataType <- ff::vmode(img$data[[1]])
  offMat$pixel_size_um <- img$pixel_size_um

  #Map ff data type to imzML encoding bytes
  if(offMat$int_dataType == "single")
  {
    EncodingByteSize <- 4
    offMat$int_dataType <- "float"
  }
  else if( offMat$int_dataType == "integer")
  {
    EncodingByteSize <- 4
    offMat$int_dataType <- "int"
  }
  else if( offMat$int_dataType == "double")
  {
    EncodingByteSize <- 8
    offMat$int_dataType <- "double"
  }
  else
  {
    stop("Error: The used ff::vmode is not supported by imzML\n");  
  }

  offMat$run_data <- data.frame( x =  img$pos[,"x"],
                                 y =  img$pos[,"y"],
                                 mzLength = rep(length(img$mass), nrow(img$pos)),
                                 mzOffset = rep(16, nrow(img$pos)),
                                 intLength = rep(length(img$mass), nrow(img$pos)),
                                 intOffset =  16+(length(img$mass)*8)+(EncodingByteSize*length(img$mass)*(0:(nrow(img$pos) - 1)))
                                 )
  

  #1- Store the binary data
  cat("Writing the .ibd file(binary)...\n")
  lapply(img$data, function(x){ ff::open.ff(x) })
  pb <- txtProgressBar(min = 0, max = nrow(img$pos), initial = 0, style = 3)
  iPixel <- 0
  ibd_fname <- file.path(path.expand(save_path), paste0(baseFileName, ".ibd"))
  ibd_conn <- file(ibd_fname, "w+b")
  intUUID <- strtoi(substring(offMat$UUID, seq(1,nchar(offMat$UUID),2), seq(2,nchar(offMat$UUID),2)), base = 16)
  writeBin(intUUID, ibd_conn, size = 1, endian="little") #Write the UUID
  writeBin(img$mass, ibd_conn, size = 8, endian="little", useBytes = T) #Store mass axis at offset 16 (mass aixis is always a double)
  for( i in 1:length(img$data))
  {
    dm <- img$data[[i]][,]
    for( j in 1:nrow(dm))
    {
      iPixel <- iPixel + 1
      setTxtProgressBar(pb, iPixel)
      writeBin(dm[j, ], ibd_conn, size = EncodingByteSize, endian="little", useBytes = T)
    }
    rm(dm)
  }
  close(pb)
  close(ibd_conn)
  
  cat("Calculating MD5 checksum...\n")
  offMat$SHA <- ""
  offMat$MD5 <- toupper(digest::digest( ibd_fname, algo = "md5", file = T))
  
  #2- Store the xml part
  cat("Writing the .imzML file(XML)...\n")
  imzML_fname <- file.path(path.expand(save_path), paste0(baseFileName, ".imzML"))
  bResult <- CimzMLStore( imzML_fname ,offMat)
  if(bResult)
  {
    cat("imzML exported successfully.\n")
  }
  else
  {
    cat("imzML exported with errors.\n")
  }
}
