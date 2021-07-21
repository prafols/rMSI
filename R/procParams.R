#########################################################################
#     rMSI - R package for MSI data handling and visualization
#     Copyright (C) 2021 Pere Rafols Soler
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

###############################################################
###            DATA DESCRIPTION R5 OBJECTS                  ###
###############################################################

DataInfo <- setRefClass("DataInfo", 
                         fields = list(
                          version = "character",
                          raw_data_path = "data.frame",
                          out_path = "character",
                          roi_list = "list"
                          ),
                         
                         #Constructor
                         method = list(
                           initialize = function(..., version = as.character(packageVersion("rMSI")))
                           {
                             raw_data_path <<- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("imzML", "subimage_roi_xml")) 
                             callSuper(..., raw_data_path = raw_data_path,  version = version)
                           })
                        )
DataInfo$methods(
  
  #' getNumberOfImages.
  #'
  #' get the number of subimages listed in a DataInfo object.
  #' The parseROIs() method must be run befor in order to fill the subimage list.
  #' 
  #' @return the number of images.
  #' 
  #'
  getNumberOfImages = function()
  {
    if(length(roi_list) == 0)
    {
      stop("ERROR: empty roi_list, call parseROIs() first!")
    }
    
    return (length(roi_list))
  },
  
  
  #' getImgPathPos.
  #' 
  #' Get a single image information: name, path and roi positions.
  #'
  #' @param index image index.
  #'
  #' @return a list with name, path and roi position.
  #'
  getImgPathPos = function(index)
  {
    if(index < 1 || index > length(roi_list))
    {
      stop("ERROR: index out of roi_list bounds")
    }
    
    if(length(roi_list) == 0)
    {
      stop("ERROR: empty roi_list, call parseROIs() first!")
    }
    
    return(roi_list[[index]])
  },
  
  #' appendImzMLDataPath.
  #' 
  #' Append a imzML raw data path to the list of images to process. 
  #' A given imzML file (specified in the path argument) may contain multiple MS images (in example: for an acqusition with multiple tissue sections).
  #' The subimage_roi_xml XML file must contain the pixel positions of each pixel in each image ROI. Then, the imzML will be splitted to multiple images.
  #'
  #' @param path_imzML complete path to the data imzML file.
  #' @param subimage_roi_xml optional path to the XML file containing the ROIs to split the imzML file in multiple images.
  #'
  #'
  appendImzMLDataPath = function(path_imzML, subimage_roi_xml = NA)
  {
    aux_df <- setNames( data.frame( path_imzML, subimage_roi_xml), names(raw_data_path))
    raw_data_path <<- rbind(raw_data_path, aux_df)
  },
  
  #' parseROIs.
  #' 
  #' Parses the ROI information previously loaded in the raw_data_path data.frame with the appendImzMLDataPath() method.
  #' The results are sotred in the roi_lst list. Each item in the roi_lst corresponds to a single MS image to process.
  #'
  parseROIs = function()
  {
    roi_list <<- list() #Clear previous roi list
    for(i in nrow(raw_data_path))
    {
      cat(paste0("Parsing ROI info of imzML ", i ,  " of ", nrow(raw_data_path)))
      
      if(!is.na(raw_data_path$subimage_roi_xml[i]))
      {
        rois <- ParseBrukerXML(raw_data_path$subimage_roi_xml[i])  
        for( roi in rois)
        {
          roi_list[[length(roi_list)+1]] <<- list( name = paste0(unlist(strsplit(basename(raw_data_path$imzML[i]), split = "\\."))[1] , "_ROI_", roi$name ), 
                                                imzML = raw_data_path$imzML[i], 
                                                ROIpos = roi$pos)
        }
        
      }
      else
      {
        #No ROI file provided so just a single image
        roi_list[[length(roi_list)+1]] <<- list( name = unlist(strsplit(basename(raw_data_path$imzML[i]), split = "\\."))[1], 
                                               imzML = raw_data_path$imzML[i], 
                                               ROIpos = NULL) #NULL pos means all pixels in the imzML
      }
    }
  }
)

#TODO document me!
#' ImzMLDataDescription.
#' 
#' Create an empty DataInfo object to latter fill it.
#'
#' @return
#' @export
#'
#' @examples
ImzMLDataDescription <- function()
{
  data_info <- DataInfo()
  return (data_info)
}

###############################################################
###            PROCESSING PARAMETERS R5 OBJECTS             ###
###############################################################

SmoothingParams <- setRefClass("SmoothingParams", 
                               fields = list(
                                 enable = "logical",
                                 kernelSize = "integer"
                               ),
                               
                               #Constructor
                               method = list(
                                 initialize = function(...,
                                                       enable = T,
                                                       kernelSize = as.integer(7)
                                 )
                                 {
                                   callSuper(..., enable = enable, kernelSize = kernelSize)
                                 })
)

AlignmentParams <- setRefClass("AlignmentParams", 
                               fields = list(
                                 enable = "logical",
                                 bilinear = "logical",
                                 iterations = "integer",
                                 maxShiftppm = "numeric",
                                 refLow = "numeric",
                                 refMid = "numeric",
                                 refHigh = "numeric",
                                 overSampling = "integer",
                                 winSizeRelative = "numeric"
                               ),
                               
                               #Constructor
                               method = list(
                                 initialize = function(...,
                                                       enable = T,
                                                       bilinear = F,
                                                       iterations = as.integer(1),
                                                       maxShiftppm = 200,
                                                       refLow = 0.1,
                                                       refMid = 0.5,
                                                       refHigh = 0.9,
                                                       overSampling = as.integer(2),
                                                       winSizeRelative = 0.6
                                                       )
                                 {
                                   callSuper(..., enable = enable, bilinear = bilinear, iterations = iterations, maxShiftppm = maxShiftppm, 
                                             refLow = refLow, refMid = refMid, refHigh = refHigh, overSampling = overSampling, winSizeRelative = winSizeRelative)
                                 })
)



PreProcParams <- setRefClass("PreProcParams", 
                             fields = list(
                              merge = "logical", #TRUE to process multiple images using a common mass axis
                              smoothing = "SmoothingParams",
                              alignment = "AlignmentParams" 
                              #TODO compplete with all params, each stage its own class
                              )
                            )

AnnotationParams <- setRefClass("AnnotationParams", 
                                fields = list(
                                  ILSThreshold = "numeric"
                                  #TODO complete me!
                                  )
                                )

ProcParams <- setRefClass("ProcParams",
                          fields = list(
                           version = "character",
                           preprocessing = "PreProcParams",
                           annotations = "AnnotationParams",
                           outputpath = "character"
                           ),
                          
                          #Constructor
                          method = list(
                            initialize = function(..., version = as.character(packageVersion("rMSI")))
                                          {
                                            callSuper(..., version = version)
                                          })
                         )

#Set read-only fields
ref <- getRefClass("ProcParams")
ref$lock("version")

#Adding methods to the ProcParams class
ProcParams$methods(
                    getMergedProcessing = function()
                    {
                      return(preprocessing$merge)
                    },
                    
                    setMergedProcessing = function(bMerge)
                    {
                      preprocessing$merge <<- bMerge
                    },
                    
                    setOutputPath = function(datapath)
                    {
                      outputpath <<- path.expand(datapath)
                    }
                    
                   )

#TODO document me!
#' ProcessingParameters.
#' 
#' Create an empty ProcParams object to latter fill it.
#' 
#' @return
#' @export
#'
#' @examples
ProcessingParameters <- function()
{
  proc_params <- ProcParams()
  return (proc_params)
}
