% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/libimgramdisk.R
\name{insertRasterImageAtCols}
\alias{insertRasterImageAtCols}
\title{Inserts a image at specified Cols of a rMSI object.}
\usage{
insertRasterImageAtCols(Img, Cols, raster_matrix)
}
\arguments{
\item{Img}{the rMSI object where the data is stored (ramdisk).}

\item{Cols}{the columns indexes from which data will be inserted}

\item{raster_matrix}{a raster image represented as a matrix with pixel values.}
}
\description{
A raster image provided as a matrix is inserted at given Cols with a gaussian shape.
The raster_matrix has nrows as X direction.
}
