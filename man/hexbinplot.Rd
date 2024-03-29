\name{Formula methods}
\alias{hexbinplot}
\alias{hexbinplot,formula,Raster-method}
\alias{xyplot,formula,Raster-method}
\alias{hexbinplot,formula,SpatRaster-method}
\alias{xyplot,formula,SpatRaster-method}
\title{Formula methods}
\description{Formula methods
}
\usage{
\S4method{xyplot}{formula,Raster}(x, data, dirXY, maxpixels=1e5,
            alpha=0.05,
            xscale.components=xscale.raster, yscale.components=yscale.raster,
            par.settings=rasterTheme(),...)
\S4method{hexbinplot}{formula,Raster}(x, data, dirXY, maxpixels=1e6,
            xscale.components=xscale.raster, yscale.components=yscale.raster,
            par.settings=rasterTheme(),...)
\S4method{xyplot}{formula,SpatRaster}(x, data, dirXY, maxpixels=1e5,
            alpha=0.05,
            xscale.components=xscale.raster, yscale.components=yscale.raster,
            par.settings=rasterTheme(),...)
\S4method{hexbinplot}{formula,SpatRaster}(x, data, dirXY, maxpixels=1e6,
            xscale.components=xscale.raster, yscale.components=yscale.raster,
            par.settings=rasterTheme(),...)
}
\arguments{
\item{x}{A \code{formula} describing the variables to be related. It may
  include the layer names (which are internally converted to valid ones
  with \code{\link{make.names}}) and the \code{x}, \code{y} variables
  representing the coordinates of the Raster object. Besides, if dirXY
  is not missing, the variable \code{dirXY} can also be included in the
  \code{formula}.
  }
\item{data}{A \code{Raster} object or a \code{SpatRaster} object.}
\item{dirXY}{A direction as a function of the coordinates (see
  \code{\link{xyLayer}}).}
\item{maxpixels}{A numeric, for \code{\link{sampleRegular}} or \code{spatSample}.}
\item{alpha}{A numeric, transparency of the points.}
\item{xscale.components, yscale.components,par.settings}{Customization of \code{lattice}. See \code{xyplot}
for details.}
\item{\dots}{Additional arguments for the \code{xyplot} and
  \code{\link{hexbinplot}} functions.}
}

\author{Oscar Perpiñán Lamigueiro}

\examples{
library(raster)
library(terra)

f <- system.file("external/test.grd", package="raster")
r <- raster(f)
names(r)

xyplot(test~y, data=r, alpha=0.5)

\dontrun{
dataURL <- "https://raw.github.com/oscarperpinan/bookvis/master/data/"

##Solar irradiation data from CMSAF http://dx.doi.org/10.5676/EUM_SAF_CM/RAD_MVIRI/V001
old <- setwd(tempdir())
download.file(paste0(dataURL, "SISmm2008_CMSAF.zip"),
   "SISmm2008_CMSAF.zip", method='wget')
unzip("SISmm2008_CMSAF.zip")

listFich <- dir(pattern='\\\\.nc')
stackSIS <- stack(listFich)
SISmm <- stackSIS*24 ##from irradiance (W/m2) to irradiation Wh/m2
setwd(old)

names(SISmm) <- month.abb

##Relation between the January & February versus July radiation for four
##differents longitude regions.
xyplot(Jan+Feb~Jul|cut(x, 4), data=SISmm, auto.key=list(space='right'))
##Faster with hexbinplot
hexbinplot(Jan~Jul|cut(x, 6), data=SISmm)
}
}
\keyword{methods}
\keyword{spatial}
