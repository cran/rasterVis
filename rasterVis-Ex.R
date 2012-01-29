pkgname <- "rasterVis"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('rasterVis')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("bwplot-methods")
### * bwplot-methods

flush(stderr()); flush(stdout())

### Name: bwplot-methods
### Title: Box and whisker plots of Raster objects.
### Aliases: bwplot bwplot,RasterStackBrick,missing-method
###   bwplot,formula,Raster-method
### Keywords: methods spatial

### ** Examples

f <- system.file("external/test.grd", package="raster")
r <- raster(f)
s <- stack(r, r-500, r+500)
bwplot(s)

## Not run: 
##D ##Solar irradiation data from CMSAF
##D ##Data available from http://www.box.net/shared/rl51y1t9sldxk54ogd44
##D 
##D old <- getwd()
##D ##change to your folder...
##D setwd('CMSAF')
##D listFich <- dir(pattern='2008')
##D stackSIS <- stack(listFich)
##D stackSIS <- stackSIS*24 ##from irradiance (W/m2) to irradiation Wh/m2
##D setwd(old)
##D 
##D idx <- seq(as.Date('2008-01-15'), as.Date('2008-12-15'), 'month')
##D 
##D SISmm <- setZ(stackSIS, idx)
##D layerNames(SISmm) <- month.abb
##D 
##D bwplot(SISmm)
##D bwplot(SISmm, FUN=as.yearqtr)##FUN applies to z if not NULL
## End(Not run)
## Not run: 
##D ##http://neo.sci.gsfc.nasa.gov/Search.html?group=64
##D pop <- raster('875430rgb-167772161.0.FLOAT.TIFF')
##D pop[pop==99999] <- NA
##D levelplot(pop, zscaleLog=10, par.settings=BTCTheme,
##D           panel=panel.levelplot.raster, interpolate=TRUE)
##D 
##D ##http://neo.sci.gsfc.nasa.gov/Search.html?group=20
##D landClass <- raster('241243rgb-167772161.0.TIFF')
##D landClass[landClass==254] <- NA
##D 
##D 
##D s <- stack(pop, landClass)
##D layerNames(s) <- c('pop', 'landClass')
##D 
##D bwplot(asinh(pop) ~ landClass|cut(y, 3), data=s,
##D        layout=c(3, 1), violin=FALSE)
##D 
##D bwplot(asinh(pop) ~ cut(y, 5)|landClass, data=s,
##D        scales=list(x=list(rot=45)), layout=c(4, 5),
##D        strip=strip.custom(strip.levels=TRUE))
## End(Not run)



cleanEx()
nameEx("chooseRegion")
### * chooseRegion

flush(stderr()); flush(stdout())

### Name: Interaction
### Title: Interaction with trellis objects.
### Aliases: identifyRaster chooseRegion identifyRaster
###   identifyRaster,Raster-method
### Keywords: methods spatial

### ** Examples

f <- system.file("external/test.grd", package="raster")
r <- raster(f)
levelplot(r)
##Do not close the last graphical window
##Use the left button of the mouse to identify points and the right button to finish
chosen_r <- identifyRaster(r, values=TRUE)
chosen_r
s <- stack(r, r-500, r+500)
levelplot(s)
chosen_s <- identifyRaster(s, values=TRUE)
chosen_s

## Not run: 
##D ##The package mgcv is needed for the next example
##D ##Use the left button of the mouse to build a border with points, and the right button to finish.
##D ##The points enclosed by the border will be highlighted and returned as a SpatialPoints object.
##D levelplot(s)
##D reg <- chooseRegion()
##D summary(reg)
## End(Not run)

## Not run: 
##D ##Solar irradiation data from CMSAF
##D ##Data available from http://www.box.net/shared/rl51y1t9sldxk54ogd44
##D 
##D old <- getwd()
##D ##change to your folder...
##D setwd('CMSAF')
##D listFich <- dir(pattern='2008')
##D stackSIS <- stack(listFich)
##D stackSIS <- stackSIS*24 ##from irradiance (W/m2) to irradiation Wh/m2
##D setwd(old)
##D 
##D idx <- seq(as.Date('2008-01-15'), as.Date('2008-12-15'), 'month')
##D 
##D SISmm <- setZ(stackSIS, idx)
##D layerNames(SISmm) <- month.abb
##D 
##D levelplot(SISmm)
##D 
##D ##Do not close the last graphical window
##D ##Interaction
##D ##Use the left button of the mouse to identify points and the right button to finish
##D chosen <- identifyRaster(SISmm, layer=3, values=TRUE)
##D chosen
##D ##Use the left button of the mouse to build a border with points, and the right button to finish.
##D ##The points enclosed by the border will be highlighted and returned as a SpatialPoints object.
##D reg <- chooseRegion()
##D summary(reg)
## End(Not run)



cleanEx()
nameEx("densityplot-methods")
### * densityplot-methods

flush(stderr()); flush(stdout())

### Name: densityplot-methods
### Title: Density plots for Raster objects.
### Aliases: densityplot densityplot,RasterLayer,missing-method
###   densityplot,RasterStackBrick,missing-method
###   densityplot,formula,Raster-method
### Keywords: methods spatial

### ** Examples

f <- system.file("external/test.grd", package="raster")
r <- raster(f)
densityplot(r)
s <- stack(r, r+500, r-500)
densityplot(s)

## Not run: 
##D ##Solar irradiation data from CMSAF
##D ##Data available from http://www.box.net/shared/rl51y1t9sldxk54ogd44
##D 
##D old <- getwd()
##D ##change to your folder...
##D setwd('CMSAF')
##D listFich <- dir(pattern='2008')
##D stackSIS <- stack(listFich)
##D stackSIS <- stackSIS*24 ##from irradiance (W/m2) to irradiation Wh/m2
##D setwd(old)
##D 
##D idx <- seq(as.Date('2008-01-15'), as.Date('2008-12-15'), 'month')
##D 
##D SISmm <- setZ(stackSIS, idx)
##D layerNames(SISmm) <- month.abb
##D 
##D densityplot(SISmm)
##D densityplot(SISmm, FUN=as.yearqtr)##FUN applies to z if not NULL
## End(Not run)
## Not run: 
##D ##http://neo.sci.gsfc.nasa.gov/Search.html?group=64
##D pop <- raster('875430rgb-167772161.0.FLOAT.TIFF')
##D pop[pop==99999] <- NA
##D levelplot(pop, zscaleLog=10, par.settings=BTCTheme,
##D           panel=panel.levelplot.raster, interpolate=TRUE)
##D 
##D ##http://neo.sci.gsfc.nasa.gov/Search.html?group=20
##D landClass <- raster('241243rgb-167772161.0.TIFF')
##D landClass[landClass==254] <- NA
##D 
##D 
##D s <- stack(pop, landClass)
##D layerNames(s) <- c('pop', 'landClass')
##D 
##D densityplot(~asinh(pop)|landClass, data=s,
##D             scales=list(relation='free'),
##D             strip=strip.custom(strip.levels=TRUE))
## End(Not run)



cleanEx()
nameEx("gplot-methods")
### * gplot-methods

flush(stderr()); flush(stdout())

### Name: gplot-methods
### Title: Use ggplot to plot a Raster* object
### Aliases: gplot gplot,Raster-method
### Keywords: methods spatial

### ** Examples
 
## Not run: 
##D r <- raster(system.file("external/test.grd", package="raster"))
##D s <- stack(r, r*2)
##D layerNames(s) <- c('meuse', 'meuse x 2')
##D 
##D library(ggplot2)
##D 
##D theme_set(theme_bw())
##D gplot(s) + geom_tile(aes(fill = value)) +
##D           facet_wrap(~ variable) +
##D           scale_fill_gradient(low = 'white', high = 'blue') +
##D           coord_equal()
## End(Not run)



cleanEx()
nameEx("hexbinplot")
### * hexbinplot

flush(stderr()); flush(stdout())

### Name: Formula methods
### Title: Formula methods
### Aliases: hexbinplot hexbinplot,formula,Raster-method
###   xyplot,formula,Raster-method
### Keywords: methods spatial

### ** Examples

f <- system.file("external/test.grd", package="raster")
r <- raster(f)
layerNames(r)

xyplot(test~y, data=r, alpha=0.5)

## Not run: 
##D ##Solar irradiation data from CMSAF
##D ##Data available from http://www.box.net/shared/rl51y1t9sldxk54ogd44
##D 
##D old <- getwd()
##D ##change to your folder...
##D setwd('CMSAF')
##D listFich <- dir(pattern='2008')
##D stackSIS <- stack(listFich)SISmm <- SISmm*24 ##from irradiance (W/m2) to irradiation Wh/m2
##D setwd(old)
##D 
##D layerNames(SISmm) <- month.abb
##D 
##D ##Relation between the January & February versus July radiation for four
##D ##differents longitude regions.
##D xyplot(Jan+Feb~Jul|cut(x, 4), data=SISmm, auto.key=list(space='right'))
##D ##Faster with hexbinplot
##D hexbinplot(Jan~Jul|cut(x, 6), data=SISmm)
## End(Not run)



cleanEx()
nameEx("histogram-methods")
### * histogram-methods

flush(stderr()); flush(stdout())

### Name: histogram-methods
### Title: Histogram of Raster objects.
### Aliases: histogram histogram,RasterLayer,missing-method
###   histogram,RasterStackBrick,missing-method
###   histogram,formula,Raster-method
### Keywords: methods spatial

### ** Examples

f <- system.file("external/test.grd", package="raster")
r <- raster(f)
histogram(r)
s <- stack(r, r+500, r-500)
histogram(s)

## Not run: 
##D ##Solar irradiation data from CMSAF
##D ##Data available from http://www.box.net/shared/rl51y1t9sldxk54ogd44
##D 
##D old <- getwd()
##D ##change to your folder...
##D setwd('CMSAF')
##D listFich <- dir(pattern='2008')
##D stackSIS <- stack(listFich)
##D stackSIS <- stackSIS*24 ##from irradiance (W/m2) to irradiation Wh/m2
##D setwd(old)
##D 
##D idx <- seq(as.Date('2008-01-15'), as.Date('2008-12-15'), 'month')
##D 
##D SISmm <- setZ(stackSIS, idx)
##D layerNames(SISmm) <- month.abb
##D 
##D histogram(SISmm)
##D histogram(SISmm, FUN=as.yearqtr)
## End(Not run)

## Not run: 
##D ##http://neo.sci.gsfc.nasa.gov/Search.html?group=64
##D pop <- raster('875430rgb-167772161.0.FLOAT.TIFF')
##D pop[pop==99999] <- NA
##D levelplot(pop, zscaleLog=10, par.settings=BTCTheme,
##D           panel=panel.levelplot.raster, interpolate=TRUE)
##D 
##D ##http://neo.sci.gsfc.nasa.gov/Search.html?group=20
##D landClass <- raster('241243rgb-167772161.0.TIFF')
##D landClass[landClass==254] <- NA
##D 
##D 
##D s <- stack(pop, landClass)
##D layerNames(s) <- c('pop', 'landClass')
##D 
##D histogram(~asinh(pop)|landClass, data=s,
##D             scales=list(relation='free'),
##D             strip=strip.custom(strip.levels=TRUE))
## End(Not run)




cleanEx()
nameEx("horizonplot-methods")
### * horizonplot-methods

flush(stderr()); flush(stdout())

### Name: horizonplot-methods
### Title: Horizon plots of Raster objects.
### Aliases: horizonplot horizonplot,RasterStackBrick-method
###   horizonplot,RasterStackBrick,missing-method
### Keywords: methods spatial

### ** Examples

## Not run: 
##D ##Solar irradiation data from CMSAF
##D ##Data available from http://www.box.net/shared/rl51y1t9sldxk54ogd44
##D 
##D old <- getwd()
##D ##change to your folder...
##D setwd('CMSAF')
##D listFich <- dir(pattern='2008')
##D stackSIS <- stack(listFich)
##D stackSIS <- stackSIS*24 ##from irradiance (W/m2) to irradiation Wh/m2
##D setwd(old)
##D 
##D idx <- seq(as.Date('2008-01-15'), as.Date('2008-12-15'), 'month')
##D 
##D SISmm <- setZ(stackSIS, idx)
##D layerNames(SISmm) <- month.abb
##D 
##D horizonplot(SISmm)
## End(Not run)

## Not run: 
##D library(zoo)
##D 
##D url <- "ftp://ftp.wiley.com/public/sci_tech_med/spatio_temporal_data/"
##D sst.dat = read.table(paste(url, "SST011970_032003.dat", sep=''), header = FALSE) 
##D sst.ll = read.table(paste(url, "SSTlonlat.dat", sep=''), header = FALSE)
##D 
##D spSST <- SpatialPointsDataFrame(sst.ll, sst.dat)
##D gridded(spSST) <- TRUE
##D proj4string(spSST) = "+proj=longlat +datum=WGS84"
##D SST <- brick(spSST)
##D 
##D idx <- seq(as.Date('1970-01-01'), as.Date('2003-03-01'), by='month')
##D idx <- as.yearmon(idx)
##D SST <- setZ(SST, idx)
##D layerNames(SST) <- as.character(idx)
##D horizonplot(SST)
## End(Not run)




cleanEx()
nameEx("hovmoller-methods")
### * hovmoller-methods

flush(stderr()); flush(stdout())

### Name: hovmoller-methods
### Title: Hovmoller plots
### Aliases: hovmoller hovmoller,RasterStackBrick-method
### Keywords: methods spatial methods

### ** Examples

## Not run: 
##D ##Solar irradiation data from CMSAF
##D ##Data available from http://www.box.net/shared/rl51y1t9sldxk54ogd44
##D 
##D old <- getwd()
##D ##change to your folder...
##D setwd('CMSAF')
##D listFich <- dir(pattern='2008')
##D stackSIS <- stack(listFich)
##D stackSIS <- stackSIS*24 ##from irradiance (W/m2) to irradiation Wh/m2
##D setwd(old)
##D 
##D idx <- seq(as.Date('2008-01-15'), as.Date('2008-12-15'), 'month')
##D 
##D SISmm <- setZ(stackSIS, idx)
##D layerNames(SISmm) <- month.abb
##D 
##D hovmoller(SISmm, dirXY=y, xlab='Latitude')##latitude
##D hovmoller(SISmm, dirXY=y, labels=FALSE, xlab='Latitude')##without labels
##D hovmoller(SISmm, dirXY=y, add.contour=FALSE, xlab='Latitude')##without contours
##D 
##D hovmoller(SISmm, dirXY=sqrt(x^2+y^2))##a function of coordinates...
## End(Not run)

## Not run: 
##D library(zoo)
##D 
##D url <- "ftp://ftp.wiley.com/public/sci_tech_med/spatio_temporal_data/"
##D sst.dat = read.table(paste(url, "SST011970_032003.dat", sep=''), header = FALSE) 
##D sst.ll = read.table(paste(url, "SSTlonlat.dat", sep=''), header = FALSE)
##D 
##D spSST <- SpatialPointsDataFrame(sst.ll, sst.dat)
##D gridded(spSST) <- TRUE
##D proj4string(spSST) = "+proj=longlat +datum=WGS84"
##D SST <- brick(spSST)
##D 
##D idx <- seq(as.Date('1970-01-01'), as.Date('2003-03-01'), by='month')
##D idx <- as.yearmon(idx)
##D SST <- setZ(SST, idx)
##D layerNames(SST) <- as.character(idx)
##D hovmoller(SST, panel=panel.levelplot.raster,
##D           interpolate=TRUE, par.settings=RdBuTheme)
## End(Not run)



cleanEx()
nameEx("levelplot-methods")
### * levelplot-methods

flush(stderr()); flush(stdout())

### Name: levelplot-methods
### Title: Level and contour plots of Raster objects.
### Aliases: levelplot contourplot levelplot,Raster,missing-method
###   contourplot,Raster,missing-method
### Keywords: methods spatial

### ** Examples

f <- system.file("external/test.grd", package="raster")
r <- raster(f)
levelplot(r)

## defining the scales for the marginal plot
levelplot(r, scales.margin=list(x=c(100, 600), y=c(100, 1000)))
## if a component of the list is null, it is internally calculated
levelplot(r, scales.margin=list(x=c(100, 1000)))

## log-scale
levelplot(r^2, zscaleLog=TRUE, contour=TRUE)

s <- stack(r, r+500, r-500)
levelplot(s, contour=TRUE)
contourplot(s, labels=list(cex=0.4), cuts=12)

##Add a layer of sampling points
##and change the theme
pts <- sampleRandom(r, size=20, sp=TRUE)
levelplot(r, par.settings=BTCTheme) + layer(sp.points(pts, col='red'))
contourplot(r, labels=FALSE) + layer(sp.points(pts, col='red'))

## Not run: 
##D ##Solar irradiation data from CMSAF
##D ##Data available from http://www.box.net/shared/rl51y1t9sldxk54ogd44
##D 
##D old <- getwd()
##D ##change to your folder...
##D setwd('CMSAF')
##D listFich <- dir(pattern='2008')
##D stackSIS <- stack(listFich)
##D stackSIS <- stackSIS*24 ##from irradiance (W/m2) to irradiation Wh/m2
##D setwd(old)
##D 
##D idx <- seq(as.Date('2008-01-15'), as.Date('2008-12-15'), 'month')
##D 
##D SISmm <- setZ(stackSIS, idx)
##D layerNames(SISmm) <- month.abb
##D 
##D levelplot(SISmm)
##D 
##D levelplot(SISmm, layers=1, FUN.margin=median, contour=TRUE)
## End(Not run)


cleanEx()
nameEx("plot3d")
### * plot3d

flush(stderr()); flush(stdout())

### Name: plot3D
### Title: Interactive 3D plot of a RasterLayer
### Aliases: plot3D plot3D,RasterLayer-method
### Keywords: methods spatial

### ** Examples

if (require(rgl)) {
data(volcano)
r <- raster(volcano)
drape <- cut(r, 5)
plot3D(r, drape=drape, zfac=4)
decorate3d(xlab = "x", ylab = "y", zlab = "z", axes=TRUE)
}



cleanEx()
nameEx("splom-methods")
### * splom-methods

flush(stderr()); flush(stdout())

### Name: splom-methods
### Title: Scatter plot matrices of Raster objects.
### Aliases: splom splom,RasterStackBrick,missing-method
### Keywords: methods spatial

### ** Examples

## Not run: 
##D ##Solar irradiation data from CMSAF
##D ##Data available from http://www.box.net/shared/rl51y1t9sldxk54ogd44
##D 
##D old <- getwd()
##D ##change to your folder...
##D setwd('CMSAF')
##D listFich <- dir(pattern='2008')
##D stackSIS <- stack(listFich)
##D stackSIS <- stackSIS*24 ##from irradiance (W/m2) to irradiation Wh/m2
##D setwd(old)
##D 
##D idx <- seq(as.Date('2008-01-15'), as.Date('2008-12-15'), 'month')
##D 
##D SISmm <- setZ(stackSIS, idx)
##D layerNames(SISmm) <- month.abb
##D 
##D splom(SISmm)
## End(Not run)


cleanEx()
nameEx("vectorplot")
### * vectorplot

flush(stderr()); flush(stdout())

### Name: vectorplot-methods
### Title: Vector plots of Raster objects.
### Aliases: vectorplot vectorplot,Raster-method
### Keywords: methods spatial

### ** Examples

df <- expand.grid(x=seq(-2, 2, .1), y=seq(-2, 2, .1))
df$z <- with(df, (3*x^2 + y)*exp(-x^2-y^2))
r1 <- rasterFromXYZ(df)
df$z <- with(df, x*exp(-x^2-y^2))
r2 <- rasterFromXYZ(df)
df$z <- with(df, y*exp(-x^2-y^2))
r3 <- rasterFromXYZ(df)

projection(r1) <- projection(r2) <- projection(r3) <- CRS("+proj=longlat +datum=WGS84")

vectorplot(r1, par.settings=RdBuTheme)
vectorplot(r2, par.settings=RdBuTheme)
vectorplot(r3, par.settings=RdBuTheme)



cleanEx()
nameEx("xyLayer")
### * xyLayer

flush(stderr()); flush(stdout())

### Name: xyLayer
### Title: xyLayer
### Aliases: xyLayer
### Keywords: spatial

### ** Examples

f <- system.file("external/test.grd", package="raster")
r <- raster(f)
dirX <- xyLayer(r, x)
dirXY <-xyLayer(r, sqrt(x^2 + y^2))
levelplot(dirXY, margin=FALSE)



cleanEx()
nameEx("xyplot-methods")
### * xyplot-methods

flush(stderr()); flush(stdout())

### Name: xyplot-methods
### Title: xyplot for Raster objects
### Aliases: xyplot xyplot,RasterStackBrick,missing-method
### Keywords: methods spatial

### ** Examples

## Not run: 
##D ##Solar irradiation data from CMSAF
##D ##Data available from http://www.box.net/shared/rl51y1t9sldxk54ogd44
##D 
##D old <- getwd()
##D ##change to your folder...
##D setwd('CMSAF')
##D listFich <- dir(pattern='2008')
##D stackSIS <- stack(listFich)
##D stackSIS <- stackSIS*24 ##from irradiance (W/m2) to irradiation Wh/m2
##D setwd(old)
##D 
##D idx <- seq(as.Date('2008-01-15'), as.Date('2008-12-15'), 'month')
##D 
##D SISmm <- setZ(stackSIS, idx)
##D layerNames(SISmm) <- month.abb
##D 
##D xyplot(SISmm)
## End(Not run)



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
