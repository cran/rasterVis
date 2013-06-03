pkgname <- "rasterVis"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('rasterVis')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
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
##D names(SISmm) <- month.abb
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
##D names(s) <- c('pop', 'landClass')
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
##D names(SISmm) <- month.abb
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
##D names(SISmm) <- month.abb
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
##D names(s) <- c('pop', 'landClass')
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
##D names(s) <- c('meuse', 'meuse x 2')
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
names(r)

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
##D names(SISmm) <- month.abb
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
##D names(SISmm) <- month.abb
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
##D names(s) <- c('pop', 'landClass')
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
##D names(SISmm) <- month.abb
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
##D names(SST) <- as.character(idx)
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
##D names(SISmm) <- month.abb
##D 
##D ## Latitude as default
##D hovmoller(SISmm, xlab='Latitude')
##D 
##D ## With contour lines and labels
##D hovmoller(SISmm, labels=TRUE, add.contour=TRUE,
##D           xlab='Latitude')
##D 
##D ## Smooth color regions with latticeExtra::panel.2dsmoother
##D hovmoller(SISmm, panel=panel.2dsmoother, n=1000,
##D           labels=FALSE, add.contour=TRUE,
##D           xlab='Latitude')
##D 
##D ## Using a function of coordinates
##D hovmoller(SISmm, dirXY=sqrt(x^2+y^2))
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
##D names(SST) <- as.character(idx)
##D hovmoller(SST, panel=panel.levelplot.raster,
##D           xscale.components=xscale.raster.subticks,
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

## Change the color theme
levelplot(r, par.settings=GrTheme())
levelplot(r, par.settings=PuOrTheme())

myTheme=rasterTheme(region=brewer.pal('Blues', n=9))
levelplot(r, par.settings=myTheme)

## Define the legend breaks
my.at <- seq(100, 1850, 500)
levelplot(r, at=my.at)

myColorkey <- list(at=my.at, ## where the colors change
                   labels=list(
                     at=my.at ## where to print labels
                     ))
levelplot(r, at=my.at, colorkey=myColorkey)

myColorkey <- list(at=my.at, ## where the colors change
                   labels=list(
                     labels=letters[seq_along(my.at)], ## labels
                     at=my.at ## where to print labels
                     ))
levelplot(r, at=my.at, colorkey=myColorkey)

## shrink and border color
rCenter <- (maxValue(r) + minValue(r)) / 2
levelplot(r - rCenter, par.settings=RdBuTheme(), shrink=c(.8, 15), border='black')

## With subticks
levelplot(r, xscale.components=xscale.raster.subticks,
             yscale.components=yscale.raster.subticks)

levelplot(r, xscale.components=xscale.raster.subticks,
             yscale.components=yscale.raster.subticks,
             scales=list(x=list(rot=30, cex=0.8)))

## log-scale
levelplot(r^2, zscaleLog=TRUE, contour=TRUE)

## Customizing axis and title labels
levelplot(r, margin=FALSE,
          main=list('My plot', col='red'),
          xlab=c('This is the', 'X-Axis'),
          ylab=list('Y-Axis', rot=30, fontface='bold')
          )

## xlim and ylim to display a smaller region
levelplot(r, xlim=c(179000, 181000), ylim=c(329500, 334000))

## RasterStacks
s <- stack(r, r+500, r-500)
levelplot(s, contour=TRUE)
contourplot(s, labels=list(cex=0.4), cuts=12)

## Use of layout
levelplot(s, layout=c(1, 3))
levelplot(s, layout=c(1, 1))

## names.attr to change the labels of each panel
levelplot(s, names.attr=c('R', 'R + 500', 'R - 500'))

## defining the scales for the marginal plot
levelplot(r, scales.margin=list(x=c(100, 600), y=c(100, 1000)))
## if a component of the list is null, it is internally calculated
levelplot(r, scales.margin=list(x=c(100, 1000)))

## Add a layer of sampling points
## and change the theme
pts <- sampleRandom(r, size=20, sp=TRUE)

## Using +.trellis and layer from latticeExtra
levelplot(r, par.settings = BTCTheme) + layer(sp.points(pts, col = 'red'))
contourplot(r, labels = FALSE) + layer(sp.points(pts, col = 'red'))

## or with a custom panel function
levelplot(r, par.settings=BTCTheme,
          panel=function(...){
            panel.levelplot(...)
            sp.points(pts, col='red')
            })


## Categorical data
r <- raster(nrow=10, ncol=10)
r[] = 1
r[51:100] = 3
r[3:6, 1:5] = 5
r <- ratify(r)
     
rat <- levels(r)[[1]]
rat$landcover <- c('Pine', 'Oak', 'Meadow')
rat$class <- c('A1', 'B2', 'C3')
levels(r) <- rat
r

levelplot(r, col.regions=c('palegreen', 'midnightblue', 'indianred1'))

## with 'att' you can choose another variable from the RAT
levelplot(r, att=2, col.regions=c('palegreen', 'midnightblue', 'indianred1'))
levelplot(r, att='class', col.regions=c('palegreen', 'midnightblue', 'indianred1'))

r2 <- raster(r)
r2[] = 3
r2[51:100] = 1
r2[3:6, 1:5] = 5

r3 <- init(r, function(n)sample(c(1, 3, 5), n, replace=TRUE))

## Multilayer categorical Raster* are displayed only if their RATs are the same
levels(r2) <- levels(r3) <- levels(r)

s <- stack(r, r2, r3)
names(s) <- c('A', 'B', 'C')

levelplot(s)
levelplot(s, att=2)

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
##D names(SISmm) <- month.abb
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
extent(r) <- c(0, 610, 0, 870)
drape <- cut(r, 5)
plot3D(r, drape=drape, zfac=2)
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
##D names(SISmm) <- month.abb
##D 
##D splom(SISmm)
## End(Not run)


cleanEx()
nameEx("vectorplot")
### * vectorplot

flush(stderr()); flush(stdout())

### Name: vectorplot-methods
### Title: Vector plots of Raster objects.
### Aliases: vectorplot vectorplot,Raster-method streamplot
###   streamplot,Raster-method streamplot,RasterStack-method
### Keywords: methods spatial

### ** Examples

## Not run: 
##D proj <- CRS('+proj=longlat +datum=WGS84')
##D 
##D df <- expand.grid(x=seq(-2, 2, .01), y=seq(-2, 2, .01))
##D df$z <- with(df, (3*x^2 + y)*exp(-x^2-y^2))
##D r1 <- rasterFromXYZ(df, crs=proj)
##D df$z <- with(df, x*exp(-x^2-y^2))
##D r2 <- rasterFromXYZ(df, crs=proj)
##D df$z <- with(df, y*exp(-x^2-y^2))
##D r3 <- rasterFromXYZ(df, crs=proj)
##D s <- stack(r1, r2, r3)
##D names(s) <- c('R1', 'R2', 'R3')
##D 
##D vectorplot(r1)
##D vectorplot(r2, par.settings=RdBuTheme())
##D vectorplot(r3, par.settings=PuOrTheme())
##D 
##D 
##D ## If no cluster is provided, streamplot uses parallel::mclapply except
##D ## with Windows. Therefore, next code could spend a long time under
##D ## Windows.
##D streamplot(r1)
##D 
##D ## With a cluster
##D hosts <- rep('localhost', 4)
##D cl <- makeCluster(hosts)
##D streamplot(r2, cl=cl,
##D            par.settings=streamTheme(symbol=brewer.pal(n=5,
##D                                                       name='Reds')))
##D stopCluster(cl)
##D 
##D ## Without parallel
##D streamplot(r3, parallel=FALSE,
##D            par.settings=streamTheme(symbol=brewer.pal(n=5,
##D                                                       name='Greens')))
##D 
##D ## Configuration of droplets and streamlets
##D streamplot(s, layout=c(1, 3), droplet=list(pc=.2), streamlet=list(L=20),
##D            par.settings=streamTheme(cex=.6))
## End(Not run)




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
##D names(SISmm) <- month.abb
##D 
##D xyplot(SISmm)
## End(Not run)



### * <FOOTER>
###
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
