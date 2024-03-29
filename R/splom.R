setGeneric('splom')

.splom <- function(dat, plot.loess, colramp, varname.cex, ...)
{
    diag.panel = function(x,...)
    {
        yrng <- current.panel.limits()$ylim
        d <- density(x)
        d$y <- with(d, yrng[1] + 0.95 * diff(yrng) * y / max(y) )
        panel.lines(d)
        diag.panel.splom(x,...)
    }
    
    lower.panel = function(x, y,
                           plot.loess = plot.loess,
                           ...)
    {
        panel.hexbinplot(x, y, ...)
        if (plot.loess) panel.loess(x, y, ..., col = 'red')
    }
    splom(~ dat,
          colramp = colramp,
          varname.cex = varname.cex,
          plot.loess = plot.loess,
          panel = panel.hexbinplot,
          diag.panel = diag.panel,
          lower.panel = lower.panel,
          pscale = 0, ...)
}

setMethod('splom',
          signature(x='Raster', data = 'missing'),
          definition=function(x, data = NULL, maxpixels = 1e5,
                              plot.loess = FALSE,
                              colramp = BTC,
                              varname.cex = 0.6,
                              ...)
          {
              nms <- names(x)
              
              if (maxpixels < raster::ncell(x)) 
                  dat <- sampleRandom(x, maxpixels)
              else 
                  dat <- raster::values(x)
              
              colnames(dat) <- nms
              
              .splom(dat,
                     plot.loess,
                     colramp,
                     varname.cex,
                     ...)
          })


setMethod('splom',
          signature(x='SpatRaster', data = 'missing'),
          definition=function(x, data = NULL, maxpixels = 1e5,
                              plot.loess = FALSE,
                              colramp = BTC,
                              varname.cex = 0.6,
                              ...)
          {
              
              nms <- names(x)
              
              if (maxpixels < terra::ncell(x)) 
                  dat <- spatSample(x, maxpixels, method = "random")
              else 
                  dat <- terra::values(x)
              
              colnames(dat) <- nms
              
              .splom(dat,
                     plot.loess,
                     colramp,
                     varname.cex,
                     ...)
          })

