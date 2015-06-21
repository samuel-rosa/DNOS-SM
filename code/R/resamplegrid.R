# DESCRIPTION ##################################################################
# Evaluation of grid resampling methods and window size to calculate DEM
# derivatives. I evaluate three methods for resampling grids (SRTM DEM - 90 m) 
# to a smaller pixel size. The methods are:
# - nearest neighbour (1 cell),
# - bilinear (4 cells), and
# - bicubic (16 cells).
# All there methods are implemented in GRASS GIS (function r.resample.interp).
# I also evaluate how the window size used to calculate DEM derivatives affects
# the correlation between DEM derivatives and soil properties. Four window sizes
# are evaluated: 9, 29, 49, and 69 pixels. These values were choosen based on 
# the window size limit imposed by GRASS 6.4 (69 pixels) and the availability of
# data beyond the boundary of the study area (1000 m).
# Four DEM derivatives are evaluated:
# - relative elevation,
# - slope, and
# - topographic position index.
# SETTINGS #####################################################################
rm(list=ls())
gc()
# Load packages
require(sp)
require(rgdal)
require(maptools)
require(raster)
require(spgrass6)
require(xtable)
# Load data
load("sm-dnos-general.RData")
load("sm-dnos-grid-resample.RData")
ls()

# GRASS GIS definitions (mapset = "predictions")
initGRASS(gisBase = "/usr/lib/grass64/",
          gisDbase = GRASSgisDbase, location = "dnos-sm-rs",
          mapset = "predictions", pid = Sys.getpid(), override = TRUE)
system("g.region rast=dnos.raster")

# EVALUATE GRID RESAMPLING METHODS (SRTM DEM) ##################################
# A pixel size of 15 m is used to maximize the window size (69*15 = 1035 m) used
# to calculate DEM derivatives according to the availability of elevation data
# beyond the boundary of the study areas (1000 m).
system(paste("g.list rast"))

# nearest neighbour
system(paste("r.resamp.interp method=nearest input=srtm.dem90 output=srtm.dem15.near"))
system(paste("r.info srtm.dem15.near"))

# bilinear
system(paste("r.resamp.interp --o method=bilinear input=srtm.dem90 output=srtm.dem15.bili"))
system(paste("r.info srtm.dem15.bili"))

# bicubic
system(paste("r.resamp.interp method=bicubic input=srtm.dem90 output=srtm.dem15.bicu"))
system(paste("r.info srtm.dem15.bicu"))

# DEM derivative (slope)
# Here I evaluate how the resampling method affects the correlation
# between DEM derivatives and soil properties

window <- c(9, 29, 49, 69)

# near DEM
lapply(X = window, FUN = function(X){
  system(paste("r.param.scale input=srtm.dem15.near output=slope.near.",
               X, " size=", X, " param=slope", sep = ""))})
slope.near <- lapply(X = window, FUN = function(X){
  readRAST6(paste("slope.near.", X, sep = ""))})

# bilinear DEM
lapply(X = window, FUN = function(X){
  system(paste("r.param.scale input=srtm.dem15.bili output=slope.bili.",
               X, " size=", X, " param=slope", sep = ""))})
slope.bili <- lapply(X = window, FUN = function(X){
  readRAST6(paste("slope.bili.", X, sep = ""))})

# bicubic DEM
lapply(X = window, FUN = function(X){
  system(paste("r.param.scale input=srtm.dem15.bicu output=slope.bicu.",
               X, " size=", X, " param=slope", sep = ""))})
slope.bicu <- lapply(X = window, FUN = function(X){
  readRAST6(paste("slope.bicu.", X, sep = ""))})

# plot
dev.off()
pdf(paste(resample.dir, "resample-maps.pdf", sep = ""), height = 12, width = 16,
    pointsize = 30)
par(mfrow = c(3, 4), plt = c(0.10, 1, 0, 0.9), mgp = c(0, 0, 0))
image(slope.near[[1]], asp = 1, col = gray(seq(1:64)/64))
title(ylab = "nearest")
mtext(paste("window size = ", window[1], sep = ""), cex = 0.6)
image(slope.near[[2]], asp = 1, col = gray(seq(1:64)/64))
mtext(paste("window size = ", window[2], sep = ""), cex = 0.6)
image(slope.near[[3]], asp = 1, col = gray(seq(1:64)/64))
mtext(paste("window size = ", window[3], sep = ""), cex = 0.6)
image(slope.near[[4]], asp = 1, col = gray(seq(1:64)/64))
mtext(paste("window size = ", window[4], sep = ""), cex = 0.6)
image(slope.bili[[1]], asp = 1, col = gray(seq(1:64)/64))
title(ylab = "bilinear")
image(slope.bili[[2]], asp = 1, col = gray(seq(1:64)/64))
image(slope.bili[[3]], asp = 1, col = gray(seq(1:64)/64))
image(slope.bili[[4]], asp = 1, col = gray(seq(1:64)/64))
image(slope.bicu[[1]], asp = 1, col = gray(seq(1:64)/64))
title(ylab = "bicubic")
image(slope.bicu[[2]], asp = 1, col = gray(seq(1:64)/64))
image(slope.bicu[[3]], asp = 1, col = gray(seq(1:64)/64))
image(slope.bicu[[4]], asp = 1, col = gray(seq(1:64)/64))
dev.off()

# CONCLUSION: both resampling methods 'nearest neighbor'and 'bilinear' are
# inapropriate because calculated DEM derivatives (slope) present spurious
# artifacts that do not comply with reality.

# sample slope values to observation points
slope.near.pts <- lapply(X = slope.near,
                     FUN = function(X){
                       slope = extract(x = as(X, "RasterLayer"),
                                       y = as(cal.data, "SpatialPoints"),
                                       method = "simple")
                       slope
                       })
slope.bili.pts <- lapply(X = slope.bili,
                         FUN = function(X){
                           slope = extract(x = as(X, "RasterLayer"),
                                           y = as(cal.data, "SpatialPoints"),
                                           method = "simple")
                           slope
                         })
slope.bicu.pts <- lapply(X = slope.bicu,
                         FUN = function(X){
                           slope = extract(x = as(X, "RasterLayer"),
                                           y = as(cal.data, "SpatialPoints"),
                                           method = "simple")
                           slope
                         })
tmp <- data.frame(window = window,
                  near = round(as.numeric(lapply(slope.near.pts, mean))),
                  near.sd = round(as.numeric(lapply(slope.near.pts, sd))),
                  bili = round(as.numeric(lapply(slope.bili.pts, mean))),
                  bili.sd = round(as.numeric(lapply(slope.bili.pts, sd))),       
                  bicu = round(as.numeric(lapply(slope.bicu.pts, mean))),
                  bicu = round(as.numeric(lapply(slope.bicu.pts, sd))))
tmp
xtable(tmp, digits = 0)
rm(tmp)

# CONCLUSION: resampling methods cause a minor change in the mean and standard
# deviation of DEM derivatives. This indicates that different resampling do not
# cause significant changes in raster values. The main difference in on the
# creation (or removal) of spurious artifacts.

# correlation with soil properties
a <- lapply(X = slope.near.pts, FUN = function(X){round(cor(X, cal.data$clay), 2)})
b <- lapply(X = slope.bili.pts, FUN = function(X){round(cor(X, cal.data$clay), 2)})
c <- lapply(X = slope.bicu.pts, FUN = function(X){round(cor(X, cal.data$clay), 2)})

slope.cor <- data.frame(window, as.numeric(a), as.numeric(b), as.numeric(c))
colnames(slope.cor) <- c("window", "nearest", "bilinear", "bicubic")
slope.cor
xtable(slope.cor)
rm(a, b, c)

# CONCLUSION: because resampling methods do not significantly alter raster
# values, the correlation strength is maintained. On the other hand, increasing
# the window size used to calculate DEM derivatives increases the correlation
# with soil properties (clay content).

# FINAL CONCLUSION: the bicubic resampling method is the most appropriate to
# resample the SRTM DEM (90 m) to 15 m because it does not produce spurious
# artifacts when DEM derivatives are calculated. Besides, raster values are not
# strongly affected, indicating that it produces an apropriate representation
# of the reality.

# EFFECT OF WINDOW SIZE ON THE CORRELATION BETWEEN DEM DERIVATIVES AND SOIL ####
# PROPERTIES


# Relative elevation (ELEV_90)
lapply(X = window, FUN = function(X){
  system(paste("r.param.scale --o input=ELEV_90 output=rel.elev90.",
               X, " size=", X, " param=elev", sep = ""))})
rel.elev90 <- lapply(X = window, FUN = function(X){
  readRAST6(paste("rel.elev90.", X, sep = ""))})

# Relative elevation (ELEV_10)
lapply(X = window, FUN = function(X){
  system(paste("r.param.scale --o input=ELEV_10 output=rel.elev10.",
               X, " size=", X, " param=elev", sep = ""))})
rel.elev10 <- lapply(X = window, FUN = function(X){
  readRAST6(paste("rel.elev10.", X, sep = ""))})

# Slope (ELEV_90)
lapply(X = window, FUN = function(X){
  system(paste("r.param.scale input=ELEV_90 output=slope90.",
               X, " size=", X, " param=slope", sep = ""))})
slope90 <- lapply(X = window, FUN = function(X){
  readRAST6(paste("slope90.", X, sep = ""))})

# Slope (ELEV_10)
lapply(X = window, FUN = function(X){
  system(paste("r.param.scale input=ELEV_10 output=slope10.",
               X, " size=", X, " param=slope", sep = ""))})
slope10 <- lapply(X = window, FUN = function(X){
  readRAST6(paste("slope10.", X, sep = ""))})

# Topographic Position Index (ELEV_90)
lapply(X = c(window*pred.cell), function(X){
  system(paste("saga_cmd ta_morphometry 18 -DEM=", dem.dir,
               "ELEV_90.sgrd -TPI=", dem.dir, "TPI90-", X,
               ".sgrd -RADIUS_MIN=3 -RADIUS_MAX=", X, sep = ""))})
tpi90 <- lapply(X = c(window*pred.cell), FUN = function(X){
  readGDAL(paste(dem.dir, "TPI90-", X, ".sdat", sep = ""))})

# Topographic Position Index (ELEV_10)
lapply(X = c(window*pred.cell), function(X){
  system(paste("saga_cmd ta_morphometry 18 -DEM=", dem.dir,
               "ELEV_10.sgrd -TPI=", dem.dir, "TPI10-", X,
               ".sgrd -RADIUS_MIN=3 -RADIUS_MAX=", X, sep = ""))})
tpi10 <- lapply(X = c(window*pred.cell), FUN = function(X){
  readGDAL(paste(dem.dir, "TPI10-", X, ".sdat", sep = ""))})


# plot (ELEV_90)
dev.off()
pdf(paste(resample.dir, "dem90-derivatives.pdf", sep = ""), height = 12, width = 16,
    pointsize = 30)
par(mfrow = c(3, 4), plt = c(0.10, 1, 0, 0.9), mgp = c(0, 0, 0))
image(rel.elev90[[1]], asp = 1, col = gray(seq(1:64)/64))
title(ylab = "relative elevation (ELEV_90)")
mtext(paste("window size = ", window[1], sep = ""), cex = 0.6)
image(rel.elev90[[2]], asp = 1, col = gray(seq(1:64)/64))
mtext(paste("window size = ", window[2], sep = ""), cex = 0.6)
image(rel.elev90[[3]], asp = 1, col = gray(seq(1:64)/64))
mtext(paste("window size = ", window[3], sep = ""), cex = 0.6)
image(rel.elev90[[4]], asp = 1, col = gray(seq(1:64)/64))
mtext(paste("window size = ", window[4], sep = ""), cex = 0.6)
image(slope90[[1]], asp = 1, col = gray(seq(1:64)/64))
title(ylab = "slope (ELEV_90)")
image(slope90[[2]], asp = 1, col = gray(seq(1:64)/64))
image(slope90[[3]], asp = 1, col = gray(seq(1:64)/64))
image(slope90[[4]], asp = 1, col = gray(seq(1:64)/64))
image(tpi90[[1]], asp = 1, col = gray(seq(1:64)/64))
title(ylab = "topographic position index (ELEV_90)")
image(tpi90[[2]], asp = 1, col = gray(seq(1:64)/64))
image(tpi90[[3]], asp = 1, col = gray(seq(1:64)/64))
image(tpi90[[4]], asp = 1, col = gray(seq(1:64)/64))
dev.off()

# plot (ELEV_10)
dev.off()
pdf(paste(resample.dir, "dem10-derivatives.pdf", sep = ""), height = 12, width = 16,
    pointsize = 30)
par(mfrow = c(3, 4), plt = c(0.10, 1, 0, 0.9), mgp = c(0, 0, 0))
image(rel.elev10[[1]], asp = 1, col = gray(seq(1:64)/64))
title(ylab = "relative elevation (ELEV_10)")
mtext(paste("window size = ", window[1], sep = ""), cex = 0.6)
image(rel.elev10[[2]], asp = 1, col = gray(seq(1:64)/64))
mtext(paste("window size = ", window[2], sep = ""), cex = 0.6)
image(rel.elev10[[3]], asp = 1, col = gray(seq(1:64)/64))
mtext(paste("window size = ", window[3], sep = ""), cex = 0.6)
image(rel.elev10[[4]], asp = 1, col = gray(seq(1:64)/64))
mtext(paste("window size = ", window[4], sep = ""), cex = 0.6)
image(slope10[[1]], asp = 1, col = gray(seq(1:64)/64))
title(ylab = "slope (ELEV_10)")
image(slope10[[2]], asp = 1, col = gray(seq(1:64)/64))
image(slope10[[3]], asp = 1, col = gray(seq(1:64)/64))
image(slope10[[4]], asp = 1, col = gray(seq(1:64)/64))
image(tpi10[[1]], asp = 1, col = gray(seq(1:64)/64))
title(ylab = "topographic position index (ELEV_10)")
image(tpi10[[2]], asp = 1, col = gray(seq(1:64)/64))
image(tpi10[[3]], asp = 1, col = gray(seq(1:64)/64))
image(tpi10[[4]], asp = 1, col = gray(seq(1:64)/64))
dev.off()

# CONCLUSION: small window sizes yield DEM derivatives with a lot of detail
# which describe small terrain features. Large window sizes smooth DEM
# derivatives. A final conclusion cannot be derived solely from the plots.

# Sample raster data to observation points
rel.elev90.pts <- lapply(X = rel.elev90,
                         FUN = function(X){
                           elev90 = extract(x = as(X, "RasterLayer"),
                                          y = as(cal.data, "SpatialPoints"),
                                          method = "simple")
                           elev90})
slope90.pts <- lapply(X = slope90,
                         FUN = function(X){
                           slope90 = extract(x = as(X, "RasterLayer"),
                                            y = as(cal.data, "SpatialPoints"),
                                            method = "simple")
                           slope90})
tpi90.pts <- lapply(X = tpi90,
                    FUN = function(X){
                      tpi90 = extract(x = as(X, "RasterLayer"),
                                      y = as(cal.data, "SpatialPoints"),
                                      method = "simple")
                      tpi90})
round(as.numeric(lapply(rel.elev90.pts, mean)))
round(as.numeric(lapply(rel.elev90.pts, sd)))
round(as.numeric(lapply(slope90.pts, mean)))
round(as.numeric(lapply(slope90.pts, sd)))
round(as.numeric(lapply(tpi90.pts, mean)))
round(as.numeric(lapply(tpi90.pts, sd)))


rel.elev10.pts <- lapply(X = rel.elev10,
                         FUN = function(X){
                           elev10 = extract(x = as(X, "RasterLayer"),
                                            y = as(cal.data, "SpatialPoints"),
                                            method = "simple")
                           elev10})
slope10.pts <- lapply(X = slope10,
                    FUN = function(X){
                      slope10 = extract(x = as(X, "RasterLayer"),
                                      y = as(cal.data, "SpatialPoints"),
                                      method = "simple")
                      slope10})
tpi10.pts <- lapply(X = tpi10,
                    FUN = function(X){
                      tpi10 = extract(x = as(X, "RasterLayer"),
                                      y = as(cal.data, "SpatialPoints"),
                                      method = "simple")
                      tpi10})
round(as.numeric(lapply(rel.elev10.pts, mean)))
round(as.numeric(lapply(rel.elev10.pts, sd)))
round(as.numeric(lapply(slope10.pts, mean)))
round(as.numeric(lapply(slope10.pts, sd)))
round(as.numeric(lapply(tpi10.pts, mean)))
round(as.numeric(lapply(tpi10.pts, sd)))

tmp <- data.frame(window = as.integer(window),
                  rel.elev90 = round(as.numeric(lapply(rel.elev90.pts, mean))),
                  rel.elev10 = round(as.numeric(lapply(rel.elev10.pts, mean))),
                  slope90 = round(as.numeric(lapply(slope90.pts, mean))),
                  slope10 = round(as.numeric(lapply(slope10.pts, mean))),
                  tpi90 = round(as.numeric(lapply(tpi90.pts, mean))),
                  tpi10 = round(as.numeric(lapply(tpi10.pts, mean))))
tmp
xtable(tmp, digits = 0)
rm(tmp)

# CONCLUSION: mean relative elevation at observation points do not change when
# the window size is increased. But the smoothing effect is observed in
# the decrease of the standard deviation. the same occurs with the slope.
# On the other hand, mean topographic position index and slope at observation
# points decrease with increasing window size. The standard deviation of TPI
# increases with the window size. The same result was observed for both DEMs.

# Correlation with clay content
tmp <- data.frame(window = as.integer(window),
                  rel.elev90 = as.numeric(lapply(rel.elev90.pts,
                                                 function(X){round(cor(X, cal.data$clay), 4)})),
                  rel.elev10 = as.numeric(lapply(rel.elev10.pts,
                                                 function(X){round(cor(X, cal.data$clay), 4)})), 
                  slope90 = as.numeric(lapply(slope90.pts,
                                              function(X){round(cor(X, cal.data$clay), 4)})),
                  slope10 = as.numeric(lapply(slope10.pts,
                                              function(X){round(cor(X, cal.data$clay), 4)})),
                  tpi90 = as.numeric(lapply(tpi90.pts,
                                            function(X){round(cor(X, cal.data$clay), 4)})),
                  tpi10 = as.numeric(lapply(tpi10.pts,
                                            function(X){round(cor(X, cal.data$clay), 4)})))
tmp
xtable(tmp, digits = 3)
rm(tmp)

# CONCLUSION: increasing the window size to calculate the relative elevation
# did not alter the correlation between this DEM derivative and soil clay
# content. However, the correlation between soil clay content and topographic
# position index and slope had a large increase with increasing the window size.
# The same result was observed for both DEMs.

# FINAL CONCLUSION: increasing the window size to calculate DEM derivatives
# caused the correlation between DEM derivatives and soil clay content to remain
# the same (relative elevation) or increase (slope and topographic position
# index). Therefore, the largest possible window size should be used to
# calculate DEM derivatives. This size is limited to the availablity of data
# beyond the boundary of the study area.

# INTRODUCE SPATIAL FILTER TO SIMULATE MULTI-SCALE #############################
system("g.list rast")

# slope
size <- c(3, 7, 15, 31, 69, 139)
lapply(size, function(X){
  system(paste("r.neighbors -c --o in=SLP_90_3 output=tmp", X, " method=average size=", X, sep = ""))
})

system(paste("r.multi.scale input=ELEV_90 output=SLP_90_", 139, " size=", 139, " param=slope", sep = ""))

lapply(size, function(X){
  system(paste("r.null map=SLP_90_", X, " setnull=0", sep = ""))
})

system("r.univar map=tmp139 zones=BASIN_90")
system("r.univar map=SLP_90_139 zones=BASIN_90")
system("r.colors map=SLP_90_139 color=slope")
system("r.mask --o in=BASIN_90")
system("d.mon x3");system("d.rast.leg SLP_90_139")
system("r.colors map=tmp139 raster=SLP_90_139")
system("d.mon x2");system("d.rast.leg tmp139")
system("d.mon x1");system("d.histogram SLP_90_139")
system("d.mon x0");system("d.histogram tmp139")

system("r.mapcalc 'tmp=tmp139-SLP_90_139'")
system("d.mon x2");system("d.rast.leg tmp")
system("d.mon x1");system("d.histogram tmp")

system("g.remove MASK")

# read rasters
filter <- lapply(X = size, FUN = function(X){readRAST6(paste("tmp", X, sep = ""))})

grass <- lapply(X = size, FUN = function(X){readRAST6(paste("SLP_90_", X, sep = ""))})

# sample data
filter.pts <- lapply(X = filter,
                    FUN = function(X){
                      filter = extract(x = as(X, "RasterLayer"),
                                       y = as(cal.data, "SpatialPoints"),
                                       method = "simple")
                      filter})

grass.pts <- lapply(X = grass,
                    FUN = function(X){
                      grass = extract(x = as(X, "RasterLayer"),
                                      y = as(cal.data, "SpatialPoints"),
                                      method = "simple")
                      grass})

# correlations
data.frame(size = size,
           quadratic = as.numeric(lapply(grass.pts, function(X){round(cor(X, cal.data$clay), 2)})),
           filter = as.numeric(lapply(filter.pts, function(X){round(cor(X, cal.data$clay), 2)}))
           )




# SAVE DATA ####################################################################
ls()
save(pred.cell, rel.elev10, rel.elev90, rel.elev10.pts, rel.elev90.pts,
     resample.dir, slope.bicu, slope.bicu.pts, slope10, slope10.pts, slope90,
     slope90.pts,
     slope.bili, slope.bili.pts, slope.cor, slope.near, slope.near.pts, tpi10,
     tpi90, tpi10.pts, tpi90.pts, window,
     file = "sm-dnos-grid-resample.RData")

# End!