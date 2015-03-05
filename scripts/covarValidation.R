# DESCRIPTION ##################################################################
# Validation of the environmental covariates available, including:
# - three digital elevation models (SRTM, TOPODATA, and DSG),
# - two land use maps (DSG (1980, 1992) and Samuel-Rosa et al. (2011)),
# - two area-class soil maps (Azolin & Mutti (1988) and Miguel et al. (2012)),
# - two geological maps (Gasparetto et al. (1988) and Maciel Filho (1990)), and
# - a 2009 Google Earth image.
# SETTINGS #####################################################################
rm(list = ls())
gc()
# Load Packages
require(sp)
require(raster)
require(rgdal)
require(VecStatGraphs2D)
require(vec2dtransf)
require(pedometrics)
require(spgrass6)
# Load data
load("sm-dnos-covar-validation.RData")
load("sm-dnos-general.RData")
ls()

# user defined functions
class.match <- function(observed, predicted){
  data = data.frame(observed, predicted, stringsAsFactors = FALSE)
  s = seq(1, length(observed), 1)
  ind = lapply(seq(1, length(data[, 1]), 1), function(X){
    agrep(data[X, 1], data[X, 2], max.distance = 0)
  })
  ind = !is.na(as.numeric(ind))
  return(ind)
}

# GRASS GIS definitions
initGRASS(gisBase = "/usr/lib/grass64/", gisDbase = GRASSgisDbase,
          location = "dnos-sm-rs", mapset = "predictions",
          pid = Sys.getpid(), override = TRUE)

# LOAD ENVIRONMENTAL COVARIATES ################################################

# ELEV_90 (SRTM v41 DEM)
ELEV_90 <- raster(readRAST6("ELEV_90"))
ELEV_90
image(ELEV_90, col = topo.colors(100), asp = 1)

# ELEV_30 (TOPODATA DEM)
ELEV_30 <- raster(readRAST6("ELEV_30"))
ELEV_30
image(ELEV_30, col = topo.colors(100), asp = 1)

# ELEV_10 (CONTOUR LINES)
ELEV_10 <- raster(readRAST6("ELEV_10"))
ELEV_10
image(ELEV_10, col = topo.colors(100), asp = 1)

# SOIL_25
SOIL_25 <- raster(readRAST6("SOIL_25"))
SOIL_25
image(SOIL_25, col = topo.colors(100), asp = 1)
plot(val_data, col = "red", pch = 20, add = TRUE)

# SOIL_100
SOIL_100 <- raster(readRAST6("SOIL_100"))
SOIL_100
image(SOIL_100, col = topo.colors(100), asp = 1)
plot(val.data, col = "red", pch = 20, add = TRUE)

# LU2009
LU2009 <- raster(readRAST6("LU2009"))
LU2009
image(LU2009, col = topo.colors(100), asp = 1)
plot(val.data, col = "red", pch = 20, add = TRUE)

# LU1980
LU1980 <- raster(readRAST6("LU1980"))
LU1980
image(LU1980, col = topo.colors(100), asp = 1)
plot(val.data, col = "red", pch = 20, add = TRUE)

# GEO_25
GEO_25 <- raster(readRAST6("predGEO_25"))
GEO_25
image(GEO_25, col = topo.colors(100), asp = 1)
plot(val.data, col = "red", pch = 20, add = TRUE)

# GEO_50
GEO_50 <- raster(readRAST6("GEO_50"))
GEO_50
image(GEO_50, col = topo.colors(100), asp = 1)
plot(val.data, col = "red", pch = 20, add = TRUE)

# LOAD POINT DATA ##############################################################

# Ground Control Points (post-processed) ---------------------------------------
gcp.pos <- shapefile(paste(point.dir, "gcp_points.shp", sep = ""))
colnames(gcp.pos@data)[colnames(gcp.pos@data) == "id"] <- "siteID"
gcp.pos <- gcp.pos[order(gcp.pos$siteID),]
gcp.pos <- spTransform(gcp.pos, wgs1984utm22s)
plot(gcp.pos)

# Ground Control Points (Topo map 1:25.000) ------------------------------------
gcp.car25 <- shapefile(paste(point.dir, "gcp_carta25.shp", sep = ""))
colnames(gcp.car25@data)[colnames(gcp.car25@data) == "id"] <- "siteID"
gcp.car25 <- gcp.car25[order(gcp.car25$siteID),]
gcp.car25 <- spTransform(gcp.car25, wgs1984utm22s)
plot(gcp.car25)

# Ground Control Points (Google Earth) -----------------------------------------
gcp.goo <- shapefile(paste(point.dir, "gcp_google.shp", sep = ""))
colnames(gcp.goo@data)[colnames(gcp.goo@data) == "Name"] <- "siteID"
gcp.goo <- gcp.goo[order(as.numeric(gcp.goo$siteID)),]
gcp.goo <- spTransform(gcp.goo, wgs1984utm22s)
plot(gcp.goo)

# Ground Control Points (RapidEye) ---------------------------------------------
gcp_eye <- shapefile(paste(point.dir, "gcp_rapideye.shp", sep = ""))
colnames(gcp_eye@data)[colnames(gcp_eye@data) == "id"] <- "siteID"
gcp_eye <- gcp_eye[order(gcp_eye$siteID),]
gcp_eye <- spTransform(gcp_eye, wgs1984utm22s)
plot(gcp_eye)

# Ground Control Points (Geological map 1:25K) ---------------------------------
gcp.geo25 <- shapefile(paste(point.dir, "gcp_geology25.shp", sep = ""))
colnames(gcp.geo25@data)[colnames(gcp.geo25@data) == "id"] <- "siteID"
gcp.geo25 <- gcp.geo25[order(gcp.geo25$siteID),]
gcp.geo25 <- spTransform(gcp.geo25, wgs1984utm22s)
plot(gcp.geo25)

# Ground Control Points (Geological map 1:50K) ---------------------------------
gcp.geo50 <- shapefile(paste(point.dir, "gcp_geology50.shp", sep = ""))
colnames(gcp.geo50@data)[colnames(gcp.geo50@data) == "id"] <- "siteID"
gcp.geo50 <- gcp.geo50[order(gcp.geo50$siteID),]
gcp.geo50 <- spTransform(gcp.geo50, wgs1984utm22s)
plot(gcp.geo50)

# Ground Control Points (Soil map 1:100K) --------------------------------------
gcp.soil100 <- shapefile(paste(point.dir, "gcp_soil1988.shp", sep = ""))
colnames(gcp.soil100@data)[colnames(gcp.soil100@data) == "id"] <- "siteID"
gcp.soil100 <- gcp.soil100[order(gcp.soil100$siteID),]
gcp.soil100 <- spTransform(gcp.soil100, wgs1984utm22s)
plot(gcp.soil100)

# Ground Control Points (Landsat 5 TM 2010dez26) -------------------------------
gcp_landsat <- shapefile(paste(point.dir, "gcp_2010dez26.shp", sep = ""))
colnames(gcp_landsat@data) <- "siteID"
gcp_landsat@data$siteID <- as.numeric(gcp_landsat$siteID)
gcp_landsat <- gcp_landsat[order(gcp_landsat$siteID),]
gcp_landsat <- spTransform(gcp_landsat, wgs1984utm22s)
plot(gcp_landsat)

# Check georeferencing
proj4string(ELEV_90)
proj4string(ELEV_30)
proj4string(ELEV_10)
proj4string(val.data)
proj4string(gcp.car25)
proj4string(gcp.goo)
proj4string(gcp_eye)
proj4string(gcp.geo25)
proj4string(gcp.geo50)
proj4string(gcp.soil100)
proj4string(gcp_landsat)

# DATA ANALYSIS (LOCATIONAL VALIDATION) ########################################
ls()

# Topographic map (1:25K) ======================================================

# Prepare spsurvey.object
spsurvey.car25 <- 
  spsurvey.analysis(design = coordenadas(gcp.car25),
                    data.cont = deltagcp(gcp.pos, gcp.car25),
                    popcorrect = TRUE, pcfsize = length(gcp.car25$siteID),
                    support = rep(1, length(gcp.car25$siteID)),
                    wgt = rep(1, length(gcp.car25$siteID)), vartype = "SRS")
# CDF
cdf.car25 <- cont.analysis(spsurvey.obj = spsurvey.car25)
levels(cdf.car25$Pct$Indicator)
cdfStats(cdf.car25, "dx", all = TRUE)
cdfTable(cdf.car25)
# Plot azimuth
dev.off()
cairo_pdf(paste(covar.validation.dir, "azim-car25.pdf", sep = ""))
par(cex = 0.5, lwd = 0.1, ps=25)
DrawHistogram(deltagcp(gcp.pos, gcp.car25)$azimuth, Direction = 2)
title(xlab = "Topographic maps (1:25k)")
dev.off()

# Google Earth  ================================================================
# Prepare spsurvey.object 
spsurvey.goo <- 
  spsurvey.analysis(design = coordenadas(gcp.goo),
                    data.cont = deltagcp(gcp.pos, gcp.goo),
                    popcorrect = TRUE, pcfsize = length(gcp.goo$siteID),
                    support = rep(1, length(gcp.goo$siteID)),
                    wgt = rep(1, length(gcp.goo$siteID)), vartype = "SRS")
# CDF 
cdf.goo <- cont.analysis(spsurvey.obj = spsurvey.goo)
levels(cdf.goo$Pct$Indicator)
cdfstats(cdf.goo, "dx", all = TRUE)
cdftable(cdf.goo)
# Plot azimuth 
dev.off()
cairo_pdf(paste(covar.validation.dir, "azim-goo.pdf", sep = ""))
par(cex = 0.5, lwd = 0.1, ps=25)
DrawHistogram(deltagcp(gcp.pos, gcp.goo)$azimuth, Direction = 2)
title(xlab = "Google Earth imagery")
dev.off()

# Geological map (1:25K) =======================================================
# Prepare spsurvey.object 
spsurvey.geo25 <- 
  spsurvey.analysis(design = coordenadas(gcp.geo25),
                    data.cont = deltagcp(gcp.pos, gcp.geo25),
                    popcorrect = TRUE, pcfsize = length(gcp.geo25$siteID),
                    support = rep(1, length(gcp.geo25$siteID)),
                    wgt = rep(1, length(gcp.geo25$siteID)), vartype = "SRS")

# CDF 
cdf.geo25 <- cont.analysis(spsurvey.obj = spsurvey.geo25)
levels(cdf.geo25$Pct$Indicator)
cdfstats(cdf.geo25, "dy", all = TRUE)
cdftable(cdf.geo25)
# Plot azimuth 
dev.off()
cairo_pdf(paste(covar.validation.dir, "azim-geo25.pdf", sep = ""))
par(cex = 0.5, lwd = 0.1, ps=25)
DrawHistogram(deltagcp(gcp.pos, gcp.geo25)$azimuth, Direction = 2)
title(xlab = "Geological map (1:25k)")
dev.off()

# Geological map (1:50K) =======================================================
# Prepare spsurvey.object 
spsurvey.geo50 <- 
  spsurvey.analysis(design = coordenadas(gcp.geo50),
                    data.cont = deltagcp(gcp.pos, gcp.geo50),
                    popcorrect = TRUE, pcfsize = length(gcp.geo50$siteID),
                    support = rep(1, length(gcp.geo50$siteID)),
                    wgt = rep(1, length(gcp.geo50$siteID)), vartype = "SRS")
# CDF 
cdf.geo50 <- cont.analysis(spsurvey.obj = spsurvey.geo50)
levels(cdf.geo50$Pct$Indicator)
cdfstats(cdf.geo50, "dx", all = TRUE)
cdftable(cdf.geo50)
# Plot azimuth 
dev.off()
cairo_pdf(paste(covar.validation.dir, "azim-geo50.pdf", sep = ""))
par(cex = 0.5, lwd = 0.1, ps=25)
DrawHistogram(deltagcp(gcp.pos, gcp.geo50)$azimuth, Direction = 2)
title(xlab = "Geological map (1:50k)")
dev.off()

# Soil map (1:100K) ============================================================
# Prepare spsurvey.object 
spsurvey.soil100 <- 
  spsurvey.analysis(design = coordenadas(gcp.soil100),
                    data.cont = deltagcp(gcp.pos, gcp.soil100),
                    popcorrect = TRUE, pcfsize = length(gcp.soil100$siteID),
                    support = rep(1, length(gcp.soil100$siteID)),
                    wgt = rep(1, length(gcp.soil100$siteID)), vartype = "SRS")
# CDF 
cdf.soil100 <- cont.analysis(spsurvey.obj = spsurvey.soil100)
levels(cdf.soil100$Pct$Indicator)
cdfstats(cdf.soil100, "module", all = TRUE)
cdftable(cdf.soil100)
# Plot azimuth 
dev.off()
cairo_pdf(paste(covar.validation.dir, "azim-soil100.pdf", sep = ""))
par(cex = 0.5, lwd = 0.1, ps=25)
DrawHistogram(deltagcp(gcp.pos, gcp.soil100)$azimuth, Direction = 2)
title(xlab = "Soil map (1:100k)")
dev.off()

# Landsat 5 TM (2010dez26) =====================================================
# Prepare spsurvey.object
spsurvey_landsat <-
 spsurvey.analysis(design = coordenadas(gcp_landsat),
                   data.cont = deltagcp(gcp.pos, gcp_landsat),
                   popcorrect = FALSE,
                   pcfsize = length(gcp_landsat$siteID),
                   support = rep(1, length(gcp_landsat$siteID)),
                   wgt = rep(1, length(gcp_landsat$siteID)),
                   vartype = "SRS")
# CDF
cdf_landsat <- cont.analysis(spsurvey.obj = spsurvey_landsat)
levels(cdf_landsat$Pct$Indicator)
cdfstats(cdf_landsat, "azimuth", all = TRUE)
cdftable(cdf_landsat)
# Plot azimuth
dev.off()
cairo_pdf(paste(covar.validation.dir, "azim-2010dez26.pdf", sep = ""))
par(cex = 0.5, lwd = 0.1, ps=25)
DrawHistogram(deltagcp(gcp.pos, gcp_landsat)$azimuth, Direction = 2)
title(xlab = "Landsat 5 TM")
dev.off()

# RapidEye =====================================================================
# Prepare spsurvey.object
spsurvey_eye <- 
  spsurvey.analysis(design = coordenadas(gcp_eye),
                    data.cont = deltagcp(gcp.pos, gcp_eye),
                    popcorrect = FALSE, 
                    pcfsize = length(gcp_eye$siteID),
                    support = rep(1, length(gcp_eye$siteID)),
                    wgt = rep(1, length(gcp_eye$siteID)), vartype = "SRS")
# CDF
cdf_eye <- cont.analysis(spsurvey.obj = spsurvey_eye)
levels(cdf_eye$Pct$Indicator)
cdfstats(cdf_eye, "azimuth", all = TRUE)
cdftable(cdf_eye)
# Plot azimuth
dev.off()
cairo_pdf(paste(covar.validation.dir, "azim-eye.pdf", sep = ""))
par(cex = 0.5, lwd = 0.1, ps=25)
DrawHistogram(deltagcp(gcp.pos, gcp_eye)$azimuth, Direction = 2)
title(xlab = "RapidEye imagery")
dev.off()

# DATA ANALYSIS (ATTRIBUTE VALIDATION) #########################################

# ELEV_90 ======================================================================
# Extract elevation data
ptsELEV_90 <- extract(x = as(ELEV_90, "RasterLayer"),
                      y = as(val.data, "SpatialPoints"),
                      method = "simple", sp = TRUE)
ptsELEV_90$ELEV_90 <- round(ptsELEV_90$ELEV_90)
colnames(ptsELEV_90@data) <- c("z")
ptsELEV_90$siteID <- rep(1:12, rep(5, 12))
ptsELEV_90
# Prepare spsurvey.object
dELEV_90 <- deltagcp(val.data, ptsELEV_90, type = "z", aggregate = TRUE)
# central transect coordinates
a <- rep(x = c(F, F, T, F, F), times = 12)
xcoord <- coordinates(val.data)[,1]
xcoord <- xcoord[which(a)]
ycoord <- coordinates(val.data)[,2]
ycoord <- ycoord[which(a)]
rm(a)
# spsurvey.object
spsurvey.ELEV_90 <- 
  spsurvey.analysis(design = data.frame(siteID = dELEV_90$siteID, xcoord, ycoord),
                    data.cont = dELEV_90, wgt = rep(1, length(dELEV_90$siteID)), 
                    support = rep(1, 12), vartype = "SRS")
rm(xcoord, ycoord)
# CDF
cdf.ELEV_90 <- cont.analysis(spsurvey.obj = spsurvey.ELEV_90)
levels(cdf.ELEV_90$Pct$Indicator)
cdfstats(cdf.ELEV_90, "dz", all = TRUE)
cdftable(cdf.ELEV_90, type = "z")
# CDF plot
dev.off()
pdf(paste(covar.validation.dir, "cdf-ELEV-90.pdf", sep = ""),
    width = 7, height = 2, pointsize = 8)
par(mfrow = c(1, 3), cex = 0.8)
cdfplot(obj = cdf.ELEV_90, ind = "dz", xlbl = "Error (m)", xlim = c(-40, -5),
        legloc = "TL")
cdfplot(obj = cdf.ELEV_90, ind = "abs.dz", xlbl = "Absolute error (m)",
        xlim = c(5, 40))
cdfplot(obj = cdf.ELEV_90, ind = "sq.dz", xlim = c(0, 1500),
        xlbl = expression(paste("Squared error (", m^2, ")", sep = "")))
dev.off()

# ELEV_30 ======================================================================
# Extract elevation data
ptsELEV_30 <- extract(x = as(ELEV_30, "RasterLayer"),
                      y = as(val.data, "SpatialPoints"),
                      method = "simple", sp = TRUE)
ptsELEV_30$ELEV_30 <- round(ptsELEV_30$ELEV_30)
colnames(ptsELEV_30@data) <- c("z")
ptsELEV_30$siteID <- rep(1:12, rep(5, 12))
ptsELEV_30
# Prepare spsurvey.object
dELEV_30 <- deltagcp(val.data, ptsELEV_30, type = "z", aggregate = TRUE)
# central transect coordinates
a <- rep(x = c(F, F, T, F, F), times = 12)
xcoord <- coordinates(val.data)[,1]
xcoord <- xcoord[which(a)]
ycoord <- coordinates(val.data)[,2]
ycoord <- ycoord[which(a)]
rm(a)
# spsurvey.object
spsurvey.ELEV_30 <-
  spsurvey.analysis(design = data.frame(siteID = dELEV_30$siteID, xcoord, ycoord),
                    data.cont = dELEV_30, wgt = rep(1, length(dELEV_30$siteID)), 
                    support = rep(1, 12), vartype = "SRS")
rm(xcoord, ycoord)
# CDF
cdf.ELEV_30 <- cont.analysis(spsurvey.obj = spsurvey.ELEV_30)
levels(cdf.ELEV_30$Pct$Indicator)
cdfstats(cdf.ELEV_30, "dz", all = TRUE)
cdftable(cdf.ELEV_30, type = "z")
# CDF plot
dev.off()
pdf(paste(covar.validation.dir, "cdf-ELEV-30.pdf", sep = ""),
    width = 7, height = 2, pointsize = 8)
par(mfrow = c(1, 3), cex = 0.8)
cdfplot(obj = cdf.ELEV_30, ind = "dz", xlbl = "Error (m)", xlim = c(-40, -5), 
        legloc = "TL")
cdfplot(obj = cdf.ELEV_30, ind = "abs.dz", xlbl = "Absolute error (m)",
        xlim = c(5, 40))
cdfplot(obj = cdf.ELEV_30, ind = "sq.dz", xlim = c(0, 1500),
        xlbl = expression(paste("Squared error (", m^2, ")", sep = "")))
dev.off()

# ELEV_10 ======================================================================
# Extract elevation data
ptsELEV_10 <- extract(x = as(ELEV_10, "RasterLayer"),
                      y = as(val.data, "SpatialPoints"), 
                      method = "simple", sp = TRUE)
ptsELEV_10$ELEV_10 <- round(ptsELEV_10$ELEV_10)
colnames(ptsELEV_10@data) <- c("z")
ptsELEV_10$siteID <- rep(1:12, rep(5, 12))
ptsELEV_10
# Prepare spsurvey.object
dELEV_10 <- gcpDiff(val.data, ptsELEV_10, type = "z", aggregate = TRUE)
# central transect coordinates
a <- rep(x = c(F, F, T, F, F), times = 12)
xcoord <- coordinates(val.data)[,1]
xcoord <- xcoord[which(a)]
ycoord <- coordinates(val.data)[,2]
ycoord <- ycoord[which(a)]
rm(a)
# spsurvey.object
spsurvey.ELEV_10 <- 
  spsurvey.analysis(design = data.frame(siteID = dELEV_10$siteID, xcoord, ycoord),
                    data.cont = dELEV_10, wgt = rep(1, length(dELEV_10$siteID)), 
                    support = rep(1, 12), vartype = "SRS")
rm(xcoord, ycoord)
# CDF
cdf.ELEV_10 <- cont.analysis(spsurvey.obj = spsurvey.ELEV_10)
levels(cdf.ELEV_10$Pct$Indicator)
cdfstats(cdf.ELEV_10, "dz", all = TRUE)
cdftable(cdf.ELEV_10, type = "z")
# CDF plot
dev.off()
pdf(paste(covar.validation.dir, "cdf-ELEV-10.pdf", sep = ""),
    width = 7, height = 2, pointsize = 8)
par(mfrow = c(1, 3), cex = 0.8)
cdfplot(obj = cdf.ELEV_10, ind = "dz", xlbl = "Error (m)", xlim = c(-40, -5),
        legloc = "TL")
cdfplot(obj = cdf.ELEV_10, ind = "abs.dz", xlbl = "Absolute error (m)",
        xlim = c(5, 40))
cdfplot(obj = cdf.ELEV_10, ind = "sq.dz", xlim = c(0, 1500),
        xlbl = expression(paste("Squared error (", m^2, ")", sep = "")))
dev.off()

# SOIL_25 ======================================================================

# create objects to store the data
ptsSOIL_25 <- as.list(c(NA, NA))
names(ptsSOIL_25) <- c("val_data", "cal_data")
cdfTAXA_25 <- as.list(c(NA, NA))
names(cdfTAXA_25) <- c("val_data", "cal_data")

# validation data
ptsSOIL_25$val_data <- extract(x = SOIL_25, y = as(val_data, "SpatialPoints"), 
                               method = "simple", sp = TRUE)
indx <- match(ptsSOIL_25$val_data$SOIL_25, SOIL_25.id$code, nomatch = 0)
ptsSOIL_25$val_data$pred <- rep(NA, length(ptsSOIL_25$val_data$SOIL_25))
ptsSOIL_25$val_data$pred[indx != 0] <- as.character(SOIL_25.id[indx, 2])
rm(indx)
ptsSOIL_25$val_data$pred <- sub(c("C-R"), c("CX RL"), ptsSOIL_25$val_data$pred)
ptsSOIL_25$val_data$ind <- class.match(val_data$TAXON, ptsSOIL_25$val_data$pred)
cdfTAXA_25$val_data <- category.est(catvar = ptsSOIL_25$val_data$ind,
                                    wgt = rep(1, length(ptsSOIL_25$val_data$ind)),
                                    x = coordinates(ptsSOIL_25$val_data)[,1],
                                    y = coordinates(ptsSOIL_25$val_data)[,2])
cdfTAXA_25

# calibration data
ptsSOIL_25$cal_data <- extract(x = SOIL_25, y = as(cal_data, "SpatialPoints"), 
                               method = "simple", sp = TRUE)
indx <- match(ptsSOIL_25$cal_data$SOIL_25, SOIL_25.id$code, nomatch = 0)
ptsSOIL_25$cal_data$pred <- rep(NA, length(ptsSOIL_25$cal_data$SOIL_25))
ptsSOIL_25$cal_data$pred[indx != 0] <- as.character(SOIL_25.id[indx, 2])
rm(indx)
ptsSOIL_25$cal_data$pred <- sub(c("C-R"), c("CX RL"), ptsSOIL_25$cal_data$pred)
ptsSOIL_25$cal_data$ind <- class.match(cal_data$TAXON, ptsSOIL_25$cal_data$pred)
cdfTAXA_25$cal_data <- category.est(catvar = ptsSOIL_25$cal_data$ind,
                                    wgt = rep(1, length(ptsSOIL_25$cal_data$ind)),
                                    x = coordinates(ptsSOIL_25$cal_data)[,1],
                                    y = coordinates(ptsSOIL_25$cal_data)[,2])
cdfTAXA_25
rm(SOIL_25)

# SOIL_100 =====================================================================

# create objects to store the data
ptsSOIL_100 <- as.list(c(NA, NA))
names(ptsSOIL_100) <- c("val_data", "cal_data")
cdfTAXA_100 <- as.list(c(NA, NA))
names(cdfTAXA_100) <- c("val_data", "cal_data")

# validation data
ptsSOIL_100$val_data <- extract(x = SOIL_100, y = as(val_data, "SpatialPoints"),
                                method = "simple", sp = TRUE)
indx <- match(ptsSOIL_100$val_data$SOIL_100, SOIL_100.id$code, nomatch = 0)
ptsSOIL_100$val_data$pred <- rep(NA, length(ptsSOIL_100$val_data$SOIL_100))
ptsSOIL_100$val_data$pred[indx != 0] <- as.character(SOIL_100.id[indx, 2])
rm(indx)
# conversion table
SOIL_100.id$SiBCS <- c("CX", "GX", "PBAC", "PBAC", "PV RL", "PV", "SX",
                       "PVA", "RL", "RL", "RL CX", "RL", "RL", "MT RL")
indx <- match(ptsSOIL_100$val_data$pred, SOIL_100.id$MU, nomatch = 0)
ptsSOIL_100$val_data$predSiBCS <- rep(NA, length(ptsSOIL_100$val_data$SOIL_100))
ptsSOIL_100$val_data$predSiBCS[indx != 0] <- as.character(SOIL_100.id[indx, "SiBCS"])
rm(indx)
ptsSOIL_100$val_data$predSiBCS
ptsSOIL_100$val_data$ind <- class.match(val_data$TAXON, 
                                        ptsSOIL_100$val_data$predSiBCS)
cdfTAXA_100$val_data <- category.est(catvar = ptsSOIL_100$val_data$ind,
                                     wgt = rep(1, length(ptsSOIL_100$val_data$ind)),
                                     x = coordinates(ptsSOIL_100$val_data)[,1], 
                                     y = coordinates(ptsSOIL_100$val_data)[,2])
cdfTAXA_100

# calibration data
ptsSOIL_100$cal_data <- extract(x = SOIL_100, y = as(cal_data, "SpatialPoints"),
                                method = "simple", sp = TRUE)
indx <- match(ptsSOIL_100$cal_data$SOIL_100, SOIL_100.id$code, nomatch = 0)
ptsSOIL_100$cal_data$pred <- rep(NA, length(ptsSOIL_100$cal_data$SOIL_100))
ptsSOIL_100$cal_data$pred[indx != 0] <- as.character(SOIL_100.id[indx, 2])
rm(indx)
# conversion table
SOIL_100.id$SiBCS <- c("CX", "GX", "PBAC", "PBAC", "PV RL", "PV", "SX",
                       "PVA", "RL", "RL", "RL CX", "RL", "RL", "MT RL")
indx <- match(ptsSOIL_100$cal_data$pred, SOIL_100.id$MU, nomatch = 0)
ptsSOIL_100$cal_data$predSiBCS <- rep(NA, length(ptsSOIL_100$cal_data$SOIL_100))
ptsSOIL_100$cal_data$predSiBCS[indx != 0] <- as.character(SOIL_100.id[indx, "SiBCS"])
rm(indx)
ptsSOIL_100$cal_data$predSiBCS
ptsSOIL_100$cal_data$ind <- class.match(cal_data$TAXON, 
                                        ptsSOIL_100$cal_data$predSiBCS)
cdfTAXA_100$cal_data <- category.est(catvar = ptsSOIL_100$cal_data$ind,
                                     wgt = rep(1, length(ptsSOIL_100$cal_data$ind)),
                                     x = coordinates(ptsSOIL_100$cal_data)[,1], 
                                     y = coordinates(ptsSOIL_100$cal_data)[,2])
cdfTAXA_100
rm(SOIL_100)

# LU1980 =====================================================================
# Extract attribute data
ptsLU1980 <- extract(x = LU1980, y = as(val.data, "SpatialPoints"), 
                     method = "simple", sp = TRUE)
indx <- match(ptsLU1980$LU1980, LU1980.id$code, nomatch = 0)
ptsLU1980$pred <- rep(NA, length(ptsLU1980$LU1980))
ptsLU1980$pred[indx != 0] <- as.character(LU1980.id[indx, "LU"])
rm(indx)
ptsLU1980
ptsLU1980$ind <- class.match(val.data$LAND, ptsLU1980$pred)
cdf.LU1980 <- category.est(catvar = ptsLU1980$ind, 
                           wgt = rep(1, length(ptsLU1980$ind)),
                           x = coordinates(ptsLU1980)[,1], 
                           y = coordinates(ptsLU1980)[,2])
cdf.LU1980[2,1:6][,-2]
rm(LU1980)

# LU2009 =====================================================================
# Extract attribute data
ptsLU2009 <- extract(x = LU2009, y = as(val.data, "SpatialPoints"), 
                     method = "simple", sp = TRUE)
indx <- match(ptsLU2009$LU2009, LU2009.id$code, nomatch = 0)
ptsLU2009$pred <- rep(NA, length(ptsLU2009$LU2009))
ptsLU2009$pred[indx != 0] <- as.character(LU2009.id[indx, "LU"])
rm(indx)
ptsLU2009
ptsLU2009$ind <- class.match(val.data$LAND, ptsLU2009$pred)
cdf.LU2009 <- category.est(catvar = ptsLU2009$ind, 
                           wgt = rep(1, length(ptsLU2009$ind)),
                           x = coordinates(ptsLU2009)[,1], 
                           y = coordinates(ptsLU2009)[,2])
cdf.LU2009[2,1:6][,-2]
rm(LU2009)

# GEO_25 =======================================================================
# Extract attribute data
ptsGEO_25 <- extract(x = GEO_25, y = as(val.data, "SpatialPoints"),
                     method = "simple", sp = TRUE)
indx <- match(ptsGEO_25$predGEO_25, geo25.id$code, nomatch = 0)
ptsGEO_25$pred <- rep(NA, length(ptsGEO_25$predGEO_25))
ptsGEO_25$pred[indx != 0] <- as.character(geo25.id[indx, "GEO"])
rm(indx)
# conversion table
geo25.id
geo25.id$PARENT <- c("igneous", "sedimentary", "sedimentary", "igneous",
                     "sedimentary", "sedimentary", "sedimentary")
indx <- match(ptsGEO_25$pred, geo25.id$GEO, nomatch = 0)
ptsGEO_25$PARENT <- rep(NA, length(ptsGEO_25$predGEO_25))
ptsGEO_25$PARENT[indx != 0] <- as.character(geo25.id[indx, "PARENT"])
rm(indx)
ptsGEO_25$ind <- class.match(val.data$PARENT, ptsGEO_25$PARENT)
cdf.GEO_25 <- category.est(catvar = ptsGEO_25$ind, wgt = rep(1, length(ptsGEO_25$ind)),
                           x = coordinates(ptsGEO_25)[,1], y = coordinates(ptsGEO_25)[,2])
cdf.GEO_25[2,1:6][,-2]
rm(GEO_25)

# GEO_50 =======================================================================
# Extract attribute data
ptsGEO_50 <- extract(x = GEO_50, y = as(val.data, "SpatialPoints"),
                     method = "simple", sp = TRUE)
# conversion table
system("r.category GEO_50")
GEO_50.id <- data.frame(code = c(1:5),
                        GEO = c("Serra Geral Formation - Inferior Sequence",
                                "Botucatu Formation", "Caturrita Formation",
                                "Serra Geral Formation - Superior Sequence",
                                "Santa Maria Formation - Alemoa Member"),
                        PARENT = c("igneous", "sedimentary", "sedimentary",
                                   "igneous", "sedimentary"),
                        stringsAsFactors = FALSE)
indx <- match(ptsGEO_50$GEO_50, GEO_50.id$code, nomatch = 0)
ptsGEO_50$pred <- rep(NA, length(ptsGEO_50$GEO_50))
ptsGEO_50$pred[indx != 0] <- as.character(GEO_50.id[indx, "GEO"])
rm(indx)
indx <- match(ptsGEO_50$pred, as.character(GEO_50.id$GEO), nomatch = 0)
ptsGEO_50$PARENT <- rep(NA, length(ptsGEO_50$GEO_50))
ptsGEO_50$PARENT[indx != 0] <- as.character(GEO_50.id[indx, "PARENT"])
rm(indx)
ptsGEO_50$ind <- class.match(val.data$PARENT, ptsGEO_50$PARENT)
cdf.GEO_50 <- category.est(catvar = ptsGEO_50$ind, wgt = rep(1, length(ptsGEO_50$ind)),
                           x = coordinates(ptsGEO_50)[,1], y = coordinates(ptsGEO_50)[,2])
cdf.GEO_50[2,1:6][,-2]
rm(GEO_50)

# DATA CORRECTION - FIT AFFINE TRANSFORMATION MODELS ###########################
# CAUTION: Corrections are applied to data in their respective working scripts!
ls()

# Landsat ======================================================================
affine_landsat <- 
  AffineTransformation(controlPoints = data.frame(gcp_landsat$coords.x1,
                                                  gcp_landsat$coords.x2,
                                                  gcp.pos$coords.x1,
                                                  gcp.pos$coords.x2))
calculateParameters(affine_landsat)
getParameters(affine_landsat)
plotGridTransformation(affine_landsat, bbox(gcp.pos), 10)
# applyTransformation(affine_landsat, obj)
getRMSE(affine_landsat)

# RapidEye =====================================================================
affine_eye <- 
  AffineTransformation(controlPoints = data.frame(gcp_eye$coords.x1,
                                                  gcp_eye$coords.x2,
                                                  gcp.pos$coords.x1,
                                                  gcp.pos$coords.x2))
calculateParameters(affine_eye)
getParameters(affine_eye)
plotGridTransformation(affine_eye, bbox(gcp.pos), 10)
# applyTransformation(affine_eye, obj)
getRMSE(affine_eye)

# Topographic Maps (1:25k) =====================================================
affine.car25 <- 
  AffineTransformation(controlPoints = data.frame(gcp.car25$coords.x1,
                                                  gcp.car25$coords.x2,
                                                  gcp.pos$coords.x1,
                                                  gcp.pos$coords.x2))
calculateParameters(affine.car25)
getParameters(affine.car25)
plotGridTransformation(affine.car25, bbox(gcp.pos), 10)
# applyTransformation(affine.car25, obj)
getRMSE(affine.car25)

# Geological map (1:25k) =======================================================
tmp <- gcp.pos[!is.na(match(gcp.pos$siteID, gcp.geo25$siteID)) == TRUE,]
affine.geo25 <-
  AffineTransformation(controlPoints = data.frame(gcp.geo25$coords.x1,
                                                  gcp.geo25$coords.x2,
                                                  tmp$coords.x1,
                                                  tmp$coords.x2))
rm(tmp)
calculateParameters(affine.geo25)
getParameters(affine.geo25)
plotGridTransformation(affine.geo25, bbox(gcp.pos), 10)
# applyTransformation(affine.geo25, land)

# Geological map (1:50k) =======================================================
tmp <- gcp.pos[!is.na(match(gcp.pos$siteID, gcp.geo50$siteID)) == TRUE,]
affine.geo50 <-
  AffineTransformation(controlPoints = data.frame(gcp.geo50$coords.x1,
                                                  gcp.geo50$coords.x2,
                                                  tmp$coords.x1,
                                                  tmp$coords.x2))
rm(tmp)
calculateParameters(affine.geo50)
getParameters(affine.geo50)
plotGridTransformation(affine.geo50, bbox(gcp.pos), 10)
#applyTransformation(affine.geo50, land)

# Save data ####################################################################
ls()
save(affine.car25, affine.geo25, affine.geo50, affine_eye, affine_landsat,
     # cummulative distribution functions
     cdf.car25, cdf_eye, cdf.geo25, cdf.geo50, cdf.goo, cdf.soil100, cdf.GEO_50,
     cdf.GEO_25, cdf_landsat,
     cdf.ELEV_10, cdf.ELEV_30, cdf.ELEV_90, cdfTAXA_25, cdfTAXA_100,
     cdf.LU1980, cdf.LU2009,
     dELEV_10, dELEV_30, dELEV_90,
     # ground control points
     gcp.car25, gcp_eye, gcp.geo25, gcp.geo50, gcp.goo, gcp.pos, gcp.soil100,
     gcp_landsat,
     # point data
     ptsELEV_90, ptsELEV_30, ptsELEV_10, ptsSOIL_25, ptsSOIL_100, ptsLU2009,
     ptsLU1980, ptsGEO_25, ptsGEO_50,
     spsurvey_landsat, spsurvey.car25, spsurvey.ELEV_10, spsurvey.ELEV_30, 
     spsurvey.ELEV_30,
     spsurvey_eye, spsurvey.geo25, spsurvey.geo50, spsurvey.goo, spsurvey.soil100,
     covar.validation.dir, val.data,  class.match,
     SOIL_25.id, SOIL_100.id, LU1980.id, LU2009.id, geo25.id, geo50.id, GEO_50.id,
     file = "sm-dnos-covar-validation.RData")
# End!