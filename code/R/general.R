# DESCRIPTION ##################################################################
# General configurations and pre-processing. We have some exploratory data
# analysis, the preparation of figures of the study area, etc.

# SETTINGS #####################################################################
rm(list = ls())
gc()
require(sp)
require(rgdal)
require(maptools)
require(raster)
require(spgrass6)
require(plyr)
library(plotKML)
library(car)
library(stringr)
library(fitdistrplus)
library(timeDate)
library(lattice)
library(latticeExtra)
require(grid)
library(maptools)
library(rgeos)
require(pedometrics)
load("sm-dnos-general.RData")
ls()

# SET DATA DIRECTORIES AND DEFINITIONS #########################################

# Directories ------------------------------------------------------------------
# Exploratory data analysis
explora_dir <- path.expand("~/Dropbox/dnos-sm-rs/explora/")
# Validation of the covariates
cov_valid_dir <- path.expand("~/Dropbox/dnos-sm-rs/covar-valid/")
# R data
rdata_dir <- path.expand("~/Dropbox/dnos-sm-rs/rdata/")
# Hydrological data
hydro_dir <- path.expand("~/Dropbox/dnos-sm-rs/hydrology/")
# DEM and its derivatives
dem_dir <- path.expand("~/Dropbox/dnos-sm-rs/dem/")
# Geological maps
geo_dir <- path.expand("~/Dropbox/dnos-sm-rs/geology/")
# Land use maps
land_dir <- path.expand("~/Dropbox/dnos-sm-rs/land-use/")
# Area-class soil maps
soil_dir <- path.expand("~/Dropbox/dnos-sm-rs/soil-class/")
# Point data
point_dir <- path.expand("~/Dropbox/dnos-sm-rs/point-data/")
# Resampling methors for grids
resample_dir <- path.expand("~/Dropbox/dnos-sm-rs/grid-resample/")
# Boudaries
bound_dir <- path.expand("~/Dropbox/dnos-sm-rs/boundaries/")
# 1st thesis article
firstArticle_dir <- path.expand("~/Dropbox/dnos-sm-rs/firstArticle/")
# Point pattern analysys
ppa_dir <- path.expand("~/Dropbox/dnos-sm-rs/ppa/")

# Cell size of the prediction grid ---------------------------------------------
# We choose the cellsize to be 5 metres because this is the smallest
# cellsize among the environmental covariates (RapidEye images). This way
# we do not lose any information due to the aggregation of spatial data.
cellsize <- 5

# GRASS GIS addons
# Sys.getenv("GRASS_ADDON_PATH")
# Sys.setenv("GRASS_ADDON_PATH" = path.expand("~/softwares/GRASS_ADDON_PATH"))
# grass.addons <- path.expand("~/softwares/GRASS_ADDON_PATH/")
# system(paste("g.extension r.gauss")) # install
# system(paste(grass.addons, "r.gauss -help", sep = ""))

# GRASS GIS DBase --------------------------------------------------------------
dbGRASS <- path.expand("~/dbGRASS")
initGRASS(gisBase = "/usr/lib/grass64/", gisDbase = dbGRASS, 
          location = "dnos-sm-rs", mapset = "predictions", pid = Sys.getpid(),
          override = TRUE)
#writeRAST6(dnos.raster, "dnos.raster", overwrite = TRUE)
system("g.region rast=dnos.raster")
gmeta6()

# Coordinate reference systems -------------------------------------------------
sirgas2000 <- CRS("+init=epsg:4674")
sirgas2000utm22s <- CRS("+init=epsg:31982")
wgs1984utm22s <- CRS("+init=epsg:32722")
wgs1984 <- CRS("+init=epsg:4326")
ca_utm22s <- CRS("+init=epsg:22522")

# SAVE DATA ####################################################################
ls()
save(bound_dir, ca_utm22s, cellsize, cov_valid_dir, dbGRASS, dem_dir, 
     explora_dir, firstArticle_dir, geo_dir, hydro_dir, land_dir, point_dir,
     ppa_dir, rdata_dir, resample_dir, sirgas2000, sirgas2000utm22s, 
     soil_dir, wgs1984, wgs1984utm22s, 
     file = paste(rdata_dir, "sm-dnos-general.RData", sep = ""))






# LOAD AND PROCESS DATA ########################################################

# Laboratory data ==============================================================
labData <- read.table("data/labData.csv", dec = ".", head = TRUE, sep = ";",
                       stringsAsFactors = FALSE, na.strings = "na")
head(labData)
coordinates(labData) <- ~ longitude + latitude
proj4string(labData) <- sirgas2000
labData <- spTransform(labData, wgs1984utm22s)
str(labData)
plot(labData)




# Validation Data ==============================================================
val_data <- read.table(paste(point.dir, "validation-data.csv", sep = ""),
                       sep = ";", head = TRUE, dec = ".", na.strings = "na")
str(val_data)
val_data$sampleid <- as.character(val_data$sampleid)
coordinates(val_data) <- ~ longitude + latitude
proj4string(val_data) <- sirgas2000
val_data <- spTransform(val_data, wgs1984utm22s)
plot(val_data, pch = 20, cex = 0.5)
plotKML(val_data, points_names = c(341:400))
# geoidal heights
write.table(data.frame(coordinates(spTransform(val_data, wgs1984))[,2],
                       coordinates(spTransform(val_data, wgs1984))[,1]),
            paste(covar.validation.dir, "INPUT.DAT", sep = ""),
            col.names = FALSE, row.names = FALSE)
setwd(covar.validation.dir)
system("wine F477")
setwd(rdata.dir)
N.egm <- read.table(paste(covar.validation.dir, "OUTF477.DAT", sep = ""))[,3]
N.egm
# H = h - N
# H:   Ellipsoid height
# h:   Orthometric height (mean-sea-level height)
# N:   Geoidal undulation
val_data$z <- as.numeric(val_data$altDGPS + N.egm)
rm(N.egm)
val_data$TAXON <- as.character(val_data$TAXON)

# Calibration Data =============================================================
cal_data <- read.table(paste(point_dir, "calibration-data.csv", sep = ""),
                       dec = ".", head = TRUE, sep = ";",
                       stringsAsFactors = FALSE, na.strings = "na")
head(cal_data)
coordinates(cal_data) <- ~ longitude + latitude
proj4string(cal_data) <- sirgas2000
cal_data <- spTransform(cal_data, wgs1984utm22s)
str(cal_data)
plot(cal_data)

# Calibration data - fill gaps (bulk density) ----------------------------------
density <- read.table(paste(point_dir, "calibration-data.csv", sep = ""),
                      dec = ".", head = TRUE, sep = ";",
                      stringsAsFactors = FALSE, na.strings = "na")
density <- density[, c("density_n", "density_sd", "density")]
cal_data@data <- data.frame(cal_data@data, density)
colnames(cal_data@data)
# fit model
density_model <- lm(density ~ TAXON + PARENT + LAND + carbon + clay + sand,
                    data = cal_data@data, na.action = na.exclude)
anova(density_model)
summary(density_model)
plot(cal_data$density ~ round(fitted(density_model), 2),
     xlim = c(0, 2), ylim = c(0, 2))
abline(a = 0, b = 1)
# predict
newdata <- data.frame(sampleid = cal_data$sampleid, density = cal_data$density, 
                      TAXON = cal_data$TAXON, PARENT = cal_data$PARENT,
                      LAND = cal_data$LAND, carbon = cal_data$carbon,
                      clay = cal_data$clay, sand = cal_data$sand)
newdata <- newdata[is.na(newdata$density), ]
newdata$fit <- round(predict(density_model, newdata, se.fit = TRUE)$fit, 2)
newdata$se.fit <- round(predict(density_model, newdata, se.fit = TRUE)$se.fit, 2)
newdata
cal_data$density[match(newdata$sampleid, cal_data$sampleid)] <- newdata$fit
cal_data$density_se.fit <- NA
cal_data$density_se.fit[match(newdata$sampleid, cal_data$sampleid)] <- newdata$fit
cal_data$density_se.fit
rm(density, newdata)

# Extent =======================================================================

# extent of the region in which there is data available
data.extent <- extent(c(220817, 237524, 6711252, 6724461))

# extent of the prediction grid (15 m) -----------------------------------------

# set bounding box
grow <- 1500
min.x <- floor((bbox(cal_data)["longitude", "min"] - grow)/cellsize)*cellsize
max.x <- ceiling((bbox(cal_data)["longitude", "max"] + grow)/cellsize)*cellsize
min.y <- floor((bbox(cal_data)["latitude", "min"] - grow)/cellsize)*cellsize
max.y <- ceiling((bbox(cal_data)["latitude", "max"] + grow)/cellsize)*cellsize

# calculate the number of grid cells
cells.x <- (max.x - min.x)/cellsize
cells.y <- (max.y - min.y)/cellsize
cells <- cells.x*cells.y

# create prediction grid (10 m)
dnos.raster <- SpatialGrid(GridTopology(cellcentre.offset = c(min.x+(cellsize/2),
                                                              min.y+(cellsize/2)),
               cellsize = c(cellsize, cellsize), cells.dim = c(cells.x, cells.y)),
               proj4string = wgs1984utm22s)
str(dnos.raster)
plot(extent(dnos.raster), asp = 1)
plot(cal_data, add = TRUE)
dnos.raster <- SpatialGridDataFrame(dnos.raster, data.frame(tmp = rep(1, cells)))
str(dnos.raster)
rm(min.x, max.x, min.y, max.y, cells.x, cells.y, cells, grow)
tmp <- as(extent(dnos.raster), "SpatialPolygons")
proj4string(tmp) <- wgs1984utm22s
shapefile(tmp, paste(bounding.dir, "dnos-raster-bbox.shp", sep = ""), overwrite=TRUE)
rm(tmp)

# extent of the entire catchment area
# dnos.extent <- extent(bbox(dnos.lim))
# extent of the study area
# sub.extent <- extent(bbox(dnos.sub.lim))
# #----------------------------------------------------
# # Limite da área de estudo
# #----------------------------------------------------
# dnos.limite <-
#   shapefile(filename="/home/alessandro/rdata/dnos_sm/shapefiles/limite_sub_bacia.shx",
#             verbose=TRUE)
# str(dnos.limite)
# (proj4string(dnos.limite) <- proj4string(dnos.coords))
# plot(dnos.limite)
# (tmp.lim <- as(dnos.limite, "SpatialPolygons"))
# class(dnos.raster)
# (tmp.ras <- as(dnos.raster, "RasterLayer"))
# (tmp.ras <- rasterize(tmp.lim, tmp.ras))
# str(dnos.raster.limite <- as(tmp.ras, "SpatialPointsDataFrame"))
# gridded(dnos.raster.limite) <- TRUE
# str(dnos.raster.limite)
# ls()
# rm(tmp.lim, tmp.ras)

#----------------------------------------------------
# Localização da área de estudo
#----------------------------------------------------
# brasil <-
#   shapefile(filename="/home/lgcs-mds/rdata/dnos_sm/shapefiles/brasil.shp",
#             verbose=TRUE)
# proj4string(brasil) <- proj4string(dnos.coords)
# summary(brasil@data)
# str(brasil@data$UF_05)
# class(brasil@data$UF_05)
# 
# sm.N = char2dms("29d41'00\"S")
# sm.E = char2dms("53d48'00\"W")
# sm = data.frame(name='sm', x=as.numeric(sm.E), y=as.numeric(sm.N))
# coordinates(sm) = c('x', 'y')
# proj4string(sm) = proj4string(brasil)
# class(sm)

# Plot----------------------------------------------
# dev.off()
# pts <- list("sp.points", sm, pch=20, col = "black")
# arrow <- list("SpatialPolygonsRescale",
#               layout.north.arrow(type=1),
#               offset = c(-32, 0), scale = 3)
# png(filename = "/home/lgcs-mds/rdata/dnos_sm/figures/location.png")
# print(spplot(brasil, zcol = "UF_05", aspect = "iso",
#        colorkey = FALSE, scales = list(draw = TRUE),
#        sp.layout = list(pts, arrow),
#        xlab = "E", ylab = "N",
#        panel = function(x,y, ...) {
#          panel.polygonsplot(x, y, ...);
#          panel.grid(h=-1,v=-1,
#                     col="lightgray",
#                     lty=2)}))
# dev.off()
# rm(pts, arrow, sm.E, sm.N, sm)

#####################################################
### Visualização geral dos dados
#####################################################
#----------------------------------------------------
# Plotar as observações de calibração
#----------------------------------------------------
# dev.off()
# pdf(file="/home/alessandro/rdata/dnos_sm/figures/sample_points.pdf",
#     width=8, height=8)
# plot(dnos.limite, axes=TRUE, asp=1,
#      xlim=bbox(dnos.coords)[1,],
#      ylim=bbox(dnos.coords)[2,])
# title(main="Location of the sample points",
#       sub="Black: calibration (340); Red: validation (60)",
#       xlab="E (m)", ylab="N (m)")
# points(dnos.coords, col="black", pch=20)
# points(dnos.validation, col="red", pch=20)
# grid(col="gray", lty=2)
# dev.off()

# SAVE DATA ####################################################################
ls()
save(bc_lambda, density_model,
     cal_data, val_data,
     cellsize, data.extent, grass.addons,
     GRASSgisDbase, 
     # directories
     bounding.dir, covar.validation.dir, explora.dir, dem.dir, geo.dir,
     hydro.dir, land.dir, point_dir, rdata.dir, resample.dir, update_covars_dir,
     soil.dir, ppp_dir,
     # coordinate reference systems
     ca_utm22s, sirgas2000, sirgas2000utm22s, wgs1984, wgs1984utm22s,
     file = "sm-dnos-general.RData")
# End!