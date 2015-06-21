# DESCRIPTION ##################################################################
# General configurations and pre-processing. Here we define variables that will
# be used throughout the entire project. They identify data directories and 
# coordinate reference systems. The prediction grid is also defined here, as
# well as the GRASS GIS database.

# SETTINGS #####################################################################
rm(list = ls())
gc()
require(sp)
require(rgdal)
require(raster)
require(spgrass6)

# SET DATA DIRECTORIES AND DEFINITIONS #########################################

# Data directories  ------------------------------------------------------------
boundary_data <- path.expand("~/projects/dnos-sm-rs/data/boundaries/")
geology_data <- path.expand("~/projects/dnos-sm-rs/data/geology/")
hydrology_data <- path.expand("~/projects/dnos-sm-rs/data/hydrology/")
landuse_data <- path.expand("~/projects/dnos-sm-rs/data/landuse/")
point_data <- path.expand("~/projects/dnos-sm-rs/data/point/")
r_data <- path.expand("~/projects/dnos-sm-rs/data/R/")
soilclass_data <- path.expand("~/projects/dnos-sm-rs/data/soilclass/")
terrain_data <- path.expand("~/projects/dnos-sm-rs/data/terrain/")

# Prediction grid --------------------------------------------------------------

# Cell size of the prediction grid
# We choose the cellsize to be 5 metres because this is the smallest
# cellsize among the covariates (RapidEye imagery). This way we do not lose any
# information due to the aggregation of spatial data.
cellsize <- 5

# Extent of the region in which there is data available
data_extent <- extent(c(220817, 237524, 6711252, 6724461))

# Load point soil data (field data)
file <- paste(point_data, "fieldData.csv", sep = "")
field_data <- read.table(file, dec = ".", head = TRUE, sep = ";",
                         stringsAsFactors = FALSE, na.strings = "na")
colnames(field_data)
field_data <- field_data[1:340, 3:4] # do not load compiled and validation data
coordinates(field_data) <- ~ longitude + latitude
proj4string(field_data) <- sirgas2000
field_data <- spTransform(field_data, wgs1984utm22s)
str(field_data)
plot(field_data)

# Extent of the prediction grid
# set bounding box
grow <- 1500
bb <- bbox(field_data)
min_x <- floor((bb["longitude", "min"] - grow) / cellsize) * cellsize
max_x <- ceiling((bb["longitude", "max"] + grow) / cellsize) * cellsize
min_y <- floor((bb["latitude", "min"] - grow) / cellsize) * cellsize
max_y <- ceiling((bb["latitude", "max"] + grow) / cellsize) * cellsize

# calculate the number of grid cells
cells_x <- (max_x - min_x) / cellsize
cells_y <- (max_y - min_y) / cellsize
cells <- cells_x * cells_y

# create prediction grid
dnos_raster <- GridTopology(cellcentre.offset = c(min_x + (cellsize / 2),
                                                  min_y + (cellsize / 2)),
                            cellsize = c(cellsize, cellsize), 
                            cells.dim = c(cells_x, cells_y))
dnos_raster <- SpatialGrid(dnos_raster, proj4string = wgs1984utm22s)
str(dnos_raster)
plot(extent(dnos_raster), asp = 1)
plot(field_data, add = TRUE)
dnos_raster <- SpatialGridDataFrame(dnos_raster, 
                                    data.frame(tmp = rep(1, cells)))
str(dnos_raster)
rm(min_x, max_x, min_y, max_y, cells_x, cells_y, cells, grow)

# bounding box of the dnos_raster
tmp <- as(extent(dnos_raster), "SpatialPolygons")
proj4string(tmp) <- wgs1984utm22s
shapefile(tmp, paste(boundary_data, "dnos-raster-bbox2.shp", sep = ""))
rm(tmp)

# # Create GRASS GIS DBASE structure -------------------------------------------
dbGRASS <- path.expand("~/dbGRASS")
initGRASS(gisBase = "/usr/lib/grass64/", gisDbase = dbGRASS,
          location = "dnos-sm-rs", mapset = "predictions", pid = Sys.getpid())
writeRAST6(dnos_raster, "dnos.raster")
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
save(boundary_data, cellsize, dbGRASS, terrain_data, geology_data,
     hydrology_data, landuse_data, point_data, r_data, soilclass_data,
     sirgas2000, sirgas2000utm22s, wgs1984, wgs1984utm22s, ca_utm22s,
     file = paste(r_data, "general.RData", sep = ""))






# LOAD AND PROCESS DATA ########################################################
# This is going to be moved to another script!

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

