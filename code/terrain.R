# DESCRIPTION ##################################################################
# Preparation of the digital elevation models (DEM) and calculation of their
# derivatives.
# Two DEMs are used, the SRTM, and the product of the interpolation of contour 
# lines and stream data.
# The references are:
# - Jarvis, A.; Reuter, H. I.; Nelson, A. & Guevara, E. Hole-filled SRTM for 
#   the globe version 4. 2008.
# - DSG. Santa Maria - SE. Sheet SH.22-V-C-IV/1-SE. Brasília: Ministério do
#   Exército, Departamento de Engenharia e Comunicações, Diretoria do Serviço
#   Geográfico do Exército, 1992.
# - DSG. Santa Maria - NE. Sheet SH.22-V-C-IV-1-NE. Brasília: Ministério do
#   Exército, Departamento de Engenharia e Comunicações, Diretoria do Serviço
#   Geográfico do Exército, 1992.
# - DSG. Camobi - SO. Sheet SH.22-V-C-IV/2-SO. Brasília: Ministério do Exército,
#   Departamento de Engenharia e Comunicações, Diretoria do Serviço Geográfico
#   do Exército, 1980.
# SETTINGS #####################################################################
rm(list = ls())
gc()
# Load packages
require(lattice)
require(sp)
require(rgdal)
require(maptools)
require(raster)
require(RSAGA)
require(spgrass6)
require(vec2dtransf)
require(latticeExtra)
# Load data
load("sm-dnos-general.RData")
load("sm-dnos-covar-validation.RData")
load("sm-dnos-terrain.RData")
ls()
# GRASS GIS definitions
initGRASS(gisBase = "/usr/lib/grass64/", gisDbase = GRASSgisDbase, 
          location = "dnos-sm-rs", mapset = "predictions", 
          pid = Sys.getpid(), override = TRUE)
system("g.region rast=dnos.raster")

# DIGITAL ELEVATION MODELS #####################################################

# Topographic map - 1:25,000 ===================================================

# load buffer, contours, points, rivers, and lakes
buffer25 <- shapefile(paste(bounding.dir, "buffer-bacia-carta.shp", sep = ""))
contours25 <- shapefile(paste(dem.dir, "contour-lines25.shp", sep = ""))
contours25$level <- as.numeric(contours25$level)
points25 <- shapefile(paste(dem.dir, "points25.shp", sep = ""))
points25$level <- as.numeric(points25$level)
lakes25 <- shapefile(paste(hydro.dir, "lakes25.shp", sep = ""))
rivers25 <- shapefile(paste(hydro.dir, "drainage25.shp", sep = ""))

# transform to WGS84
buffer25 <- spTransform(buffer25, wgs1984utm22s)
contours25 <- spTransform(contours25, wgs1984utm22s)
points25 <- spTransform(points25, wgs1984utm22s)
lakes25 <- spTransform(lakes25, wgs1984utm22s)
rivers25 <- spTransform(rivers25, wgs1984utm22s)

# apply affine transformation
buffer25 <- applyTransformation(affine.car25, buffer25)
contours25 <- applyTransformation(affine.car25, contours25)
points25 <- applyTransformation(affine.car25, points25)
lakes25 <- applyTransformation(affine.car25, lakes25)
rivers25 <- applyTransformation(affine.car25, rivers25)

# create reference grid (based on prediction grid - 5 m)
# this is more usefull when processing in GRASS GIS
grow <- 150
min.x <- floor((bbox(buffer25)["x", "min"] - grow)/cellsize)*cellsize
max.x <- ceiling((bbox(buffer25)["x", "max"] + grow)/cellsize)*cellsize
min.y <- floor((bbox(buffer25)["y", "min"] - grow)/cellsize)*cellsize
max.y <- ceiling((bbox(buffer25)["y", "max"] + grow)/cellsize)*cellsize
cells.x <- (max.x - min.x)/cellsize
cells.y <- (max.y - min.y)/cellsize
cells <- cells.x*cells.y
grid.carta25 <- SpatialGrid(GridTopology(cellcentre.offset = c(min.x+(cellsize/2),
                                                               min.y+(cellsize/2)),
                                         cellsize = c(cellsize, cellsize),
                                         cells.dim = c(cells.x, cells.y)),
                            proj4string = wgs1984utm22s)
str(grid.carta25)
plot(extent(grid.carta25), asp = 1)
plot(buffer25, add = TRUE)
grid.carta25 <- SpatialGridDataFrame(grid.carta25, data.frame(tmp = rep(1, cells)))
str(grid.carta25)
rm(min.x, max.x, min.y, max.y, cells.x, cells.y, cells, grow)

# export shapefiles 
shapefile(contours25, paste(dem.dir, "contours25-affine.shp", sep = ""), overwrite = TRUE)
shapefile(points25, paste(dem.dir, "points25-affine.shp", sep = ""))
shapefile(lakes25, paste(dem.dir, "lakes25-affine.shp", sep = ""))
shapefile(rivers25, paste(dem.dir, "rivers25-affine.shp", sep = ""))
tmp <- as(extent(grid.carta25), "SpatialPolygons")
shapefile(tmp, paste(dem.dir, "extent25-affine.shp", sep = ""))

# interpolate DEM using ANUDEM in ArcGIS
# this is to use the full potential of contour lines, points, lakes and steam data
# ...

# import into GRASS GIS mapset "predicitions"
initGRASS(gisBase = "/usr/lib/grass64/", gisDbase = GRASSgisDbase, override = TRUE,
          location = "dnos-sm-rs", mapset = "predictions", pid = Sys.getpid())
system("g.region rast=dnos.raster")
system(paste("r.in.gdal --o in=", dem.dir, "dem-carta25.tif out=dem.carta25", sep = ""))
system("r.info -r dem.carta25")

# remove stair-like artifacts
system("r.neighbors --o in=dem.carta25 out=dem.carta25.smooth method=average size=7")
system(" r.mapcalc 'tmp=dem.carta25-dem.carta25.smooth'")
system("r.info -r tmp")
system("d.mon x0");system("d.rast.leg tmp")
system("d.mon x1");system("d.histogram tmp")

# calculate geoidal heights (MAPGEO2010 and EGM69) (Hegm = Hmapgeo + Nmapgeo - Negm)
# coords <- coordinates(spTransform(cal.data, wgs1984))
# write.table(data.frame(coords[,2], coords[,1]), paste(dem.dir, "INPUT.DAT", sep = ""),
#             col.names = FALSE, row.names = FALSE)
# rm(coords)
# setwd(dem.dir)
# system("wine F477")
# setwd(rdata.dir)
# mapgeo <- read.table(paste(dem.dir, "INPUT_saida.txt", sep = ""))[,1]
# egm <- read.table(paste(dem.dir, "OUTF477.DAT", sep = ""))[,3]
# quantile(mapgeo);quantile(egm);quantile(c(mapgeo-egm))
# rm(mapgeo, egm)
system("r.mapcalc dem.carta25.egm=dem.carta25.smooth+1")
system("r.info -r dem.carta25.egm")
system("nviz elev=dem.carta25.egm")

# round values to 2 digits
system("r.mapcalc 'ELEV_10=(double(round(dem.carta25.egm*100))/100)'")
system("r.info -r ELEV_10")

# export as SAGA grid (sgrd)
system(paste("r.out.gdal input=ELEV_10 format=SAGA output=", dem.dir, "ELEV_10.sdat", sep = ""))
system(paste("saga_cmd pj_proj4 0 -CRS_METHOD=0 -PRECISE -CRS_PROJ4='", wgs1984utm22s,
             "' -GRIDS=", dem.dir, "ELEV_10.sgrd", sep = ""))

# SRTM v41 (90 m) ==============================================================

# Warp (WGS 1984 / UTM  zone 22 S) and clip
system(paste("gdalwarp -s_srs epsg:4326 -t_srs epsg:32722 -te ",
             data.extent@xmin, " ", data.extent@ymin, " ", data.extent@xmax, " ",
             data.extent@ymax, " -r cubic -of GTiff -overwrite ", dem.dir,
             "srtm-v41.tif ", dem.dir, "srtm-dem90.tif", sep = ""))
system(paste("gdalinfo ", dem.dir, "srtm-dem90.tif", sep = ""))

# import to GRASS GIS (mapset = "predictions")
initGRASS(gisBase = "/usr/lib/grass64/", gisDbase = GRASSgisDbase, location = "dnos-sm-rs",
          mapset = "predictions", pid = Sys.getpid(), override = TRUE)
system("g.region rast=dnos.raster")
system(paste("r.in.gdal --o input=", dem.dir, "srtm-dem90.tif output=srtm.dem90", sep = ""))
system("r.info srtm.dem90")

# resample to 5 m (dnos.raster)
system(paste("r.resamp.interp --o method=bicubic input=srtm.dem90 output=srtm.dem15"))
system(paste("r.info srtm.dem15"))

# fill sinks
system("r.fill.dir --o input=srtm.dem15 elevation=ELEV_90 direction=FDIR_90 type=answers")
system("r.info ELEV_90")
system("r.mapcalc DIFF_90=ELEV_90-srtm.dem15");system("r.info DIFF_90")
system("d.mon x0"); system("d.rast DIFF_90")
system("d.rast FDIR_90")
system("nviz elev=ELEV_90")

# round values to 2 digits
system("r.mapcalc 'ELEV_90=(double(round(ELEV_90*100))/100)'")
system("r.info -r ELEV_90")

# export as SAGA grid (sgrd)
system("r.info ELEV_90")
system(paste("r.out.gdal input=ELEV_90 format=SAGA output=", dem.dir, "ELEV_90.sdat", sep = ""))
system(paste("saga_cmd pj_proj4 0 -CRS_METHOD=0 -PRECISE -CRS_PROJ4='", wgs1984utm22s,
             "' -GRIDS=", dem.dir, "ELEV_90.sgrd", sep = ""))

# TOPODATA (30 m) =====================================================================

# Warp (to WGS 1984 / UTM zone 22 S) and clip
system(paste("gdalwarp -s_srs epsg:4326 -t_srs epsg:32722 -te ",
             data.extent@xmin, " ", data.extent@ymin, " ", data.extent@xmax,
             " ", data.extent@ymax, " -r cubic -of GTiff -overwrite ", dem.dir,
             "topodata.tif ", dem.dir, "topodata-dem30.tif", sep = ""))
system(paste("gdalinfo ", dem.dir, "topodata-dem30.tif", sep = ""))

# import to GRASS GIS (mapset = "predictions")
initGRASS(gisBase = "/usr/lib/grass64/", gisDbase = GRASSgisDbase, location = "dnos-sm-rs",
          mapset = "predictions", pid = Sys.getpid(), override = TRUE)
system("g.region rast=dnos.raster")
system(paste("r.in.gdal input=", dem.dir, "topodata-dem30.tif output=topodata.dem30", sep = ""))
system(paste("r.info topodata.dem30"))

# resample to 5 m (dnos.raster)
system("r.resamp.interp --o method=bicubic input=topodata.dem30 output=topodata.dem15")
system("r.info topodata.dem15")

# fill sinks
system("r.fill.dir --o input=topodata.dem15 elevation=ELEV_30 direction=FDIR_30 type=answers")
system("r.info ELEV_30")
system("r.mapcalc DIFF_30=ELEV_30-topodata.dem15");system("r.info DIFF_30")
system("d.mon x0");system("d.rast DIFF_30")
system("d.rast FDIR_30")
system("nviz elev=ELEV_30")

# round values to 2 digits
system("r.mapcalc 'ELEV_30=(double(round(ELEV_30*100))/100)'")
system("r.info -r ELEV_30")

# Export as SAGA grid (sgrd)
system("r.info ELEV_30")
system(paste("r.out.gdal input=ELEV_30 format=SAGA output=", dem.dir, "ELEV_30.sdat", sep = ""))
system(paste("saga_cmd pj_proj4 0 -CRS_METHOD=0 -PRECISE -CRS_PROJ4='", wgs1984utm22s,
             "' -GRIDS=", dem.dir, "ELEV_30.sgrd", sep = ""))

# DERIVE TERRAIN ATTRIBUTES ####################################################

# ELEV_90 ======================================================================

# ELEV_90 - ELEVATION ----------------------------------------------------------
# meters, two decimal places, used to calculate DEM derivatives and modelling
system("r.info -r ELEV_90")
system("r.mapcalc 'ELEV_90=(double(round(ELEV_90 * 100)) * 0.01)'")
system("r.colors map=ELEV_90 color=haxby")
system("d.mon x3")
system("d.rast.leg ELEV_90")
system("d.mon x4")
system("d.histogram ELEV_90")

# ELEV_90 - SLOPE --------------------------------------------------------------
# degrees, two decimal places, used to calculate DEM derivatives and modelling
# system("r.param.scale -help")
size <- c(3, 7, 15, 31, 63, 127, 255)
cmd <- paste("r.multi.scale --o input=ELEV_90 output=SLP_90_", size, 
             " s_tol=0.1 size=", size, " param=slope", sep = "")
lapply(cmd, system)
cmd <- paste("r.info -r SLP_90_", size, sep = "")
lapply(cmd, system)
# correct flat surfaces
cmd <- paste("r.mapcalc 'SLP_90_", size, "=SLP_90_", size, "+0.1'", sep = "")
lapply(cmd, system)
cmd <- paste("r.colors map=SLP_90_", size, " color=slope", sep = "")
lapply(cmd, system)
system("d.mon x0");system("d.rast.leg SLP_90_3")
system("r.info -r SLP_90_3")
# round to 2 decimal places
cmd <- paste("r.mapcalc 'SLP_90_", size, "=(double(round(SLP_90_", size, 
             "*100))/100)'", sep = "")
lapply(cmd, system)
cmd <- paste("r.info -r SLP_90_", size, sep = "")
lapply(cmd, system)

# ELEV_90 - ASPECT -------------------------------------------------------------
# degrees, two decimal places, used to calculate DEM derivatives and modelling
# system("r.param.scale -help")
system("r.multi.scale -help")
size <- c(3, 7, 15, 31, 63, 127, 255)
cmd <- paste("r.multi.scale --o input=ELEV_90 output=ASP_90_", size,
             " s_tol=0.1 size=", size, " param=aspect", sep = "")
lapply(cmd, system)
system("d.mon x0");system("d.histogram ASP_90_255")
system("d.mon x1");system("d.rast.leg ASP_90_255")
system("r.info -r ASP_90_3")
# convert aspect values to 0-360: if(asp < 0, asp + 360, asp)
cmd <- paste("r.mapcalc 'ASP_90_", size, "=if(ASP_90_", size, "<0,ASP_90_", size,
             "+360,ASP_90_", size, ")'", sep = "")
lapply(cmd, system)
system("r.info -r ASP_90_3")
# rotate aspect values: if(asp < 90, asp + 270, asp - 90)
cmd <- paste("r.mapcalc 'ASP_90_", size, "=if(ASP_90_", size, "<90,ASP_90_", size,
             "+270,ASP_90_", size, "-90)'", sep = "")
lapply(cmd, system)
system("r.info -r ASP_90_3")

# ELEV_90 - NORTHERNNESS -------------------------------------------------------
# abs(180-asp)
# degrees, integers, used for modelling
cmd <- paste("r.mapcalc 'NOR_90_", 
             size, "=abs(180-ASP_90_", size, ")'", sep = "")
lapply(cmd, system)
# round to 2 decimal places
cmd <- paste("r.mapcalc 'NOR_90_", size,
             "=double(round(NOR_90_", size, "*100))*0.01'", sep = "")
lapply(cmd, system)
cmd <- paste("r.info -r NOR_90_", size, sep = "")
lapply(cmd, system)
# set color
cmd <- paste("r.colors map=NOR_90_", size, " color=grey", sep = "")
lapply(cmd, system)
system("d.mon x0");system("d.histogram NOR_90_3")
system("d.mon x1");system("d.rast.leg NOR_90_3")

# ELEV_90 - FLOW ACCUMULATION --------------------------------------------------
# square meters, integers, used to calculate DEM derivatives
system("r.watershed -help")
system("r.watershed --o -a -f --v elevation=ELEV_90 accumulation=ACC_90")
system(paste("r.mapcalc 'ACC_90=ACC_90*(", cellsize*cellsize, ")'", sep = ""))
# round
system("r.mapcalc 'ACC_90=round(ACC_90)'")
system("r.info ACC_90")
system("d.mon x2")
system("d.rast.leg ACC_90")

# ELEV_90 - TOPOGRAPHIC WETNESS INDEX ------------------------------------------
# ln(a/tan(beta))
# unitless, two decimal places, used for modelling
size <- c(3, 7, 15, 31, 63, 127, 255)
cmd <- paste("r.mapcalc 'TWI_90_", size, "=log((ACC_90/", cellsize, 
             ")/tan(SLP_90_", size, "))'", sep = "")
lapply(cmd, system)
system("r.info -r TWI_90_3")
# round to 2 decimal places
cmd <- paste("r.mapcalc 'TWI_90_", size, "=double(round(TWI_90_", size, 
             "*100))*0.01'", sep = "")
lapply(cmd, system)
cmd <- paste("r.info -r TWI_90_", size, sep = "")
lapply(cmd, system)
system("d.mon x0");system("d.rast.leg TWI_90_255")

# ELEV_90 - STREAM POWER INDEX -------------------------------------------------
# ln(a*tan(beta))
# unitless, two decimal places, used for modelling
size <- c(3, 7, 15, 31, 63, 127, 255)
cmd <- paste("r.mapcalc 'SPI_90_", size, "=log((ACC_90/", cellsize, 
             ")*tan(SLP_90_", size, "))'", sep = "")
lapply(cmd, system)
# round to 2 decimal places
cmd <- paste("r.mapcalc 'SPI_90_", size, "=double(round(SPI_90_", size,
             "*100))*0.01'", sep = "")
lapply(cmd, system)
cmd <- paste("r.info -r SPI_90_", size, sep = "")
lapply(cmd, system)
system("d.mon x0")
system("d.rast.leg SPI_90_3")
system("d.mon x1")
system("d.histogram SPI_90_3")

# ELEV_90 - TOPOGRAPHIC POSITION INDEX -----------------------------------------
# meters, two decimal places, used for modelling
size <- c(3, 7, 15, 31, 63, 127, 255)
cmd <- paste("saga_cmd ta_morphometry 18 -DEM=", dem.dir, "ELEV_90.sgrd -TPI=",
             dem.dir, "TPI_90_", size, ".sgrd -RADIUS_MIN=3 -RADIUS_MAX=",
             (size*cellsize)/2, sep = "")
lapply(cmd, system)
# export to tif with CRS data
cmd <- paste("gdal_translate -a_srs epsg:32722 ", dem.dir, "TPI_90_", size,
             ".sdat ", dem.dir, "TPI_90_", size, ".tif", sep = "")
lapply(cmd, system)
# import into GRASS
cmd <- paste("r.in.gdal --o in=", dem.dir, "TPI_90_", size, ".tif out=TPI_90_",
             size, sep = "")
lapply(cmd, system)
# round to 2 decimal places
cmd <- paste("r.mapcalc 'TPI_90_", size, "=(double(round(TPI_90_", size,
             "*100))/100)'", sep = "")
lapply(cmd, system)
cmd <- paste("r.info -r TPI_90_", size, sep = "")
lapply(cmd, system)
system("d.mon x0")
system("d.rast.leg TPI_90_127")
system("d.mon x1")
system("d.histogram TPI_90_3")

# ELEV_90 - MORPHOMETRIC PROTECTION INDEX --------------------------------------
size <- c(3, 7, 15, 31, 63, 127, 255)
system("saga_cmd ta_morphometry 7")
cmd <- paste("saga_cmd ta_morphometry 7 -DEM=", dem.dir, 
             "ELEV_90.sgrd -PROTECTION=", dem.dir, "MPI_90_", size, 
             ".sgrd -RADIUS=", size, sep = "")
lapply(cmd, system)
cmd <- paste("gdal_translate -a_srs epsg:32722 ", dem.dir, "MPI_90_", size,
             ".sdat ", dem.dir, "MPI_90_", size, ".tif", sep = "")
lapply(cmd, system)
cmd <- paste("r.in.gdal --o in=", dem.dir, "MPI_90_", size, ".tif out=MPI_90_", 
             size, sep = "")
lapply(cmd, system)
cmd <- paste("r.info -r MPI_90_", size, sep = "")
lapply(cmd, system)
cmd <- paste("r.mapcalc 'MPI_90_", size, "=(double(round(MPI_90_", size,
             "*100))/100)'", sep = "")
lapply(cmd, system)
system("d.mon x0")
system("d.erase")
system("d.rast.leg MPI_90_255")
system("d.erase")
system("d.histogram MPI_90_255")

# ELEV_90 - PHYSIOGRAPHIC REGIONS ----------------------------------------------
# Decision tree rules
# tpi > 30: 'hillslope'
# tpi < 30: 'other'
### elev < 200 : 'other'
##### mpi < 0.1: 'depression'
##### mpi > 0.1: 'hillslope'
### elev > 200 : 'other'
##### elev > 400: 'other'       
####### mpi > 0.1: 'hillslope'
####### mpi < 0.1: 'plateau'
##### elev < 400: 'hillslope'

size <- c(3, 7, 15, 31, 63, 127, 255)
system("g.list group")
cmd <- paste("i.group group=tmpgroup", size, " subgroup=tmpgroup", size, 
             " input=MPI_90_", size, ",ELEV_90,TPI_90_", size, sep = "")
lapply(cmd, system)
cmd <- paste("i.cluster tmpgroup", size, " subgroup=tmpgroup", size, 
             " sigfile=cluster", size, " classes=3", sep = "")
lapply(cmd, system)
cmd <- paste("i.maxlik --o group=tmpgroup", size, " subgroup=tmpgroup", size, 
             " sigfile=cluster", size, " class=tmprast", size, sep = "")
lapply(cmd, system)
system("d.mon x0")
system("d.erase")
system("d.rast.leg tmprast255")
cmd <- paste("g.remove group=tmpgroup", size, sep = "")
lapply(cmd, system)
# majority filter
system("r.neighbors --o input=tmprast3 output=STR_90_3 method=mode size=15")
system("r.neighbors --o input=tmprast7 output=STR_90_7 method=mode size=15")
system("r.neighbors --o input=tmprast15 output=STR_90_15 method=mode size=15")
system("r.neighbors --o input=tmprast31 output=STR_90_31 method=mode size=15")
system("r.neighbors --o input=tmprast63 output=STR_90_63 method=mode size=15")
system("r.neighbors --o input=tmprast127 output=STR_90_127 method=mode size=15")
system("r.neighbors --o input=tmprast255 output=STR_90_255 method=mode size=15")

system("d.mon x0")
system("d.erase")
system("d.rast.leg STR_90_255")

cmd <- paste("g.remove rast=tmprast", size, sep = "")
lapply(cmd, system)

# ELEV_10 ======================================================================

# ELEV_10 - ELEVATION ----------------------------------------------------------
# meters, two decimal places, used to calculate DEM derivatives and modelling
system("r.info -r ELEV_10")
system("r.colors map=ELEV_10 color=haxby")
system("d.mon x0")
system("d.rast.leg ELEV_10")
system("d.mon x1")
system("d.histogram ELEV_10")

# ELEV_10 - SLOPE --------------------------------------------------------------
# degrees, two decimal places, used to calculate DEM derivatives and modelling
# system("r.param.scale -help")
system("r.multi.scale -help")
size <- c(3, 7, 15, 31, 63, 127, 255)
cmd <- paste("r.multi.scale --o input=ELEV_10 output=SLP_10_", size, 
             " s_tol=0.1 size=", size, " param=slope", sep = "")
lapply(cmd, system)
# correct flat surfaces
cmd <- paste("r.mapcalc 'SLP_10_", size, "=SLP_10_", size, "+0.1'", sep = "")
lapply(cmd, system)
system("d.mon x0");system("d.rast.leg SLP_10_3")
system("r.info -r SLP_10_3")
# round to 2 decimal places
cmd <- paste("r.mapcalc 'SLP_10_", size, "=(double(round(SLP_90_", size, 
             "*100))/100)'", sep = "")
lapply(cmd, system)
cmd <- paste("r.info -r SLP_10_", size, sep = "")
lapply(cmd, system)

# ELEV_10 - ASPECT -------------------------------------------------------------
# degrees, two decimal places, used to calculate DEM derivatives and modelling
# system("r.param.scale -help")
system("r.multi.scale -help")
size <- c(3, 7, 15, 31, 63, 127, 255)
cmd <- paste("r.multi.scale --o input=ELEV_10 output=ASP_10_", size, 
             " s_tol=0.1 size=", size, " param=aspect", sep = "")
lapply(cmd, system)
system("d.mon x1");system("d.rast.leg ASP_10_15")
system("r.info -r ASP_10_3")
# convert aspect values to 0-360: if(asp < 0, asp + 360, asp)
cmd <- paste("r.mapcalc 'ASP_10_", size, "=if(ASP_10_", size, "<0,ASP_10_",
             size, "+360,ASP_10_", size, ")'", sep = "")
lapply(cmd, system)
system("r.info -r ASP_10_3")
# rotate aspect values: if(asp < 90, asp + 270, asp - 90)
cmd <- paste("r.mapcalc 'ASP_10_", size, "=if(ASP_10_", size, "<90,ASP_10_",
             size, "+270,ASP_10_", size, "-90)'", sep = "")
lapply(cmd, system)
system("r.info -r ASP_10_3")

# ELEV_10 - NORTHERNNESS -------------------------------------------------------
# abs(180-asp)
# degrees, integers, used for modelling
cmd <- paste("r.mapcalc 'NOR_10_", size, "=abs(180-ASP_10_", size, ")'", 
             sep = "")
lapply(cmd, system)
# round to 2 decimal places
cmd <- paste("r.mapcalc 'NOR_10_", size, "=double(round(NOR_10_", size, 
             "*100))*0.01'", sep = "")
lapply(cmd, system)
cmd <- paste("r.info -r NOR_10_", size, sep = "")
lapply(cmd, system)
# set color
cmd <- paste("r.colors map=NOR_10_", size, " color=grey", sep = "")
lapply(cmd, system)
system("d.mon x0")
system("d.histogram NOR_10_3")
system("d.mon x1")
system("d.rast.leg NOR_10_3")

# ELEV_10 - FLOW ACCUMULATION --------------------------------------------------
# square meters, integers, used to calculate DEM derivatives
# system("r.watershed -help")
system("r.watershed --o -a -f --v elevation=ELEV_10 depression=LAKE_25 accumulation=ACC_10")
system(paste("r.mapcalc 'ACC_10=ACC_10*(", cellsize*cellsize, ")'", sep = ""))
# round
system("r.mapcalc 'ACC_10=round(ACC_10)'")
system("r.info ACC_10")
system("d.mon x1")
system("d.rast.leg ACC_10")

# ELEV_10 - TOPOGRAPHIC WETNESS INDEX ------------------------------------------
# ln(a/tan(beta))
# unitless, two decimal places, used for modelling
size <- c(3, 7, 15, 31, 63, 127, 255)
cmd <- paste("r.mapcalc 'TWI_10_", size, "=log((ACC_10/", cellsize,
             ")/tan(SLP_10_", size, "))'", sep = "")
lapply(cmd, system)
system("r.info -r TWI_10_3")
# round to 2 decimal places
cmd <- paste("r.mapcalc 'TWI_10_", size, "=double(round(TWI_10_", size, 
             "*100))*0.01'", sep = "")
lapply(cmd, system)
cmd <- paste("r.info -r TWI_10_", size, sep = "")
lapply(cmd, system)
system("d.mon x0")
system("d.rast.leg TWI_10_3")

# ELEV_10 - STREAM POWER INDEX -------------------------------------------------
# ln(a*tan(beta))
# unitless, two decimal places, used for modelling
size <- c(3, 7, 15, 31, 63, 127, 255)
cmd <- paste("r.mapcalc 'SPI_10_", size, "=log((ACC_10/", cellsize, 
             ")*tan(SLP_10_", size, "))'", sep = "")
lapply(cmd, system)
# round to 2 decimal places
cmd <- paste("r.mapcalc 'SPI_10_", size, "=double(round(SPI_10_", size, 
             "*100))*0.01'", sep = "")
lapply(cmd, system)
cmd <- paste("r.info -r SPI_10_", size, sep = "")
lapply(cmd, system)
system("d.mon x0")
system("d.rast.leg SPI_10_3")
system("d.mon x1")
system("d.histogram SPI_10_3")

# ELEV_10 - TOPOGRAPHIC POSITION INDEX -----------------------------------------
# meters, two decimal places, used for modelling
size <- c(3, 7, 15, 31, 63, 127, 255)
cmd <- paste("saga_cmd ta_morphometry 18 -DEM=", dem.dir, "ELEV_10.sgrd -TPI=",
             dem.dir, "TPI_10_", size, ".sgrd -RADIUS_MIN=3 -RADIUS_MAX=",
             (size*cellsize)/2, sep = "")
lapply(cmd, system)
# export to tif with CRS data
cmd <- paste("gdal_translate -a_srs epsg:32722 ", dem.dir, "TPI_10_", size,
             ".sdat ", dem.dir, "TPI_10_", size, ".tif", sep = "")
lapply(cmd, system)
# import into GRASS
cmd <- paste("r.in.gdal --o in=", dem.dir, "TPI_10_", size, ".tif out=TPI_10_", 
             size, sep = "")
lapply(cmd, system)
# round to 2 decimal places
cmd <- paste("r.mapcalc 'TPI_10_", size, "=(double(round(TPI_10_", size,
             "*100))/100)'", sep = "")
lapply(cmd, system)
cmd <- paste("r.info -r TPI_10_", size, sep = "")
lapply(cmd, system)
system("d.mon x0")
system("d.rast.leg TPI_10_3")
system("d.mon x0")
system("d.histogram TPI_10_3")

# ELEV_10 - MORPHOMETRIC PROTECTION INDEX --------------------------------------
size <- c(3, 7, 15, 31, 63, 127, 255)
system("saga_cmd ta_morphometry 7")
cmd <- paste("saga_cmd ta_morphometry 7 -DEM=", dem.dir, 
             "ELEV_10.sgrd -PROTECTION=", dem.dir, "MPI_10_", size, 
             ".sgrd -RADIUS=", size, sep = "")
lapply(cmd, system)
cmd <- paste("gdal_translate -a_srs epsg:32722 ", dem.dir, "MPI_10_", size,
             ".sdat ", dem.dir, "MPI_10_", size, ".tif", sep = "")
lapply(cmd, system)
cmd <- paste("r.in.gdal --o in=", dem.dir, "MPI_10_", size, ".tif out=MPI_10_", 
             size, sep = "")
lapply(cmd, system)
cmd <- paste("r.info -r MPI_10_", size, sep = "")
lapply(cmd, system)
cmd <- paste("r.mapcalc 'MPI_10_", size, "=(double(round(MPI_10_", size,
             "*100))/100)'", sep = "")
lapply(cmd, system)
system("d.mon x0")
system("d.erase")
system("d.rast.leg MPI_10_255")
system("d.erase")
system("d.histogram MPI_10_255")

# ELEV_10 - PHYSIOGRAPHIC REGIONS ----------------------------------------------

# Decision tree rules
# tpi > 30: 'hillslope'
# tpi < 30: 'other'
### elev < 200 : 'other'
##### mpi < 0.1: 'depression'
##### mpi > 0.1: 'hillslope'
### elev > 200 : 'other'
##### elev > 400: 'other'       
####### mpi > 0.1: 'hillslope'
####### mpi < 0.1: 'plateau'
##### elev < 400: 'hillslope'

# Reclassify values ==================================================

#strata90 <- raster("STRATA_90.sdat")
#strata90tab <- read.table("STRATA_90.txt", header = TRUE)
#strata90tab
#levels(factor(values(tmp)))
#plot(strata90)
#rcl <- cbind(strata90tab$MINIMUM, c(1, 2, 2, 3, 2, 2))
#strata90 <- reclassify(strata90, rcl, filename = "tmp",
#                       format = "SAGA", overwrite = TRUE)
#plot(strata90)
#plot(raster("tmp.tif"))
#rm(rcl)
# Filter strata ======================================================
#rsaga.get.usage(lib = "grid_filter", module = "Majority Filter",
#                env = sagaenv)

strata90 <- raster("STRATA_90.sdat")
plot(strata90)

# DELINEATE CATCHMENTS AND STREAM NETWORK ######################################

# ELEV_90 ======================================================================

# streams
system("r.watershed")
system("r.watershed --o -f -a elev=ELEV_90 stream=STREAM_90 threshold=500")
system("d.mon x0");system("d.rast STREAM_90")

# the stream raster usually requires thinning
system("r.thin --o in=STREAM_90 out=STREAM_90thin iterations=500")
system("r.to.vect --o in=STREAM_90thin out=STREAM_90 feature=line")
system("d.mon x0");system("d.vect STREAM_90")
system("v.clean --o in=STREAM_90 out=STREAM_90clean type=line error=STREAM_90error tool=rmdangle thresh=250")
system("d.mon x0");system("d.vect STREAM_90error")
system("d.mon x1");system("d.vect STREAM_90clean")
system("d.mon x2");system("d.vect STREAM_90")

# basins
#threshold: minimum number of cells
system("r.watershed --o -f -a elev=ELEV_90 drainage=DRAINAGE_90 basin=BASIN_90 threshold=100000")
system("d.mon x0")
system("d.rast BASIN_90")

# convert to vectors
system("r.to.vect --o in=BASIN_90 out=BASIN_90 feature=area")
system("d.mon x0")
system("d.vect BASIN_90")

# Make a hillshade raster for displaying "3D"
system("r.shaded.relief --o map=ELEV_90 shade=SHADE_90 zmult=3")

# display layers
system("d.mon x0")
system("d.his h=ELEV_90 i=SHADE_90")
system("d.vect map=STREAM_90clean color=blue")
system("d.vect map=BASIN_90 type=boundary color=red")

# delineate basin using water outlet (identifyed using QGIS)6715781.97
system("r.water.outlet --o drainage=DRAINAGE_90 basin=BASIN_90 easting=229822.497 northing=6715837.518")
system("d.mon x0")
system("d.rast BASIN_90")
system("r.to.vect --o in=BASIN_90 out=BASIN_90 feature=area")

# display layers
system("d.mon x0")
system("d.his h=ELEV_90 i=SHADE_90")
system("d.vect map=STREAM_90clean color=blue width=2")
system("d.vect map=BASIN_90 type=boundary color=red width=2")

# ELEV_10 ======================================================================

# convert vector lakes to raster
system("v.to.rast --o in=lakes25 out=LAKE_25 use=val value=1")
system("d.mon x0")
system("d.rast LAKE_25")

# streams
writeVECT6(rivers25, "STREAM_10")
system("d.mon x0")
system("d.vect STREAM_10")
system("v.to.rast in=STREAM_10 out=STREAM_10 use=val value=1")

# basins
system("r.watershed --o -f -a elev=ELEV_10 depression=LAKE_25 drainage=DRAINAGE_10 basin=BASIN_10 threshold=50000")
system("d.mon x0")
system("d.rast BASIN_10")

# convert to vectors
system("r.to.vect --o in=BASIN_10 out=BASIN_10 feature=area")
system("v.build BASIN_10")
system("d.mon x0")
system("d.vect BASIN_10")

# create buffer using the positional accuracy (RMSE)
system("v.buffer --o in=BASIN_10 out=buffer_BASIN_10 distance=30")
system("v.to.rast in=buffer_BASIN_10 out=buffer_BASIN_10 use=val value=1")

# Make a hillshade raster for displaying "3D"
system("r.shaded.relief --o map=ELEV_10 shade=SHADE_10 zmult=3")

# display layers
system("d.mon x1")
system("d.his h=ELEV_10 i=SHADE_10")
system("d.vect map=STREAM_10 color=blue")
system("d.vect map=BASIN_10 type=boundary color=red")

# delineate basin using water outlet identifyed using QGIS
# point located on the bridge (-53.78969,-29.65868,229969.610,6715780.680)
system("r.water.outlet --o drainage=DRAINAGE_10 basin=BASIN_10 easting=229969.610 northing=6715780.680")
system("d.mon x0")
system("d.rast BASIN_10")
system("r.to.vect --o in=BASIN_10 out=BASIN_10 feature=area")

# display layers
system("d.erase")
system("d.his h=ELEV_10 i=SHADE_10")
system("d.vect map=STREAM_10 color=blue width=2")
system("d.vect map=BASIN_10 type=boundary color=red width=2")

# SAVE FIGURES WITH DIGITAL ELEVATION MODELS ###################################
# Use as boundary the basin estimated with the digital elevation models derived
# from the contour lines plus the 30 m buffer (buffer_BASIN_10).
system("r.mask -o buffer_BASIN_10")
boundary <- readVECT6("buffer_BASIN_10")@bbox
breaks <- seq(50, 550, 5)
colors <- terrain.colors(length(breaks) - 1)
# ELEV_10
map1 <- readRAST6("ELEV_10")
map1@bbox <- boundary
map1@data[, 1] <- as.integer(map1@data[, 1])
map1 <- spplot(map1, main = "", at = breaks, col.regions = colors)
# ELEV_90
map2 <- readRAST6("ELEV_90")
map2@bbox <- boundary
map2@data[, 1] <- as.integer(map2@data[, 1])
map2 <- spplot(map2, main = "", at = breaks, col.regions = colors)
# save
dev.off()
pdf(file = paste(dem.dir, "ELEV_10.pdf", sep = ""),
    width = 6.3/cm(1), height = 6.3/cm(1))
trellis.par.set(fontsize = list(text = 7, points = 5),
                plot.line = list(lwd = 0.001),
                axis.line = list(lwd = 0.01),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
print(map1)
dev.off()
rm(map1)
gc()
dev.off()
pdf(file = paste(dem.dir, "ELEV_90.pdf", sep = ""),
    width = 6.3/cm(1), height = 6.3/cm(1))
trellis.par.set(fontsize = list(text = 7, points = 5),
                plot.line = list(lwd = 0.001),
                axis.line = list(lwd = 0.01),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
print(map2)
dev.off()
rm(map2)
gc()
rm(boundary)
gc()
system("g.remove MASK")







# Convergence Index
window = 69
system(paste(grass.addons, "r.convergence input= output= window=", window, " weights=inverse", sep = ""))
rm(window)

# Distance from Urban Areas
system(paste("r.grow.distance -m input= output= value= metric=euclidean"))

# Multiresolution Index of Valley Bottom Flatness
system(paste("saga_cmd ta_morphometry 8"))
rsaga.get.usage(lib = "ta_morphometry", module = 8, env = sagaenv)
rsaga.geoprocessor(lib = "ta_morphometry", module = 8, env = sagaenv,
                   param = list(DEM = "SRTM_90.sgrd",
                                MRVBF = "MIVBF_90.sgrd",
                                MRRTF = "MIRTF_90.sgrd"))
spplot(readGDAL("MIRTF_90.sdat"), scales = list(draw = TRUE), col.regions = terrain.colors(64))

# Downslope Distance Gradient
# Method: gradient (degree)
system(paste("saga_cmd ta_morphometry 9"))
rsaga.geoprocessor(lib = "ta_morphometry", module = 9, env = sagaenv,
                   param = list(DEM = "SRTM_90.sgrd",
                                GRADIENT = "DDG_90.sgrd",
                                DIFFERENCE = "DDGd_90.sgrd",
                                DISTANCE = 10,
                                OUTPUT = 2))
spplot(readGDAL("DDG_90.sdat"), scales = list(draw = TRUE), col.regions = terrain.colors(64))

# Mass Balance Index
# Needs 'vertical distance to channel network'
system(paste("saga_cmd ta_morphometry 10"))
rsaga.geoprocessor(lib = "ta_morphometry", module = 10, env = sagaenv,
                   param = list(DEM = "SRTM_90.sgrd",
                                HREL = "VDCN_90.sgrd",
                                MBI = "MBI_90.sgrd",
                                TSLOPE = 15,
                                TCURVE = 0.01,
                                THREL = 15))
spplot(readGDAL("MBI_90.sdat"), scales = list(draw = TRUE), col.regions = terrain.colors(64))

# Relative Heights and Slope Positions
# Slope Height, Valley Depth, Normalized Height,
# Standardized Height, Mid-Slope Position
system(paste("saga_cmd ta_morphometry 14"))
rsaga.geoprocessor(lib = "ta_morphometry", module = 14, env = sagaenv,
                   param = list(DEM = "SRTM_90.sgrd",
                                HO = "SLOPH_90.sgrd",
                                HU = "VALL_90.sgrd",
                                NH = "NORMH_90.sgrd",
                                SH = "STANH_90.sgrd",
                                MS = "MID_90.sgrd",
                                W = 0.5, T = 10, E = 2))
spplot(readGDAL("SLOPH_90.sdat"), scales = list(draw = TRUE), col.regions = terrain.colors(64))

# Terrain Ruggedness Index (TRI)
# RADIUS: 540m / 90m = 6 cells --> follows ten Caten
# et al. Spatial resolution of a digital elevation
# model defined by the wavelet function. PAB, v.47, p.449-457, 2012.
system(paste("saga_cmd ta_morphometry 16"))
rsaga.geoprocessor(lib = "ta_morphometry", module = 16, env = sagaenv,
                   param = list(DEM = "SRTM_90.sgrd",
                                TRI = "TRI_90.sgrd",
                                RADIUS = 6,
                                DISTANCE_WEIGHTING_DW_WEIGHTING = 0))
spplot(readGDAL("TRI_90.sdat"), scales = list(draw = TRUE), col.regions = terrain.colors(64))

# Vector Ruggedness Measure (VRM)
# RADIUS: 540m / 90m = 6 cells --> follows ten Caten
# et al. Spatial resolution of a digital elevation
# model defined by the wavelet function. PAB, v.47, p.449-457, 2012.
system(paste("saga_cmd ta_morphometry 17"))
rsaga.geoprocessor(lib = "ta_morphometry", module = 17, env = sagaenv,
                   param = list(DEM = "SRTM_90.sgrd",
                                VRM = "VRM_90.sgrd",
                                RADIUS = 6,
                                DISTANCE_WEIGHTING_DW_WEIGHTING = 0))
spplot(readGDAL("VRM_90.sdat"), scales = list(draw = TRUE), col.regions = terrain.colors(64))

# Topographic Position Index (TPI)
# RADIUS_MAX: 540m --> follows ten Caten
# et al. Spatial resolution of a digital elevation
# model defined by the wavelet function. PAB, v.47, p.449-457, 2012.
system(paste("saga_cmd ta_morphometry 18"))
rsaga.geoprocessor(lib = "ta_morphometry", module = 18, env = sagaenv,
                   param = list(DEM = "SRTM_90.sgrd",
                                TPI = "TPI_90.sgrd",
                                RADIUS_MIN = 0,
                                RADIUS_MAX = 540))
spplot(readGDAL("TPI_90.sdat"), scales = list(draw = TRUE), col.regions = terrain.colors(64))

# Terrain Surface Texture
# RADIUS = 540m / 90m = 6 cells --> based on ten Caten # et al.
# Spatial resolution of a digital elevation model defined by the
# wavelet function. PAB, v.47, p.449-457, 2012.
system(paste("saga_cmd ta_morphometry 20"))
rsaga.geoprocessor(lib = "ta_morphometry", module = 20, env = sagaenv,
                   param = list(DEM = "SRTM_90.sgrd",
                                TEXTURE = "TST_90.sgrd",
                                RADIUS = 6))
spplot(readGDAL("TST_90.sdat"), scales = list(draw = TRUE), col.regions = terrain.colors(64))

# Terrain Surface Convexity (TSC)
# RADIUS = 540m / 90m = 6 cells --> based on ten Caten # et al.
# Spatial resolution of a digital elevation model defined by the
# wavelet function. PAB, v.47, p.449-457, 2012.
system(paste("saga_cmd ta_morphometry 21"))
rsaga.get.usage(lib = "ta_morphometry", module = 21, env = sagaenv)
rsaga.geoprocessor(lib = "ta_morphometry", module = 21, env = sagaenv,
                   param = list(DEM = "SRTM_90.sgrd",
                                CONVEX = "TSC_90.sgrd",
                                RADIUS = 6))
spplot(readGDAL("TSC_90.sdat"), scales = list(draw = TRUE), col.regions = terrain.colors(64))

# Library ta_hydrology ===============================================

# Catchment Area (Parallel) ------------------------------------------
# Catchment Area, Catchment Height, Catchment Slope,
# Catchment Aspect, Flow Path Length
# 'CSLOPE' e 'CASPECT' são dados em radianos
rsaga.get.modules(lib = "ta_hydrology")
rsaga.get.usage(lib = "ta_hydrology", module = 0, env = sagaenv)
rsaga.geoprocessor(lib = "ta_hydrology", module = 0, env = sagaenv,
                   param = list(ELEVATION = "SRTM_90.sgrd",
                                CAREA = "CAREA_90.sgrd",
                                CHEIGHT = "CHEIG_90.sgrd",
                                CSLOPE = "CSLOPrad_90.sgrd",
                                CASPECT = "CASPErad_90.sgrd",
                                FLWPATH = "FLOL_90.sgrd",
                                Method = 4))
spplot(readGDAL("CAREA_90.sdat"), scales = list(draw = TRUE), col.regions = terrain.colors(64))

# Convert 'CSLOPrad' to percentage 
formula <- paste("(tan(a))*100", sep = "")
rsaga.grid.calculus("CSLOPrad_90.sgrd", "CSLOP_90.sgrd", formula)
rm(formula)
spplot(readGDAL("CSLOP_90.sdat"), scales = list(draw = TRUE), col.regions = terrain.colors(64))

# Convert 'CASPErad to degrees
fac = round(180/pi, 4)
formula = paste(fac, "*a", sep = "")
rsaga.grid.calculus("CASPErad_90.sgrd", "CASPEdeg_90.sgrd", formula)
rm(fac, formula)

# Convert 'CATCHMENT ASPECT' to 'CATCHMENT NORTHERNESS'
fac <- 180
formula <- paste("abs(", fac, "-a)", sep = "")
rsaga.grid.calculus("CASPEdeg_90.sgrd", "CNORT_90.sgrd", formula)
rm(fac, formula)
spplot(readGDAL("CNORT_90.sdat"), scales = list(draw = TRUE), col.regions = terrain.colors(64))

# Slope Length -------------------------------------------------------
rsaga.get.usage(lib = "ta_hydrology", module = 7, env = sagaenv)
rsaga.geoprocessor(lib = "ta_hydrology", module = 7, env = sagaenv,
                   param = list(DEM = "SRTM_90.sgrd",
                                LENGTH = "SLOPL_90.sgrd"))
spplot(readGDAL("SLOPL_90.sdat"), scales = list(draw = TRUE), col.regions = terrain.colors(64))

# Flow Width and Specific Catchment Area -----------------------------
rsaga.get.usage(lib = "ta_hydrology", module = 19, env = sagaenv)
rsaga.geoprocessor(lib = "ta_hydrology", module = 19, env = sagaenv,
                   param = list(DEM = "SRTM_90.sgrd",
                                WIDTH = "FLOW_90.sgrd",
                                TCA = "CAREA_90.sgrd",
                                SCA = "SCA_90.sgrd",
                                METHOD = 1))
spplot(readGDAL("SCA_90.sdat"), scales = list(draw = TRUE), col.regions = terrain.colors(100))

# Topographic Wetness Index (TWI) ------------------------------------
rsaga.get.usage(lib = "ta_hydrology", module = 20, env = sagaenv)
rsaga.geoprocessor(lib = "ta_hydrology", module = 20, env = sagaenv,
                   param = list(SLOPE = "SLOPrad_90.sgrd",
                                AREA = "SCA_90.sgrd",
                                TWI = "TWI_90.sgrd"))
spplot(readGDAL("TWI_90.sdat"), scales = list(draw = TRUE), col.regions = terrain.colors(64))

# Stream Power Index (SPI) -------------------------------------------
rsaga.get.usage(lib = "ta_hydrology", module = 21, env = sagaenv)
rsaga.geoprocessor(lib = "ta_hydrology", module = 21, env = sagaenv,
                   param = list(SLOPE = "SLOPrad_90.sgrd",
                                AREA = "SCA_90.sgrd",
                                SPI = "SPI_90.sgrd"))
spplot(readGDAL("SPI_90.sdat"), scales = list(draw = TRUE), col.regions = terrain.colors(100))

# LS Factor ----------------------------------------------------------
rsaga.get.usage(lib = "ta_hydrology", module = 22, env = sagaenv)
rsaga.geoprocessor(lib = "ta_hydrology", module = 22, env = sagaenv,
                   param = list(SLOPE = "SLOPrad_90.sgrd",
                                AREA = "SCA_90.sgrd",
                                LS = "LS_90.sgrd"))
spplot(readGDAL("LS_90.sdat"), scales = list(draw = TRUE), col.regions = terrain.colors(100))

rsaga.get.modules("ta_hydrology")

# Channel Network --------------------------------------------------------------
system(paste("saga_cmd ta_channels 0 -ELEVATION=", dem.dir, "SRTM_90.sgrd -CHNLNTWRK=", dem.dir, "CHAN_90.sgrd -CHNLROUTE=", dem.dir, "CHAN_DIR_90.sgrd -SHAPES=", dem.dir, "CHAN_90.shp -INIT_GRID=", dem.dir, "CAREA_90.sgrd -INIT_METHOD=2 -INIT_VALUE=1000 -MINLEN=25", sep = ""))
spplot(readGDAL(paste(dem.dir, "CHAN_90.sdat", sep = "")), scales = list(draw = TRUE), col.regions = terrain.colors(64))

# Vertical Distance to Channel Network -----------------------------------------
rsaga.get.modules(lib = "ta_channels")
rsaga.get.usage(lib = "ta_channels", module = 3, env = sagaenv)
rsaga.geoprocessor(lib = "ta_channels", module = 3, env = sagaenv,
                   param=list(ELEVATION = "SRTM_90.sgrd",
                              CHANNELS = "CHAN_90.sgrd",
                              DISTANCE = "VDCN_90.sgrd",
                              BASELEVEL = "baselevel_90.sgrd"))
spplot(readGDAL("VDCN_90.sdat"), scales = list(draw = TRUE), col.regions = terrain.colors(64))

# Overland Flow Distance to Channel Network --------------------------
rsaga.get.usage(lib = "ta_channels", module = 4, env = sagaenv)
rsaga.geoprocessor(lib = "ta_channels", module = 4, env = sagaenv,
                   param = list(ELEVATION = "SRTM_90.sgrd",
                                CHANNELS = "CHAN_90.sgrd",
                                DISTANCE = "OFDCN_90.sgrd",
                                DISTVERT = "VOFDCN_90.sgrd",
                                DISTHORZ = "HOFDCN_90.sgrd",
                                METHOD = 1))
spplot(readGDAL("HOFDCN_90.sdat"), scales = list(draw = TRUE), col.regions = terrain.colors(64))

# Library ta_lighting ================================================
rsaga.get.modules(lib="ta_lighting")

# Potential Incoming Solar Radiation ---------------------------------
rsaga.get.usage(lib = "ta_lighting", module = 2, env = sagaenv)
rsaga.geoprocessor(lib = "ta_lighting", module = 2, env = sagaenv,
                   param = list(GRD_DEM = "SRTM_90.sgrd",
                                GRD_TOTAL = "SOLAR_90.sgrd",
                                LATITUDE = -29.65,
                                PERIOD = 2, DHOUR = 1, DDAYS = 1,
                                HOUR_RANGE_MIN = 5,
                                HOUR_RANGE_MAX = 22,
                                DAY_A = 0, MON_A = 0,
                                DAY_B = 30, MON_B = 11, METHOD = 0))
spplot(readGDAL("SOLAR_90.sdat"), scales = list(draw = TRUE), col.regions = terrain.colors(64))

# Watershed Basins (90 m) ------------------------------------------------------
system(paste("saga_cmd ta_channels 1 -ELEVATION=", dem.dir, "SRTM_90.sgrd -CHANNELS=", dem.dir, "CHAN_90.sgrd -BASINS=", dem.dir, "BASINS_90.sgrd -MINSIZE=25", sep = ""))
spplot(readGDAL(paste(dem.dir, "BASINS_90.sdat", sep = "")), scales = list(draw = TRUE), col.regions = terrain.colors(64))

# Vectorising Grid Classes
system(paste("saga_cmd shapes_grid 6 -GRID=", dem.dir, "BASINS_90.sgrd -POLYGONS=", dem.dir, "BASINS_90.shp", sep = ""))
spplot(shapefile(paste(dem.dir, "BASINS_90.shp", sep = "")))

# Derive terrain attributes - TOPODATA #########################################

# Library ta_morphometry =============================================
rsaga.get.modules(lib = "ta_morphometry")

# Elevation ----------------------------------------------------------
tmp <- read.sgrd("TOPO_30")
str(tmp)
write.sgrd(tmp, "ELEV_30")
rm(tmp)
spplot(readGDAL("ELEV_30.sdat"), scales = list(draw = TRUE),
       col.regions = topo.colors(64))

# Slope, Aspect, Curvature -------------------------------------------
# 'Slope' and 'Aspect' are given in radians
# Method: Fit 2.Degree Polynom (Zevenbergen & Thorne 1987)
rsaga.get.modules(lib = "ta_morphometry")
rsaga.get.usage(lib = "ta_morphometry", module = 0, env = sagaenv)
rsaga.geoprocessor(lib = "ta_morphometry", module = 0,
                   param = list(ELEVATION = "TOPO_30.sgrd",
                                SLOPE = "SLOPrad_30.sgrd",
                                ASPECT = "ASPErad_30.sgrd",
                                CURV = "CURV_30.sgrd",
                                HCURV = "HCURV_30.sgrd",
                                VCURV = "VCURV_30.sgrd",
                                METHOD = 5), env = sagaenv)
spplot(readGDAL("SLOPrad_30.sdat"), scales = list(draw = TRUE),
       col.regions = topo.colors(64))

# Convert 'ASPECT' to degrees
fac = round(180/pi, 4)
formula = paste(fac, "*a", sep = "")
rsaga.grid.calculus("ASPErad_30.sgrd", "ASPEdeg_30.sgrd", formula)
rm(fac, formula)

# Convert 'ASPECT' to 'NORTHERNESS'
fac <- 180
formula <- paste("abs(", fac, "-a)", sep = "")
rsaga.grid.calculus("ASPEdeg_30.sgrd", "NORT_30.sgrd", formula)
rm(fac, formula)
spplot(readGDAL("NORT_30.sdat"), scales = list(draw = TRUE),
       col.regions = topo.colors(64))

# Convert 'SLOPE' to percentage 
formula <- paste("(tan(a))*100", sep = "")
rsaga.grid.calculus("SLOPrad_30.sgrd", "SLOP_30.sgrd", formula)
rm(formula)
spplot(readGDAL("SLOP_30.sdat"), scales = list(draw = TRUE),
       col.regions = topo.colors(64))

# Convergence index --------------------------------------------------
# Method: Gradient with 3x3 neighbours
rsaga.get.usage(lib = "ta_morphometry", module = 1, env = sagaenv)
rsaga.geoprocessor(lib = "ta_morphometry", module = 1,
                   param = list(ELEVATION ="TOPO_30.sgrd",
                                RESULT = "CONV_30.sgrd",
                                METHOD = 1,
                                NEIGHBOURS = 1), env = sagaenv)
spplot(readGDAL("CONV_30.sdat"), scales = list(draw = TRUE),
       col.regions = terrain.colors(64))

# Morphometric protection index --------------------------------------
# 'RADIUS' definido conforme resultados de ten Caten
# et al. Spatial resolution of a digital elevation
# model defined by the wavelet function. PAB, v.47,
# p.449-457, 2012.
rsaga.get.usage(lib = "ta_morphometry", module = 7, env = sagaenv)
rsaga.geoprocessor(lib = "ta_morphometry", module = 7, env = sagaenv,
                   param = list(DEM = "TOPO_30.sgrd",
                                PROTECTION = "MPI_30.sgrd",
                                RADIUS = 120))
spplot(readGDAL("MPI_30.sdat"), scales = list(draw = TRUE),
       col.regions = terrain.colors(64))

# Multiresolution Index of Valley Bottom Flatness --------------------
rsaga.get.usage(lib = "ta_morphometry", module = 8, env = sagaenv)
rsaga.geoprocessor(lib = "ta_morphometry", module = 8, env = sagaenv,
                   param = list(DEM = "TOPO_30.sgrd",
                                MRVBF = "MIVBF_30.sgrd",
                                MRRTF = "MIRTF_30.sgrd"))
spplot(readGDAL("MIVBF_30.sdat"), scales = list(draw = TRUE),
       col.regions = terrain.colors(64))

# Downslope Distance Gradient ----------------------------------------
# Method: gradient (degree)
rsaga.get.usage(lib = "ta_morphometry", module = 9, env = sagaenv)
rsaga.geoprocessor(lib = "ta_morphometry", module = 9, env = sagaenv,
                   param = list(DEM = "TOPO_30.sgrd",
                                GRADIENT = "DDG_30.sgrd",
                                DIFFERENCE = "DDGd_30.sgrd",
                                DISTANCE = 10,
                                OUTPUT = 2))
spplot(readGDAL("DDG_30.sdat"), scales = list(draw = TRUE),
       col.regions = terrain.colors(64))

# Mass Balance Index -------------------------------------------------
# Needs 'vertical distance to channel network'
rsaga.get.usage(lib = "ta_morphometry", module = 10, env = sagaenv)
rsaga.geoprocessor(lib = "ta_morphometry", module = 10, env = sagaenv,
                   param = list(DEM = "TOPO_30.sgrd",
                                HREL = "VDCN_30.sgrd",
                                MBI = "MBI_30.sgrd",
                                TSLOPE = 15,
                                TCURVE = 0.01,
                                THREL = 15))
spplot(readGDAL("MBI_30.sdat"), scales = list(draw = TRUE),
       col.regions = terrain.colors(64))

# Relative Heights and Slope Positions -------------------------------
# Slope Height, Valley Depth, Normalized Height,
# Standardized Height, Mid-Slope Position
rsaga.get.usage(lib = "ta_morphometry", module = 14, env = sagaenv)
rsaga.geoprocessor(lib = "ta_morphometry", module = 14, env = sagaenv,
                   param = list(DEM = "TOPO_30.sgrd",
                                HO = "SLOPH_30.sgrd",
                                HU = "VALL_30.sgrd",
                                NH = "NORMH_30.sgrd",
                                SH = "STANH_30.sgrd",
                                MS = "MID_30.sgrd",
                                W = 0.5, T = 10, E = 2))
spplot(readGDAL("SLOPH_30.sdat"), scales = list(draw = TRUE),
       col.regions = terrain.colors(64))

# Terrain Ruggedness Index (TRI) -------------------------------------
# RADIUS: 540m / 30m = 18 cells --> follows ten Caten
# et al. Spatial resolution of a digital elevation
# model defined by the wavelet function. PAB, v.47, p.449-457, 2012.
rsaga.get.usage(lib = "ta_morphometry", module = 16, env = sagaenv)
rsaga.geoprocessor(lib = "ta_morphometry", module = 16, env = sagaenv,
                   param = list(DEM = "TOPO_30.sgrd",
                                TRI = "TRI_30.sgrd",
                                RADIUS = 18,
                                DISTANCE_WEIGHTING_DW_WEIGHTING = 0))
spplot(readGDAL("TRI_30.sdat"), scales = list(draw = TRUE),
       col.regions = terrain.colors(64))

# Vector Ruggedness Measure (VRM) ------------------------------------
# RADIUS: 540m / 30m = 18 cells --> follows ten Caten
# et al. Spatial resolution of a digital elevation
# model defined by the wavelet function. PAB, v.47, p.449-457, 2012.
rsaga.get.usage(lib = "ta_morphometry", module = 17, env = sagaenv)
rsaga.geoprocessor(lib = "ta_morphometry", module = 17, env = sagaenv,
                   param = list(DEM = "TOPO_30.sgrd",
                                VRM = "VRM_30.sgrd",
                                RADIUS = 18,
                                DISTANCE_WEIGHTING_DW_WEIGHTING = 0))
spplot(readGDAL("VRM_30.sdat"), scales = list(draw = TRUE),
       col.regions = terrain.colors(64))

# Topographic Position Index (TPI) -----------------------------------
# RADIUS_MAX: 540m --> follows ten Caten
# et al. Spatial resolution of a digital elevation
# model defined by the wavelet function. PAB, v.47, p.449-457, 2012.
rsaga.get.usage(lib = "ta_morphometry", module = 18, env = sagaenv)
rsaga.geoprocessor(lib = "ta_morphometry", module = 18, env = sagaenv,
                   param = list(DEM = "TOPO_30.sgrd",
                                TPI = "TPI_30.sgrd",
                                RADIUS_MIN = 0,
                                RADIUS_MAX = 540))
spplot(readGDAL("TPI_30.sdat"), scales = list(draw = TRUE),
       col.regions = terrain.colors(64))

# Terrain Surface Texture --------------------------------------------
# RADIUS = 540m / 30m = 18 cells --> based on ten Caten # et al.
# Spatial resolution of a digital elevation model defined by the
# wavelet function. PAB, v.47, p.449-457, 2012.
rsaga.get.usage(lib = "ta_morphometry", module = 20, env = sagaenv)
rsaga.geoprocessor(lib = "ta_morphometry", module = 20, env = sagaenv,
                   param = list(DEM = "TOPO_30.sgrd",
                                TEXTURE = "TST_30.sgrd",
                                RADIUS = 18))
spplot(readGDAL("TST_30.sdat"), scales = list(draw = TRUE),
       col.regions = terrain.colors(64))

# Terrain Surface Convexity (TSC) ------------------------------------
# RADIUS = 540m / 30m = 18 cells --> based on ten Caten # et al.
# Spatial resolution of a digital elevation model defined by the
# wavelet function. PAB, v.47, p.449-457, 2012.
rsaga.get.usage(lib = "ta_morphometry", module = 21, env = sagaenv)
rsaga.geoprocessor(lib = "ta_morphometry", module = 21, env = sagaenv,
                   param = list(DEM = "TOPO_30.sgrd",
                                CONVEX = "TSC_30.sgrd",
                                RADIUS = 18))
spplot(readGDAL("TSC_30.sdat"), scales = list(draw = TRUE),
       col.regions = terrain.colors(64))

# Library ta_hydrology ===============================================

# Catchment Area (Parallel) ------------------------------------------
# Catchment Area, Catchment Height, Catchment Slope,
# Catchment Aspect, Flow Path Length
# 'CSLOPE' e 'CASPECT' são dados em radianos
rsaga.get.modules(lib = "ta_hydrology")
rsaga.get.usage(lib = "ta_hydrology", module = 0, env = sagaenv)
rsaga.geoprocessor(lib = "ta_hydrology", module = 0, env = sagaenv,
                   param = list(ELEVATION = "TOPO_30.sgrd",
                                CAREA = "CAREA_30.sgrd",
                                CHEIGHT = "CHEIG_30.sgrd",
                                CSLOPE = "CSLOPrad_30.sgrd",
                                CASPECT = "CASPErad_30.sgrd",
                                FLWPATH = "FLOL_30.sgrd",
                                Method = 4))
spplot(readGDAL("CAREA_30.sdat"), scales = list(draw = TRUE),
       col.regions = terrain.colors(64))

# Convert 'CSLOPrad' to percentage 
formula <- paste("(tan(a))*100", sep = "")
rsaga.grid.calculus("CSLOPrad_30.sgrd", "CSLOP_30.sgrd", formula)
rm(formula)
spplot(readGDAL("CSLOP_30.sdat"), scales = list(draw = TRUE),
       col.regions = terrain.colors(64))

# Convert 'CASPErad to degrees
fac = round(180/pi, 4)
formula = paste(fac, "*a", sep = "")
rsaga.grid.calculus("CASPErad_30.sgrd", "CASPEdeg_30.sgrd", formula)
rm(fac, formula)

# Convert 'CATCHMENT ASPECT' to 'CATCHMENT NORTHERNESS'
fac <- 180
formula <- paste("abs(", fac, "-a)", sep = "")
rsaga.grid.calculus("CASPEdeg_30.sgrd", "CNORT_30.sgrd", formula)
rm(fac, formula)
spplot(readGDAL("CNORT_30.sdat"), scales = list(draw = TRUE),
       col.regions = terrain.colors(64))

# Slope Length -------------------------------------------------------
rsaga.get.usage(lib = "ta_hydrology", module = 7, env = sagaenv)
rsaga.geoprocessor(lib = "ta_hydrology", module = 7, env = sagaenv,
                   param = list(DEM = "TOPO_30.sgrd",
                                LENGTH = "SLOPL_30.sgrd"))
spplot(readGDAL("SLOPL_30.sdat"), scales = list(draw = TRUE),
       col.regions = terrain.colors(64))

# Flow Width and Specific Catchment Area -----------------------------
rsaga.get.usage(lib = "ta_hydrology", module = 19, env = sagaenv)
rsaga.geoprocessor(lib = "ta_hydrology", module = 19, env = sagaenv,
                   param = list(DEM = "TOPO_30.sgrd",
                                WIDTH = "FLOW_30.sgrd",
                                TCA = "CAREA_30.sgrd",
                                SCA = "SCA_30.sgrd",
                                METHOD = 1))
spplot(readGDAL("SCA_30.sdat"), scales = list(draw = TRUE),
       col.regions = terrain.colors(100))

# Topographic Wetness Index (TWI) ------------------------------------
rsaga.get.usage(lib = "ta_hydrology", module = 20, env = sagaenv)
rsaga.geoprocessor(lib = "ta_hydrology", module = 20, env = sagaenv,
                   param = list(SLOPE = "SLOPrad_30.sgrd",
                                AREA = "SCA_30.sgrd",
                                TWI = "TWI_30.sgrd"))
spplot(readGDAL("TWI_30.sdat"), scales = list(draw = TRUE),
       col.regions = terrain.colors(64))

# Stream Power Index (SPI) -------------------------------------------
rsaga.get.usage(lib = "ta_hydrology", module = 21, env = sagaenv)
rsaga.geoprocessor(lib = "ta_hydrology", module = 21, env = sagaenv,
                   param = list(SLOPE = "SLOPrad_30.sgrd",
                                AREA = "SCA_30.sgrd",
                                SPI = "SPI_30.sgrd"))
spplot(readGDAL("SPI_30.sdat"), scales = list(draw = TRUE),
       col.regions = terrain.colors(100))

# LS Factor ----------------------------------------------------------
rsaga.get.usage(lib = "ta_hydrology", module = 22, env = sagaenv)
rsaga.geoprocessor(lib = "ta_hydrology", module = 22, env = sagaenv,
                   param = list(SLOPE = "SLOPrad_30.sgrd",
                                AREA = "SCA_30.sgrd",
                                LS = "LS_30.sgrd"))
spplot(readGDAL("LS_30.sdat"), scales = list(draw = TRUE),
       col.regions = terrain.colors(64))

# Channel Network (30 m) -------------------------------------------------------
system(paste("saga_cmd ta_channels 0 -ELEVATION=", dem.dir, "TOPO_30.sgrd -CHNLNTWRK=", dem.dir, "CHAN_30.sgrd -CHNLROUTE=", dem.dir, "CHAN_DIR_30.sgrd -SHAPES=", dem.dir, "CHAN_30.shp -INIT_GRID=", dem.dir, "CAREA_30.sgrd -INIT_METHOD=2 -INIT_VALUE=1000 -MINLEN=75", sep = ""))
spplot(readGDAL(paste(dem.dir, "CHAN_30.sdat", sep = "")), scales = list(draw = TRUE), col.regions = terrain.colors(64))

# Vertical Distance to Channel Network -------------------------------
rsaga.get.modules(lib = "ta_channels")
rsaga.get.usage(lib = "ta_channels", module = 3, env = sagaenv)
rsaga.geoprocessor(lib = "ta_channels", module = 3, env = sagaenv,
                   param=list(ELEVATION = "TOPO_30.sgrd",
                              CHANNELS = "CHAN_30.sgrd",
                              DISTANCE = "VDCN_30.sgrd",
                              BASELEVEL = "baselevel_30.sgrd"))
spplot(readGDAL("VDCN_30.sdat"), scales = list(draw = TRUE),
       col.regions = terrain.colors(64))

# Overland Flow Distance to Channel Network --------------------------
rsaga.get.usage(lib = "ta_channels", module = 4, env = sagaenv)
rsaga.geoprocessor(lib = "ta_channels", module = 4, env = sagaenv,
                   param = list(ELEVATION = "TOPO_30.sgrd",
                                CHANNELS = "CHAN_30.sgrd",
                                DISTANCE = "OFDCN_30.sgrd",
                                DISTVERT = "VOFDCN_30.sgrd",
                                DISTHORZ = "HOFDCN_30.sgrd",
                                METHOD = 1))
spplot(readGDAL("HOFDCN_30.sdat"), scales = list(draw = TRUE),
       col.regions = terrain.colors(64))

# Watershed Basins (30 m) ------------------------------------------------------
system(paste("saga_cmd ta_channels 1 -ELEVATION=", dem.dir, "TOPO_30.sgrd -CHANNELS=", dem.dir, "CHAN_30.sgrd -BASINS=", dem.dir, "BASINS_30.sgrd -MINSIZE=75", sep = ""))
spplot(readGDAL(paste(dem.dir, "BASINS_30.sdat", sep = "")), scales = list(draw = TRUE), col.regions = terrain.colors(64))

# Vectorising Grid Classes
system(paste("saga_cmd shapes_grid 6 -GRID=", dem.dir, "BASINS_30.sgrd -POLYGONS=", dem.dir, "BASINS_30.shp", sep = ""))
spplot(shapefile(paste(dem.dir, "BASINS_30.shp", sep = "")))

# Library ta_lighting ================================================
rsaga.get.modules(lib="ta_lighting")

# Potential Incoming Solar Radiation ---------------------------------
rsaga.get.usage(lib = "ta_lighting", module = 2, env = sagaenv)
rsaga.geoprocessor(lib = "ta_lighting", module = 2, env = sagaenv,
                   param = list(GRD_DEM = "TOPO_30.sgrd",
                                GRD_TOTAL = "SOLAR_30.sgrd",
                                LATITUDE = -29.65,
                                PERIOD = 2, DHOUR = 1, DDAYS = 1,
                                HOUR_RANGE_MIN = 5,
                                HOUR_RANGE_MAX = 22,
                                DAY_A = 0, MON_A = 0,
                                DAY_B = 30, MON_B = 11, METHOD = 0))
spplot(readGDAL("SOLAR_30.sdat"), scales = list(draw = TRUE),
       col.regions = terrain.colors(64))

# Save all terrain attributes (dnos-covars.RData) ####################

# SRTM ===============================================================
tmp <- list.files(pattern = "_90.sdat")
# raster stack
terrain90 <- stack(tmp)
# set names
names(terrain90) <- str_sub(tmp, start = 1, end = -6)
rm(tmp)
proj4string(terrain90) <- CRS("+proj=utm +zone=22 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
# check
str(terrain90[[1]])
plot(terrain90$NORT_90)

# TOPODATA ===========================================================
tmp <- list.files(pattern = "_30.sdat")
terrain30 <- stack(tmp)
names(terrain30) <- str_sub(tmp, start = 1, end = -6)
rm(tmp)
proj4string(terrain30) <- CRS("+proj=utm +zone=22 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
str(terrain30[[1]])
plot(terrain30$NORT_30)

# save RData #########################################################
ls()
save(buffer25, contours25, contours25hd, contours25ld, contours25md,
     grid.carta25,
     lakes25, points25, rivers25,
     file = "sm-dnos-terrain.RData")
# save(terrain90, terrain30,
#      file="/home/alessandro/PROJETOS/DNOS-SM/data/dnos.covars.RData")

