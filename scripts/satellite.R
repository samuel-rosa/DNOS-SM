# DESCRIPTION ##################################################################
# Processing of satellite imagery
# SETTINGS #####################################################################
rm(list = ls())
gc()
require(sp)
require(rgdal)
require(maptools)
require(raster)
require(spgrass6)
require(lattice)
require(latticeExtra)
load("sm-dnos-general.RData")
load("sm-dnos-covar-validation.RData")
load("/home/alessandro/PROJECTS/pedometrics/pedometrics/cooking/color_ramps.RData")
initGRASS(gisBase = "/usr/lib/grass64/", gisDbase = GRASSgisDbase, 
          location = "dnos-sm-rs", mapset = "predictions",
          pid = Sys.getpid(), override = TRUE)
system("g.region rast=dnos.raster")

# LOAD AND PROCESS DATA ########################################################

# Landsat 5 TM 2010dez26 =======================================================
maps <- c("BLUE_30", "GREEN_30", "RED_30", "NIR_30a", "NIR_30b",
          "MIR_30", "NDVI_30", "SAVI_30")
dir <- path.expand("~/PROJETOS/RS-BRASIL/landsat5tm/2010dez26/")

# AFFINE TRANSFORMATION
# get pixel coordinates (this has to be done manually)
coo <- coordinates(gcp_landsat)
a <- seq(1, dim(coo)[1], 1)
cmd <- paste("gdallocationinfo -geoloc ", dir, "BLUE_30.tif ",
             coo[a, 1], " ", coo[a, 2], sep = "")
lapply(cmd, system)
# Location: (5171P,6302L)
# Location: (5217P,6172L)
# Location: (5209P,5990L)
# Location: (5140P,6239L)
# Location: (5065P,6141L)
# Location: (5159P,6068L)
# Location: (5097P,5972L)
# Location: (4953P,6239L)
# Location: (4872P,5998L)
# Location: (4951P,6111L)
# Location: (5248P,6190L)
# Location: (5209P,6080L)
# Location: (5131P,6088L)
# Location: (5114P,6178L)
rm(coo, a, cmd)
pixel <- rbind(
  c("5171 6302"),
  c("5217 6172"),
  c("5209 5990"),
  c("5140 6239"),
  c("5065 6141"),
  c("5159 6068"),
  c("5097 5972"),
  c("4953 6239"),
  c("4872 5998"),
  c("4951 6111"),
  c("5248 6190"),
  c("5209 6080"),
  c("5131 6088"),
  c("5114 6178")
)
# get destination coordinates
dest <- coordinates(gcp.pos)
gcp <- sapply(seq(1, dim(dest)[1], 1), function (X) {
  paste("-gcp", pixel[X, ], dest[X, 1], dest[X, 2], "", sep = " ")
}) 
gcp <- paste(gcp, sep = "", collapse = "")
# add GCP to image (gdal_translate)
cmd <- paste("gdal_translate -of GTiff ", gcp, dir, maps,
             ".tif /tmp/", maps, ".tif", sep = "")
lapply(cmd, system)
rm(cmd, dest, gcp, pixel)
# gdalwarp
cmd <- paste("gdalwarp -s_srs epsg:32722 -t_srs epsg:32722 -r near -order 1 -co COMPRESS=NONE /tmp/", maps,
             ".tif ", dir, maps, "_affine.tif", sep = "")
lapply(cmd, system)
rm(cmd)

# READ FILES INTO GRASS GIS
cmd <- paste("r.in.gdal --o ", dir, maps, "_affine.tif out=", maps, sep = "")
lapply(cmd, system)
cmd <- paste("r.info ", maps, sep = "")
lapply(cmd, system)

# resample
system("r.resample")
output <- list("resBLUE_30", "resGREEN_30", "resRED_30", "resNIR_30a", "resNIR_30b", "resMIR_30",
               "resNDVI_30", "resSAVI_30")
cmd <- paste("r.resample --o in=", maps, " out=", output, sep = "")
lapply(cmd, system)
cmd <- paste("r.info ", output, sep = "")
lapply(cmd, system)

# rename
cmd <- paste("r.mapcalc '", maps, "=", output, "'", sep = "")
lapply(cmd, system)
system(paste("g.remove ", paste(output, collapse = ","), sep = ""))
cmd <- paste("r.info ", maps, sep = "")
lapply(cmd, system)
rm(cmd, dir, maps, output)

# RapidEye 2225403_2012-11-16 ====================================================
dir <- path.expand("~/PROJETOS/RS-BRASIL/rapideye/2225403_2012-11-16/")
maps <- list("BLUE_5", "GREEN_5", "RED_5", "EDGE_5", "NIR_5",
             "NDVI_5a", "NDVI_5b", "SAVI_5a", "SAVI_5b")

# AFFINE TRANSFORMATION
# get pixel coordinates (this has to be done manually)
coo <- coordinates(gcp_eye)
a <- seq(1, dim(coo)[1], 1)
cmd <- paste("gdallocationinfo -geoloc ", dir, "BLUE_5.tif ",
             coo[a, 1], " ", coo[a, 2], sep = "")
lapply(cmd, system)
# Location: (3880P,4519L)
# Location: (4156P,3737L)
# Location: (4114P,2639L)
# Location: (3694P,4136L)
# Location: (3254P,3552L)
# Location: (3807P,3113L)
# Location: (3437P,2541L)
# Location: (2575P,4132L)
# Location: (2085P,2690L)
# Location: (2561P,3373L)
# Location: (4344P,3843L)
# Location: (4111P,3187L)
# Location: (3651P,3235L)
# Location: (3544P,3776L)
rm(coo, a, cmd)
pixel <- rbind(
  c("3880 4519"),
  c("4156 3737"),
  c("4114 2639"),
  c("3694 4136"),
  c("3254 3552"),
  c("3807 3113"),
  c("3437 2541"),
  c("2575 4132"),
  c("2085 2690"),
  c("2561 3373"),
  c("4344 3843"),
  c("4111 3187"),
  c("3651 3235"),
  c("3544 3776")
)
# get destination coordinates
dest <- coordinates(gcp.pos)
gcp <- sapply(seq(1, dim(dest)[1], 1), function (X) {
  paste("-gcp", pixel[X, ], dest[X, 1], dest[X, 2], "", sep = " ")
}) 
gcp <- paste(gcp, sep = "", collapse = "")
# add GCP to image (gdal_translate)
cmd <- paste("gdal_translate -of GTiff ", gcp, dir, maps,
             ".tif /tmp/", maps, ".tif", sep = "")
lapply(cmd, system)
rm(cmd, dest, gcp, pixel)
# gdalwarp
cmd <- paste("gdalwarp -s_srs epsg:32722 -t_srs epsg:32722 -r near -order 1 -co COMPRESS=NONE /tmp/",
             maps, ".tif ", dir, maps, "_affine.tif", sep = "")
lapply(cmd, system)
rm(cmd)

# READ FILES INTO GRASS GIS AND RESAMPLE
cmd <- paste("r.in.gdal --o ", dir, maps, "_affine.tif out=", maps, sep = "")
lapply(cmd, system)
# resample
output <- list("resBLUE_5", "resGREEN_5", "resRED_5", "resEDGE_5", "resNIR_5",
               "resNDVI_5a",
               "resNDVI_5b", "resSAVI_5a", "resSAVI_5b")
cmd <- paste("r.resample --o in=", maps, " out=", output, sep = "")
lapply(cmd, system)
# rename
cmd <- paste("r.mapcalc '", maps, "=", output, "'", sep = "")
lapply(cmd, system)
system(paste("g.remove ", paste(output, collapse = ","), sep = ""))
cmd <- paste("r.info ", maps, sep = "")
lapply(cmd, system)
rm(cmd, maps, output)

# SAVE PLOTS ###################################################################
system("r.mask -o buffer_BASIN_10")
boundary <- readVECT6("buffer_BASIN_10")@bbox
dir <- "/home/alessandro/PROJECTS/DNOS-SM/satellite/"
# NDVI_30
map1 <- readRAST6("NDVI_30")
map1@bbox <- boundary
map1 <- spplot(map1, main = "", col.regions = ndvi$colors, at = ndvi$breaks)
dev.off()
pdf(file = paste(dir, "NDVI_30.pdf", sep = ""), 
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
# NDVI_5b
map2 <- readRAST6("NDVI_5b")
map2@bbox <- boundary
map2 <- spplot(map2, main = "", col.regions = ndvi$colors, at = ndvi$breaks)
dev.off()
pdf(file = paste(dir, "NDVI_5b.pdf", sep = ""), 
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

# End!