# DESCRIPTION ##################################################################
# Analysis of the land use data available.
# SETTINGS #####################################################################
rm(list = ls())
gc()
require(sp)
require(rgdal)
require(maptools)
require(raster)
require(RSAGA)
require(spatstat)
require(lattice)
library(lattticeExtra)
require(xtable)
library(plotKML)
require(VecStatGraphs2D)
require(vec2dtransf)
require(spgrass6)
library(plyr)
library(grid)
load("sm-dnos-general.RData")
load("sm-dnos-covar-validation.RData")
load("sm-dnos-land-use.RData")
load("/home/alessandro/Dropbox/color_ramps.RData")
ls()
# GRASS GIS
initGRASS(gisBase = "/usr/lib/grass64/", gisDbase = GRASSgisDbase, 
          override = TRUE, location = "dnos-sm-rs", mapset = "predictions", 
          pid = Sys.getpid())
system("g.region rast=dnos.raster")
#
# LOAD AND PREPROCESS DATA #####################################################
#
# Land Use 2008-09 =============================================================
LU2009 <- shapefile(paste(land.dir, "land-use-2009.shp", sep = ""))
LU2009 <- spTransform(LU2009, wgs1984utm22s)
LU2009
plot(LU2009)

# land use IDs
LU2009.id <- data.frame(unique(LU2009$code), unique(LU2009$land_use))
colnames(LU2009.id) <- c("code", "LU")
LU2009.id <- LU2009.id[order(LU2009.id$code),]
LU2009.id

# export to GRASS GIS (mapset = predictions)
initGRASS(gisBase = "/usr/lib/grass64/", gisDbase = GRASSgisDbase, 
          override = TRUE, location = "dnos-sm-rs", mapset = "predictions", 
          pid = Sys.getpid())
system("g.region rast=dnos.raster")
writeVECT6(LU2009, "LU2009", v.in.ogr_flags = "overwrite")
system("v.info -c LU2009")

# shapes to grid
system("v.to.rast --o in=LU2009 out=LU2009 use=attr column=code")
system("r.info LU2009"); system("r.category LU2009")
system("d.mon x0");system("d.rast LU2009")

# add category names
LU2009.id
system("r.category LU2009 rules=- <<EOF
       1:native forest
       2:shrubland
       3:animal husbandry
       4:crop agriculture
       5:forestry
       6:urban
       7:water
       EOF")
system("r.category LU2009")
system("d.mon x0");system("d.rast.leg LU2009")

# Land Use 1980 ================================================================
LU1980 <- shapefile(paste(land.dir, "land-use-1980.shp", sep = ""))
LU1980 <- spTransform(LU1980, wgs1984utm22s)
LU1980
plot(LU1980)

# land use IDs
LU1980.id <- data.frame(unique(LU1980$code), unique(LU1980$land_use))
colnames(LU1980.id) <- c("code", "LU")
LU1980.id <- LU1980.id[order(LU1980.id$code),]
LU1980.id

# apply affine transformation based on validation results
# see 'sm-dnos-covar-validation.R' for details
is.projected(LU1980)
LU1980 <- applyTransformation(affine.car25, LU1980)
plot(LU1980)

# export to GRASS GIS (mapset = predictions)
initGRASS(gisBase = "/usr/lib/grass64/", gisDbase = GRASSgisDbase,
          override = TRUE,
          location = "dnos-sm-rs", mapset = "predictions", pid = Sys.getpid())
system("g.region rast=dnos.raster")
writeVECT6(LU1980, "LU1980")
system("v.info -c LU1980")

# shapes to grid
system("v.to.rast --o in=LU1980 out=LU1980 use=attr column=code")
system("r.info LU1980"); system("r.category LU1980")
system("d.mon x0");system("d.rast.leg LU1980")

# add category names
LU1980.id
system("r.category LU1980 rules=- <<EOF
       1:native forest
       2:animal husbandry
       3:forestry
       4:urban
       5:water
       EOF")
system("r.category LU1980")
system("d.mon x0");system("d.rast.leg LU1980")

# CREATE ENVIRONMENTAL COVARIATES ##############################################
initGRASS(gisBase = "/usr/lib/grass64/", gisDbase = GRASSgisDbase,
          override = TRUE, location = "dnos-sm-rs",
          mapset = "predictions", pid = Sys.getpid())
system("g.region rast=dnos.raster")
#
# LU2009 =======================================================================
system("g.remove MASK")
system("r.category LU2009")
system("d.mon x0")
system("d.rast.leg LU2009")
system("d.vect map=BASIN_10 fcolor=none")
# There are seven classes in the study area. Two of them are not used: water and
# urban. Thus, there are p = 5 - 1 = 4 unique classes. Other classes can be
# created by the combination of unique classes.
#
# LU2009a (Native forest)
# Classes: 1 = Native forest
#          0 = Other (Shrubland + Animal husbandry + Crop agriculture + Forestry)
system("r.category LU2009")
system("r.reclass --o in=LU2009 out=LU2009a <<EOF
       1 = 1 Native forest
       * = 0 Other
       end
       EOF")
system("r.mapcalc LU2009a=LU2009a")
system("r.category LU2009a rules=- <<EOF
       1:Native forest
       0:Other
       EOF")
system("r.category LU2009a")
system("r.info LU2009a")
system("d.mon x0")
system("d.erase")
system("d.rast.leg LU2009a")

# LU2009b (Shrubland)
# Classes: 1 = Shrubland
#          0 = Other (Native forest + Animal husbandry + Crop agriculture + Forestry)
system("r.category LU2009")
system("r.reclass --o in=LU2009 out=LU2009b <<EOF
       2 = 1 Shrubland
       * = 0 Other
       end
       EOF")
system("r.mapcalc LU2009b=LU2009b")
system("r.category LU2009b rules=- <<EOF
       1:Shrubland
       0:Other
       EOF")
system("r.category LU2009b")
system("r.info LU2009b")
system("d.mon x0")
system("d.erase")
system("d.rast.leg LU2009b")
#
# LU2009c (Animal husbandry)
# Classes: 1 = Animal husbandry
#          0 = Other (Native forest + Shrubland + Crop agriculture + Forestry)
system("r.category LU2009")
system("r.reclass --o in=LU2009 out=LU2009c <<EOF
       3 = 1 Animal husbandry
       * = 0 Other
       end
       EOF")
system("r.mapcalc LU2009c=LU2009c")
system("r.category LU2009c rules=- <<EOF
       1:Animal husbandry
       0:Other
       EOF")
system("r.category LU2009c")
system("r.info LU2009c")
system("d.mon x0")
system("d.erase")
system("d.rast.leg LU2009c")
#
# LU2009d (Crop agriculture)
# Classes: 1 = Crop agriculture
#          0 = Other (Native forest + Shrubland + Animal husbandry + Forestry)
system("r.category LU2009")
system("r.reclass --o in=LU2009 out=LU2009d <<EOF
       4 = 1 Crop agriculture
       * = 0 Other
       end
       EOF")
system("r.mapcalc LU2009d=LU2009d")
system("r.category LU2009d rules=- <<EOF
       1:Crop agriculture
       0:Other
       EOF")
system("r.category LU2009d")
system("r.info LU2009d")
system("d.mon x0")
system("d.erase")
system("d.rast.leg LU2009d")
#
# LU2009e (Distance from settlement)
system("r.category LU2009")
system("r.reclass --o in=LU2009 out=LU2009e <<EOF
       6 = 1 Urban
       * = NULL
       end
       EOF")
system("r.mapcalc LU2009e=LU2009e")
system("d.mon x0")
system("d.rast.leg LU2009e")
system("r.grow.distance --o input=LU2009e distance=tmp")
system("r.mapcalc LU2009e=tmp")
system("d.erase")
system("d.rast.leg LU2009e")
#
# LU2009f (Distance from boundaries)
system("v.info LU2009")
system("v.type --o input=LU2009 output=tmpLU2009 type=boundary,line")
system("v.category input=tmpLU2009 output=tmptmpLU2009 option=add")
system("v.to.rast --o input=tmptmpLU2009 output=tmpLU2009 use=cat")
system("g.remove vect=tmpLU2009,tmptmpLU2009")
system("r.mapcalc 'tmpLU2009=if(tmpLU2009>0,1,0)'")
system("r.reclass --o in=tmpLU2009 out=tmptmpLU2009f <<EOF
       1 = 1 Boundary
       * = NULL
       end
       EOF")
system("r.mapcalc LU2009f=tmptmpLU2009f")
system("g.remove tmptmpLU2009f,tmpLU2009")
system("d.erase")
system("d.rast LU2009f")
system("r.grow.distance --o input=LU2009f distance=tmp")
system("r.mapcalc LU2009f=tmp")
system("d.erase")
system("d.rast.leg LU2009f")
#
# LU1980 =======================================================================
system("r.category LU1980")
system("d.mon x0")
system("d.rast.leg LU1980")
system("d.vect map=BASIN_10 fcolor=none")
# There are five classes in the study area. Two of them are not used: water and
# urban. Thus, there are p = 3 - 1 = 2 unique classes. Other classes can be
# created by the combination of unique classes.

# LU1980a (Native forest)
# Classes: 1 = Native forest
#          0 = Other (Animal husbandry + Forestry)
system("r.category LU1980")
system("r.reclass --o in=LU1980 out=LU1980a <<EOF
       1 = 1 Native forest
       * = 0 Other
       end
       EOF")
system("r.mapcalc LU1980a=LU1980a")
system("r.category LU1980a rules=- <<EOF
       1:Native forest
       0:Other
       EOF")
system("r.category LU1980a")
system("r.info LU1980a")
system("d.mon x0")
system("d.erase")
system("d.rast.leg LU1980a")

# LU1980b (Animal husbandry)
# Classes: 1 = Animal husbandry
#          0 = Other
system("r.category LU1980")
system("r.reclass --o in=LU1980 out=LU1980b <<EOF
       2 = 1 Animal husbandry
       * = 0 Other
       end
       EOF")
system("r.mapcalc LU1980b=LU1980b")
system("r.category LU1980b rules=- <<EOF
       1:Animal husbandry
       0:Other
       EOF")
system("r.category LU1980b")
system("r.info LU1980b")
system("d.mon x0")
system("d.erase")
system("d.rast.leg LU1980b")

# Land use change ==============================================================
# Classes: 1 = with change, 2 = without change
# Code: LUdiff
system("r.category LU2009")
system("r.category LU1980")
system("r.reclass --o in=LU1980 out=tmpLU1980 <<EOF
       1 = 1
       2 = 3
       3 = 5
       4 = 6
       5 = 7
       end
       EOF")
system("r.mapcalc 'LUdiff = if(LU2009 == tmpLU1980, 0, 1)'")
system("r.category LUdiff rules=- <<EOF
       1:with change
       0:without change
       EOF")
system("r.category LUdiff")
system("r.info LUdiff")
system("d.mon x0")
system("d.erase")
system("d.rast.leg LUdiff")
system("g.remove rast=tmpLU1980")

# SAVE FIGURES WITH LAND USE MAPS ##############################################
# Use as boundary the basin estimated with the digital elevation models derived
# from the contour lines plus the 30 m buffer (buffer_BASIN_10).
system("r.mask -o buffer_BASIN_10")
boundary <- readVECT6("buffer_BASIN_10")@bbox
# LU1980
map <- readRAST6("LU1980")
map@bbox <- boundary
map@data[, 1] <- as.factor(map@data[, 1])
levels(map@data[, 1])
system("r.category LU1980")
data(worldgrids_pal)
names(worldgrids_pal$glc2000)
map@data[, 1] <- revalue(map@data[, 1], c("1" = "FS",  # Semi-deciduous forest
                                          "2" = "H",   # animal husbandry
                                          "3" = "PF",  # Plantation forestry
                                          "4" = "S"))  # urban (Settlement)
color <- c(worldgrids_pal$glc2000[1], worldgrids_pal$glc2000[13],
           worldgrids_pal$glc2000[4], worldgrids_pal$glc2000[22])
p1 <- spplot(map, main = "", col.regions = color, asp = 1)
names(p1$legend) <- "inside"
p1$legend$inside$x <- 0.84
p1$legend$inside$y <- 0.875
p1$legend$inside$args$key$space <- "left"
dev.off()
pdf(file = paste(land.dir, "LU1980.pdf", sep = ""),
    width = 6.3/cm(1), height = (6.3/cm(1))*1.14)
trellis.par.set(fontsize = list(text = 7, points = 5),
                plot.line = list(lwd = 0.001),
                axis.line = list(lwd = 0.01),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
plot(p1)
dev.off()

rm(map, color, p1)
gc()

# LU2009
map <- readRAST6("LU2009")
map@bbox <- boundary
map@data[, 1] <- as.factor(map@data[, 1])
levels(map@data[, 1])
system("r.category LU2009")
map@data[, 1] <- revalue(map@data[, 1], c("1" = "FS",  # Semi-deciduous forest
                                          "2" = "SS",  # Semi-deciduous shrub
                                          "3" = "H",   # animal husbandry
                                          "4" = "AA",  # Crop agriculture (Annual field cropping)
                                          "5" = "PF",  # Plantation forestry
                                          "6" = "S",  # urban (Settlement)
                                          "7" = "O"))  # Other land uses (water reservoir)
data(worldgrids_pal)
names(worldgrids_pal$glc2000)
color <- c(worldgrids_pal$glc2000[1], worldgrids_pal$glc2000[12],
           worldgrids_pal$glc2000[13], worldgrids_pal$glc2000[16],
           worldgrids_pal$glc2000[4], worldgrids_pal$glc2000[22],
           worldgrids_pal$glc2000[20])
p1 <- spplot(map, main = "", col.regions = color, asp = 1)
names(p1$legend) <- "inside"
p1$legend$inside$x <- 0.84
p1$legend$inside$y <- 0.80
p1$legend$inside$args$key$space <- "left"
dev.off()
pdf(file = paste(land.dir, "LU2009.pdf", sep = ""),
    width = 6.3/cm(1), height = (6.3/cm(1))*1.14)
trellis.par.set(fontsize = list(text = 7, points = 5),
                plot.line = list(lwd = 0.001),
                axis.line = list(lwd = 0.01),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
plot(p1)
dev.off()

rm(map, color, p1)
gc()
rm(boundary)
gc()
system("g.remove MASK")

# SAVE DATA ####################################################################
ls()
save(LU1980.id, LU2009.id, land.dir,
     file = "sm-dnos-land-use.RData")

# End!