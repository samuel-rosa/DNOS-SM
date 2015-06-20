# DESCRIPTION ##################################################################
# Analysis of the geological data available.
# SETTINGS #####################################################################
rm(list = ls())
gc()
require(sp)
require(rgdal)
require(maptools)
require(raster)
require(vec2dtransf)
require(spgrass6)
library(plyr)
library(lattice)
library(latticeExtra)
library(grid)
load("sm-dnos-geology.RData")
load("sm-dnos-general.RData")
ls()
# GRASS GIS
initGRASS(gisBase = "/usr/lib/grass64/", gisDbase = GRASSgisDbase,
          override = TRUE, location = "dnos-sm-rs", mapset = "predictions", 
          pid = Sys.getpid())
# Geology 1:50,000 =============================================================
GEO_50 <- shapefile(paste(geo.dir, "geology50.shp", sep = ""))
GEO_50$formation <- as.factor(GEO_50$formation)
GEO_50$name <- as.factor(GEO_50$name)
GEO_50 <- spTransform(GEO_50, wgs1984utm22s)
spplot(GEO_50, "formation")

# geology IDs
geo50.id <- data.frame(unique(GEO_50$code), unique(GEO_50$formation))
colnames(geo50.id) <- c("code", "GEO")
(geo50.id <- geo50.id[order(geo50.id$code),])

# apply affine transformation based on validation results
# see 'sm-dnos-covar-validation.R' for details
is.projected(GEO_50)
GEO_50 <- applyTransformation(affine.geo50, GEO_50)
plot(GEO_50)
getRMSE(affine.geo50)

# export to GRASS GIS (mapset = predictions)
initGRASS(gisBase = "/usr/lib/grass64/", gisDbase = GRASSgisDbase,
          override = TRUE, location = "dnos-sm-rs", mapset = "predictions", 
          pid = Sys.getpid())
system("g.region rast=dnos.raster")
writeVECT6(GEO_50, "GEO_50", v.in.ogr_flags = "overwrite")
system("v.info -c GEO_50")

# shapes to grid
system("v.to.rast --o in=GEO_50 out=GEO_50 use=attr column=code")
system("r.info GEO_50"); system("r.category GEO_50")
system("d.mon x0");system("d.rast.leg GEO_50")

# add category names
geo50.id
system("r.category GEO_50 rules=- <<EOF
       1:Serra Geral Formation - Inferior Sequence
       2:Botucatu Formation
       3:Serra Geral Formation - Superior Sequence
       4:Caturrita Formation
       5:Santa Maria Formation - Alemoa Member
       EOF")
system("r.category GEO_50")
system("d.mon x0");system("d.rast.leg GEO_50")

# match category names and codes with GEO_25 (bellow)
system("r.reclass --o in=GEO_50 out=recGEO_50 <<EOF
       1 = 1 Serra Geral Formation - Inferior Sequence
       2 = 2 Botucatu Formation
       3 = 4 Serra Geral Formation - Superior Sequence
       4 = 3 Caturrita Formation
       5 = 5 Santa Maria Formation - Alemoa Member
       end
       EOF")
system("r.category GEO_25")
system("r.category recGEO_50")
system("d.mon x0")
system("d.rast.leg recGEO_50")
system("d.mon x1")
system("d.rast.leg GEO_25")

# Geology 1:25,000 =============================================================
GEO_25 <- shapefile(paste(geo.dir, "geologia25.shp", sep = ""))
GEO_25 <- spTransform(GEO_25, wgs1984utm22s)
GEO_25
names(GEO_25)
GEO_25$formation <- as.factor(GEO_25$formation)
GEO_25$names <- as.factor(GEO_25$names)
spplot(GEO_25, "formation")

# geology IDs
geo25.id <- data.frame(unique(GEO_25$code), unique(GEO_25$formation))
colnames(geo25.id) <- c("code", "GEO")
(geo25.id <- geo25.id[order(geo25.id$code),])

# apply affine transformation based on validation results
# see 'sm-dnos-covar-validation.R' for details
is.projected(GEO_25)
GEO_25 <- applyTransformation(affine.geo25, GEO_25)
plot(GEO_25)
getRMSE(affine.geo25)

# export to GRASS GIS (mapset = predictions)
initGRASS(gisBase = "/usr/lib/grass64/", gisDbase = GRASSgisDbase,
          override = TRUE, location = "dnos-sm-rs", mapset = "predictions", 
          pid = Sys.getpid())
system("g.region rast=dnos.raster")
writeVECT6(GEO_25, "GEO_25", v.in.ogr_flags = "overwrite")
system("v.info -c GEO_25")

# shapes to grid
system("v.to.rast --o in=GEO_25 out=GEO_25 use=attr column=code")
system("r.info GEO_25"); system("r.category GEO_25")
system("d.mon x0");system("d.rast.leg GEO_25")

# add category names
geo25.id
system("r.category GEO_25 rules=- <<EOF
       1:Serra Geral Formation - Inferior Sequence
       2:Botucatu Formation
       3:Caturrita Formation
       4:Serra Geral Formation - Superior Sequence
       5:Santa Maria Formation - Alemoa Member
       EOF")
system("r.category GEO_25")
system("d.mon x0");system("d.rast.leg GEO_25")

# update GEO_25 using recGEO_50
system("r.patch --o input=GEO_25,recGEO_50 out=predGEO_25")
system("d.mon x0")
system("d.rast.leg predGEO_25")

# Quaternary Deposits 1:25,000 =================================================
DEP_25 <- shapefile(paste(geo.dir, "quaternary25.shp", sep = ""))
DEP_25 <- spTransform(DEP_25, wgs1984utm22s)
DEP_25
names(DEP_25)
DEP_25$geology <- as.factor(DEP_25$geology)
DEP_25$name <- as.factor(DEP_25$name)
spplot(DEP_25, "geology")

# geology IDs
dep25.id <- data.frame(unique(DEP_25$code), unique(DEP_25$geology))
colnames(dep25.id) <- c("code", "GEO")
(dep25.id <- dep25.id[order(dep25.id$code),])

# apply affine transformation based on validation results
# see 'sm-dnos-covar-validation.R' for details
is.projected(DEP_25)
DEP_25 <- applyTransformation(affine.geo25, DEP_25)
plot(DEP_25)

# export to GRASS GIS (mapset = predictions)
initGRASS(gisBase = "/usr/lib/grass64/", gisDbase = GRASSgisDbase,
          override = TRUE, location = "dnos-sm-rs", mapset = "predictions", 
          pid = Sys.getpid())
system("g.region rast=dnos.raster")
writeVECT6(DEP_25, "DEP_25", v.in.ogr_flags = "overwrite")
system("v.info -c DEP_25")

# shapes to grid
system("v.to.rast --o in=DEP_25 out=DEP_25 use=attr column=code")
system("r.info DEP_25"); system("r.category DEP_25")
system("d.mon x0");system("d.rast.leg DEP_25")

# add category names
system("r.category DEP_25");dep25.id
system("r.category DEP_25 rules=- <<EOF
       1:colluvial deposits
       2:fluvial flood-plain deposits
       3:fluvial terrace deposits
       4:consolidated rocks
       EOF")
system("r.category DEP_25")
system("d.mon x0");system("d.rast.leg DEP_25")

# update DEP_25 with GEO_50
system("r.mapcalc 'DEP_50=recGEO_50'")
system("r.reclass --o in=DEP_50 out=recDEP_50 <<EOF
       * = 4 consolidated rocks
       end
       EOF")
system("d.mon x0");system("d.rast.leg recDEP_50")
system("r.patch --o input=DEP_25,recDEP_50 out=predDEP_25")
system("d.mon x0");system("d.rast.leg predDEP_25")
system("g.remove rast=DEP_50,recDEP_50")

# Geological Faults 1:50,000 ===================================================
FAU_50 <- shapefile(paste(geo.dir, "falhas.shp", sep = ""))
FAU_50 <- spTransform(FAU_50, wgs1984utm22s)
FAU_50
names(FAU_50)
spplot(FAU_50, "length")

# apply affine transformation based on validation results
# see 'sm-dnos-covar-validation.R' for details
is.projected(FAU_50)
FAU_50 <- applyTransformation(affine.geo50, FAU_50)
plot(FAU_50)

# export to GRASS GIS (mapset = predictions)
initGRASS(gisBase = "/usr/lib/grass64/", gisDbase = GRASSgisDbase,
          override = TRUE,
          location = "dnos-sm-rs", mapset = "predictions", pid = Sys.getpid())
system("g.region rast=dnos.raster")
writeVECT6(FAU_50, "FAU_50", v.in.ogr_flags = "overwrite")
system("v.info -c FAU_50")
system("d.mon x0");system("d.vect FAU_50")

# CREATE ENVIRONMENTAL COVARIATES ##############################################
initGRASS(gisBase = "/usr/lib/grass64/", gisDbase = GRASSgisDbase,
          override = TRUE, location = "dnos-sm-rs",
          mapset = "predictions", pid = Sys.getpid())
system("g.region rast=dnos.raster")
system("g.remove MASK")

# GEO_50 =======================================================================
system("r.category recGEO_50")
system("d.mon x0")
system("d.rast.leg recGEO_50")
system("d.vect map=BASIN_10 fcolor=none")
# Santa Maria Formation - Alemoa Member does not occur in the study area. Thus,
# the number of unique classes is p - 1 - 1 = 3. Other classes can be produced
# by the combination of unique classes.

# GEO_50a (Serra Geral Formation - Inferior Sequence)
# Classes: 1 = SG Inf
#          0 = Other (SG Sup + Bot + Cat)
system("r.category recGEO_50")
system("r.reclass --o in=recGEO_50 out=GEO_50a <<EOF
       1 = 1 Serra Geral Inferior
       * = 0 Other
       end
       EOF")
system("r.mapcalc 'GEO_50a=GEO_50a'")
system("r.category GEO_50a rules=- <<EOF
       1:Serra Geral Inferior
       0:Other
       EOF")
system("r.category GEO_50a")
system("r.info GEO_50a")
system("d.mon x0")
system("d.erase")
system("d.rast.leg GEO_50a")

# GEO_50b (Serra Geral Formation - Superior Sequence)
# Classes: 1 = SG Sup
#          0 = Other (SG Inf + Bot + Cat)
system("r.category recGEO_50")
system("r.reclass --o in=recGEO_50 out=GEO_50b <<EOF
       4 = 1 SG Sup
       * = 0 Other
       end
       EOF")
system("r.mapcalc 'GEO_50b=GEO_50b'")
system("r.category GEO_50b rules=- <<EOF
       1:Serra Geral Superior
       0:Other
       EOF")
system("r.category GEO_50b")
system("r.info GEO_50b")
system("d.mon x0")
system("d.erase")
system("d.rast.leg GEO_50b")

# GEO_50c (Botucatu Formation)
# Classes: 1 = Bot
#          0 = Other (SG Sup + SG Inf + Cat)
system("r.category recGEO_50")
system("r.reclass --o in=recGEO_50 out=GEO_50c <<EOF
       2 = 1 Bot
       * = 0 Other
       end
       EOF")
system("r.mapcalc 'GEO_50c=GEO_50c'")
system("r.category GEO_50c rules=- <<EOF
       1:Botucatu
       0:Other
       EOF")
system("r.category GEO_50c")
system("r.info GEO_50c")
system("d.mon x0")
system("d.erase")
system("d.rast.leg GEO_50c")

# GEO_25 (predGEO_25) ==========================================================
system("r.category predGEO_25")
system("d.mon x0")
system("d.rast.leg predGEO_25")
system("d.vect map=BASIN_10 fcolor=none")
# Santa Maria Formation - Alemoa Member does not occur in the study area. Thus,
# the number of unique classes is p - 1 - 1 = 3. Other classes can be produced
# by the combination of unique classes.

# GEO_25a (Serra Geral Formation - Inferior Sequence)
# Classes: 1 = SG Inf
#          0 = Other (SG Sup + Bot + Cat)
system("r.category predGEO_25")
system("r.reclass --o in=predGEO_25 out=GEO_25a <<EOF
       1 = 1 SG Inf
       * = 0 Other
       end
       EOF")
system("r.mapcalc 'GEO_25a=GEO_25a'")
system("r.category GEO_25a rules=- <<EOF
       1:Serra Geral Inferior
       0:Other
       EOF")
system("r.category GEO_25a")
system("r.info GEO_25a")
system("d.mon x0")
system("d.erase")
system("d.rast.leg GEO_25a")

# GEO_25b (Serra Geral Formation - Superior Sequence)
# Classes: 1 = SG Sup
#          0 = Other (SG Inf + Bot + Cat)
system("r.category predGEO_25")
system("r.reclass --o in=predGEO_25 out=GEO_25b <<EOF
       4 = 1 SG Sup
       * = 0 Other
       end
       EOF")
system("r.mapcalc 'GEO_25b=GEO_25b'")
system("r.category GEO_25b rules=- <<EOF
       1:Serra Geral Superior
       0:Other
       EOF")
system("r.category GEO_25b")
system("r.info GEO_25b")
system("d.mon x0")
system("d.erase")
system("d.rast.leg GEO_25b")

# GEO_25c (Botucatu Formation)
# Classes: 1 = Bot
#          0 = Other (SG Sup + SG Inf + Cat)
system("r.category predGEO_25")
system("r.reclass --o in=predGEO_25 out=GEO_25c <<EOF
       2 = 1 Bot
       * = 0 Other
       end
       EOF")
system("r.mapcalc 'GEO_25c=GEO_25c'")
system("r.category GEO_25c rules=- <<EOF
       1:Botucatu
       0:Other
       EOF")
system("r.category GEO_25c")
system("r.info GEO_25c")
system("d.mon x0")
system("d.erase")
system("d.rast.leg GEO_25c")

# DEP_25 (predDEP_25) ==========================================================

# GEO_25d (Quaternary Deposits)
# Classes: 1 = Dep
#          0 = Other
system("r.category predDEP_25")
system("r.reclass --o in=predDEP_25 out=GEO_25d <<EOF
       1 2 = 1 Dep
       * = 0 Other
       end
       EOF")
system("r.mapcalc 'GEO_25d=GEO_25d'")
system("r.category GEO_25d rules=- <<EOF
       1:Quaternary Deposits
       0:Other
       EOF")
system("r.category GEO_25d")
system("r.info GEO_25d")
system("d.mon x0")
system("d.erase")
system("d.rast.leg GEO_25d")

# SAVE FIGURES WITH  GEOLOGICAL MAPS ###########################################
# Use as boundary the basin estimated with the digital elevation models derived
# from the contour lines plus the 30 m buffer (buffer_BASIN_10).
system("r.mask -o buffer_BASIN_10")
boundary <- readVECT6("buffer_BASIN_10")@bbox
# GEO_50
map <- readRAST6("recGEO_50")
map@bbox <- boundary
map@data[, 1] <- as.factor(map@data[, 1])
levels(map@data[, 1])
system("r.category recGEO_50")
map@data[, 1] <- revalue(map@data[, 1], c("1" = "SG-I",  # Serra Geral (inferior)
                                          "2" = "BT",    # Botucatu
                                          "3" = "CT",    # Caturrita
                                          "4" = "SG-S")) # Serra Geral (superior)
p1 <- spplot(map, main = "", col.regions = topo.colors(4), asp = 1)
names(p1$legend) <- "inside"
p1$legend$inside$x <- 0.795
p1$legend$inside$y <- 0.875
p1$legend$inside$args$key$space <- "left"
dev.off()
pdf(file = paste(geo.dir, "GEO_50.pdf", sep = ""),
    width = 6.3/cm(1), height = (6.3/cm(1))*1.14)
trellis.par.set(fontsize = list(text = 7, points = 5),
                plot.line = list(lwd = 0.001),
                axis.line = list(lwd = 0.01),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
plot(p1)
dev.off()
rm(map, p1)
gc()
# GEO_25
map <- readRAST6("predGEO_25")
map@bbox <- boundary
map@data[, 1] <- as.factor(map@data[, 1])
levels(map@data[, 1])
system("r.category predGEO_25")
map@data[, 1] <- revalue(map@data[, 1], c("1" = "SG-I",  # Serra Geral (inferior)
                                          "2" = "BT",    # Botucatu
                                          "3" = "CT",    # Caturrita
                                          "4" = "SG-S")) # Serra Geral (superior)
map <- as(map, "SpatialPixelsDataFrame")
system("r.to.vect --o in=GEO_25d out=poly feature=area")
poly <- readVECT6("poly")
poly@bbox <- boundary
col <- ifelse(poly$value == 0, "transparent", " black")
fill <- ifelse(poly$value == 0, "transparent", "lightgray")
p1 <- spplot(map, main = "", col.regions = topo.colors(4), asp = 1, 
             sp.layout = list("sp.polygons", poly, first = FALSE, 
                              col = col, fill = fill, alpha = 0.5))
names(p1$legend) <- "inside"
p1$legend$inside$x <- 0.795
p1$legend$inside$y <- 0.875
p1$legend$inside$args$key$space <- "left"
p1$legend$inside$args$key$at <- c(p1$legend$inside$args$key$at, 5.5)
p1$legend$inside$args$key$col <- c(p1$legend$inside$args$key$col, "lightgray")
p1$legend$inside$args$key$labels$at <- c(p1$legend$inside$args$key$labels$at, 5)
p1$legend$inside$args$key$labels$labels <- c(p1$legend$inside$args$key$labels$labels, "QD")
dev.off()
pdf(file = paste(geo.dir, "GEO_25.pdf", sep = ""),
    width = 6.3/cm(1), height = (6.3/cm(1))*1.14)
trellis.par.set(fontsize = list(text = 7, points = 5),
                plot.line = list(lwd = 0.001),
                axis.line = list(lwd = 0.01),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
plot(p1)
dev.off()

rm(cropped, poly, p1, col, fill)
gc()
rm(boundary)
system("g.remove MASK")

# SAVE DATA ####################################################################
ls()
save(affine.geo25, affine.geo50, DEP_25, dep25.id,
     FAU_50, GEO_25, geo25.id, GEO_50, geo50.id,
     file = "sm-dnos-geology.RData")
# End!