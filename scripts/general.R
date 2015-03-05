////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//               UNIVERSIDADE FEDERAL RURAL DO RIO DE JANEIRO                 //
//                        INSTITUTO DE AGRONOMIA                              //
//          CURSO DE PÓS-GRADUAÇÃO EM AGRONOMIA-CIÊNCIA DO SOLO               //
//                                                                            //
//            CONTRIBUIÇÃO À CONSTRUÇÃO DE MODELOS DE PREDIÇÃO                //
//                         DE PROPRIEDADES DO SOLO                            //
//                                                                            //
//                           PROJETO DE PESQUISA                              //
//                                                                            //
//                          ALESSANDRO SAMUEL-ROSA                            //
//                                                                            //
//                        Seropédica, março de 2012.                          //
//                                                                            //
//----------------------------------------------------------------------------//
//                                                                            //
//Descrição:                                                                  //
//Projeto de pesquisa apresentada ao Curso de Pós-Graduação                   //
//em Agronomia-Ciência do Solo, da Universidade Federal                       //
//Rural do Rio de Janeiro (UFRRJ), Rio de Janeiro, como                       //
//requisito parcial para a obtenção do grau de Doutor em                      //
//Agronomia-Ciência do Solo.                                                  //
//                                                                            //
//----------------------------------------------------------------------------//
//                                                                            //
//Comitê de orientação:                                                       //
//Dra. Lúcia Helena Cunha dos Anjos - Orientador                              //
//Dr. Gustavo de Matos Vasques (Embrapa) - co-orientador                      //
//Dr. Gerard Heuvelink (ISRIC) - co-orientador                                //
//                                                                            //
//----------------------------------------------------------------------------//
//                                                                            //
//                           Script de trabalho:                              //
//                      Configuração geral dos dados                          //
//                                                                            //
//----------------------------------------------------------------------------//
//                                                                            //
//Descrição:                                                                  //
//                                                                            //
//                                                                            //
//----------------------------------------------------------------------------//
//                                                                            //
//e-mail: alessandrosamuel@yahoo.com.br                                       //
//homepage: soil-scientist.net                                                //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

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
library(raster)
library(rgeos)
load("sm-dnos-general.RData")
ls()
# SET DATA DIRECTORIES AND DEFINITIONS #########################################

# exploratory data analysis
explora.dir <- path.expand("~/PROJECTS/DNOS-SM/exploratory/")
# covar-validation directory
covar.validation.dir <- path.expand("~/PROJECTS/DNOS-SM/covar-validation/")
# rdata directory
rdata.dir <- path.expand("~/PROJECTS/DNOS-SM/rdata/")
# hydrological data
hydro.dir <- path.expand("~/PROJECTS/DNOS-SM/hydrology/")
# DEM and DEM derivatives
dem.dir <- path.expand("~/PROJECTS/DNOS-SM/terrain-attributes/")
# Geological maps
geo.dir <- path.expand("~/PROJECTS/DNOS-SM/geology/")
# Land use maps
land.dir <- path.expand("~/PROJECTS/DNOS-SM/land-use/")
# Soil maps
soil.dir <- path.expand("~/PROJECTS/DNOS-SM/soil-class/")
# Point data
point_dir <- path.expand("~/PROJECTS/DNOS-SM/point-data/")
# Resampling methors for grids
resample.dir <- path.expand("~/PROJECTS/DNOS-SM/grid-resample/")
# Boudaries
bounding.dir <- path.expand("~/PROJECTS/DNOS-SM/boundaries/")
# Updating environmental covariates (chapter/article 1)
update_covars_dir <- path.expand("~/PROJECTS/DNOS-SM/update-covars/")
# Point pattern analysys
ppp_dir <- path.expand("~/PROJECTS/DNOS-SM/point-pattern-analysis/")

# Prediction grid cell size
# cellsize <- 15
cellsize <- 5

# GRASS GIS addons
Sys.getenv("GRASS_ADDON_PATH")
Sys.setenv("GRASS_ADDON_PATH" = path.expand("~/softwares/GRASS_ADDON_PATH"))
grass.addons <- path.expand("~/softwares/GRASS_ADDON_PATH/")
system(paste("g.extension r.gauss")) # install
system(paste(grass.addons, "r.gauss -help", sep = ""))

# GRASS GIS DBase
GRASSgisDbase <- path.expand("~/GRASSgisDbase")
initGRASS(gisBase = "/usr/lib/grass64/", gisDbase = GRASSgisDbase, 
          location = "dnos-sm-rs",
          mapset = "predictions", pid = Sys.getpid(), override = TRUE)
writeRAST6(dnos.raster, "dnos.raster", overwrite = TRUE)
system("g.region rast=dnos.raster")
gmeta6()

# Coordinate reference systems
sirgas2000 <- CRS("+init=epsg:4674")
sirgas2000utm22s <- CRS("+init=epsg:31982")
wgs1984utm22s <- CRS("+init=epsg:32722")
wgs1984 <- CRS("+init=epsg:4326")
ca_utm22s <- CRS("+init=epsg:22522")

# LOAD AND PROCESS DATA ########################################################

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

# EXPLORATORY DATA ANALYSIS
# Dependent variables should have Gaussian distribution - this is one of the
# requirements of linear models. It is expected that, if the dependent varible
# has a Gaussian distribution, then the residuals also will have a Gaussian
# distribution. Transformations are performed using the power family of
# Box-Cox transformations. Only non-negative lambda values are used. This is
# a requirement fo performing the back transformation using Monte Carlo 
# simulations. When lambda is negative, the integral over the entire space of
# the probability distribution can be infinite, and thus the mean and variance 
# cannot be computed. When a negative lambda is estimated, it is replaced by
# zero (0).

bc_lambda <- list(clay = NA, carbon = NA, ecec = NA)

source("/home/alessandro/PROJECTS/pedometrics/pedometrics/cooking/plotHD.R")

# Calibration Data - clay content ----------------------------------------------
var <- cal_data$clay
# histogram with original variable
xlab  <-  expression(paste('CLAY (g ',kg^-1,')', sep = ''))
tmp <- plotHD(var, HD = "over", xlab = xlab, BoxCox = FALSE, stats = FALSE,
              scales = list(cex = c(1, 1)))
dev.off()
pdf(file = paste(explora.dir, "clay-dist-original.pdf", sep = ""),
    width = 6.3/cm(1), height = 6.3/cm(1))
trellis.par.set(fontsize = list(text = 7, points = 5),
                plot.line = list(lwd = 0.001),
                axis.line = list(lwd = 0.01),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
print(tmp)
dev.off()
rm(tmp, xlab)
gc()
# histogram with transformed variable
xlab  <-  expression(paste('Box-Cox CLAY (g ',kg^-1,')', sep = ''))
tmp <- plotHD(var, HD = "over", xlab = xlab, BoxCox = TRUE, stats = FALSE,
              scales = list(cex = c(1, 1)))
dev.off()
pdf(file = paste(explora.dir, "clay-dist-trans.pdf", sep = ""),
    width = 6.3/cm(1), height = 6.3/cm(1))
trellis.par.set(fontsize = list(text = 7, points = 5),
                plot.line = list(lwd = 0.001),
                axis.line = list(lwd = 0.01),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
print(tmp)
dev.off()
rm(tmp, var, xlab)
gc()
# transformation
bc_lambda$clay <- lambda
cal_data$clay_bc <- bcPower(cal_data$clay, lambda)
dev.off()
pdf(file = paste(explora.dir, "clay_bc.pdf", sep = ""))
xyplot(clay_bc ~ clay, data = cal_data@data, xlab = "original scale",
       ylab = "Box-Cox transformed", main = "Clay content")
dev.off()
rm(var, par, lambda, leg)
gc()

# Calibration Data - carbon content --------------------------------------------
var <- cal_data$carbon
# histogram with original variable
xlab  <-  expression(paste('SOC (g ',kg^-1,')', sep = ''))
tmp <- plotHD(var, HD = "over", xlab = xlab, BoxCox = FALSE, stats = FALSE,
              scales = list(cex = c(1, 1)))
dev.off()
pdf(file = paste(explora.dir, "carbon-dist-original.pdf", sep = ""),
    width = 6.3/cm(1), height = 6.3/cm(1))
trellis.par.set(fontsize = list(text = 7, points = 5),
                plot.line = list(lwd = 0.001),
                axis.line = list(lwd = 0.01),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
print(tmp)
dev.off()
rm(tmp, xlab)
gc()
# histogram with transformed variable
xlab  <-  expression(paste('Box-Cox SOC (g ',kg^-1,')', sep = ''))
tmp <- plotHD(var, HD = "over", xlab = xlab, BoxCox = TRUE, stats = FALSE,
              scales = list(cex = c(1, 1)))
dev.off()
pdf(file = paste(explora.dir, "carbon-dist-trans.pdf", sep = ""),
    width = 6.3/cm(1), height = 6.3/cm(1))
trellis.par.set(fontsize = list(text = 7, points = 5),
                plot.line = list(lwd = 0.001),
                axis.line = list(lwd = 0.01),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
print(tmp)
dev.off()
rm(var, tmp, xlab)
gc()

# transformation
bc_lambda$carbon <- lambda
cal_data$carbon_bc <- bcPower(cal_data$carbon, lambda)
dev.off()
pdf(file = paste(explora.dir, "carbon_bc.pdf", sep = ""))
xyplot(carbon_bc ~ carbon, data = cal_data@data, xlab = "original scale",
       ylab = "Box-Cox transformed", main = "Carbon content")
dev.off()
rm(var, par, lambda, leg)
gc()

# Calibration Data - ECEC ------------------------------------------------------
var <- cal_data$ecec
# histogram with original variable
xlab  <-  expression(paste('ECEC (',mmol, ' ', kg^-1,')', sep = ''))
tmp <- plotHD(var, HD = "over", xlab = xlab, BoxCox = FALSE, stats = FALSE,
              scales = list(cex = c(1, 1)))
dev.off()
pdf(file = paste(explora.dir, "ecec-dist-original.pdf", sep = ""),
    width = 6.3/cm(1), height = 6.3/cm(1))
trellis.par.set(fontsize = list(text = 7, points = 5),
                plot.line = list(lwd = 0.001),
                axis.line = list(lwd = 0.01),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
print(tmp)
dev.off()
rm(tmp, xlab)
gc()
# Histogram with transformed variable
xlab  <-  expression(paste('Box-Cox ECEC (',mmol, ' ', kg^-1,')', sep = ''))
tmp <- plotHD(var, HD = "over", xlab = xlab, BoxCox = TRUE, stats = FALSE,
              scales = list(cex = c(1, 1)))
dev.off()
pdf(file = paste(explora.dir, "ecec-dist-trans.pdf", sep = ""),
    width = 6.3/cm(1), height = 6.3/cm(1))
trellis.par.set(fontsize = list(text = 7, points = 5),
                plot.line = list(lwd = 0.001),
                axis.line = list(lwd = 0.01),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
print(tmp)
dev.off()
rm(var, tmp, xlab)
gc()
# transformation
bc_lambda$ecec <- lambda
cal_data$ecec_bc <- bcPower(cal_data$ecec, lambda)
dev.off()
pdf(file = paste(explora.dir, "ecec_bc.pdf", sep = ""))
xyplot(ecec_bc ~ ecec, data = cal_data@data, xlab = "original scale",
       ylab = "Box-Cox transformed", main = "ECEC")
dev.off()
rm(var, par, lambda, leg)

# Venn diagram
install.packages("VennDiagram")
require(VennDiagram)
venn.diagram(list(CLAY = c("DEM", "GEO"), SOC = c("DEM", "LAND"), ECEC = c("LAND", "GEO")),
             fill = c("red", "green", "blue"),
             filename = "tmp.tiff")

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

# Calibration Data - sample rasters --------------------------------------------
# settings
initGRASS(gisBase = "/usr/lib/grass64/", gisDbase = GRASSgisDbase,
          location = "dnos-sm-rs", mapset = "predictions",
          pid = Sys.getpid(), override = TRUE)
system("g.region rast=dnos.raster")
system("g.remove MASK")
system("r.mask dnos.raster")
# import calibration points into GRASS
system("g.remove vect=calibration")
pts <- data.frame(coordinates(cal_data), cal_data$sampleid)
coordinates(pts) <- ~ longitude + latitude
proj4string(pts) <- wgs1984utm22s
colnames(pts@data) <- "sampleid"
writeVECT6(pts, "calibration", v.in.ogr_flags = "overwrite")
rm(pts)
# setup database of calibration points
system("v.info -c calibration")
cols_int <- paste(soil1, land1, geo1, soil2, land2, geo2, sep = " + ")
cols_int <- str_replace_all(cols_int, "[+]", "INT,")
cols_int <- paste(cols_int, " INT", sep = "")
cols_double <- paste(sat1, dem1,  sat2, dem2, sep = " + ")
cols_double <- str_replace_all(cols_double, "[+]", " DOUBLE PRECISION,")
cols_double <- paste(cols_double, " DOUBLE PRECISION", sep = "")
cols <- paste(cols_int, cols_double, sep = ", ")
cmd <- paste("v.db.addcol map=calibration columns='", cols, "'", sep = "")
system(cmd)
system("v.info -c calibration")
rm(cols_int, cols_double, cols)
# sample rasters
# do it in GRASS because loading all rasters in R consumes all the memory
maps <- paste(soil1, land1, geo1, soil2, land2, geo2, sat1, dem1, sat2, dem2, sep = " + ")
maps <- str_replace_all(maps, "[ ]", "")
maps <- c(unlist(str_split(maps, "[+]")))
column <- maps
cmd <- paste("v.what.rast vect=calibration raster=", maps, " column=", column, sep = "")
lapply(cmd, system)
rm(maps, column)
# read calibration data
tmp <- readVECT6(vname = "calibration")
str(tmp)
tmp@data <- tmp@data[, - 1]
tmp$sampleid <- as.character(tmp$sampleid)
str(tmp)
which(c(cal_data$sampleid == tmp$sampleid) == FALSE)
proj4string(tmp) <- wgs1984utm22s
tmp@data <- join(cal_data@data, tmp@data, by = "sampleid")
cal_data <- tmp
str(cal_data)
rm(tmp)

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

# LOCATION OF THE STUDY AREA ###################################################

# Brazil -----------------------------------------------------------------------
brazil <- shapefile("/home/alessandro/PROJECTS/DNOS-SM/boundaries/brasil.shp")
bb <- bbox(brazil)
bb[1, 2] <- -34.7
brazil@bbox <- bb
pts <- data.frame(rbind(c(-53.790215, -29.657668), c(-47.887768, -15.788838),
                        c(-46.665337, -23.536445), c(-51.211134, -30.032484),
                        c(-59.992370, -3.080757)))
colnames(pts) <- c("long", "lat")
coordinates(pts) <- ~ long + lat
proj4string(pts) <- proj4string(brazil)
pts <- list("sp.points", pts, pch = 20, cex = 0.5, col = "black")
# scl <- list("SpatialPolygonsRescale", layout.scale.bar(), 
#             offset = c(-47, -31.5), scale = 10, 
#             fill = c("transparent", "black"))
# text1   <- list("sp.text", c(-47, -32.5), "0 deg", cex = 1.5)
# text2   <- list("sp.text", c(-37, -32.5), "10 deg", cex = 1.5)
# arrow <- list("SpatialPolygonsRescale", layout.north.arrow(type = 1), 
#               offset = c(-37.7, 1.5), scale = 2)
# set map colors
brazil$UF_05 <- as.factor(brazil$UF_05)
rs <- rep("lightgray", length(brazil$UF_05))
rs[which(brazil$UF_05 == "RS")] = "darkgray"
# prepare spplot
p <- spplot(brazil, zcol = "UF_05", aspect = "iso", col = "gray", 
            scales = list(draw = TRUE, 
                          x = list(at = seq(-70, -35, 5)),
                          y = list(at = seq(-30, 5, 5))),
            col.regions = colorRampPalette(rs)(27), colorkey = FALSE, 
            cex = 0.3, sp.layout = list(pts),
            par.settings = list(fontsize = list(text = 7, points = 5),
                                layout.widths = list(left.padding = 0, 
                                                     right.padding = 0), 
                                layout.heights = list(top.padding = 0,
                                                      bottom.padding = 0)),
            panel = function(x, y, ...) {
              panel.polygonsplot(x, y, ...)
              panel.abline(h = seq(-30, 0, 5), v = seq(-70, -40, 5),
                           col = "gray", lty = "dashed", lwd = 0.5)
              panel.text(x = -53.790215, y = -29.657668, "Santa Maria", pos = 2)
              panel.text(x = -47.887768, y = -15.788838, "Brasília", pos = 4)
              panel.text(x = -46.665337, y = -23.536445, "São Paulo", pos = 4)
              panel.text(x = -51.211134, y = -30.032484, "Porto Alegre", pos = 4)
              panel.text(x = -59.992370, y = -3.080757, "Manaus", pos = 4)
              }
            )
p
# save image
dev.off()
pdf(file = paste(explora.dir, "brazil.pdf", sep = ""), width = 9/cm(1),
    height = 9/cm(1))
# png(file = paste(explora.dir, "brazil.png", sep = ""), width = 255,
#     height = 255, res = 300, type = "cairo", antialias = "gray")
print(p)
dev.off()

rm(brazil, bb, pts, scl, text1, text2, rs, p, arrow)
gc()

# Calibration observations -----------------------------------------------------
pol <- readVECT6("buffer_BASIN_10")
drain <- readVECT6("STREAM_10")
drain <- gIntersection(drain, pol, byid = TRUE)
# scl <- list("SpatialPolygonsRescale", layout.scale.bar(), 
#             offset = c(230500, 6715600), scale = 1000, 
#             fill = c("transparent", "black"))
# text1 = list("sp.text", c(230500, 6715400), "0 m", cex = 1.5)
# text2 = list("sp.text", c(231500, 6715400), "500 m", cex = 1.5)
p <- spplot(pol, zcol = "cat", col = "gray", fill = "lightgray",
            scales = list(draw = TRUE),
            colorkey = FALSE, aspect = "iso",
            xlim = c(bbox(cal_data)[1, 1] * 0.9998, bbox(cal_data)[1, 2] * 1.0005), 
            ylim = c(bbox(cal_data)[2, 1] * 0.99999, bbox(cal_data)[2, 2] * 1.00001),
            par.settings = list(fontsize = list(text = 7, points = 5),
                                layout.widths = list(left.padding = 0, 
                                                     right.padding = 0), 
                                layout.heights = list(top.padding = 0,
                                                      bottom.padding = 0)),
            panel = function(x, y, ...) {
              panel.polygonsplot(x, y, ...)
              panel.points(x = coordinates(cal_data)[, 1],
                           y = coordinates(cal_data)[, 2],
                           pch = 20, cex = 0.5, col = "black")
              panel.abline(v = seq(227000, 232000, 1000), 
                           h = seq(6712000, 6722000, 1000),
                           col = "gray", lty = "dashed", lwd = 0.5)
              }) + layer(sp.lines(drain, col = "black", lty = "dashed", 
                                  lwd = 0.3))
p
# save image
dev.off()
pdf(file = paste(explora.dir, "cal-data.pdf", sep = ""), width = 9/cm(1),
    height = 9/cm(1))
# png(file = paste(explora.dir, "cal-data.png", sep = ""), width = 255,
#     height = 255, type = "cairo", antialias = "gray")
print(p)
dev.off()
rm(pol, scl, text1, text2, p, drain)
gc()

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