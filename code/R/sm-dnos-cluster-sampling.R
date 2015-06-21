//////////////////////////////////////////////////////////////////////
//                                                                  //
//           UNIVERSIDADE FEDERAL RURAL DO RIO DE JANEIRO           //
//                    INSTITUTO DE AGRONOMIA                        //
//        CURSO DE PÓS-GRADUAÇÃO EM AGRONOMIA-CIÊNCIA DO SOLO       //
//                                                                  //
//          CONTRIBUIÇÃO À CONSTRUÇÃO DE MODELOS DE PREDIÇÃO        //
//                      DE PROPRIEDADES DO SOLO                     //
//                                                                  //
//                      PROJETO DE DOUTORAMENTO                     //
//                                                                  //
//                       ALESSANDRO SAMUEL-ROSA                     //
//                                                                  //
//                     Seropédica, julho de 2013.                   //
//                                                                  //
//------------------------------------------------------------------//  
//                                                                  //
//Descrição:                                                        //
//Projeto de doutoramento apresentado ao Curso de Pós-Graduação em  //
//Agronomia-Ciência do Solo, da Universidade Federal Rural do Rio de//
//Janeiro (UFRRJ), Rio de Janeiro, como requisito parcial para a    //
//obtenção do grau de Doutor em Agronomia-Ciência do Solo.          //
//                                                                  //
//------------------------------------------------------------------//
//                                                                  //
//Comitê de orientação:                                             //
//Dra. Lúcia Helena Cunha dos Anjos - Orientador                    //
//Dr. Gustavo de Matos Vasques (Embrapa) - co-orientador            //
//Dr. Gerard Heuvelink (ISRIC) - co-orientador                      //
//                                                                  //
//------------------------------------------------------------------//
//                                                                  //
//                      Script de trabalho:                         //
//                       Cluster sampling                           //
//                                                                  //
//------------------------------------------------------------------//                                                         //                                                                  //
//Descrição: definição dos locais de amostragem para obtenção de    //
//observações de validação utilizando o método de amostragem  em    //
//cluster dentro de três extratos.                                  //
//                                                                  //
//------------------------------------------------------------------//
//                                                                  //
//e-mail: alessandrosamuel@yahoo.com.br                             //
//homepage: soil-scientist.net                                      //
//                                                                  //
//////////////////////////////////////////////////////////////////////
  
  
# Configuração inicial ###############################################
options()
options(prompt="> ", continue="+ ", digits = 5,
        width=70, show.signif.stars = T,
        verbose = TRUE)
par(mfrow=c(1,1))
ls()
rm(list=ls())
getwd()
setwd("/home/alessandro/PROJETOS/DNOS-SM/cluster-sampling")
list.files()

# Carregar os pacotes necessários ====================================
require(foreign)
require(sp)
require(raster)
require(rgdal)
require(RSAGA)
require(MASS)
require(plotKML)

# Configuração do ambiente de trabalho do RSAGA ======================
rsaga.env()
sagaenv <- rsaga.env()
sagaenv

# Jump to sampling

# Data processing ####################################################
ls()

# RSAGA usage ========================================================
rsaga.get.libraries()
rsaga.get.modules("grid_analysis", env = rsaga.env())
rsaga.search.modules("gdal")
rsaga.get.usage(lib = "io_gdal", module = "rsaga.import.gdal",
                env = sagaenv)

# DEM processing =====================================================

# Info----------------------------------------------------------------
system(paste("gdalinfo", "srtm-v41.tif", sep = " "))

# Warp and clip ------------------------------------------------------
system("gdalwarp")
# -s_srs: source srs
# -t_srs: target srs
# -te xmin ymin xmax ymax: extents of output file to be created
# -r: resampling method
# -of: output format
system(paste("gdalwarp",
             "-s_srs epsg:4326",
             "-t_srs epsg:31982",
             "-te 220817 6711252 237524 6724461",
             "-r cubic",
             "-of GTiff",
             "-overwrite",
             "srtm-v41.tif",
             "srtm-dem.tif", sep = " "))

# Save as SAGA format ------------------------------------------------
rsaga.import.gdal("srtm-dem.tif", "srtm-dem", env = sagaenv)

# Fill sinks ---------------------------------------------------------
rsaga.get.libraries()
rsaga.get.modules(lib = "ta_preprocessor", env = rsaga.env())
rsaga.get.usage(lib = "ta_preprocessor", module = 4,
                env = rsaga.env())
rsaga.geoprocessor(lib = "ta_preprocessor", module = 4,
                   param = list(ELEV = "srtm-dem.sgrd",
                                FILLED = "srtm-dem-fill.sgrd",
                                FDIR = "srtm-dem-fdir.sgrd",
                                WSHED = "srtm-dem-wshed.sgrd",
                                MINSLOPE = 0.01), env = sagaenv)
spplot(readGDAL("srtm-dem-fill.sdat"))

# Derive terrain attributes ==========================================

# Topographic position index -----------------------------------------
RADIUS = 540
spplot(readGDAL("tpi.sdat"), col.regions = terrain.colors(30))

# Morphometric Protection Index --------------------------------------
rsaga.get.libraries()
rsaga.get.modules(lib = "ta_morphometry")
rsaga.get.usage(lib = "ta_morphometry", module = 7, env = sagaenv)
rsaga.geoprocessor(lib = "ta_morphometry", module = 7,
                   param = list(DEM = "srtm-dem-fill.sgrd",
                                PROTECTION = "mpi.sgrd",
                                RADIUS = 180),
                   env = sagaenv)
spplot(readGDAL("mpi.sdat"), col.regions = terrain.colors(30))

# Slope, Aspect, Curvature -------------------------------------------
# 'Slope' and 'Aspect' are in radians
rsaga.get.modules(lib = "ta_morphometry")
rsaga.get.usage(lib = "ta_morphometry", module = 0, env = sagaenv)
rsaga.geoprocessor(lib = "ta_morphometry", module = 0,
                   param = list(ELEVATION = "srtm-dem-fill.sgrd",
                                SLOPE = "slope-rad.sgrd",
                                ASPECT = "aspect.sgrd",
                                METHOD = 5),
                   env = sagaenv)
spplot(readGDAL("slope-rad.sdat"), col.regions = terrain.colors(30))

# Convert slope-rad to slope-deg
fac = round(180/pi, 4)
formula = paste(fac, "*a", sep = "")
rsaga.grid.calculus("slope-rad.sgrd", "slope-deg.sgrd", formula)
rm(fac, formula)
spplot(readGDAL("slope-deg.sdat"), col.regions = terrain.colors(30))
summary(readGDAL("slope-deg.sdat"))

# Delineate strata using decision tree ===============================
# Not available in SAGA GIS 2.0.7
rsaga.get.libraries()
rsaga.get.modules(lib = "imagery_classification")

# Decision tree rules ------------------------------------------------
# tpi > 30: 'encosta'
# tpi < 30: 'outros'
### elev < 200 : 'outros'
##### mpi < 0.1: 'depressão'
##### mpi > 0.1: 'encosta'
### elev > 200 : 'outros'
##### elev > 400: 'outros'
####### mpi > 0.1: 'encosta'
####### mpi < 0.1: 'planalto'
##### elev < 400: 'encosta'

# aggregate strata ---------------------------------------------------
summary(as.factor(readGDAL("strata-tree.sdat")$band1))
spplot(readGDAL("strata-tree.sdat"), col.regions = terrain.colors(30))

# lookup table (how to make it in R?)
#change.tab <- matrix(NA, nrow = 7, ncol = 3)
#change.tab[,1] <- c("Low Value", 0, 4, 2, 14, 1, 6)
#change.tab[,2] <- c("High Value", 0, 4, 2, 14, 1, 6)
#change.tab[,3] <- c("Replace with", 1, 2, 2,  2, 2, 3)
#write(change.tab, "change-tab.txt", ncolumns = 3, sep = "\t")

# change values
rsaga.get.libraries()
rsaga.get.modules(lib = "grid_tools")
rsaga.get.usage(lib = "grid_tools", module = 12, env = sagaenv)
rsaga.geoprocessor(lib = "grid_tools", module = 12,
                   param = list(GRID_IN = "strata-tree.sgrd",
                                GRID_OUT = "strata.sgrd",
                                METHOD = 0,
                                LOOKUP = "change-tab.txt"),
                   env = sagaenv)
summary(as.factor(readGDAL("strata.sdat")$band1))
spplot(readGDAL("strata.sdat"), col.regions = terrain.colors(30))

# Majority filter ----------------------------------------------------
rsaga.get.libraries()
rsaga.get.modules(lib = "grid_filter")
rsaga.get.usage(lib = "grid_filter", module = 6, env = sagaenv)
rsaga.geoprocessor(lib = "grid_filter", module = 6,
                   param = list(INPUT = "strata.sgrd",
                                RESULT = "strata-filter.sgrd",
                                MODE = 2,
                                RADIUS = 3),
                   env = sagaenv)
spplot(readGDAL("strata-filter.sdat"),
       col.regions = terrain.colors(30))










# Sampling ###########################################################

# calculate the number of clusters per stratum =======================
strata <- as(readGDAL("strata-450m-masked.sdat"), "RasterLayer")
plot(strata)
(cells <- summary(as.factor(na.omit(strata@data@values))))
sum(cells)
(prop <- cells / sum(cells))
(ncluster <- as.vector(round(12 * prop)))

# process stratum no. 1 ==============================================
strata1 <- strata
strata1@data@values <- ifelse(strata@data@values == 1, 1, NA)

# define clusters ----------------------------------------------------
clus1 <- sampleRandom(x = strata1, size = ncluster[1],
                      xy = TRUE, asRaster = TRUE)
plot(strata1, colNA = "black")
plot(clus1, add = TRUE, col = "red")

# change cell values -------------------------------------------------
(cell <- which(clus1@data@values == 1))
clus1@data@values[cell[1]] <- 1
clus1@data@values[cell[2]] <- 2
which(clus1@data@values == 1)
which(clus1@data@values == 2)

# disaggregate clusters ----------------------------------------------
clus1.dis <- disaggregate(x = clus1, fact = 45, method = "")
  

# samples ------------------------------------------------------------
sample1 <- sampleStratified(clus1.dis, size = 5, xy = TRUE,
                            na.rm = TRUE, sp = TRUE)
plot(strata1)
plot(clus1.dis, add = TRUE, col = "red")
points(sample1, cex = 0.5, pch = 20)
proj4string(sample1) <- CRS("+init=epsg:31982")

# save shapefile and kml ---------------------------------------------
writeOGR(obj = sample1, dsn = getwd(), driver = "ESRI Shapefile",
         layer = "sample1")
kml(obj = sample1, folder.name = getwd(), file.name = "sample1.kml",
        colour = "yellow", alpha = 0.75,
        shape = "http://maps.google.com/mapfiles/kml/pal2/icon18.png")

# save coordinates ---------------------------------------------------
write(x = t(sample1@coords), file = "sample1.txt", ncolumns = 2,
      sep = "\t")

# process stratum no. 1 ==============================================
stratum2 <- strata
stratum2@data@values <- ifelse(strata@data@values == 2, 1, NA)

# define clusters ----------------------------------------------------
clus2 <- sampleRandom(x = stratum2, size = ncluster[2],
                      xy = TRUE, asRaster = TRUE)
plot(stratum2, colNA = "black")
plot(clus2, add = TRUE, col = "red")

# change cell values -------------------------------------------------
(cell <- which(clus2@data@values == 1))
clus2@data@values[cell[1]] <- 1
clus2@data@values[cell[2]] <- 2
clus2@data@values[cell[3]] <- 3
clus2@data@values[cell[4]] <- 4
clus2@data@values[cell[5]] <- 5
clus2@data@values[cell[6]] <- 6

# disaggregate clusters ----------------------------------------------
clus2.dis <- disaggregate(x = clus2, fact = 45, method = "")

# samples ------------------------------------------------------------
sample2 <- sampleStratified(clus2.dis, size = 5, xy = TRUE,
                            na.rm = TRUE, sp = TRUE)
plot(stratum2)
plot(clus2.dis, add = TRUE, col = "red")
points(sample2, cex = 0.5, pch = 20)
proj4string(sample2) <- CRS("+init=epsg:31982")

# save shapefile and kml ---------------------------------------------
writeOGR(obj = sample2, dsn = getwd(), driver = "ESRI Shapefile",
         layer = "sample2")
kml(obj = sample2, folder.name = getwd(), file.name = "sample2.kml",
    colour = "green", alpha = 0.75,
    shape = "http://maps.google.com/mapfiles/kml/pal2/icon18.png")

# save coordinates ---------------------------------------------------
write(x = t(sample2@coords), file = "sample2.txt", ncolumns = 2,
      sep = "\t")













# stratum 3 ==========================================================
strata.buffer3 <- strata.buffer
strata.buffer3@data@values <- ifelse(strata.buffer@data@values == 3,
                                     1, NA)
points3 <- sampleStratified(strata.buffer3, size = ncluster[3],
                            xy = TRUE, na.rm = TRUE, sp = TRUE)
plot(strata.buffer3)
points(points3[,2], points3[,3], pch = 20, cex = 0.6)

shapefile("points3", points3, overwrite = TRUE)

# load clusters ------------------------------------------------------
pts3 <- as(readGDAL("points3-buffer.sdat"), "RasterLayer")
plot(pts3, add = TRUE)

# sample clusters ----------------------------------------------------
pts3 <- sampleStratified(pts3, size = 5,
                         xy = TRUE, na.rm = TRUE, sp = TRUE)
points(pts3[,2], pts3[,3], pch = 20, cex = 0.6)
shapefile("pts3", pts3, overwrite = TRUE)

# plot all strata and points =========================================
strata <- as(readGDAL("strata.sdat"), "RasterLayer")
plot(strata)
plot(pts1, add = TRUE, pch = 1, cex = 0.5)
plot(pts2, add = TRUE, pch = 2, cex = 0.5)
plot(pts3, add = TRUE, pch = 3, cex = 0.5)

points(points1[,2], points1[,3], pch = 20, cex = 0.6, col = "red")
points(points2[,2], points2[,3], pch = 20, cex = 0.6, col = "red")
points(points3[,2], points3[,3], pch = 20, cex = 0.6, col = "red")


