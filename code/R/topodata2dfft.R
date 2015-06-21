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
//                2D-FFT do MDE do projeto TOPODATA                 //
//                                                                  //
//------------------------------------------------------------------//
//                                                                  //
//                                                                  //
//Descrição: avaliação da qualidade do modelo digital de elevação   //
//do projeto Topodata por meio do uso da transformada de Fourier em //
//duas dimensões. O trabalho original foi publicado no Congresso    //
//Brasileiro de Ciência do Solo de 2013, realizado em Florianópolis,//
//Santa Catarina.                                                   //
//                                                                  //
//------------------------------------------------------------------//
//                                                                  //
//e-mail: alessandrosamuel@yahoo.com.br                             //
//homepage: soil-scientist.net                                      //
//                                                                  //
//////////////////////////////////////////////////////////////////////


# Configuração inicial da área de trabalho ###########################
options()
options(prompt="> ", continue="+ ", digits = 5,
        width=70, show.signif.stars = T,
        verbose = TRUE)
par(mfrow=c(1,1))
ls()
rm(list=ls())

# Diretório de trabalho ==============================================
getwd()
setwd("/home/alessandro/PROJETOS/DNOS-SM/2d-fast-fourier")
list.files()
ls()

# Carregar os pacotes necessários ====================================
require(sp)
require(raster)
require(rgdal)
require(spgrass6)
require(RSAGA)
require(MASS)
require(lattice)

# Definições do SAGA GIS =============================================
sagaenv <- rsaga.env()

# Definições do GRASS GIS ============================================
initGRASS(gisBase = "/usr/lib/grass64",
          home = tempdir(), override = TRUE)

# Reprojetar MDEs (TOPODATA & SRTM) ##################################
ls()
list.files()

# Área de corte
bbox(dnos.raster)

# Reprojetar MDE TOPODATA ============================================

# Info----------------------------------------------------------------
system(paste("gdalinfo", "topodata_original.tif",
             sep = " "))
system("gdalwarp")

# Nearest Neighbour---------------------------------------------------
system(paste("gdalwarp", "-s_srs epsg:4326",
             "-t_srs epsg:31982",
             "-te 226000 6715000 232000 6721990",
             "-r near", "-of GTiff",
             "topodata_original.tif",
             "topodata_near.tif", sep = " "))

# Bilinear------------------------------------------------------------
system(paste("gdalwarp", "-s_srs epsg:4326",
             "-t_srs epsg:31982",
             "-te 226000 6715000 232000 6721990",
             "-r bilinear", "-of GTiff",
             "topodata_original.tif",
             "topodata_bilinear.tif", sep = " "))

# Cubic---------------------------------------------------------------
system(paste("gdalwarp", "-s_srs epsg:4326",
             "-t_srs epsg:31982",
             "-te 226000 6715000 232000 6721990",
             "-r cubic", "-of GTiff",
             "topodata_original.tif",
             "topodata_cubic.tif", sep = " "))

# Cubicspline---------------------------------------------------------
system(paste("gdalwarp", "-s_srs epsg:4326",
             "-t_srs epsg:31982",
             "-te 226000 6715000 232000 6721990",
             "-r cubicspline", "-of GTiff",
             "topodata_original.tif",
             "topodata_cubicspline.tif", sep = " "))

# Lanczos-------------------------------------------------------------
system(paste("gdalwarp", "-s_srs epsg:4326",
             "-t_srs epsg:31982",
             "-te 226000 6715000 232000 6721990",
             "-r lanczos", "-of GTiff",
             "topodata_original.tif",
             "topodata_lanczos.tif", sep = " "))

# Reprojetar MDE SRTM ================================================

# Info----------------------------------------------------------------
system(paste("gdalinfo", "srtm_original.tif",
             sep = " "))
system("gdalwarp")

# Nearest Neighbour---------------------------------------------------
system(paste("gdalwarp", "-s_srs epsg:4326",
             "-t_srs epsg:31982",
             "-te 226000 6715000 232000 6721990",
             "-r near", "-of GTiff",
             "srtm_original.tif",
             "srtm_near.tif", sep = " "))

# Bilinear------------------------------------------------------------
system(paste("gdalwarp", "-s_srs epsg:4326",
             "-t_srs epsg:31982",
             "-te 226000 6715000 232000 6721990",
             "-r bilinear", "-of GTiff",
             "srtm_original.tif",
             "srtm_bilinear.tif", sep = " "))

# Cubic---------------------------------------------------------------
system(paste("gdalwarp", "-s_srs epsg:4326",
             "-t_srs epsg:31982",
             "-te 226000 6715000 232000 6721990",
             "-r cubic", "-of GTiff",
             "srtm_original.tif",
             "srtm_cubic.tif", sep = " "))

# Cubicspline---------------------------------------------------------
system(paste("gdalwarp", "-s_srs epsg:4326",
             "-t_srs epsg:31982",
             "-te 226000 6715000 232000 6721990",
             "-r cubicspline", "-of GTiff",
             "srtm_original.tif",
             "srtm_cubicspline.tif", sep = " "))

# Lanczos-------------------------------------------------------------
system(paste("gdalwarp", "-s_srs epsg:4326",
             "-t_srs epsg:31982",
             "-te 226000 6715000 232000 6721990",
             "-r lanczos", "-of GTiff",
             "srtm_original.tif",
             "srtm_lanczos.tif", sep = " "))


# 2D Fast Fourier Transform - TOPODATA ###############################
ls()

# Nearest neighbour resampling =======================================

# Ler arquivo de dados------------------------------------------------
execGRASS("r.in.gdal", flags = "o",
          parameters = list(input = "topodata_near.tif",
                            output = "DEM"))

# Definir região de trabalho------------------------------------------
execGRASS("g.region", parameters = list(rast="DEM"))

# 2D FFT--------------------------------------------------------------
execGRASS("i.fft", flags = "verbose",
          parameters = list(input_image = "DEM",
                            real_image = "2Dfft_r",
                            imaginary_image = "2Dfft_i"))

# Definir região de trabalho------------------------------------------
execGRASS("g.region", parameters = list(rast="2Dfft_r"))

# Ler arquivo de dados------------------------------------------------
topodata.near.fft <- readRAST6(vname = "2Dfft_r",
                           plugin = FALSE)

# Informações sobre o arquivo de dados--------------------------------
summary(topodata.near.fft)
str(topodata.near.fft)
hist(topodata.near.fft$X2Dfft_r)

# Plotar--------------------------------------------------------------
topodata.near.fft$X2Dfft_r01 <- abs(topodata.near.fft$X2Dfft_r)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
topodata.near.fft$X2Dfft_r01 <- range01(topodata.near.fft$X2Dfft_r01)
p.topo.near <-
  spplot(topodata.near.fft, zcol="X2Dfft_r01",
         col.regions = gray(seq(0, 1, by=0.01)),
         aspect="iso", colorkey=TRUE,
         main = "1 arc-second TOPODATA DEM",
         sub = "Nearest neighbour resampling")
print(p.topo.near)
rm(range01)

# Plotar reverse scale------------------------------------------------
tmp <-
  spplot(topodata.near.fft, zcol="X2Dfft_r01",
         col.regions = rev(gray(seq(0, 1, by=0.01))),
         aspect="iso", colorkey=TRUE,
         main = "1 arc-second TOPODATA DEM",
         sub = "Nearest neighbour resampling")
print(tmp)
dev.off()
png(filename="tmp.png")
print(p.topo.near, split = c(1,1,2,1), more = T)
print(tmp, split = c(2,1,2,1), more = F)
dev.off()
rm(tmp)

# Bilinear resampling ================================================

# Configurar o GRASS GIS-----------------------------
initGRASS(gisBase = "/usr/lib/grass64",
          home = tempdir(), override=TRUE)

# Ler arquivo de dados-------------------------------
execGRASS("r.in.gdal", flags = "o",
          parameters = list(input = "topodata_bilinear.tif",
                            output = "DEM"))

# Definir região de trabalho-------------------------
execGRASS("g.region", parameters = list(rast="DEM"))

# 2D FFT---------------------------------------------
execGRASS("i.fft", flags = "verbose",
          parameters = list(input_image = "DEM",
                            real_image = "2Dfft_r",
                            imaginary_image = "2Dfft_i"))

# Definir região de trabalho-------------------------
execGRASS("g.region", parameters = list(rast="2Dfft_r"))

# Ler arquivo de dados-------------------------------
topodata.bilinear.fft <- readRAST6(vname = "2Dfft_r",
                           plugin = FALSE)

# Informações sobre o arquivo de dados---------------
summary(topodata.bilinear.fft)
str(topodata.bilinear.fft)
hist(topodata.bilinear.fft$X2Dfft_r)

# Plotar---------------------------------------------
topodata.bilinear.fft$X2Dfft_r01 <- abs(topodata.bilinear.fft$X2Dfft_r)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
topodata.bilinear.fft$X2Dfft_r01 <- range01(topodata.bilinear.fft$X2Dfft_r01)
p.topo.bi <-
  spplot(topodata.bilinear.fft, zcol="X2Dfft_r01",
         col.regions=gray(seq(0, 1, by=0.01)),
         aspect="iso", colorkey=TRUE,
         main = "1 arc-second TOPODATA DEM",
         sub = "Bilinear resampling")
print(p.topo.bi)
rm(range01)

#----------------------------------------------------
# Cubic resampling
#----------------------------------------------------

# Configurar o GRASS GIS-----------------------------
initGRASS(gisBase = "/usr/lib/grass64",
          home = tempdir(), override=TRUE)

# Ler arquivo de dados-------------------------------
execGRASS("r.in.gdal", flags = "o",
          parameters = list(input = "topodata_cubic.tif",
                            output = "DEM"))

# Definir região de trabalho-------------------------
execGRASS("g.region", parameters = list(rast="DEM"))

# 2D FFT---------------------------------------------
execGRASS("i.fft", flags = "verbose",
          parameters = list(input_image = "DEM",
                            real_image = "2Dfft_r",
                            imaginary_image = "2Dfft_i"))

# Definir região de trabalho-------------------------
execGRASS("g.region", parameters = list(rast="2Dfft_r"))

# Ler arquivo de dados-------------------------------
topodata.cubic.fft <- readRAST6(vname = "2Dfft_r",
                               plugin = FALSE)

# Informações sobre o arquivo de dados---------------
summary(topodata.cubic.fft)
str(topodata.cubic.fft)
hist(topodata.cubic.fft$X2Dfft_r)

# Plotar---------------------------------------------
print(spplot(topodata.cubic.fft, zcol="X2Dfft_r",
             col.regions=gray(seq(0, 1, by=0.1)),
             aspect="iso", colorkey=FALSE,
             main="Cubic resampling\n 1 arc-second SRTM DEM (TOPODATA)",
             sub = "2D Fast Fourier Transform - real component"))

#----------------------------------------------------
# Cubicspline resampling
#----------------------------------------------------

# Configurar o GRASS GIS-----------------------------
initGRASS(gisBase = "/usr/lib/grass64",
          home = tempdir(), override=TRUE)

# Ler arquivo de dados-------------------------------
execGRASS("r.in.gdal", flags = "o",
          parameters = list(input = "topodata_cubicspline.tif",
                            output = "DEM"))

# Definir região de trabalho-------------------------
execGRASS("g.region", parameters = list(rast="DEM"))

# 2D FFT---------------------------------------------
execGRASS("i.fft", flags = "verbose",
          parameters = list(input_image = "DEM",
                            real_image = "2Dfft_r",
                            imaginary_image = "2Dfft_i"))

# Definir região de trabalho-------------------------
execGRASS("g.region", parameters = list(rast="2Dfft_r"))

# Ler arquivo de dados-------------------------------
topodata.cubicspline.fft <- readRAST6(vname = "2Dfft_r",
                            plugin = FALSE)

# Informações sobre o arquivo de dados---------------
summary(topodata.cubicspline.fft)
str(topodata.cubicspline.fft)
hist(topodata.cubicspline.fft$X2Dfft_r)

# Plotar---------------------------------------------
print(spplot(topodata.cubicspline.fft, zcol="X2Dfft_r",
             col.regions=gray(seq(0, 1, by=0.1)),
             aspect="iso", colorkey=FALSE,
             main="Cubicspline resampling\n 1 arc-second SRTM DEM (TOPODATA)",
             sub = "2D Fast Fourier Transform - real component"))

#----------------------------------------------------
# Lanczoz resampling
#----------------------------------------------------
system("gdalwarp")

# Configurar o GRASS GIS-----------------------------
initGRASS(gisBase = "/usr/lib/grass64",
          home = tempdir(), override=TRUE)

# Ler arquivo de dados-------------------------------
execGRASS("r.in.gdal", flags = "o",
          parameters = list(input = "topodata_lanczos.tif",
                            output = "DEM"))

# Definir região de trabalho-------------------------
execGRASS("g.region", parameters = list(rast="DEM"))

# 2D FFT---------------------------------------------
execGRASS("i.fft", flags = "verbose",
          parameters = list(input_image = "DEM",
                            real_image = "2Dfft_r",
                            imaginary_image = "2Dfft_i"))

# Definir região de trabalho-------------------------
execGRASS("g.region", parameters = list(rast="2Dfft_r"))

# Ler arquivo de dados-------------------------------
topodata.lanczos.fft <- readRAST6(vname = "2Dfft_r",
                                  plugin = FALSE)

# Informações sobre o arquivo de dados---------------
summary(topodata.lanczos.fft)
str(topodata.lanczos.fft)
hist(topodata.lanczos.fft$X2Dfft_r)

# Plotar---------------------------------------------
print(spplot(topodata.lanczos.fft, zcol="X2Dfft_r",
             col.regions=gray(seq(0, 1, by=0.1)),
             aspect="iso", colorkey=FALSE,
             main="Lanczos resampling\n 1 arc-second SRTM DEM (TOPODATA)",
             sub = "2D Fast Fourier Transform - real component"))

#####################################################
### 2D Fast Fourier Transform - SRTM
#####################################################
ls()

#----------------------------------------------------
# Nearest neighbour resampling
#----------------------------------------------------

# Configurar o GRASS GIS-----------------------------
initGRASS(gisBase = "/usr/lib/grass64",
          home = tempdir(), override=TRUE)

# Ler arquivo de dados-------------------------------
execGRASS("r.in.gdal", flags = "o",
          parameters = list(input = "srtm_near.tif",
                            output = "DEM"))

# Definir região de trabalho-------------------------
execGRASS("g.region", parameters = list(rast="DEM"))

# 2D FFT---------------------------------------------
execGRASS("i.fft", flags = "verbose",
          parameters = list(input_image = "DEM",
                            real_image = "2Dfft_r",
                            imaginary_image = "2Dfft_i"))

# Definir região de trabalho-------------------------
execGRASS("g.region", parameters = list(rast="2Dfft_r"))

# Ler arquivo de dados-------------------------------
srtm.near.fft <- readRAST6(vname = "2Dfft_r",
                           plugin = FALSE)

# Informações sobre o arquivo de dados---------------
summary(srtm.near.fft)
str(srtm.near.fft)
hist(srtm.near.fft$X2Dfft_r)

# Plotar---------------------------------------------
srtm.near.fft$X2Dfft_r01 <- abs(srtm.near.fft$X2Dfft_r)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
srtm.near.fft$X2Dfft_r01 <- range01(srtm.near.fft$X2Dfft_r01)
p.srtm.near <-
  spplot(srtm.near.fft, zcol="X2Dfft_r01",
         col.regions=gray(seq(0, 1, by=0.01)),
         aspect="iso", colorkey=TRUE,
         main = "3 arc-seconds SRTM DEM",
         sub = "Nearest neighbour resampling")
print(p.srtm.near)
rm(range01)


#----------------------------------------------------
# Bilinear resampling
#----------------------------------------------------

# Configurar o GRASS GIS-----------------------------
initGRASS(gisBase = "/usr/lib/grass64",
          home = tempdir(), override=TRUE)

# Ler arquivo de dados-------------------------------
execGRASS("r.in.gdal", flags = "o",
          parameters = list(input = "srtm_bilinear.tif",
                            output = "DEM"))

# Definir região de trabalho-------------------------
execGRASS("g.region", parameters = list(rast="DEM"))

# 2D FFT---------------------------------------------
execGRASS("i.fft", flags = "verbose",
          parameters = list(input_image = "DEM",
                            real_image = "2Dfft_r",
                            imaginary_image = "2Dfft_i"))

# Definir região de trabalho-------------------------
execGRASS("g.region", parameters = list(rast="2Dfft_r"))

# Ler arquivo de dados-------------------------------
srtm.bilinear.fft <- readRAST6(vname = "2Dfft_r",
                               plugin = FALSE)

# Informações sobre o arquivo de dados---------------
summary(srtm.bilinear.fft)
str(srtm.bilinear.fft)
hist(srtm.bilinear.fft$X2Dfft_r)

# Plotar---------------------------------------------
srtm.bilinear.fft$X2Dfft_r01 <- abs(srtm.bilinear.fft$X2Dfft_r)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
srtm.bilinear.fft$X2Dfft_r01 <- range01(srtm.bilinear.fft$X2Dfft_r01)
p.srtm.bi <-
  spplot(srtm.bilinear.fft, zcol="X2Dfft_r01",
         col.regions=gray(seq(0, 1, by=0.01)),
         aspect="iso", colorkey = TRUE,
         main = "3 arc-seconds SRTM DEM",
         sub = "Bilinear resampling")
print(p.srtm.bi)
rm(range01)

#----------------------------------------------------
# Cubic resampling
#----------------------------------------------------

# Configurar o GRASS GIS-----------------------------
initGRASS(gisBase = "/usr/lib/grass64",
          home = tempdir(), override=TRUE)

# Ler arquivo de dados-------------------------------
execGRASS("r.in.gdal", flags = "o",
          parameters = list(input = "srtm_cubic.tif",
                            output = "DEM"))

# Definir região de trabalho-------------------------
execGRASS("g.region", parameters = list(rast="DEM"))

# 2D FFT---------------------------------------------
execGRASS("i.fft", flags = "verbose",
          parameters = list(input_image = "DEM",
                            real_image = "2Dfft_r",
                            imaginary_image = "2Dfft_i"))

# Definir região de trabalho-------------------------
execGRASS("g.region", parameters = list(rast="2Dfft_r"))

# Ler arquivo de dados-------------------------------
srtm.cubic.fft <- readRAST6(vname = "2Dfft_r",
                            plugin = FALSE)

# Informações sobre o arquivo de dados---------------
summary(srtm.cubic.fft)
str(srtm.cubic.fft)
hist(srtm.cubic.fft$X2Dfft_r)

# Plotar---------------------------------------------
print(spplot(srtm.cubic.fft, zcol="X2Dfft_r",
             col.regions=gray(seq(0, 1, by=0.1)),
             aspect="iso", colorkey=FALSE,
             main="Cubic resampling\n 3 arc-seconds SRTM DEM",
             sub = "2D Fast Fourier Transform - real component"))

#----------------------------------------------------
# Cubicspline resampling
#----------------------------------------------------

# Configurar o GRASS GIS-----------------------------
initGRASS(gisBase = "/usr/lib/grass64",
          home = getwd(), override=TRUE)

# Ler arquivo de dados-------------------------------
execGRASS("r.in.gdal", flags = "o",
          parameters = list(input = "srtm_cubicspline.tif",
                            output = "DEM"))

# Definir região de trabalho-------------------------
execGRASS("g.region", parameters = list(rast="DEM"))

# 2D FFT---------------------------------------------
execGRASS("i.fft", flags = "verbose",
          parameters = list(input_image = "DEM",
                            real_image = "2Dfft_r",
                            imaginary_image = "2Dfft_i"))

# Definir região de trabalho-------------------------
execGRASS("g.region", parameters = list(rast="2Dfft_r"))

# Ler arquivo de dados-------------------------------
srtm.cubicspline.fft <- readRAST6(vname = "2Dfft_r",
                                  plugin = FALSE)

# Informações sobre o arquivo de dados---------------
summary(srtm.cubicspline.fft)
str(srtm.cubicspline.fft)
hist(srtm.cubicspline.fft$X2Dfft_r)

# Plotar---------------------------------------------
print(spplot(srtm.cubicspline.fft, zcol="X2Dfft_r",
             col.regions=gray(seq(0, 1, by=0.1)),
             aspect="iso", colorkey=FALSE,
             main="Cubicspline resampling\n 3 arc-seconds SRTM DEM",
             sub = "2D Fast Fourier Transform - real component"))

#----------------------------------------------------
# Lanczoz resampling
#----------------------------------------------------
system("gdalwarp")

# Configurar o GRASS GIS-----------------------------
initGRASS(gisBase = "/usr/lib/grass64",
          home = tempdir(), override=TRUE)

# Ler arquivo de dados-------------------------------
execGRASS("r.in.gdal", flags = "o",
          parameters = list(input = "srtm_lanczos.tif",
                            output = "DEM"))

# Definir região de trabalho-------------------------
execGRASS("g.region", parameters = list(rast="DEM"))

# 2D FFT---------------------------------------------
execGRASS("i.fft", flags = "verbose",
          parameters = list(input_image = "DEM",
                            real_image = "2Dfft_r",
                            imaginary_image = "2Dfft_i"))

# Definir região de trabalho-------------------------
execGRASS("g.region", parameters = list(rast="2Dfft_r"))

# Ler arquivo de dados-------------------------------
srtm.lanczos.fft <- readRAST6(vname = "2Dfft_r",
                              plugin = FALSE)

# Informações sobre o arquivo de dados---------------
summary(srtm.lanczos.fft)
str(srtm.lanczos.fft)
hist(srtm.lanczos.fft$X2Dfft_r)

# Plotar---------------------------------------------
print(spplot(srtm.lanczos.fft, zcol="X2Dfft_r",
             col.regions=gray(seq(0, 1, by=0.1)),
             aspect="iso", colorkey=FALSE,
             main="Lanczos resampling\n 3 arc-seconds SRTM DEM",
             sub = "2D Fast Fourier Transform - real component"))

#####################################################
### Produtos gráficos
#####################################################

#----------------------------------------------------
# Nearest neighbour (SRTM vs. TOPODATA)
#----------------------------------------------------
dev.off()
png(filename="srtm_topo_near.png")
print(p.srtm.near, split = c(1,1,2,1), more = T)
print(p.topo.near, split = c(2,1,2,1), more = F)
dev.off()

#----------------------------------------------------
# Bilinear (SRTM vs. TOPODATA)
#----------------------------------------------------
dev.off()
png(filename="srtm_topo_bili.png")
print(p.srtm.bi, split = c(1,1,2,1), more = T)
print(p.topo.bi, split = c(2,1,2,1), more = F)
dev.off()

#----------------------------------------------------
# Near - Bilinear (SRTM vs. TOPODATA)
#----------------------------------------------------

# Topodata-------------------------------------------
topodata.near.fft$near_bili <-
  topodata.near.fft$X2Dfft_r - topodata.bilinear.fft$X2Dfft_r
topodata.near.fft$near_bili <- abs(topodata.near.fft$near_bili)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
topodata.near.fft$near_bili <- range01(topodata.near.fft$near_bili)

# SRTM----------------------------------------------------------------
srtm.near.fft$near_bili <-
  srtm.near.fft$X2Dfft_r - srtm.bilinear.fft$X2Dfft_r
srtm.near.fft$near_bili <- abs(srtm.near.fft$near_bili)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
srtm.near.fft$near_bili <- range01(srtm.near.fft$near_bili)

# Plot----------------------------------------------------------------

# configuração da janela de plotagem
dev.off()
x11(width = 14, height = 14)

# trellis settings
tps <- list(fontsize=list(text=26))
trellis.par.set(tps)

# plot
p1 <-
  spplot(srtm.near.fft, zcol="near_bili",
         col.regions = rev(gray(seq(0, 1, by=0.01))),
         aspect="iso", colorkey=TRUE,
         main = "3 arc-seconds SRTM DEM",
         sub = "Nearest neighbour - bilinear resampling")
p2 <-
  spplot(topodata.near.fft, zcol="near_bili",
         col.regions = rev(gray(seq(0, 1, by=0.01))),
         aspect="iso", colorkey=TRUE,
         main = "1 arc-second TOPODATA DEM",
         sub = "Nearest neighbour - bilinear resampling")
print(p1, split = c(1,1,2,1), more = T)
print(p2, split = c(2,1,2,1), more = F)
savePlot(filename = "near_bili_poster.tiff", type = "tiff")
dev.off()
rm(p1, p2)

#----------------------------------------------------
# Cubic and Cubicspline and Lanczos - Bilinear
# (TOPODATA)
#----------------------------------------------------

# Cubic----------------------------------------------
topodata.cubic.fft$cubic_bili <-
  topodata.cubic.fft$X2Dfft_r - topodata.bilinear.fft$X2Dfft_r
topodata.cubic.fft$cubic_bili <- abs(topodata.cubic.fft$cubic_bili)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
topodata.cubic.fft$cubic_bili <- range01(topodata.cubic.fft$cubic_bili)
rm(range01)

# Cubicspline----------------------------------------
topodata.cubicspline.fft$cubicsp_bili <-
  topodata.cubicspline.fft$X2Dfft_r - topodata.bilinear.fft$X2Dfft_r
topodata.cubicspline.fft$cubicsp_bili <- abs(topodata.cubicspline.fft$cubicsp_bili)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
topodata.cubicspline.fft$cubicsp_bili <- range01(topodata.cubicspline.fft$cubicsp_bili)
rm(range01)

# Lanczos--------------------------------------------
topodata.lanczos.fft$lanczos_bili <-
  topodata.lanczos.fft$X2Dfft_r - topodata.bilinear.fft$X2Dfft_r
topodata.lanczos.fft$lanczos_bili <- abs(topodata.lanczos.fft$lanczos_bili)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
topodata.lanczos.fft$lanczos_bili <- range01(topodata.lanczos.fft$lanczos_bili)
rm(range01)

# Plot-----------------------------------------------

# configuração da janela de plotagem
dev.off()
x11(width = 14, height = 14)

# trellis settings
tps <- list(fontsize=list(text=26))
trellis.par.set(tps)

# plot
p1 <-
  spplot(topodata.cubic.fft, zcol="cubic_bili",
         col.regions = rev(gray(seq(0, 1, by=0.01))),
         aspect="iso", colorkey=TRUE,
         main = "Cubic - bilinear",
         sub ="Topodata")
p2 <-
  spplot(topodata.cubicspline.fft, zcol="cubicsp_bili",
         col.regions = rev(gray(seq(0, 1, by=0.01))),
         aspect="iso", colorkey=TRUE,
         main = "Cubicspline - bilinear",
         sub ="Topodata")
p3 <-
  spplot(topodata.lanczos.fft, zcol="lanczos_bili",
         col.regions = rev(gray(seq(0, 1, by=0.01))),
         aspect="iso", colorkey=TRUE,
         main = "Lanczos - bilinear",
         sub ="Topodata")

# Salvar arquivo-------------------------------------
print(p1, split = c(1,1,3,1), more = T)
print(p2, split = c(2,1,3,1), more = T)
print(p3, split = c(3,1,3,1), more = T)
savePlot(filename = "all_bili_poster.tiff", type = "tiff")
dev.off()
rm(p1, p2, p3)

#####################################################
### Salvar objetos 2D-FFT
#####################################################
ls()
save(srtm.bilinear.fft, srtm.cubic.fft,
     srtm.cubicspline.fft, srtm.lanczos.fft,
     srtm.near.fft, topodata.bilinear.fft,
     topodata.cubic.fft, topodata.cubicspline.fft,
     topodata.lanczos.fft, topodata.near.fft,
     file="2dfft.RData")

#####################################################
### Configuração do ambiente de trabalho do RSAGA
#####################################################
rsaga.env()
sagaenv <- rsaga.env(workspace=".", cmd="saga_cmd",
                     version="2.0.7",
                     path="/usr/bin",
                     modules="/usr/lib/saga")
sagaenv

#####################################################
### MDE em 3D
#####################################################
ls()
#--------------------------------------------------
# Importar DEM e salvar no formato SAGA grid (sgrd)
#--------------------------------------------------
rsaga.import.gdal("topodata_bilinear.tif", "dem", env=sagaenv)
rsaga.get.libraries()
rsaga.get.modules(lib="grid_visualisation")
rsaga.get.usage(lib="garden_3d_viewer",
                module = "3D Shapes Viewer",
                env = sagaenv)

rsaga.geoprocessor(lib="grid_visualisation",
                   module="Create 3D Image",
                   param=list(DEM = "dem.sgrd",
                              IMAGE = "dem.sgrd",
                              ZEXAGG = 3,
                              Z_ROTATE = 20,
                              X_ROTATE = 0,
                              X_ROTATE_LEVEL = "Zero",
                              RGB = "3d_teste"),
                   env=sagaenv)

rsaga.get.libraries()
rsaga.get.modules(lib="io_grid_image")
rsaga.get.usage(lib="io_grid_image", 0)
rsaga.geoprocessor(lib="io_grid_image",
                   module = 0,
                   param=list(GRID = "3d_teste.sgrd",
                              FILE = "/home/alessandro/rdata/dnos_sm/terrain/2dfft/teste.png"),
                   env=sagaenv)


# 2D-FFT filtering ###################################################

# Read file and define region ========================================
# TOPODATA warped using cubic resampling
execGRASS("r.in.gdal", flags = "o",
          parameters = list(input = "/home/alessandro/PROJETOS/DNOS-SM/data/GIS/terrain-attributes/topodata-dem.tif",
                            output = "DEM"))
execGRASS("g.region", rast = "DEM")

# 2D-FFT =============================================================
execGRASS("i.fft", input_image = "DEM", real_image = "DEM_r",
          imaginary_image = "DEM_i", flags = "overwrite")

# Ampitude image -----------------------------------------------------
execGRASS("r.mapcalculator", flags = "overwrite",
          parameters = list(amap = "DEM_r", bmap = "DEM_i",
                        formula = "sqrt((DEM_r*DEM_r)+(DEM_i*DEM_i))",
                        outfile = "DEM_a"))
image(readRAST6("DEM_a"), col = gray((0:32)/32), axes = TRUE)
locator()

# Create mask ========================================================
par(mfrow = c(2,2))

# low pass filter ----------------------------------------------------
execGRASS("r.circle", flags = c("b", "overwrite"), output = "lp_mask",
          coordinate = "235282,6717229", max = 5000)
image(readRAST6("lp_mask"), col = "gray")

# high pass mask -----------------------------------------------------
execGRASS("r.mapcalculator", flags = "overwrite",
          parameters = list(outfile = "hp_mask",
                            amap = "lp_mask",
                            formula = "if(isnull(lp_mask),1,null())"))
image(readRAST6("hp_mask"), col = "gray")

# Filter =============================================================
execGRASS("g.rename", parameters = list(rast = "lp_mask,MASK"))
execGRASS("i.ifft", flags = "overwrite",
          parameters = list(real_image = "DEM_r",
                            imaginary_image = "DEM_i",
                            output_image = "DEM_fft"))
execGRASS("g.remove",
          parameters = list(rast = "MASK"))
image(readRAST6("DEM_fft"), col = topo.colors(64))

# difference ---------------------------------------------------------
execGRASS("r.mapcalculator", flags = c("overwrite"),
          parameters = list(outfile = "fft_diff",
                            amap = "DEM", bmap = "DEM_fft",
                            formula = "DEM-DEM_fft"))
image(readRAST6("fft_diff"), col = gray((0:16)/16), axes = TRUE)
contour(readRAST6("fft_diff"), add = TRUE, drawlabels = FALSE,
        nlevels = 5)

# Fill depressions ===================================================
execGRASS("r.fill.dir", flags = c("overwrite"),
          parameters = list(input = "DEM_fft",
                            elevation = "DEM_fft_fill",
                            direction = "DEM_dir"))
execGRASS("r.mapcalculator", flags = c("overwrite"),
          parameters = list(outfile = "fill_diff",
                            amap = "DEM_fft", bmap = "DEM_fft_fill",
                            formula = "DEM_fft-DEM_fft_fill"))
image(readRAST6("fill_diff"), col = heat.colors(12))


