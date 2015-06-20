# DESCRIPTION ##################################################################
# Analysis of the soil taxa data available.
# The following area-class soil maps are used:
# - Azolin, M. A. D. & Mutti, L. S. M. Solos da bacia hidrográfica do 
#   Vacacaí-Mirim. Porto Alegre: DNOS-UFSM, p. 20, 1988.
# - Miguel, P.; Dalmolin, R. S. D.; Pedron, F. A.; Samuel-Rosa, A.; Medeiros, P.
#   S. C.; Moura-Bueno, J. M. & Balbinot, A. Soil and land use dynamics in
#   Plateau Border areas of Rio Grande do Sul. Revista Brasileira de 
#   Agrociência. v. 17, p. 347-455, 2011.
# SETTINGS #####################################################################
rm(list = ls())
gc()
# Load Packages
require(sp)
require(raster)
require(rgdal)
require(spgrass6)
require(rpart)
library(plyr)
library(lattice)
library(latticeExtra)
library(grid)
# Load data
load("sm-dnos-general.RData")
load("sm-dnos-soil-class.RData")
load("/home/alessandro/PROJECTS/pedometrics/pedometrics/cooking/color_ramps.RData")
ls()

# GRASS GIS
initGRASS(gisBase = "/usr/lib/grass64/", gisDbase = GRASSgisDbase,
          override = TRUE, location = "dnos-sm-rs",
          mapset = "predictions", pid = Sys.getpid())
system("g.region rast=dnos.raster")

# LOAD AND PREPROCESS DATA #####################################################

# Soil Map 1:100k ==============================================================
SOIL_100 <- shapefile(paste(soil.dir, "soil-map-1988.shp", sep = ""))
SOIL_100 <- spTransform(SOIL_100, wgs1984utm22s)
SOIL_100
plot(SOIL_100)

# soil class IDs
SOIL_100.id <- data.frame(unique(SOIL_100$code), unique(SOIL_100$symbol))
colnames(SOIL_100.id) <- c("code", "MU")
SOIL_100.id <- SOIL_100.id[order(SOIL_100.id$code),]
SOIL_100.id

# export to GRASS GIS (mapset = predictions)

writeVECT6(SOIL_100, "SOIL_100", v.in.ogr_flags = c("overwrite"))
system("v.info -c SOIL_100")

# shapes to grid
system("v.to.rast --o in=SOIL_100 out=SOIL_100 use=attr column=code")
system("r.info SOIL_100"); system("r.category SOIL_100")
system("d.mon x0");system("d.rast.leg SOIL_100")

# add category names
SOIL_100.id
system("r.category SOIL_100 rules=- <<EOF
       1:C1
       2:HC1
       3:PB1-Co
       4:PBa1
       5:PEa-Rd
       6:PEd2
       7:PL1
       8:PVd2
       9:Rd1
       10:Rd2
       11:Re-C-Co
       12:Re3
       13:Re4
       14:TBa-Rd
       EOF")
system("r.category SOIL_100")
system("d.mon x0");system("d.rast.leg SOIL_100")

# Soil Map 1:25k ===============================================================
SOIL_25 <- shapefile(paste(soil.dir, "soil-map-miguel", sep = ""))
SOIL_25 <- spTransform(SOIL_25, wgs1984utm22s)
SOIL_25
plot(SOIL_25)

# soil class IDs
SOIL_25.id <- data.frame(unique(SOIL_25$code), unique(SOIL_25$symbol))
colnames(SOIL_25.id) <- c("code", "MU")
SOIL_25.id <- SOIL_25.id[order(SOIL_25.id$code),]
SOIL_25.id

# export to GRASS GIS (mapset = predictions)
initGRASS(gisBase = "/usr/lib/grass64/", gisDbase = GRASSgisDbase, override = TRUE,
          location = "dnos-sm-rs", mapset = "predictions", pid = Sys.getpid())
system("g.region rast=dnos.raster")
writeVECT6(SOIL_25, "SOIL_25", v.in.ogr_flags = "overwrite")
system("v.info -c SOIL_25")

# shapes to grid
system("v.to.rast --o in=SOIL_25 out=SOIL_25 use=attr column=code")
system("r.info SOIL_25"); system("r.category SOIL_25")
system("d.mon x0");system("d.rast.leg SOIL_25")

# add category names
SOIL_25.id
system("r.category SOIL_25 rules=- <<EOF
       1:PBAC
       2:PV
       3:C-R
       4:RY
       5:RL
       6:RL-RR
       7:RR
       8:SX
       EOF")
system("r.category SOIL_25")
system("d.mon x0");system("d.rast.leg SOIL_25")

# CREATE ENVIRONMENTAL COVARIATES ##############################################
initGRASS(gisBase = "/usr/lib/grass64/", gisDbase = GRASSgisDbase,
          override = TRUE, location = "dnos-sm-rs",
          mapset = "predictions", pid = Sys.getpid())
system("g.region rast=dnos.raster")
system("g.remove MASK")

# SOIL_100 =====================================================================
system("r.category SOIL_100")
system("d.mon x0")
system("d.rast.leg SOIL_100")
system("d.vect map=BASIN_10 fcolor=none")
# There are only five classes occurring in the study area. Thus, the total
# number of unique classes is p = 5 - 1 = 4. Other classes can be produced
# by the combination of unique classes.
# Rd1, Re4, Re-C-Co, TBa-Rd, C1

# SOIL_100a (Rd1 - Solo Litólico distrófico e eutrófico)
# Classes: 1 = Rd1
#          0 = Other
system("r.category SOIL_100")
system("r.reclass --o in=SOIL_100 out=SOIL_100a <<EOF
       9 = 1 Rd1
       * = 0 Other
       end
       EOF")
system("r.mapcalc SOIL_100a=SOIL_100a")
system("r.category SOIL_100a rules=- <<EOF
       1:Rd1
       0:Other
       EOF")
system("r.category SOIL_100a")
system("r.info SOIL_100a")
system("d.mon x0")
system("d.erase")
system("d.rast.leg SOIL_100a")

# SOIL_100b (Re4 - Solo Litólico eutrófico e distrófico relevo montanhoso)
# Classes: 1 = Re4
#          0 = Other
system("r.category SOIL_100")
system("r.reclass --o in=SOIL_100 out=SOIL_100b <<EOF
       13 = 1 Re4
       * = 0 Other
       end
       EOF")
system("r.mapcalc SOIL_100b=SOIL_100b")
system("r.category SOIL_100b rules=- <<EOF
       1:Re4
       0:Other
       EOF")
system("r.category SOIL_100b")
system("r.info SOIL_100b")
system("d.mon x0")
system("d.erase")
system("d.rast.leg SOIL_100b")

# SOIL_100c (Re-C-Co - Solo Litólico eutrófico relevo forte ondulado + 
#                      Cambissolo eutrófico + Colúvios)
# Classes: 1 = Re-C-Co
#          0 = Other
system("r.category SOIL_100")
system("r.reclass --o in=SOIL_100 out=SOIL_100c <<EOF
       11 = 1 Re-C-Co
       * = 0 Other
       end
       EOF")
system("r.mapcalc SOIL_100c=SOIL_100c")
system("r.category SOIL_100c rules=- <<EOF
       1:Re-C-Co
       0:Other
       EOF")
system("r.category SOIL_100c")
system("r.info SOIL_100c")
system("d.mon x0")
system("d.erase")
system("d.rast.leg SOIL_100c")

# SOIL_100d (TBa-Rd - Terra Bruna Estruturada álica + Solo Litólico)
# Classes: 1 = TBa-Rd
#          0 = Other
system("r.category SOIL_100")
system("r.reclass --o in=SOIL_100 out=SOIL_100d <<EOF
       14 = 1 TBa-Rd
       * = 0 Other
       end
       EOF")
system("r.mapcalc SOIL_100d=SOIL_100d")
system("r.category SOIL_100d rules=- <<EOF
       1:TBa-Rd
       0:Other
       EOF")
system("r.category SOIL_100d")
system("r.info SOIL_100d")
system("d.mon x0")
system("d.erase -f")
system("d.rast.leg SOIL_100d")

# SOIL_100e (Rd1 - Solo Litólico distrófico e eutrófico + 
#            Re4 - Solo Litólico eutrófico e distrófico relevo montanhoso)
# Classes: 1 = Rd1 + Re4
#          0 = Other
system("r.category SOIL_100")
system("r.reclass --o in=SOIL_100 out=SOIL_100e <<EOF
       9 13 = 1 Rd1 + Re4
       * = 0 Other
       end
       EOF")
system("r.mapcalc SOIL_100e=SOIL_100e")
system("r.category SOIL_100e rules=- <<EOF
       1:Rd1 + Re4
       0:Other
       EOF")
system("r.category SOIL_100e")
system("r.info SOIL_100e")
system("d.mon x0")
system("d.erase -f")
system("d.rast.leg SOIL_100e")

# SOIL_100f (TBa-Rd - Terra Bruna Estruturada álica + Solo Litólico +
#            C1 - Cambissolo eutrófico)
# Classes: 1 = TBa-Rd + C1
#          0 = Other
system("r.category SOIL_100")
system("r.reclass --o in=SOIL_100 out=SOIL_100f <<EOF
       14 1 = 1 TBa-Rd + C1
       * = 0 Other
       end
       EOF")
system("r.mapcalc SOIL_100f=SOIL_100f")
system("r.category SOIL_100f rules=- <<EOF
       1:TBa-Rd + C1
       0:Other
       EOF")
system("r.category SOIL_100f")
system("r.info SOIL_100f")
system("d.mon x0")
system("d.erase -f")
system("d.rast.leg SOIL_100f")

# SOIL_25 =====================================================================
system("r.category SOIL_25")
system("d.mon x0")
system("d.rast.leg SOIL_25")
system("d.vect map=BASIN_10 fcolor=none")
# There are eight classes in the study area. Thus, the number of unique classes
# that can be obtained is p = 8 - 1 = 7. Other classes can be obtained by the
# combination of unique classes.

# SOIL_25a (PBAC - Argissolo Bruno-Acinzentado)
# Classes: 1 = PBAC
#          0 = Other
system("r.category SOIL_25")
system("r.reclass --o in=SOIL_25 out=SOIL_25a <<EOF
       1 = 1 PBAC
       * = 0 Other
       end
       EOF")
system("r.mapcalc SOIL_25a=SOIL_25a")
system("r.category SOIL_25a rules=- <<EOF
       1:PBAC
       0:Other
       EOF")
system("r.category SOIL_25a")
system("r.info SOIL_25a")
system("d.mon x0")
system("d.erase")
system("d.rast.leg SOIL_25a")

# SOIL_25b (PV - Argissolo Vermelho)
# Classes: 1 = PV
#          0 = Other
system("r.category SOIL_25")
system("r.reclass --o in=SOIL_25 out=SOIL_25b <<EOF
       2 = 1 PV
       * = 0 Other
       end
       EOF")
system("r.mapcalc SOIL_25b=SOIL_25b")
system("r.category SOIL_25b rules=- <<EOF
       1:PV
       0:Other
       EOF")
system("r.category SOIL_25b")
system("r.info SOIL_25b")
system("d.mon x0")
system("d.erase")
system("d.rast.leg SOIL_25b")

# SOIL_25c (C-R - Cambissolo - Neossolo)
# Classes: 1 = C-R
#          0 = Other
system("r.category SOIL_25")
system("r.reclass --o in=SOIL_25 out=SOIL_25c <<EOF
       3 = 1 C-R
       * = 0 Other
       end
       EOF")
system("r.mapcalc SOIL_25c=SOIL_25c")
system("r.category SOIL_25c rules=- <<EOF
       1:C-R
       0:Other
       EOF")
system("r.category SOIL_25c")
system("r.info SOIL_25c")
system("d.mon x0")
system("d.erase")
system("d.rast.leg SOIL_25c")

# SOIL_25d (RL - Neossolo Litólico)
# Classes: 1 = RL
#          0 = Other
system("r.category SOIL_25")
system("r.reclass --o in=SOIL_25 out=SOIL_25d <<EOF
       5 = 1 RL
       * = 0 Other
       end
       EOF")
system("r.mapcalc SOIL_25d=SOIL_25d")
system("r.category SOIL_25d rules=- <<EOF
       1:RL
       0:Other
       EOF")
system("r.category SOIL_25d")
system("r.info SOIL_25d")
system("d.mon x0")
system("d.erase")
system("d.rast.leg SOIL_25d")

# SOIL_25h (PBAC + PV + SX - Argissolo Bruno-Acinzentado + Argissolo Vermelho + Planossolo Háplico)
# Classes: 1 = PBAC + PV + SX
#          0 = Other
system("r.category SOIL_25")
system("r.reclass --o in=SOIL_25 out=SOIL_25h <<EOF
       1 2 8 = 1 PBAC + PV + SX
       * = 0 Other
       end
       EOF")
system("r.mapcalc SOIL_25h=SOIL_25h")
system("r.category SOIL_25h rules=- <<EOF
       1:PBAC + PV + SX
       0:Other
       EOF")
system("r.category SOIL_25h")
system("r.info SOIL_25h")
system("d.mon x0")
system("d.erase")
system("d.rast.leg SOIL_25h")

# SOIL_25i (RL + RL-RR + RR - Neossolo Litólico + Neossolo Regolítico)
# Classes: 1 = RL + RL-RR + RR
#          0 = Other
system("r.category SOIL_25")
system("r.reclass --o in=SOIL_25 out=SOIL_25i <<EOF
       5 6 7 = 1 RL + RL-RR + RR
       * = 0 Other
       end
       EOF")
system("r.mapcalc SOIL_25i=SOIL_25i")
system("r.category SOIL_25i rules=- <<EOF
       1:RL + RL-RR + RR
       0:Other
       EOF")
system("r.category SOIL_25i")
system("r.info SOIL_25i")
system("d.mon x0")
system("d.erase")
system("d.rast.leg SOIL_25i")

# SOIL_25j (PV - Argissolo Vermelho +
#           RL - Neossolo Litólico +
#           RL-RR - Neossolo Litólico - Neossolo Regolítico +
#           C-R - Cambissolo + Neossolo)
# Classes: 1 = PV + RL + RL-RR + C-R
#          0 = Other
system("r.category SOIL_25")
system("r.reclass --o in=SOIL_25 out=SOIL_25j <<EOF
       2 3 5 6 = 1 PV + RL + RL-RR + C-R
       * = 0 Other
       end
       EOF")
system("r.mapcalc SOIL_25j=SOIL_25j")
system("r.category SOIL_25j rules=- <<EOF
       1:PV + RL + RL-RR + C-R
       0:Other
       EOF")
system("r.category SOIL_25j")
system("r.info SOIL_25j")
system("d.mon x0")
system("d.erase")
system("d.rast.leg SOIL_25j")

# SAVE FIGURES WITH AREA-CLASS SOIL MAPS #######################################
# Use as boundary the basin estimated with the digital elevation models derived
# from the contour lines plus the 30 m buffer (buffer_BASIN_10).
system("r.mask -o buffer_BASIN_10")
boundary <- readVECT6("buffer_BASIN_10")@bbox
# SOIL_100
map <- readRAST6("SOIL_100")
map@bbox <- boundary
map@data[, 1] <- as.factor(map@data[, 1])
levels(map@data[, 1])
system("r.category SOIL_100")
map@data[, 1] <- revalue(map@data[, 1], c("1" = "C1",
                                          "9" = "Rd1",
                                          "11" = "Re-C-Co",
                                          "13" = "Re4",
                                          "14" = "TBa-Rd"))
color <- c(sibcs$colors[2], sibcs$colors[41], sibcs$colors[21],
            sibcs$colors[39], sibcs$colors[43])
p1 <- spplot(map, main = "", col.regions = color, asp = 1)
names(p1$legend) <- "inside"
p1$legend$inside$x <- 0.73
p1$legend$inside$y <- 0.85
p1$legend$inside$args$key$space <- "left"
dev.off()
pdf(file = paste(soil.dir, "SOIL_100.pdf", sep = ""),
    width = 6.3/cm(1), height = (6.3/cm(1))*1.14)
trellis.par.set(fontsize = list(text = 7, points = 5),
                plot.line = list(lwd = 0.001),
                axis.line = list(lwd = 0.01),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
plot(p1)
dev.off()
rm(map, p1, color)
gc()
# SOIL_25
map <- readRAST6("SOIL_25")
map@bbox <- boundary
map@data[, 1] <- as.factor(map@data[, 1])
levels(map@data[, 1])
system("r.category SOIL_25")
map@data[, 1] <- revalue(map@data[, 1], c("1" = "PBAC",
                                          "2" = "PV",
                                          "3" = "C-R",
                                          "4" = "RY",
                                          "5" = "RL",
                                          "6" = "RL-RR",
                                          "7" = "RR",
                                          "8" = "SX"))
color <- c(sibcs$colors[15], sibcs$colors[17], sibcs$colors[2],
            sibcs$colors[40], sibcs$colors[39], sibcs$colors[13],
            sibcs$colors[41], sibcs$colors[52])
p1 <- spplot(map, main = "", col.regions = color, asp = 1)

names(p1$legend) <- "inside"
p1$legend$inside$args$key$space <- "left"
p1$legend$inside$x <- 0.77
p1$legend$inside$y <- 0.775
p1$legend$inside$args$key$col <- rev(p1$legend$inside$args$key$col)
p1$legend$inside$args$key$labels$labels <- rev(p1$legend$inside$args$key$labels$labels)
p1$legend
dev.off()
pdf(file = paste(soil.dir, "SOIL_25.pdf", sep = ""),
    width = 6.3/cm(1), height = (6.3/cm(1))*1.14)
trellis.par.set(fontsize = list(text = 7, points = 5),
                plot.line = list(lwd = 0.001),
                axis.line = list(lwd = 0.01),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
plot(p1)
dev.off()

rm(map, boundary, colors)
gc()
system("g.remove MASK")

# SAVE DATA ####################################################################
ls()
save(SOIL_100, SOIL_100.id, SOIL_25, SOIL_25.id, soil.dir,
     file = "sm-dnos-soil-class.RData")










#####################################################
### Configuração inicial das variáveis
#####################################################
ls()
str(dnos.taxon)
summary(dnos.taxon$Taxon)
dnos.taxon$CX <- dnos.taxon$Taxon == "CX"
dnos.taxon$PA <- dnos.taxon$Taxon == "PA"
dnos.taxon$PBAC <- dnos.taxon$Taxon == "PBAC"
dnos.taxon$ALISOLS <- dnos.taxon$Taxon == "PBAC"
dnos.taxon$PV <- dnos.taxon$Taxon == "PV"
dnos.taxon$PVA <- dnos.taxon$Taxon == "PVA"
dnos.taxon$RF <- dnos.taxon$Taxon == "RF"
dnos.taxon$RL <- dnos.taxon$Taxon == "RL"
dnos.taxon$RQ <- dnos.taxon$Taxon == "RQ"
dnos.taxon$RR <- dnos.taxon$Taxon == "RR"
dnos.taxon$SX <- dnos.taxon$Taxon == "SX"
dnos.taxon$ARG <-
  is.element(dnos.taxon$Taxon,
             c("PV", "PVA", "PA", "PBAC", "SX"))
dnos.taxon$NEO <-
  is.element(dnos.taxon$Taxon,
             c("RL", "RR", "RQ", "RF", "CX"))
dnos.taxon$ACRISOLS <-
  is.element(dnos.taxon$Taxon,
             c("PV", "PVA", "PA"))
dnos.taxon$WATER.SOIL <-
  is.element(dnos.taxon$Taxon,
             c("SX", "RF"))
dnos.taxon$LOW.WATER.SOIL <-
  is.element(dnos.taxon$Taxon,
             c("RQ", "RL"))
dnos.taxon$YOUNG.SOIL <-
  is.element(dnos.taxon$Taxon,
             c("RR", "CX"))
str(dnos.taxon)
ls()


# Conversão em objeto espacial
str(dnos.taxon)
coordinates(dnos.taxon) <- c("X", "Y")
str(dnos.taxon)
class(dnos.taxon)


# Salvar objetos no arquivo de dados
save(dnos.coords, dnos.psd, dnos.raster,
     dnos.rocha, dnos.taxon, dnos.uso,
     file="/home/alessandro/rdata/dnos_sm/dnos_dados.RData")
ls()
rm(list=ls())


###################################################
### Análise exploratória
###################################################


# Estatísticas descritivas-------------------------
str(dnos.taxon)
su.1 <- (summary(dnos.taxon$Taxon))
su.2 <- round(summary(
  dnos.taxon$Taxon)/length(dnos.taxon$Taxon)*100, 2)
(desc.taxon <- data.frame(
  rbind(su.1, su.2), row.names=c("n", "%")))
rm(su.1, su.2)
ls()


# Barplot dos taxa--------------------------------
pdf(file="/home/alessandro/rdata/dnos_sm/figures/cal_soil_taxa.pdf",
    width=8, height=8)
plot(dnos.taxon$Taxon,
     main="Count of inferred soil taxa",
     ylim=c(0, 175), col="lightgray",
     xlab="Taxon", ylab="Count",
     sub="Braziliam System of Soil Classification, 2006")
abline(h=seq(0, 150, by=25), lty=2, col="lightgray")
dev.off()


# Posição dos taxa inferidos-----------------------
pdf(file="/home/alessandro/rdata/dnos_sm/figures/cal_soil_taxa_pos.pdf",
    width=8, height=8, pointsize=20)
print(spplot(
  dnos.taxon, zcol="Taxon",
  main="Location of inferred soil taxa",
  sub="Braziliam System of Soil Classification, 2006",
  key.space="right",
  xlab="E (m)", ylab="N (m)",
  xlim=bbox(dnos.raster)[1,],
  ylim=bbox(dnos.raster)[2,],
  aspect="iso", pch=c(1:10),
  scales = list(draw = TRUE),
  col.regions=bpy.colors(10),
  panel=function(x,y, ...) {
    panel.xyplot(x, y, ...);
    panel.grid(h=-1,v=-1,
               col="lightgray",
               lty=2)}))
dev.off()
ls()

# Posição de cada táxon inferido (all at once)-------------------
length(dnos.taxon@data)
(taxon <- colnames(dnos.taxon@data))
for(i in 3:length(dnos.taxon@data)){
  pdf(file=paste("/home/alessandro/rdata/dnos_sm/figures/cal_soil_taxa_pos_",
                 taxon[i],".pdf", sep=""), width=8, height=8)
  print(spplot(
    dnos.taxon, zcol=taxon[i],
    main=paste(taxon[i]),
    xlab="E (m)", ylab="N (m)",
    xlim=bbox(dnos.raster)[1,],
    ylim=bbox(dnos.raster)[2,],
    aspect="iso", auto.key=FALSE,
    scales = list(draw = TRUE),
    col.regions=c("black", "red"),
    panel=function(x,y, ...) {
      panel.xyplot(x, y, ...);
      panel.grid(h=-1,v=-1,
                 col="lightgray",
                 lty=2)}))
  dev.off()
}
rm(taxon, i)
ls()

###################################################
### Análise variográfica
###################################################
ls()


# Variograma omnidirecional------------------------


# tmp
str(dnos.taxon)
colnames(dnos.taxon@data)
tmp.var.ind <- variogram(YOUNG.SOIL ~ 1,
                         loc=dnos.taxon,
                         cutoff=1000,
                         width=1000/10,
                         cloud=F)
print(plot(tmp.var.ind,
           main="tmp",
           pch=20, col="blue",
           xlim=c(0, max(tmp.var.ind$dist)*1.1),
           ylim=c(0, max(tmp.var.ind$gamma)*1.1),
           plot.numbers=T))
ls()
rm(tmp.var.ind)


# Conclusões: 10 classes iniciais resultaram em 5 classes finais
#---
# RL: ACE positiva moderada a forte, com duas estruturas
# RQ: ACE negativa fraca. Poucos pontos. Juntar com RL.
### LOW.WATER.SOIL: ACE positiva moderada a forte, com duas estruturas
#---
# CX: ACE negativa
# RR: ACE positiva fraca
### YOUNG.SOIL: ACE positiva muito fraca no longo alcance
#---
# SX: ACE negativa, depois positiva, fracas.
# RF: ACE negativa fraca. Juntar com SX.
### WATER.SOIL: ACE negativa fraca no curto alcance; ACE potiviva fraca no longo alcance. Usar WATER.SOIL
#---
### PBAC: ACE negativa no curto alcance; ACE positiva no longo alcance. Será chamado ALISOLS
#---
# PV: ACE positiva
# PVA: poucas observações, localizadas na borda. Melhor juntar com PV e PA
# PA: ausência de ACE
### ACRISOLS: ACE positiva fraca, mas melhor. Usar ACRISOLS
#---
# Conclusões: 10 classes iniciais resultaram em 5 classes finais. São elas: ALISOLS, ACRISOLS, WATER.SOIL, LOW.WATER.SOIL, YOUNG.SOIL.


# ALISOLS
# Wave model; bin de 100 m
summary(dnos.taxon$ALISOLS)
ALISOLS.var.ind <- variogram(ALISOLS ~ 1,
                          loc=dnos.taxon)
ALISOLS.var.ind.plot <-
  plot(ALISOLS.var.ind,
       main="Indicator variogram",
       sub="ALISOLS",
       pch=20, col="blue",
       xlim=c(0, max(ALISOLS.var.ind$dist)*1.1),
       ylim=c(0, max(ALISOLS.var.ind$gamma)*1.1),
       plot.numbers=T)
print(ALISOLS.var.ind.plot)



# ACRISOLS
summary(dnos.taxon$ACRISOLS)
ACRISOLS.var.ind <- variogram(ACRISOLS ~ 1,
                        loc=dnos.taxon),
                        cutoff=960,
                        width=120)
ACRISOLS.var.ind.plot <-
  plot(ACRISOLS.var.ind,
       main="Indicator variogram",
       sub="ACRISOLS",
       pch=20, col="blue",
       xlim=c(0, max(ACRISOLS.var.ind$dist)*1.1),
       ylim=c(0, max(ACRISOLS.var.ind$gamma)*1.1),
       plot.numbers=T)
print(ACRISOLS.var.ind.plot)


# WATER.SOIL (bin de 150 m)
summary(dnos.taxon$WATER.SOIL)
WATER.SOIL.var.ind <- variogram(WATER.SOIL ~ 1,
                        loc=dnos.taxon),
                        cutoff=1200,
                        width=120)
WATER.SOIL.var.ind.plot <-
  plot(WATER.SOIL.var.ind,
       main="Indicator variogram",
       sub="WATER.SOIL",
       pch=20, col="blue",
       xlim=c(0, max(WATER.SOIL.var.ind$dist)*1.1),
       ylim=c(0, max(WATER.SOIL.var.ind$gamma)*1.1),
       plot.numbers=T)
print(WATER.SOIL.var.ind.plot)


# LOW.WATER.SOIL
summary(dnos.taxon$LOW.WATER.SOIL)
LOW.WATER.SOIL.var.ind <- variogram(LOW.WATER.SOIL ~ 1,
                                loc=dnos.taxon),
                                cutoff=1200,
                                width=120)
LOW.WATER.SOIL.var.ind.plot <-
  plot(LOW.WATER.SOIL.var.ind,
       main="Indicator variogram",
       sub="LOW.WATER.SOIL",
       pch=20, col="blue",
       xlim=c(0, max(LOW.WATER.SOIL.var.ind$dist)*1.1),
       ylim=c(0, max(LOW.WATER.SOIL.var.ind$gamma)*1.1),
       plot.numbers=T)
print(LOW.WATER.SOIL.var.ind.plot)


# YOUNG.SOIL (bin de 100 m)
summary(dnos.taxon$YOUNG.SOIL)
YOUNG.SOIL.var.ind <- variogram(YOUNG.SOIL ~ 1,
                                    loc=dnos.taxon),
                                    cutoff=2000,
                                    width=200)
YOUNG.SOIL.var.ind.plot <-
  plot(YOUNG.SOIL.var.ind,
       main="Indicator variogram",
       sub="YOUNG.SOIL",
       pch=20, col="blue",
       xlim=c(0, max(YOUNG.SOIL.var.ind$dist)*1.1),
       ylim=c(0, max(YOUNG.SOIL.var.ind$gamma)*1.1),
       plot.numbers=T)
print(YOUNG.SOIL.var.ind.plot)

ls()

# Superfície variográfica-----------------------------------------------


# YOUNG.SOIL (bin de 100 m)
YOUNG.SOIL.vario.map <-
  plot(variogram(YOUNG.SOIL ~ 1, dnos.taxon, map=TRUE,
                 cutoff=1000, width=100),
       main="Variogram map",
       sub="YOUNG.SOIL",
       col.regions=bpy.colors(64))
print(YOUNG.SOIL.vario.map)


# WATER.SOIL (bin de 150 m)
WATER.SOIL.vario.map <-
  plot(variogram(WATER.SOIL ~ 1, dnos.taxon, map=TRUE,
                 cutoff=1000, width=150),
       main="Variogram map",
       sub="WATER.SOIL",
       col.regions=bpy.colors(64))
print(WATER.SOIL.vario.map)


# LOW.WATER.SOIL (bin de 100 m)
LOW.WATER.SOIL.vario.map <-
  plot(variogram(LOW.WATER.SOIL ~ 1, dnos.taxon, map=TRUE,
                 cutoff=1000, width=100),
       main="Variogram map",
       sub="LOW.WATER.SOIL",
       col.regions=bpy.colors(64))
print(LOW.WATER.SOIL.vario.map)


# ACRISOLS (bin de 150 m)
ACRISOLS.vario.map <-
  plot(variogram(ACRISOLS ~ 1, dnos.taxon, map=TRUE,
                 cutoff=1000, width=150),
       main="Variogram map",
       sub="ACRISOLS",
       col.regions=bpy.colors(64))
print(ACRISOLS.vario.map)


# ALISOLS (bin de 100 m)
ALISOLS.vario.map <-
  plot(variogram(ALISOLS ~ 1, dnos.taxon, map=TRUE,
                 cutoff=1000, width=100),
       main="Variogram map",
       sub="ALISOLS",
       col.regions=bpy.colors(64))
print(ALISOLS.vario.map)


# Conclusão: Não há padrão anisotrópico bem definido para qualquer das variáveis que justifique a modelagem anisotrópica.

#--------------------------------------------------------
# Ajuste do modelo
#--------------------------------------------------------
print(show.vgms())
ls()

# ALISOLS
# Ajuste visual
print(ALISOLS.var.ind.plot)
# sill = 0.23; nugget = 0.05; psill = 0.18; range = 120; mod: Cir
(ALISOLS.var.ind.mod <- vgm(psill=0.11, model="Mat",
                            range=80, nugget=0.05,
                            kappa=5))

print(plot(ALISOLS.var.ind,
           main="Indicator variogram",
           sub="ALISOLS",
           pch=20, col="blue",
           xlim=c(0, max(ALISOLS.var.ind$dist)*1.1),
           ylim=c(0, 0.3),
           plot.numbers=T, model=ALISOLS.var.ind.mod))

# TESTES COM BESSEL
plot(c(0.08, ALISOLS.var.ind$gamma)~c(0, ALISOLS.var.ind$dist),
     xlim=c(0, max(ALISOLS.var.ind$dist)*1.1),
     ylim=c(0, 0.3))
lines(besselJ(c(0.08, ALISOLS.var.ind$gamma), nu=0.75)~
        c(0,ALISOLS.var.ind$dist))
tmp <- Variogram(c(0.08, ALISOLS.var.ind$gamma), model="bessel",
                 param=c(NA, 0.17, 0.1, 0.09, 0.75))
lines(tmp~c(0,ALISOLS.var.ind$dist))

# TESTE COM geoR
data <- data.frame(dnos.taxon)
data <- as.geodata(data, coords.col=c(2, 3), data.col="ALISOLS")
str(data)
tmp <- variog(data, max.dist=2000)
plot(tmp)
eyefit(tmp)
rm(data, tmp)

# Ajuste automático
(ALISOLS.var.ind.mod.fit <-
   fit.variogram(object=ALISOLS.var.ind,
                 fit.sills=c(FALSE, TRUE),
                 model=ALISOLS.var.ind.mod))
print(plot(ALISOLS.var.ind,
           main="Indicator variogram",
           sub="ALISOLS",
           pch=20, col="blue",
           xlim=c(0, max(ALISOLS.var.ind$dist)*1.1),
           ylim=c(0, max(ALISOLS.var.ind$gamma)*1.1),
           plot.numbers=T, model=ALISOLS.var.ind.mod.fit))


# ACRISOLS
# Ajuste visual
print(ACRISOLS.var.ind.plot)
# sill = 0.035; nugget = 0.00; psill = 0.035; range = 400; mod: Gau
(ACRISOLS.var.ind.mod <- vgm(psill=0.035,
                        model="Gau",
                        range=400/sqrt(3),
                        nugget=0))
print(plot(ACRISOLS.var.ind,
           main="Indicator variogram",
           sub="ACRISOLS",
           pch=20, col="blue",
           xlim=c(0, max(ACRISOLS.var.ind$dist)*1.1),
           ylim=c(0, max(ACRISOLS.var.ind$gamma)*1.1),
           plot.numbers=T, model=ACRISOLS.var.ind.mod))
# Ajuste automático
(ACRISOLS.var.ind.mod.fit <-
   fit.variogram(object=ACRISOLS.var.ind,
                 model=ACRISOLS.var.ind.mod,
                 warn.if.neg=TRUE))
print(plot(ACRISOLS.var.ind,
           main="Indicator variogram",
           sub="ACRISOLS",
           pch=20, col="blue",
           xlim=c(0, max(ACRISOLS.var.ind$dist)*1.1),
           ylim=c(0, max(ACRISOLS.var.ind$gamma)*1.1),
           plot.numbers=T, model=ACRISOLS.var.ind.mod.fit))
vmf <- ACRISOLS.var.ind.mod.fit
1-vmf$psill[1]/sum(vmf$psill)
rm(vmf)




# WATER.SOIL???????????????????????????????????????????????????????????????????
# Ajuste visual
print(WATER.SOIL.var.ind.plot)
# sill = 0.08; nugget = 0; psill = 0.14; range = 300; mod: Pen
(WATER.SOIL.var.ind.mod <- vgm(psill=-0.12,
                            model="Pen",
                            range=300,
                            nugget=0.25))
print(plot(WATER.SOIL.var.ind,
           main="Indicator variogram",
           sub="WATER.SOIL",
           pch=20, col="blue",
           xlim=c(0, max(WATER.SOIL.var.ind$dist)*1.1),
           ylim=c(0, max(WATER.SOIL.var.ind$gamma)*1.1),
           plot.numbers=T, model=WATER.SOIL.var.ind.mod))
# Ajuste automático
(WATER.SOIL.var.ind.mod.fit <-
   fit.variogram(object=WATER.SOIL.var.ind,
                 model=WATER.SOIL.var.ind.mod))
print(plot(WATER.SOIL.var.ind,
           main="Indicator variogram",
           sub="WATER.SOIL",
           pch=20, col="blue",
           xlim=c(0, max(WATER.SOIL.var.ind$dist)*1.1),
           ylim=c(0, max(WATER.SOIL.var.ind$gamma)*1.1),
           plot.numbers=T, model=WATER.SOIL.var.ind.mod.fit))


# LOW.WATER.SOIL
# Ajuste visual
print(LOW.WATER.SOIL.var.ind.plot)
# sill = 0.18; nugget = 0.02; psill = 0.16; range = 300; mod: Pen
(LOW.WATER.SOIL.var.ind.mod <- vgm(psill=0.16,
                               model="Pen",
                               range=300,
                               nugget=0.02))
print(plot(LOW.WATER.SOIL.var.ind,
           main="Indicator variogram",
           sub="LOW.WATER.SOIL",
           pch=20, col="blue",
           xlim=c(0, max(LOW.WATER.SOIL.var.ind$dist)*1.1),
           ylim=c(0, max(LOW.WATER.SOIL.var.ind$gamma)*1.1),
           plot.numbers=T, model=LOW.WATER.SOIL.var.ind.mod))
# Ajuste automático
(LOW.WATER.SOIL.var.ind.mod.fit <-
   fit.variogram(object=LOW.WATER.SOIL.var.ind,
                 model=LOW.WATER.SOIL.var.ind.mod))
print(plot(LOW.WATER.SOIL.var.ind,
           main="Indicator variogram",
           sub="LOW.WATER.SOIL",
           pch=20, col="blue",
           xlim=c(0, max(LOW.WATER.SOIL.var.ind$dist)*1.1),
           ylim=c(0, max(LOW.WATER.SOIL.var.ind$gamma)*1.1),
           plot.numbers=T, model=LOW.WATER.SOIL.var.ind.mod.fit))
LOW.WATER.SOIL.var.ind.mod.fit
vmf <- LOW.WATER.SOIL.var.ind.mod.fit
1-vmf$psill[1]/sum(vmf$psill)
rm(vmf)

# YOUNG.SOIL????????????????????????????????????????????????????????????????????
# Ajuste visual
print(YOUNG.SOIL.var.ind.plot)
# sill = 0.18; nugget = 0.02; psill = 0.16; range = 300; mod: Pen
(YOUNG.SOIL.var.ind.mod <- vgm(psill=0.16,
                                   model="Pen",
                                   range=300,
                                   nugget=0.02))
print(plot(YOUNG.SOIL.var.ind,
           main="Indicator variogram",
           sub="YOUNG.SOIL",
           pch=20, col="blue",
           xlim=c(0, max(LOW.WATER.SOIL.var.ind$dist)*1.1),
           ylim=c(0, max(LOW.WATER.SOIL.var.ind$gamma)*1.1),
           plot.numbers=T, model=LOW.WATER.SOIL.var.ind.mod))
# Ajuste automático
(LOW.WATER.SOIL.var.ind.mod.fit <-
   fit.variogram(object=LOW.WATER.SOIL.var.ind,
                 model=LOW.WATER.SOIL.var.ind.mod))
print(plot(LOW.WATER.SOIL.var.ind,
           main="Indicator variogram",
           sub="LOW.WATER.SOIL",
           pch=20, col="blue",
           xlim=c(0, max(LOW.WATER.SOIL.var.ind$dist)*1.1),
           ylim=c(0, max(LOW.WATER.SOIL.var.ind$gamma)*1.1),
           plot.numbers=T, model=LOW.WATER.SOIL.var.ind.mod.fit))
