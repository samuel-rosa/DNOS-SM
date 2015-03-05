# DESCRIPTION ##################################################################
# Exploratory data analysis of point soil data and environmental covariates.
# SETTINGS #####################################################################
rm(list = ls())
gc()
# Load packages
require(fitdistrplus)
require(histogram)
require(stringr)
require(sp)
require(rgdal)
require(maptools)
require(raster)
require(spgrass6)
require(car)
require(TeachingDemos)
require(lattice)
# Load data
load("sm-dnos-general.RData")
load("sm-dnos-phd-chap1.RData")
ls()
# GRASS GIS
initGRASS(gisBase = "/usr/lib/grass64/", gisDbase = GRASSgisDbase,
          location = "dnos-sm-rs", mapset = "predictions",
          pid = Sys.getpid(), override = TRUE)
system("g.region rast=dnos.raster")
system("g.remove MASK")
system("r.mask dnos.raster")
system("r.mask buffer_BASIN_10")

# EXPLORATORY DATA ANALYSIS - ENVIRONMENTAL COVARIATES #########################
# Take a look at the maps and the distribution of points across them (geographic
# space).
# Take a look at the distribution of predictor variables across observations -
# (feature space).
# NOTE: There are no requirements about the probability distribution of
#       predictor variables. They are deterministic variables and can have any 
#       type of distribution. However, they can be transformed to obtain a 
#       better correlation with the dependent variable.

# Land use 1980 (LU1980 - land1) ===============================================

# geographic space
map <- str_replace_all(land1, "[ ]", "")
map <- c(unlist(str_split(map, "[+]")))
lapply(map, function(X){
  map = map
  system(paste("d.mon x", match(X, map), sep = ""))
  system(paste("d.rast.leg ", X, sep = ""))
  system("d.vect calibration")
});rm(map)
# LU1980a: this map is fine.
# LU1980b: this map is fine.
# LU1980c: this map is fine.

# feature space
covars <- str_replace_all(land1, "[ ]", "")
covars <- c(unlist(str_split(covars, "[+]")))
covars <- cal_data@data[,match(covars, colnames(cal_data@data))]
summary(covars)
rm(covars)
# LU1980b and LU1980c are very similar. They may be aliased.

# Land use 2009 (LU2009 - land2) =============================================

# geographic space
map <- str_replace_all(land2, "[ ]", "") 
map <- c(unlist(str_split(map, "[+]")))
lapply(map, function(X){
  map = map
  system(paste("d.mon x", match(X, map), sep = ""))
  system(paste("d.rast.leg ", X, sep = ""))
  system("d.vect calibration")
});rm(map)
# LU2009a: this map is fine.
# LU2009b: nor very fine.
# LU2009c: this map is fine.
# LU2009d: not very fine.
# LU2009e: fine.
# LU2009f: this map is fine.
# LUdiff: this map is fine

# feature space
covars <- str_replace_all(land2, "[ ]", "")
covars <- c(unlist(str_split(covars, "[+]")))
covars <- cal_data@data[,match(covars, colnames(cal_data@data))]
summary(covars)
rm(covars)
# LU2009e and LU2009f might be aliased.

# Area-class soil map 1:100,000 (SOIL_100 - soil1) ============================

# geographic space
map <- str_replace_all(soil1, "[ ]", "") 
map <- c(unlist(str_split(map, "[+]")))
lapply(map, function(X){
  map = map
  system(paste("d.mon x", match(X, map), sep = ""))
  system(paste("d.rast.leg ", X, sep = ""))
  system("d.vect calibration")
});rm(map)
# SOIL_100a: not sure if this is fine. Only 11 observations.
# SOIL_100b: fine.
# SOIL_100c: fine.
# SOIL_100d: fine.
# SOIL_100e: fine.
# SOIL_100f: fine.

# feature space
covars <- str_replace_all(soil1, "[ ]", "")
covars <- c(unlist(str_split(covars, "[+]")))
covars <- cal_data@data[,match(covars, colnames(cal_data@data))]
summary(covars)
rm(covars)
# SOIL_100a: not sure if this is fine. Only 11 observations. Remove from models.

# Area-class soil map 1:25,000 (SOIL_25 - soil2) ===============================

# geographic space
map <- str_replace_all(soil2, "[ ]", "") 
map <- c(unlist(str_split(map, "[+]")))
lapply(map[1:5], function(X){
  map = map
  system(paste("d.mon x", match(X, map)-1, sep = ""))
  system(paste("d.rast.leg ", X, sep = ""))
  system("d.vect calibration")
})
lapply(map[6:10], function(X){
  map = map
  system(paste("d.mon x", match(X, map)-6, sep = ""))
  system(paste("d.rast.leg ", X, sep = ""))
  system("d.vect calibration")
})
rm(map)
# SOIL_25a: fine
# SOIL_25b: this map is fine.
# SOIL_25c: this map is fine.
# SOIL_25d: this map is fine.
# SOIL_25e: not sure if this is fine. Remove from models.
# SOIL_25f: not sure if this is fine. Remove from models.
# SOIL_25g: this is not fine. Remove from models.
# SOIL_25h: fine
# SOIL_25i: fine
# SOIL_25j: fine

# feature space
covars <- str_replace_all(soil2, "[ ]", "")
covars <- c(unlist(str_split(covars, "[+]")))
covars <- cal_data@data[,match(covars, colnames(cal_data@data))]
summary(covars)
rm(covars)

system("d.mon x0")
system("d.rast.leg SOIL_25f")
system("d.vect calibration")

# Geological map 1:50000 (GEO_50 - geo1) =======================================

# geographic space
map <- str_replace_all(geo1, "[ ]", "") 
map <- c(unlist(str_split(map, "[+]")))
lapply(map, function(X){
  map = map
  system(paste("d.mon x", match(X, map), sep = ""))
  system(paste("d.rast.leg ", X, sep = ""))
  system("d.vect calibration")
});rm(map)
# GEO_25a: this map is not realy fine, but it is worth tring to use it!
# GEO_25b: this map is fine!
# GEO_25c: not sure
# GEO_25d: this map is fine!
# GEO_25e: this map is fine!

# feature space
covars <- str_replace_all(geo1, "[ ]", "")
covars <- c(unlist(str_split(covars, "[+]")))
covars <- cal_data@data[,match(covars, colnames(cal_data@data))]
summary(covars)
rm(covars)
# all fine

# Geological map 1:25000 (GEO_25 - geo2) =======================================

# geographic space
map <- str_replace_all(geo2, "[ ]", "")
map <- c(unlist(str_split(map, "[+]")))
lapply(map, function(X){
  map = map
  system(paste("d.mon x", match(X, map), sep = ""))
  system(paste("d.rast.leg ", X, sep = ""))
  system("d.vect calibration")
});rm(map)
# GEO_25a: this map is not realy fine, but it is worth tring to use it!
# GEO_25b: this map is fine!
# GEO_25c: this map is fine!
# GEO_25d: fine
# GEO_25e: fine
# GEO_25f: fine

# feature space
covars <- str_replace_all(geo2, "[ ]", "")
covars <- c(unlist(str_split(covars, "[+]")))
covars <- cal_data@data[,match(covars, colnames(cal_data@data))]
summary(covars)
rm(covars)

# Digital elevarion model with 90 meters resolution (ELEV_90 - dem1) ===========

# geographic space
map <- str_replace_all(dem1, "[ ]", "")
map <- c(unlist(str_split(map, "[+]")))

# ELEV_90: uniform distribution
system("d.mon x0")
system(paste("d.rast.leg ", map[1], sep = ""))
system("d.vect calibration")
system("d.mon x1")
system(paste("d.histogram ", map[1], sep = ""))

# log-normal distribution
lapply(map[2:8], function(X){
  map = map
  system(paste("d.mon x", match(X, map)-2, sep = ""))
  system(paste("d.rast.leg ", X, sep = ""))
  system("d.vect calibration")
})
lapply(map[2:8], function(X){
  map = map
  system(paste("d.mon x", match(X, map)-2, sep = ""))
  system(paste("d.histogram ", X, sep = ""))
})

# TPI_90: normal distribution
lapply(map[9:15], function(X){
  map = map
  system(paste("d.mon x", match(X, map)-9, sep = ""))
  system(paste("d.rast.leg ", X, sep = ""))
  system("d.vect calibration")
})
lapply(map[9:15], function(X){
  map = map
  system(paste("d.mon x", match(X, map)-9, sep = ""))
  system(paste("d.histogram ", X, sep = ""))
})

# NOR_90: uniform distribution
lapply(map[16:22], function(X){
  system(paste("r.colors map=", X, " color=grey", sep = ""))
})
lapply(map[16:22], function(X){
  map = map
  system(paste("d.mon x", match(X, map)-16, sep = ""))
  system(paste("d.rast.leg ", X, sep = ""))
  system("d.vect calibration")
})
lapply(map[16:22], function(X){
  map = map
  system(paste("d.mon x", match(X, map)-16, sep = ""))
  system(paste("d.histogram ", X, sep = ""))
})

# TWI_90: log-normal distribution
lapply(map[23:29], function(X){
  map = map
  system(paste("d.mon x", match(X, map)-23, sep = ""))
  system(paste("d.rast.leg ", X, sep = ""))
  system("d.vect calibration")
})
lapply(map[23:29], function(X){
  map = map
  system(paste("d.mon x", match(X, map)-23, sep = ""))
  system(paste("d.histogram ", X, sep = ""))
})

# SPI_90: normal distribution
lapply(map[30:36], function(X){
  map = map
  system(paste("d.mon x", match(X, map)-30, sep = ""))
  system(paste("d.rast.leg ", X, sep = ""))
  system("d.vect calibration")
})
lapply(map[30:36], function(X){
  map = map
  system(paste("d.mon x", match(X, map)-30, sep = ""))
  system(paste("d.histogram ", X, sep = ""))
})
rm(map)

# feature space
covars <- str_replace_all(dem1, "[ ]", "")
covars <- c(unlist(str_split(covars, "[+]")))
covars <- cal_data@data[,match(covars, colnames(cal_data@data))]
length(covars)
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "gray", ...)
}
pairs(covars[1:12], panel = panel.smooth, cex = 0.5, pch = 20,
      diag.panel = panel.hist, cex.labels = 1.5, font.labels = 2)
pairs(covars[13:24], panel = panel.smooth, cex = 0.5, pch = 20,
      diag.panel = panel.hist, cex.labels = 1.5, font.labels = 2)
pairs(covars[25:36], panel = panel.smooth, cex = 0.5, pch = 20,
      diag.panel = panel.hist, cex.labels = 1.5, font.labels = 2)

par(mfrow = c(1, 3))
hist(unlist(covars["SPI_10_63"]))
hist(sqrt(unlist(covars["SPI_10_63"])))
hist(log(unlist(covars["SPI_10_63"])))
dev.off()
# ELEV_90:                 uniform
# SLP_90_3 - SLP_90_127:   slight positive skew -> square root
# SLP_90_255:              normal
# TPI_90:                  normal
# NOR_90_3 - NOR_90_127:   normal/uniform
# NOR_90_255:              slight positive skew -> square root
# TWI_90:                  strong positive skew -> logarithm
# SPI_90:                  very slight positive skew -> do nothing
rm(covars)

# CORRELATION WITH SOIL ATTRIBUTES
covars <- str_replace_all(dem1, "[ ]", "")
covars <- c(unlist(str_split(covars, "[+]")))
covars <- cal_data@data[,match(covars, colnames(cal_data@data))]
length(covars)

# clay content
var <- cal_data$clay
var <- bcPower(var, clay.lambda)
par(mfrow = c(3,4))
for(i in 1:length(covars[1:12])){
  nam = names(covars[1:12][i])
  plot(var ~ covars[,nam], main = nam)
}
for(i in 1:length(covars[13:24])){
  nam = names(covars[13:24][i])
  plot(var ~ covars[,nam], main = nam)
}
for(i in 1:length(covars[25:36])){
  nam = names(covars[25:36][i])
  plot(var ~ covars[,nam], main = nam)
}
# watch the behavior of TWI

# carbon content
var <- cal_data$carbon
var <- bcPower(var, carbon.lambda)
par(mfrow = c(3,4))
for(i in 1:length(covars[1:12])){
  nam = names(covars[1:12][i])
  plot(var ~ covars[,nam], main = nam)
}
for(i in 1:length(covars[13:24])){
  nam = names(covars[13:24][i])
  plot(var ~ covars[,nam], main = nam)
}
for(i in 1:length(covars[25:36])){
  nam = names(covars[25:36][i])
  plot(var ~ covars[,nam], main = nam)
}
# watch the behavior of TWI

# ecec content
var <- cal_data$ecec
var <- bcPower(var, ecec.lambda)
par(mfrow = c(3,4))
for(i in 1:length(covars[1:12])){
  nam = names(covars[1:12][i])
  plot(var ~ covars[,nam], main = nam)
}
for(i in 1:length(covars[13:24])){
  nam = names(covars[13:24][i])
  plot(var ~ covars[,nam], main = nam)
}
for(i in 1:length(covars[25:36])){
  nam = names(covars[25:36][i])
  plot(var ~ covars[,nam], main = nam)
}
# watch the behavior of TWI
rm(covars)

# histogram of predictor variables and calibration data
par(mfrow = c(1, 2))
y <- "SAVI_5a"
map <- readRAST6(paste(y))@data[, 1]
cal <- cal_data@data[, y]
hist(map, main = y, freq = FALSE)
rug(cal)
legend(mean(cal), mean(hist(map, plot = FALSE)$density),
       legend = round(cor(cal, cal_data$carbon), 2))

# Digital elevarion model derived from contour lines (ELEV_10 - dem2) ==========

# geographic space
map <- str_replace_all(dem2, "[ ]", "")
map <- c(unlist(str_split(map, "[+]")))

# ELEV_10: uniform distribution
system("d.mon x0");system(paste("d.rast.leg ", map[1], sep = ""));system("d.vect calibration")
system("d.mon x1");system(paste("d.histogram ", map[1], sep = ""))

# SLP_10: log-normal distribution
lapply(map[2:8], function(X){
  map = map
  system(paste("d.mon x", match(X, map)-2, sep = ""))
  system(paste("d.rast.leg ", X, sep = ""))
  system("d.vect calibration")
})
lapply(map[2:8], function(X){
  map = map
  system(paste("d.mon x", match(X, map)-2, sep = ""))
  system(paste("d.histogram ", X, sep = ""))
})

# TPI_10: normal distribution
lapply(map[9:15], function(X){
  map = map
  system(paste("d.mon x", match(X, map)-9, sep = ""))
  system(paste("d.rast.leg ", X, sep = ""))
  system("d.vect calibration")
})
lapply(map[9:15], function(X){
  map = map
  system(paste("d.mon x", match(X, map)-9, sep = ""))
  system(paste("d.histogram ", X, sep = ""))
})

# NOR_10: uniform distribution
lapply(map[16:22], function(X){system(paste("r.colors map=", X, " color=grey", sep = ""))})
lapply(map[16:22], function(X){
  map = map
  system(paste("d.mon x", match(X, map)-16, sep = ""))
  system(paste("d.rast.leg ", X, sep = ""))
  system("d.vect calibration")
})
lapply(map[16:22], function(X){
  map = map
  system(paste("d.mon x", match(X, map)-16, sep = ""))
  system(paste("d.histogram ", X, sep = ""))
})

# TWI_10: log-normal distribution
lapply(map[23:29], function(X){
  map = map
  system(paste("d.mon x", match(X, map)-23, sep = ""))
  system(paste("d.rast.leg ", X, sep = ""))
  system("d.vect calibration")
})
lapply(map[23:29], function(X){
  map = map
  system(paste("d.mon x", match(X, map)-23, sep = ""))
  system(paste("d.histogram ", X, sep = ""))
})

# SPI_10: quasi-normal distribution
lapply(map[30:36], function(X){
  map = map
  system(paste("d.mon x", match(X, map)-30, sep = ""))
  system(paste("d.rast.leg ", X, sep = ""))
  system("d.vect calibration")
})
lapply(map[30:36], function(X){
  map = map
  system(paste("d.mon x", match(X, map)-30, sep = ""))
  system(paste("d.histogram ", X, sep = ""))
})
rm(map)

# feature space
covars <- str_replace_all(dem2, "[ ]", "")
covars <- c(unlist(str_split(covars, "[+]")))
covars <- calibration@data[,match(covars, colnames(calibration@data))]
length(covars)
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "gray", ...)
}
pairs(covars[1:12], panel = panel.smooth, cex = 0.5, pch = 20,
      diag.panel = panel.hist, cex.labels = 2, font.labels = 2)
pairs(covars[13:24], panel = panel.smooth, cex = 0.5, pch = 20,
      diag.panel = panel.hist, cex.labels = 2, font.labels = 2)
pairs(covars[25:36], panel = panel.smooth, cex = 0.5, pch = 20,
      diag.panel = panel.hist, cex.labels = 2, font.labels = 2)
# ELEV_10:                 uniform
# SLP_10_3 - SLP_10_127:   slight positive skew -> square root
# SLP_10_255:              normal
# TPI_10:                  normal
# NOR_10_3 - NOR_10_127:   normal/uniform
# NOR_10_255:              slight positive skew -> square root
# TWI_10:                  strong positive skew -> logarithm
# SPI_10:                  very slight positive skew -> do nothing
rm(covars)

# CORRELATION WITH SOIL ATTRIBUTES
covars <- str_replace_all(dem2, "[ ]", "")
covars <- c(unlist(str_split(covars, "[+]")))
covars <- cal_data@data[,match(covars, colnames(cal_data@data))]
length(covars)

# clay content
var <- cal_data$clay
var <- bcPower(var, bc_lambda$clay)
par(mfrow = c(3,4))
for(i in 1:length(covars[1:12])){
  nam = names(covars[1:12][i])
  plot(var ~ covars[,nam], main = nam)
}
for(i in 1:length(covars[13:24])){
  nam = names(covars[13:24][i])
  plot(var ~ covars[,nam], main = nam)
}
for(i in 1:length(covars[25:36])){
  nam = names(covars[25:36][i])
  plot(var ~ covars[,nam], main = nam)
}

# carbon content
var <- cal_data$carbon
var <- bcPower(var, bc_lambda$carbon)
par(mfrow = c(3,4))
for(i in 1:length(covars[1:12])){
  nam = names(covars[1:12][i])
  plot(var ~ covars[,nam], main = nam)
}
for(i in 1:length(covars[13:24])){
  nam = names(covars[13:24][i])
  plot(var ~ covars[,nam], main = nam)
}
for(i in 1:length(covars[25:36])){
  nam = names(covars[25:36][i])
  plot(var ~ covars[,nam], main = nam)
}
a <- vector()
nam <- vector()
for(i in 1:length(covars[9:15])){
  nam[i] <- names(covars[9:15][i])
  a[i] <- cor(var, covars[, nam[i]])
}
cbind(nam, a)

# ecec content
var <- cal_data$ecec
var <- bcPower(var, ecec.lambda)
par(mfrow = c(3,4))
for(i in 1:length(covars[1:12])){
  nam = names(covars[1:12][i])
  plot(var ~ covars[,nam], main = nam)
}
for(i in 1:length(covars[13:24])){
  nam = names(covars[13:24][i])
  plot(var ~ covars[,nam], main = nam)
}
for(i in 1:length(covars[25:36])){
  nam = names(covars[25:36][i])
  plot(var ~ covars[,nam], main = nam)
}
# watch the behavior of TWI
rm(covars, var)

# Orbital image from Landsat 5 TM (sat1) =======================================

# geographic space
map <- str_replace_all(sat1, "[ ]", "")
map <- c(unlist(str_split(map, "[+]")))
lapply(map[1:6], function(X){
  map = map
  system(paste("d.mon x", match(X, map), sep = ""))
  system(paste("d.rast.leg ", X, sep = ""))
  system("d.vect calibration")
})
lapply(map[1:6], function(X){
  map = map
  system(paste("d.mon x", match(X, map), sep = ""))
  system(paste("d.histogram ", X, sep = ""))
})
# BLUE_30:  log-normal
# GREEN_30: log-normal
# RED_30:   log-normal
# NIR_30a:  negative skew
# NIR_30b:  log-normal
# MIR_30:   negative skew
lapply(map[7:8], function(X){
  map = map
  system(paste("d.mon x", match(X, map)-7, sep = ""))
  system(paste("d.rast.leg ", X, sep = ""))
  system("d.vect calibration")
})
lapply(map[7:8], function(X){
  map = map
  system(paste("d.mon x", match(X, map)-7, sep = ""))
  system(paste("d.histogram ", X, sep = ""))
})
# NDVI_30: negative skew
# SAVI_30: negative skew
rm(map)

# feature space
covars <- str_replace_all(sat1, "[ ]", "")
covars <- c(unlist(str_split(covars, "[+]")))
covars <- calibration@data[,match(covars, colnames(calibration@data))]
length(covars)
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "gray", ...)
}
pairs(covars, panel = panel.smooth, cex = 0.5, pch = 20,
      diag.panel = panel.hist, cex.labels = 1.5, font.labels = 2)

plot.hist <- function(data){
  par(mfrow = c(1, 3))
  hist(unlist(covars[data]))
  hist(exp(unlist(covars[data])))
  hist((unlist(covars[data]))^5)
}
plot.hist("NIR_30b")
# BLUE_30:  very slight positive skew -> do nothing
# GREEN_30: slight positive skew      -> square root
# RED_30:   strong positive skew      -> logarithm
# NIR_30a:  strong negative skew      -> x^5
# NIR_30b:  very small range          -> drop
# MIR_30:   very small range          -> drop
# NDVI_30:  negative skew             -> x^5
# SAVI_30:  negative skew             -> x^5
rm(covars)

# CORRELATION WITH SOIL ATTRIBUTES
covars <- str_replace_all(sat1, "[ ]", "")
covars <- c(unlist(str_split(covars, "[+]")))
covars <- cal_data@data[,match(covars, colnames(cal_data@data))]

# clay content
var <- cal_data$clay
var <- bcPower(var, clay.lambda)
par(mfrow = c(3,4))
for(i in 1:length(covars)){
  nam = names(covars[i])
  plot(var ~ covars[,nam], main = nam)
}

# carbon content
var <- cal_data$carbon
var <- bcPower(var, carbon.lambda)
par(mfrow = c(3,4))
for(i in 1:length(covars)){
  nam = names(covars[i])
  plot(var ~ covars[,nam], main = nam)
}

# ecec content
var <- cal_data$ecec
var <- bcPower(var, ecec.lambda)
par(mfrow = c(3,4))
for(i in 1:length(covars)){
  nam = names(covars[i])
  plot(var ~ covars[,nam], main = nam)
}
rm(covars, var)

# Orbital image from RapidEye (sat2) ===========================================

# geographic space
map <- str_replace_all(sat2, "[ ]", "")
map <- c(unlist(str_split(map, "[+]")))
lapply(map[1:5], function(X){
  map = map
  system(paste("d.mon x", match(X, map), sep = ""))
  system(paste("d.rast.leg ", X, sep = ""))
  system("d.vect calibration")
})
lapply(map[1:5], function(X){
  map = map
  system(paste("d.mon x", match(X, map), sep = ""))
  system(paste("d.histogram ", X, sep = ""))
})
# BLUE_5:  log-normal
# GREEN_5: log-normal
# RED_5:   log-normal
# EDGE_5:  log-normal
# NIR_5:   log-normal
lapply(map[6:9], function(X){
  map = map
  system(paste("d.mon x", match(X, map)-6, sep = ""))
  system(paste("d.rast.leg ", X, sep = ""))
  system("d.vect calibration")
})
lapply(map[6:9], function(X){
  map = map
  system(paste("d.mon x", match(X, map)-6, sep = ""))
  system(paste("d.histogram ", X, sep = ""))
})
rm(map)
# NDVI_5a: negative skew
# NDVI_5b: normal
# SAVI_5a: negative skew
# SAVI_5b: normal

# feature space
covars <- str_replace_all(sat2, "[ ]", "")
covars <- c(unlist(str_split(covars, "[+]")))
covars <- calibration@data[,match(covars, colnames(calibration@data))]
length(covars)
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "gray", ...)
}
pairs(covars, panel = panel.smooth, cex = 0.5, pch = 20,
      diag.panel = panel.hist, cex.labels = 1.5, font.labels = 2)

plot.hist <- function(data){
  par(mfrow = c(1, 3))
  hist(unlist(covars[data]))
  hist((unlist(covars[data]))^2)
  hist(exp(unlist(covars[data])))
}
plot.hist("SAVI_5b")
# BLUE_5:   slight positive skew  -> square root
# GREEN_5:  normal
# RED_5:    positive skew         -> logarithm
# EDGE_5:   normal
# NIR_5:    normal
# NDVI_5a:  negative skew         -> x^3
# NDVI_5b:  normal
# SAVI_5a:  slight negative skew  -> x^3
# SAVI_5b:  slight negative skew  -> x^2
rm(covars)

# CORRELATION WITH SOIL ATTRIBUTES
covars <- str_replace_all(sat2, "[ ]", "")
covars <- c(unlist(str_split(covars, "[+]")))
covars <- cal_data@data[,match(covars, colnames(cal_data@data))]
length(covars)

# clay content
var <- cal_data$clay
var <- bcPower(var, clay.lambda)
par(mfrow = c(3,3))
for(i in 1:length(covars)){
  nam = names(covars[i])
  plot(var ~ covars[,nam], main = nam)
}

# carbon content
var <- cal_data$carbon
var <- bcPower(var, carbon.lambda)
par(mfrow = c(3,3))
for(i in 1:length(covars)){
  nam = names(covars[i])
  plot(var ~ covars[,nam], main = nam)
}

# ecec content
var <- cal_data$ecec
var <- bcPower(var, ecec.lambda)
par(mfrow = c(3,3))
for(i in 1:length(covars)){
  nam = names(covars[i])
  plot(var ~ covars[,nam], main = nam)
}
rm(covars, var)





















# SAVE DATA ####################################################################
ls()


# End!
