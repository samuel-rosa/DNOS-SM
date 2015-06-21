# Federal Rural University of Rio de Janeiro
# Graduate School in Agronomy - Soil Science
# Project: Contribution to the Construction of Models to Predict Soil Properties
# Author: Alessandro Samuel-Rosa
# Wageningen, August of 2014
#
# Supervisory committee:
# Lúcia Helena Cunha dos Anjos (UFRRJ)
# Gustavo de Mattos Vasques (Embrapa)
# Gerard B M Heuvelink (ISRIC)
#
# Working script: spatio-temporal point pattern analysis
#
# e-mail: alessandrosamuel@yahoo.com.br
# homepage: soil-scientist.net
#
# SETTINGS #####################################################################
options(device = x11, stringsAsFactors = FALSE)
rm(list = ls())
gc()
require(pedometrics)
require(spatstat)
require(splancs)
require(zoo)
require(maptools)
require(sp)
require(plotKML)
require(plyr)
require(raster)
require(rgdal)
require(spgrass6)
require(lattice)
require(stpp)
require(lgcp)
load("sm-dnos-general.RData")
load("sm-dnos-ppa.RData")
ls()
# GRASS GIS
initGRASS(gisBase = "/usr/lib/grass64/", gisDbase = GRASSgisDbase,
          override = TRUE, location = "dnos-sm-rs", mapset = "predictions", 
          pid = Sys.getpid())
system("g.remove MASK")
system("g.list rast")
#
# DATA PREPARATION #############################################################
# 
# Observation window (owin) ----------------------------------------------------
dnos_owin0 <- ripras(x = coordinates(cal_data)[1:340, 1],
                     y = coordinates(cal_data)[1:340, 2])
unitname(dnos_owin0) <- c("meter", "meters")
str(dnos_owin0)
plot(dnos_owin0)
plot(cal_data, add = TRUE, cex = 0.3, pch = 3)
dnos_owin0 <- as(dnos_owin0, "SpatialPolygons")
dnos_owin0 <- as(dnos_owin0, "SpatialPolygonsDataFrame")
proj4string(dnos_owin0) <- proj4string(cal_data)
shapefile(dnos_owin0, paste(ppp_dir, "dnos-owin.shp", sep = ""), overwrite = TRUE)
writeVECT6(dnos_owin0, "dnos_owin", v.in.ogr_flags = c("overwrite"))
system("d.mon x0")
system("d.vect dnos_owin")
system("v.to.rast --o in=dnos_owin out=dnos_owin use=val")
system("d.rast dnos_owin")
system("r.mask dnos_owin")
dnos_owin <- raster(readRAST6("dnos.owin"))
dnos_owin <- crop(dnos_owin, extent(dnos_owin0))
dnos_owin <- as(dnos_owin, "SpatialGridDataFrame")
dnos_owin <- as.owin(as.im(dnos_owin))
str(dnos_owin)
plot(dnos_owin)
rm(dnos_owin0)
#
# Planar point patern (ppp)
dnos_ppp <- ppp(x = coordinates(cal_data)[1:340, 1],
                y = coordinates(cal_data)[1:340, 2],
                window = dnos_owin, check = TRUE)
unitname(dnos_ppp) <- c("meter", "meters")
str(dnos_ppp)
summary(dnos_ppp)
plot(dnos_ppp, cex = 0.3, pch = 3)
rm(dnos_owin)
gc()
#
# Dirichlet Tessellation (tess) ------------------------------------------------
dnos_tess <- dirichlet(dnos_ppp) # wait!
str(dnos_tess)
plot(dnos_tess)
dnos_tess$image$v <- matrix(as.numeric(dnos_tess$image$v), dim(dnos_tess$image))
dnos_tess <- raster(dnos_tess$image)
image(dnos_tess, col = rev(heat.colors(340)), asp = 1)
# reclassify values using the 'DATE' field of 'cal_data' object
head(cal_data$DATE)
day <- as.POSIXlt(cal_data$DATE[1:340], format = "%m/%d/%Y")
day <- as.factor(round(julian(day, origin = day[1])))
levels(day) <- order(as.numeric(levels(day)))
day
rcl <- cbind(seq(1:340), day)
dnos_tess <- reclassify(dnos_tess, rcl)
plot(dnos_tess, col = rev(terrain.colors(22)), asp = 1)
dnos_tess@data@values <- as.factor(dnos_tess@data@values)
dnos_tess <- as(dnos_tess, "SpatialGridDataFrame")
proj4string(dnos_tess) <- proj4string(cal_data)
image(dnos_tess, col = rev(heat.colors(22)), aspect = "iso")
dnos_tess$layer <- as.integer(dnos_tess$layer)
writeRAST6(dnos_tess, "dnos_tess", overwrite = TRUE)
rm(day, rcl, dnos_tess)
gc()
system("d.mon x0")
system("d.erase")
system("d.rast.leg dnos_tess")
#
# SPATIAL ANALYSIS #############################################################
#
# Complete Spatial Randomness --------------------------------------------------
# Inhomogeneous G-function
nearG <- Ginhom(dnos_ppp, correction = "best", sigma = bw.diggle, adjust = 2)
dev.off()
pdf(paste(ppp_dir, "gest.pdf", sep = ""))
plot(nearG, main = "Inhomogeneous G-function",
     xlab = "r (meters)", xaxp = c(0, 475, 19), ylim = c(-0.2, 1.2),
     legendpos = "float", legendargs = list(bty = "n"), xlim = c(0, 475))
dev.off()
nsim <- 99
alpha <- 1/(1 + nsim)
sub <- paste("n = ", nsim, " Monte Carlo simulations (alpha = ", alpha, ")", 
             sep = "")
nearGsim <- envelope(dnos_ppp, nsim = nsim, fun = Ginhom, global = TRUE,
                     correction = "best", sigma = bw.diggle, adjust = 2) # WAIT!
dev.off()
pdf(paste(ppp_dir, "gest-sim.pdf", sep = ""))
plot(nearGsim, sub = sub, main = "Envelope of the inhomogeneous G-function",
     xlab = "r (meters)", ylim = c(-0.2, 1.2), legendpos = "float",
     legendargs = list(bty = "n"), xaxp = c(0, 475, 19), xlim = c(0, 475))
dev.off()
rm(nsim, alpha)
gc()
#
# First-Order Properties -------------------------------------------------------
# Kernel density
dnos_kernel <- density(dnos_ppp, kernel = "gaussian", edge = TRUE,
                       sigma = bw.diggle(dnos_ppp), diggle = TRUE) # wait!
str(dnos_kernel)
dnos_kernel <- as(dnos_kernel, "SpatialGridDataFrame")
dnos_kernel$v <- dnos_kernel$v / max(dnos_kernel$v, na.rm = TRUE)
proj4string(dnos_kernel) <- proj4string(cal_data)
pts <- list("sp.points", cal_data[1:340, ], pch = 3, cex = 0.1, col = "black")
dev.off()
pdf(paste(ppp_dir, "rel-kernel.pdf", sep = ""))
spplot(dnos_kernel, aspect = "iso", sub = "Relative kernel dentity estimate",
       col.regions = topo.colors(100), sp.layout = list(pts))
dev.off()
writeRAST6(dnos_kernel, "dnos_kernel", overwrite = TRUE)
system("d.erase")
system("d.rast.leg dnos_kernel")
rm(pts, dnos_kernel)
gc()
# Empty space distances
dnos_empty <- distmap(dnos_ppp)
dnos_empty <- as(dnos_empty, "SpatialGridDataFrame")
pts <- list("sp.points", cal_data[1:340, ], pch = 3, cex = 0.1, col = "black")
dev.off()
pdf(paste(ppp_dir, "empty-space.pdf", sep = ""))
spplot(dnos_empty, aspect = "iso", sub = "Empty space distances",
       col.regions = topo.colors(100), sp.layout = list(pts))
dev.off()
proj4string(dnos_empty) <- proj4string(cal_data)
writeRAST6(dnos_empty, "dnos_empty", overwrite = TRUE)
rm(pts, dnos_empty)
gc()
system("d.mon x0")
system("d.erase")
system("d.rast.leg dnos_empty")
# Stienen diagram
dnos_stienen <- dnos_ppp %mark% (nndist(dnos_ppp)/2)
dnos_stienen <- as(dnos_stienen, "SpatialPointsDataFrame")
dev.off()
pdf(paste(ppp_dir, "stienen.pdf", sep = ""))
bubble(dnos_stienen, fill = FALSE, col = "black", maxsize = 2, main = "",
       sub = "Stienen diagram", aspect = "iso", key.space = "n", pch = 1, 
       lty = "b")
dev.off()
proj4string(dnos_stienen) <- proj4string(cal_data)
plotKML(obj = dnos_stienen, folder.name = paste(ppp_dir, sep = ""),
        point_names = "marks", file.name = "stienen.kml", alpha = 0.75,
        balloon = TRUE,
        shape = "http://maps.google.com/mapfiles/kml/pal2/icon18.png")
gc()
#
# Second-Order Properties ------------------------------------------------------
# Inhomogeneous K-function
nearK <- Kinhom(dnos_ppp, correction = "best", sigma = bw.diggle, adjust = 2)
dev.off()
pdf(paste(ppp_dir, "kest.pdf", sep = ""))
plot(nearK, main = "Inhomogeneous K-function", xlim = c(0, 475), 
     xaxp = c(0, 475, 19), xlab = "r (meters)", legendpos = "float",
     legendargs = list(bty = "n"))
dev.off()
nsim <- 99
alpha <- 1/(1 + nsim)
sub <- paste("n = ", nsim, " Monte Carlo simulations (alpha = ", alpha, ")", 
             sep = "")
nearKsim <- envelope(dnos_ppp, nsim = nsim, fun = Kinhom, global = TRUE, # WAIT!
                     correction = "best", sigma = bw.diggle, adjust = 2)
dev.off()
pdf(paste(ppp_dir, "kest-sim.pdf", sep = ""))
plot(nearKsim, sub = sub, main = "Envelope of the inhomogeneous K-function",
     xlab = "r (meters)", legendpos = "float", xlim = c(0, 475),
     xaxp = c(0, 475, 19), legendargs = list(bty = "n"))
dev.off()
gc()
#
# Model fitting ================================================================
# Prepare covariates
system("g.list rast")
ycoord <- readRAST6("dnos.raster")
colnames(ycoord@data) <- "ycoord"
a <- coordinates(ycoord)[, 2]
ycoord$ycoord <-  (a - min(a)) / (max(a) - min(a))
rm(a)
xcoord <- readRAST6("dnos.raster")
colnames(xcoord@data) <- "xcoord"
a <- coordinates(xcoord)[, 1]
xcoord$xcoord <- (a - min(a)) / (max(a) - min(a))
rm(a)
STR_90_255 <- readRAST6("STR_90_255")
STR_90_255$STR_90_255 <- as.factor(STR_90_255$STR_90_255)
dnos_tess <- readRAST6("dnos_tess")
dnos_tess$dnos_tess <- as.factor(dnos_tess$dnos_tess)
form <- formula(~ LU2009a + LU2009b + LU2009c + LU2009d + LU2009e + LU2009f +
                  STR_90_255 + dnos_tess + ycoord + xcoord)
# Nonstationary Poisson point proccess model
ppp_fit <-
  ppm(dnos_ppp, form,
      covariates = list(LU2009a    = as.im(readRAST6("LU2009a")), # Native forest
                        LU2009b    = as.im(readRAST6("LU2009b")), # Shrubland
                        LU2009c    = as.im(readRAST6("LU2009c")), # Animal husbandry
                        LU2009d    = as.im(readRAST6("LU2009d")), # Crop agriculture
                        LU2009e    = as.im(readRAST6("LU2009e")), # Dist from settlements
                        LU2009f    = as.im(readRAST6("LU2009f")), # Dist from boundaries
                        STR_90_255 = as.im(STR_90_255),
                        dnos_tess  = as.im(dnos_tess),
                        ycoord     = as.im(ycoord),
                        xcoord     = as.im(xcoord)))
# Model diagnostics
attributes(ppp_fit)
round(ppp_fit$coef, 3)
anova.ppm(ppp_fit, test = "Chisq")
dev.off()
pdf(paste(ppp_dir, "diagn-ppm-fit.pdf", sep = ""))
diagnose.ppm(ppp_fit, cex = 0.3) # wait!
dev.off()
# Residuals
ppp_res <- residuals(ppp_fit)
ppp_res <- Smooth(ppp_res) # wait!
ppp_res <- as(ppp_res, "SpatialGridDataFrame")
ppp_res$v <- ppp_res$v / max(abs(ppp_res$v), na.rm = TRUE)
proj4string(ppp_res) <- proj4string(cal_data)
pts <- list("sp.points", cal_data[1:340, ], pch = 3, cex = 0.1, col = "black")
dev.off()
pdf(paste(ppp_dir, "rel-residuals.pdf", sep = ""))
spplot(ppp_res, aspect = "iso", sub = "Relative residuals",
       col.regions = bpy.colors(64), sp.layout = list(pts))
dev.off()
rm(pts)
writeRAST6(ppp_res, "ppp_res", overwrite = TRUE)
rm(ppp_res)
gc()
system("d.mon x0")
system("d.rast.leg ppp_res")
#
# Fitted trend
ppp_trend <-
  predict(ppp_fit, locations = dnos_ppp$window,
          covariates = list(LU2009a = as.im(readRAST6("LU2009a")),
                            LU2009b = as.im(readRAST6("LU2009b")),
                            LU2009c = as.im(readRAST6("LU2009c")),
                            LU2009d = as.im(readRAST6("LU2009d")),
                            LU2009e = as.im(readRAST6("LU2009e")),
                            LU2009f = as.im(readRAST6("LU2009f")),
                            STR_90_255 = as.im(STR_90_255),
                            dnos_tess = as.im(dnos_tess),
                            ycoord     = as.im(ycoord),
                            xcoord     = as.im(xcoord)))
rm(dnos_tess, form, xcoord, ycoord, STR_90_255)
gc()
ppp_trend <- as(ppp_trend, "SpatialGridDataFrame")
ppp_trend$v <- ppp_trend$v / max(abs(ppp_trend$v), na.rm = TRUE)
proj4string(ppp_trend) <- proj4string(cal_data)
pts <- list("sp.points", cal_data[1:340, ], pch = 3, cex = 0.1, col = "black")
dev.off()
pdf(paste(ppp_dir, "rel-trend.pdf", sep = ""))
spplot(ppp_trend, aspect = "iso", sub = "Relative fitted trend",
       col.regions = topo.colors(64), sp.layout = list(pts))
dev.off()
rm(pts)
writeRAST6(ppp_trend, "ppp_trend", overwrite = TRUE)
rm(ppp_trend)
gc()
system("d.mon x0")
system("d.rast.leg ppp_trend")
#
# Save figures with covariates
# Land use
map <- readRAST6("LU2009")
map@bbox <- cal_data@bbox
map@data[, 1] <- as.factor(map@data[, 1])
map@data[, 1] <- revalue(map@data[, 1], c("1" = "FS",
                                          "2" = "SS",
                                          "3" = "H",
                                          "4" = "AA",
                                          "5" = "PF",
                                          "6" = "S",
                                          "7" = "O"))
data(worldgrids_pal)
color <- c(worldgrids_pal$glc2000[1], worldgrids_pal$glc2000[12],
           worldgrids_pal$glc2000[13], worldgrids_pal$glc2000[16],
           worldgrids_pal$glc2000[4], worldgrids_pal$glc2000[22],
           worldgrids_pal$glc2000[20])
pts <- list("sp.points", cal_data[1:340, ], pch = 3, col = "black")
p1 <- spplot(map, main = "", col.regions = color, asp = 1, 
             colorkey = FALSE, sp.layout = list(pts), sub = "Land use")
dev.off()
pdf(file = paste(ppp_dir, "land-use.pdf", sep = ""))
plot(p1)
dev.off()
rm(map, pts)
gc()
# Physiographic strata
map <- readRAST6("STR_90_255")
map@bbox <- cal_data@bbox
map@data[, 1] <- as.factor(map@data[, 1])
pts <- list("sp.points", cal_data[1:340, ], pch = 3, col = "black")
p1 <- spplot(map, main = "", asp = 1, colorkey = FALSE, sp.layout = list(pts),
             sub = "Physiographic strata")
dev.off()
pdf(file = paste(ppp_dir, "physio-strata.pdf", sep = ""))
plot(p1)
dev.off()
rm(map, pts)
gc()
# Field campaigns
map <- readRAST6("dnos_tess")
map@bbox <- cal_data@bbox
map@data[, 1] <- as.factor(map@data[, 1])
pts <- list("sp.points", cal_data[1:340, ], pch = 3, col = "black")
p1 <- spplot(map, main = "", asp = 1, colorkey = FALSE, sp.layout = list(pts),
             sub = "Field campaigns", col.regions = rev(heat.colors(22)))
dev.off()
pdf(file = paste(ppp_dir, "field-campaigns.pdf", sep = ""))
plot(p1)
dev.off()
rm(map, pts)
gc()
# Topography
map <- readRAST6("SHADE_10")
map@bbox <- cal_data@bbox
pts <- list("sp.points", cal_data[1:340, ], pch = 3, col = "black")
p1 <- spplot(map, main = "", asp = 1, colorkey = FALSE, sp.layout = list(pts),
             sub = "Topography", col.regions = gray(seq(0, 1, 1/100)))
dev.off()
pdf(file = paste(ppp_dir, "topography.pdf", sep = ""))
plot(p1)
dev.off()
rm(map, pts)
gc()
#
# Envelope of the detrended G function
# nsim <- 99
# alpha <- 1/(1 + nsim)
# sub <- paste("n = ", nsim, " Monte Carlo simulations (alpha = ", alpha, ")", 
#              sep = "")
# fitGsim <- envelope(dnos_ppp, nsim = nsim, fun = Ginhom, global = TRUE, 
#                     lambda = ppp_fit)
# dev.off()
# pdf(paste(ppp_dir, "fit-gest-sim.pdf", sep = ""))
# plot(fitGsim, xaxp = c(0, 475, 19), legendargs = list(bty = "n"),
#      xlab = "r (meters)", ylim = c(-0.2, 1.2), sub = sub, legendpos = "float",
#      main = "Envelope of the detrended\ninhomogeneous G function")
# dev.off()
#
# Envelope of the detrended K function
# nsim <- 3
# alpha <- 1/(1 + nsim)
# sub <- paste("n = ", nsim, " Monte Carlo simulations (alpha = ", alpha, ")", 
#              sep = "")
# fitKsim <- envelope(dnos_ppp, nsim = nsim, fun = Kinhom, global = TRUE, 
#                     lambda = ppp_fit)
# dev.off()
# pdf(paste(ppp_dir, "fit-kest-sim.pdf", sep = ""))
# plot(fitKsim, xaxp = c(0, 475, 19), legendargs = list(bty = "n"),
#      xlab = "r (meters)", ylim = c(-0.2, 1.2), sub = sub, legendpos = "float",
#      main = "Envelope of the detrended\ninhomogeneous K function")
# dev.off()
#
# TEMPORAL ANALYSIS ############################################################
# 
# nearest neighbour distances
dnos_near <- nndistG(as.points(coordinates(cal_data[1:340, ])))$dists
round(quantile(dnos_near, probs = seq(0, 1, by = 0.25)), 0)
plotHD(dnos_near)
day <- as.POSIXlt(cal_data$DATE[1:340], format = "%m/%d/%Y")
day <- as.factor(round(julian(day, origin = day[1])))
levels(day) <- order(as.numeric(levels(day)))
dnos_near <- data.frame(id = c(1:length(dnos_near)), near = dnos_near, 
                        day = day)
dnos_near <- arrange(dnos_near, day)
dnos_near$seq <- c(1:dim(dnos_near)[1])
stats <- data.frame(mean_near = c(by(dnos_near$near, dnos_near$day, mean)),
                    sd_near = c(by(dnos_near$near, dnos_near$day, sd)),
                    min_near = c(by(dnos_near$near, dnos_near$day, min)),
                    max_near = c(by(dnos_near$near, dnos_near$day, max)),
                    mean_day = c(by(c(1:dim(dnos_near)[1]), dnos_near$day, 
                                    mean)))
dnos_near <- list(all = dnos_near, stats = stats)
rm(stats)
gc()
# auto-correlation of first order r(1) and polynomial fit
auto_cor <- round(cor(dnos_near$all$near[1:length(dnos_near$all$near) - 1],
                      dnos_near$all$near[2:length(dnos_near$all$near)]), 2)
poly_fit <- lm(dnos_near$all$near ~ poly(dnos_near$all$seq, degree = 3))
summary(poly_fit)
plot(dnos_near$all$near ~ dnos_near$all$seq)
lines(dnos_near$all$seq, predict(poly_fit))
text(40, 300, paste("r(1) = ", auto_cor, sep = ""))
dev.off()
pdf(paste(ppp_dir, "nndistG.pdf", sep = ""))
xyplot(dnos_near$all$near ~ c(1:length(dnos_near$all$day)),
       xlab = "Observation sequence", main = "", 
       ylab = "Nearest neighbour distance (m)", asp = 1,
       sub = "Solid line = rolling mean (k = 34); Dashed line = polynomial fit",
       panel = function(x, y) {
         panel.xyplot(x, y)
         panel.lines(rollmean(x = dnos_near$all$near, k = 34), col = "black")
         panel.lines(predict(poly_fit), col = "black", lty = 2)
         panel.text(0, max(dnos_near$all$near), pos = 4,
                    labels = paste("r(1) = ", auto_cor, sep = ""))
         panel.text(0, max(dnos_near$all$near) * 0.95, pos = 4,
                    labels = paste("Adj R2 = ", 
                                   round(summary(poly_fit)$r.squared, 2), 
                                   sep = ""))
       })
dev.off()
# daily polynomial fit (mean, minimum, and maximum)
poly_fit_mean <- lm(dnos_near$stats$mean_near ~ poly(dnos_near$stats$mean_day, 
                                                     degree = 3))
poly_fit_min <- lm(dnos_near$stats$min_near ~ poly(dnos_near$stats$mean_day, 
                                                     degree = 3))
poly_fit_max <- lm(dnos_near$stats$max_near ~ poly(dnos_near$stats$mean_day, 
                                                   degree = 3))
dev.off()
pdf(paste(ppp_dir, "nndistG-daily.pdf", sep = ""))
xyplot(dnos_near$stats$mean_near ~ c(1:length(dnos_near$stats$mean_day)),
       xlab = "Field campaign", main = "", 
       ylim = c(0, max(dnos_near$stats$max_near) * 1.05),
       ylab = "Nearest neighbour distance (m)", asp = 1,
       sub = "Black = daily average; Red = daily minimum; Blue = daily maximum",
       panel = function(x, y) {
         panel.xyplot(x, y,, col = "black", pch = 20)
         panel.lines(predict(poly_fit_mean), col = "black", lty = 1)
         panel.lines(predict(poly_fit_max), col = "blue", lty = 2)
         panel.lines(predict(poly_fit_min), col = "red", lty = 2)
         panel.points(dnos_near$stats$min_near ~ 
                        c(1:length(dnos_near$stats$mean_day)),
                      col = "red", pch = 20)
         panel.points(dnos_near$stats$max_near ~ 
                        c(1:length(dnos_near$stats$mean_day)),
                      col = "blue", pch = 20)
       })
dev.off()
pdf(paste(ppp_dir, "nndistG-days.pdf", sep = ""))
trellis.par.set(fontsize = list(text = 8, points = 4))
xyplot(dnos_near$all$near ~ dnos_near$all$seq | dnos_near$all$day,
       xlab = "Observation sequence", main = "", 
       ylab = "Nearest neighbour distance (m)", asp = 1,
       panel = function(x, y) {
         panel.xyplot(x, y)
         panel.lines(rollmean(x = dnos_near$all$near, k = 34), col = "black",
                     lwd = 0.5)
         panel.points(dnos_near$stats$mean_near ~ dnos_near$stats$mean_day,
                      col = "red", pch = 2)
       })
dev.off()
gc()
# Animation
tmp <- data.frame(coordinates(cal_data)[1:340, ], dnos_near$all$day)
colnames(tmp) <- c("x", "y", "t")
tmp <- as.3dpoints(tmp)
animation(tmp, runtime = 15, cex = 0.5)
# Sampling model (motivation plot) -----------------------------------
dev.off()
pdf("sample-model.pdf", width = 3, height = 3, pointsize = 8)
par(mgp = c(2,0,0))
plot(dnos_near$all$near, cex = 0.8, pch = 20, xlab = "Time", ylim = c(0, 290),
     ylab = "Nearest neighbour distance", bty = "n",
     xaxp = c(0, 350, 14), type = "n", yaxt="n", xaxt="n")
lines(x = c(0, 160), y = c(25, 250))
lines(x = c(160, 160), y = c(25, 250), lty = 2)
text(x = 80, y = 250, labels = "Phase I", pos = 3)
text(x = 80, y = 18, labels = "MFM", pos = 1)
arrows(x0 = 0, x1 = 160, y0 = 20, y1 = 20, code = 2, length = 0.05)
lines(x = c(160, 285), y = c(250, 250))
lines(x = c(285, 285), y = c(25, 250), lty = 2)
text(x = 227.5, y = 250, labels = "Phase II", pos = 3)
text(x = 227.5, y = 18, labels = "OFM", pos = 1)
arrows(x0 = 160, x1 = 285, y0 = 20, y1 = 20, code = 2, length = 0.05)
lines(x = c(285, 340), y = c(250, 25))
text(x = 312.5, y = 250, labels = "Phase III", pos = 3)
text(x = 312.5, y = 18, labels = "MFM", pos = 1)
arrows(x0 = 285, x1 = 340, y0 = 20, y1 = 20, code = 2, length = 0.05)
dev.off()

# Save ppp data ######################################################
save(dnos_ppp, dnos_near, dnos_stienen, nearG, nearGsim, nearK, nearKsim, 
     ppp_fit, poly_fit, poly_fit_mean, poly_fit_min, poly_fit_max,
     file = "sm-dnos-ppa.RData")
# End!