# Ground Control Points (Landsat 5 TM 2009jan05) -------------------------------
gcp.2009jan05 <- shapefile(paste(point.dir, "gcp_2009jan05.shp", sep = ""))
colnames(gcp.2009jan05@data) <- "siteID"
gcp.2009jan05@data$siteID <- as.numeric(gcp.2009jan05$siteID)
gcp.2009jan05 <- gcp.2009jan05[order(gcp.2009jan05$siteID),]
gcp.2009jan05 <- spTransform(gcp.2009jan05, wgs1984utm22s)
plot(gcp.2009jan05)

# Ground Control Points (Landsat 5 TM 2008dez20) -------------------------------
gcp.2008dez20 <- shapefile(paste(point.dir, "gcp_2008dez20.shp", sep = ""))
colnames(gcp.2008dez20@data) <- "siteID"
gcp.2008dez20@data$siteID <- as.numeric(gcp.2008dez20$siteID)
gcp.2008dez20 <- gcp.2008dez20[order(gcp.2008dez20$siteID),]
gcp.2008dez20 <- spTransform(gcp.2008dez20, wgs1984utm22s)
plot(gcp.2008dez20)

# Ground Control Points (Landsat 5 TM 2006jan29) -------------------------------
gcp.2006jan29 <- shapefile(paste(point.dir, "gcp_2006jan29.shp", sep = ""))
colnames(gcp.2006jan29@data) <- "siteID"
gcp.2006jan29@data$siteID <- as.numeric(gcp.2006jan29$siteID)
gcp.2006jan29 <- gcp.2006jan29[order(gcp.2006jan29$siteID),]
gcp.2006jan29 <- spTransform(gcp.2006jan29, wgs1984utm22s)
plot(gcp.2006jan29)

# Ground Control Points (Landsat 5 TM 2006dez31) -------------------------------
gcp.2006dez31 <- shapefile(paste(point.dir, "gcp_2006dez31.shp", sep = ""))
colnames(gcp.2006dez31@data) <- "siteID"
gcp.2006dez31@data$siteID <- as.numeric(gcp.2006dez31$siteID)
gcp.2006dez31 <- gcp.2006dez31[order(gcp.2006dez31$siteID),]
gcp.2006dez31 <- spTransform(gcp.2006dez31, wgs1984utm22s)
plot(gcp.2006dez31)

# Ground Control Points (Landsat 5 TM 2005jan26) -------------------------------
gcp.2005jan26 <- shapefile(paste(point.dir, "gcp_2005jan26.shp", sep = ""))
colnames(gcp.2005jan26@data) <- "siteID"
gcp.2005jan26@data$siteID <- as.numeric(gcp.2005jan26$siteID)
gcp.2005jan26 <- gcp.2005jan26[order(gcp.2005jan26$siteID),]
gcp.2005jan26 <- spTransform(gcp.2005jan26, wgs1984utm22s)
plot(gcp.2005jan26)

# Ground Control Points (Landsat 5 TM 1998dez25) -------------------------------
gcp.1998dez25 <- shapefile(paste(point.dir, "gcp_1998dez25.shp", sep = ""))
colnames(gcp.1998dez25@data) <- "siteID"
gcp.1998dez25@data$siteID <- as.numeric(gcp.1998dez25$siteID)
gcp.1998dez25 <- gcp.1998dez25[order(gcp.1998dez25$siteID),]
gcp.1998dez25 <- spTransform(gcp.1998dez25, wgs1984utm22s)
plot(gcp.1998dez25)

# Ground Control Points (Landsat 5 TM 1996dez03) -------------------------------
gcp.1996dez03 <- shapefile(paste(point.dir, "gcp_1996dez03.shp", sep = ""))
colnames(gcp.1996dez03@data) <- "siteID"
gcp.1996dez03@data$siteID <- as.numeric(gcp.1996dez03$siteID)
gcp.1996dez03 <- gcp.1996dez03[order(gcp.1996dez03$siteID),]
gcp.1996dez03 <- spTransform(gcp.1996dez03, wgs1984utm22s)
plot(gcp.1996dez03)

# Ground Control Points (Landsat 5 TM 1994jan12) -------------------------------
gcp.1994jan12 <- shapefile(paste(point.dir, "gcp_1994jan12.shp", sep = ""))
colnames(gcp.1994jan12@data) <- "siteID"
gcp.1994jan12@data$siteID <- as.numeric(gcp.1994jan12$siteID)
gcp.1994jan12 <- gcp.1994jan12[order(gcp.1994jan12$siteID),]
gcp.1994jan12 <- spTransform(gcp.1994jan12, wgs1984utm22s)
plot(gcp.1994jan12)

# Ground Control Points (Landsat 5 TM 1991jan04) -------------------------------
gcp.1991jan04 <- shapefile(paste(point.dir, "gcp_1991jan04.shp", sep = ""))
colnames(gcp.1991jan04@data) <- "siteID"
gcp.1991jan04@data$siteID <- as.numeric(gcp.1991jan04$siteID)
gcp.1991jan04 <- gcp.1991jan04[order(gcp.1991jan04$siteID),]
gcp.1991jan04 <- spTransform(gcp.1991jan04, wgs1984utm22s)
plot(gcp.1991jan04)

# Landsat 5 TM (2009jan05) =====================================================
# Prepare spsurvey.object
spsurvey.2009jan05 <- 
  spsurvey.analysis(design = coordenadas(gcp.2009jan05),
                    data.cont = deltagcp(gcp.pos, gcp.2009jan05),
                    popcorrect = TRUE,
                    pcfsize = length(gcp.2009jan05$siteID),
                    support = rep(1, length(gcp.2009jan05$siteID)),
                    wgt = rep(1, length(gcp.2009jan05$siteID)),
                    vartype = "SRS")
# CDF
cdf.2009jan05 <- cont.analysis(spsurvey.obj = spsurvey.2009jan05)
levels(cdf.2009jan05$Pct$Indicator)
cdfstats(cdf.2009jan05, "dx", all = TRUE)
cdftable(cdf.2009jan05)
# Plot azimuth
dev.off()
cairo_pdf(paste(covar.validation.dir, "azim-2009jan05.pdf", sep = ""))
par(cex = 0.5, lwd = 0.1, ps=25)
DrawHistogram(deltagcp(gcp.pos, gcp.2009jan05)$azimuth, Direction = 2)
title(xlab = "Landsat imagery (2009jan05)")
dev.off()

# Landsat 5 TM (2008dez20) =====================================================
# Prepare spsurvey.object
spsurvey.2008dez20 <- 
  spsurvey.analysis(design = coordenadas(gcp.2008dez20),
                    data.cont = deltagcp(gcp.pos, gcp.2008dez20),
                    popcorrect = TRUE,
                    pcfsize = length(gcp.2008dez20$siteID),
                    support = rep(1, length(gcp.2008dez20$siteID)),
                    wgt = rep(1, length(gcp.2008dez20$siteID)),
                    vartype = "SRS")
# CDF
cdf.2008dez20 <- cont.analysis(spsurvey.obj = spsurvey.2008dez20)
levels(cdf.2008dez20$Pct$Indicator)
cdfstats(cdf.2008dez20, "dx", all = TRUE)
cdftable(cdf.2008dez20)
# Plot azimuth
dev.off()
cairo_pdf(paste(covar.validation.dir, "azim-2008dez20.pdf", sep = ""))
par(cex = 0.5, lwd = 0.1, ps=25)
DrawHistogram(deltagcp(gcp.pos, gcp.2008dez20)$azimuth, Direction = 2)
title(xlab = "Landsat imagery (2008dez20)")
dev.off()

# Landsat 5 TM (2006jan29) =====================================================
# Prepare spsurvey.object
spsurvey.2006jan29 <- 
  spsurvey.analysis(design = coordenadas(gcp.2006jan29),
                    data.cont = deltagcp(gcp.pos, gcp.2006jan29),
                    popcorrect = TRUE,
                    pcfsize = length(gcp.2006jan29$siteID),
                    support = rep(1, length(gcp.2006jan29$siteID)),
                    wgt = rep(1, length(gcp.2006jan29$siteID)),
                    vartype = "SRS")
# CDF
cdf.2006jan29 <- cont.analysis(spsurvey.obj = spsurvey.2006jan29)
levels(cdf.2006jan29$Pct$Indicator)
cdfstats(cdf.2006jan29, "dx", all = TRUE)
cdftable(cdf.2006jan29)
# Plot azimuth
dev.off()
cairo_pdf(paste(covar.validation.dir, "azim-2006jan29.pdf", sep = ""))
par(cex = 0.5, lwd = 0.1, ps=25)
DrawHistogram(deltagcp(gcp.pos, gcp.2006jan29)$azimuth, Direction = 2)
title(xlab = "Landsat imagery (2006jan29)")
dev.off()

# Landsat 5 TM (2006dez31) =====================================================
# Prepare spsurvey.object
spsurvey.2006dez31 <- 
  spsurvey.analysis(design = coordenadas(gcp.2006dez31),
                    data.cont = deltagcp(gcp.pos, gcp.2006dez31),
                    popcorrect = TRUE,
                    pcfsize = length(gcp.2006dez31$siteID),
                    support = rep(1, length(gcp.2006dez31$siteID)),
                    wgt = rep(1, length(gcp.2006dez31$siteID)),
                    vartype = "SRS")
# CDF
cdf.2006dez31 <- cont.analysis(spsurvey.obj = spsurvey.2006dez31)
levels(cdf.2006dez31$Pct$Indicator)
cdfstats(cdf.2006dez31, "dx", all = TRUE)
cdftable(cdf.2006dez31)
# Plot azimuth
dev.off()
cairo_pdf(paste(covar.validation.dir, "azim-2006dez31.pdf", sep = ""))
par(cex = 0.5, lwd = 0.1, ps=25)
DrawHistogram(deltagcp(gcp.pos, gcp.2006dez31)$azimuth, Direction = 2)
title(xlab = "Landsat imagery (2006dez31)")
dev.off()

# Landsat 5 TM (2005jan26) =====================================================
# Prepare spsurvey.object
spsurvey.2005jan26 <- 
  spsurvey.analysis(design = coordenadas(gcp.2005jan26),
                    data.cont = deltagcp(gcp.pos, gcp.2005jan26),
                    popcorrect = TRUE,
                    pcfsize = length(gcp.2005jan26$siteID),
                    support = rep(1, length(gcp.2005jan26$siteID)),
                    wgt = rep(1, length(gcp.2005jan26$siteID)),
                    vartype = "SRS")
# CDF
cdf.2005jan26 <- cont.analysis(spsurvey.obj = spsurvey.2005jan26)
levels(cdf.2005jan26$Pct$Indicator)
cdfstats(cdf.2005jan26, "dx", all = TRUE)
cdftable(cdf.2005jan26)
# Plot azimuth
dev.off()
cairo_pdf(paste(covar.validation.dir, "azim-2005jan26.pdf", sep = ""))
par(cex = 0.5, lwd = 0.1, ps=25)
DrawHistogram(deltagcp(gcp.pos, gcp.2005jan26)$azimuth, Direction = 2)
title(xlab = "Landsat imagery (2005jan26)")
dev.off()

# Landsat 5 TM (1998dez25) =====================================================
# Prepare spsurvey.object
spsurvey.1998dez25 <- 
  spsurvey.analysis(design = coordenadas(gcp.1998dez25),
                    data.cont = deltagcp(gcp.pos, gcp.1998dez25),
                    popcorrect = TRUE,
                    pcfsize = length(gcp.1998dez25$siteID),
                    support = rep(1, length(gcp.1998dez25$siteID)),
                    wgt = rep(1, length(gcp.1998dez25$siteID)),
                    vartype = "SRS")
# CDF
cdf.1998dez25 <- cont.analysis(spsurvey.obj = spsurvey.1998dez25)
levels(cdf.1998dez25$Pct$Indicator)
cdfstats(cdf.1998dez25, "dx", all = TRUE)
cdftable(cdf.1998dez25)
# Plot azimuth
dev.off()
cairo_pdf(paste(covar.validation.dir, "azim-1998dez25.pdf", sep = ""))
par(cex = 0.5, lwd = 0.1, ps=25)
DrawHistogram(deltagcp(gcp.pos, gcp.1998dez25)$azimuth, Direction = 2)
title(xlab = "Landsat imagery (1998dez25)")
dev.off()

# Landsat 5 TM (1996dez03) =====================================================
# Prepare spsurvey.object
spsurvey.1996dez03 <- 
  spsurvey.analysis(design = coordenadas(gcp.1996dez03),
                    data.cont = deltagcp(gcp.pos, gcp.1996dez03),
                    popcorrect = TRUE,
                    pcfsize = length(gcp.1996dez03$siteID),
                    support = rep(1, length(gcp.1996dez03$siteID)),
                    wgt = rep(1, length(gcp.1996dez03$siteID)),
                    vartype = "SRS")
# CDF
cdf.1996dez03 <- cont.analysis(spsurvey.obj = spsurvey.1996dez03)
levels(cdf.1996dez03$Pct$Indicator)
cdfstats(cdf.1996dez03, "dx", all = TRUE)
cdftable(cdf.1996dez03)
# Plot azimuth
dev.off()
cairo_pdf(paste(covar.validation.dir, "azim-1996dez03.pdf", sep = ""))
par(cex = 0.5, lwd = 0.1, ps=25)
DrawHistogram(deltagcp(gcp.pos, gcp.1996dez03)$azimuth, Direction = 2)
title(xlab = "Landsat imagery (1996dez03)")
dev.off()

# Landsat 5 TM (1994jan12) =====================================================
# Prepare spsurvey.object
spsurvey.1994jan12 <- 
  spsurvey.analysis(design = coordenadas(gcp.1994jan12),
                    data.cont = deltagcp(gcp.pos, gcp.1994jan12),
                    popcorrect = TRUE,
                    pcfsize = length(gcp.1994jan12$siteID),
                    support = rep(1, length(gcp.1994jan12$siteID)),
                    wgt = rep(1, length(gcp.1994jan12$siteID)),
                    vartype = "SRS")
# CDF
cdf.1994jan12 <- cont.analysis(spsurvey.obj = spsurvey.1994jan12)
levels(cdf.1994jan12$Pct$Indicator)
cdfstats(cdf.1994jan12, "dx", all = TRUE)
cdftable(cdf.1994jan12)
# Plot azimuth
dev.off()
cairo_pdf(paste(covar.validation.dir, "azim-1994jan12.pdf", sep = ""))
par(cex = 0.5, lwd = 0.1, ps=25)
DrawHistogram(deltagcp(gcp.pos, gcp.1994jan12)$azimuth, Direction = 2)
title(xlab = "Landsat imagery (1994jan12)")
dev.off()

# Landsat 5 TM (1991jan04) =====================================================
# Prepare spsurvey.object
spsurvey.1991jan04 <- 
  spsurvey.analysis(design = coordenadas(gcp.1991jan04),
                    data.cont = deltagcp(gcp.pos, gcp.1991jan04),
                    popcorrect = TRUE,
                    pcfsize = length(gcp.1991jan04$siteID),
                    support = rep(1, length(gcp.1991jan04$siteID)),
                    wgt = rep(1, length(gcp.1991jan04$siteID)),
                    vartype = "SRS")
# CDF
cdf.1991jan04 <- cont.analysis(spsurvey.obj = spsurvey.1991jan04)
levels(cdf.1991jan04$Pct$Indicator)
cdfstats(cdf.1991jan04, "dx", all = TRUE)
cdftable(cdf.1991jan04)
# Plot azimuth
dev.off()
cairo_pdf(paste(covar.validation.dir, "azim-1991jan04.pdf", sep = ""))
par(cex = 0.5, lwd = 0.1, ps=25)
DrawHistogram(deltagcp(gcp.pos, gcp.1991jan04)$azimuth, Direction = 2)
title(xlab = "Landsat imagery (1991jan04)")
dev.off()