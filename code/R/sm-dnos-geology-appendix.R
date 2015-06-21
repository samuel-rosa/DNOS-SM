# FIT MODEL TO PREDICT GEO_25 ##################################################

# Large Scale Model ============================================================

# initiate GRASS GIS
initGRASS(gisBase = "/usr/lib/grass64/", gisDbase = GRASSgisDbase, override = TRUE,
          location = "dnos-sm-rs", mapset = "predictions", pid = Sys.getpid())
system("g.region rast=dnos.raster")

# generate random points
system("v.random --o out=sampleGEO_25 n=20000")
# system("d.mon x0");system("d.vect sampleGEO_25")

# create table and sample from rasters
# system("v.db.droptable -f map=sampleGEO_25")
system("v.db.addtable --o map=sampleGEO_25 layer=1 columns='cat INTEGER, ELEV_10 DOUBLE PRECISION, GEO_25 INTEGER, GEO_50 INTEGER'")
system("v.what.rast vector=sampleGEO_25 raster=ELEV_10 layer=1 column=ELEV_10")
system("v.what.rast vector=sampleGEO_25 raster=GEO_25 layer=1 column=GEO_25")
system("v.what.rast vector=sampleGEO_25 raster=GEO_50 layer=1 column=GEO_50")
system("v.info -c sampleGEO_25")

# import data to R
sampleGEO_25 <- readVECT6("sampleGEO_25", layer = 1)
str(sampleGEO_25@data)
sampleGEO_25$GEO_25 <- as.factor(sampleGEO_25$GEO_25)
sampleGEO_25$GEO_50 <- as.factor(sampleGEO_25$GEO_50)
str(sampleGEO_25@data)
covarsGEO_25 <- stack(list(raster(readRAST6("ELEV_10")), raster(readRAST6("GEO_50"))))
summary(values(covarsGEO_25))

# fit classification tree model: large scale model
fitGEO_25a <- rpart(GEO_25 ~ ELEV_10, data = sampleGEO_25, method = "class")
fitGEO_25a
printcp(fitGEO_25a)
rsq.rpart(fitGEO_25a)
plot(fitGEO_25a, uniform=TRUE);text(fitGEO_25a, use.n = TRUE, all=TRUE, cex=.8)
plotcp(fitGEO_25a)

# make predictions
predGEO_25a <- predict(object = covarsGEO_25, model = fitGEO_25a, na.rm = TRUE, type = "class")
predGEO_25a
predGEO_25a <- as(predGEO_25a, "SpatialGridDataFrame")
predGEO_25a$value <- as.integer(predGEO_25a$value)
image(predGEO_25a, asp = 1)

# export to GRASS
writeRAST6(predGEO_25a, "predGEO_25a", overwrite = TRUE)
system("d.mon x0");system("d.rast.leg predGEO_25a")

# add category names
system("r.support map=predGEO_25a raster=GEO_25")
system("r.category predGEO_25a")
system("d.mon x0");system("d.rast.leg predGEO_25a")

# Small Scale Model ============================================================

# select data
predGEO_25b.ext <- extent(shapefile(paste(bounding.dir,"geo25-model.shp", sep = "")))
covarsGEO_25b <- stack(list(raster(readRAST6("ELEV_10")), raster(readRAST6("GEO_50")),
                            raster(readRAST6("GEO_25"))))
covarsGEO_25b <- crop(covarsGEO_25b, predGEO_25b.ext)
sampleGEO_25b <- crop(covarsGEO_25b, extent(covarsGEO_25b$GEO_25))
covarsGEO_25b <- dropLayer(covarsGEO_25b, 3)
sampleGEO_25b <- sampleRandom(sampleGEO_25b, 1000, sp = TRUE)
sampleGEO_25b$GEO_25 <- as.factor(sampleGEO_25b$GEO_25)
sampleGEO_25b$GEO_50 <- as.factor(sampleGEO_25b$GEO_50)
# str(sampleGEO_25b)
summary(sampleGEO_25b$GEO_25)
pts1 <- data.frame(sampleGEO_25b)[which(sampleGEO_25b$GEO_25 == 1),]
set.seed(1964)
pts1 <- pts1[sample(nrow(pts1), 200),]
pts2 <- data.frame(sampleGEO_25b)[which(sampleGEO_25b$GEO_25 != 1),]
sampleGEO_25b <- rbind(pts1, pts2);rm(pts1, pts2)
coordinates(sampleGEO_25b) <- ~ x + y
summary(sampleGEO_25b$GEO_25)
image(covarsGEO_25b, 2, add = TRUE)
plot(sampleGEO_25b, add = TRUE)

# fit model
fitGEO_25b <- rpart(GEO_25 ~ ELEV_10 + GEO_50, data = sampleGEO_25b, method = "class")
fitGEO_25b
printcp(fitGEO_25b)
rsq.rpart(fitGEO_25b)
plot(fitGEO_25b, uniform=TRUE);text(fitGEO_25b, use.n = TRUE, all=TRUE, cex=.8)
plotcp(fitGEO_25b)

# make predictions
predGEO_25b <- predict(object = covarsGEO_25b, model = fitGEO_25b, type = "class")
predGEO_25b
predGEO_25b <- as(predGEO_25b, "SpatialGridDataFrame")
predGEO_25b$value <- as.integer(predGEO_25b$value)
image(predGEO_25b, asp = 1)
plot(sampleGEO_25b, add = TRUE)

# export to GRASS
writeRAST6(predGEO_25b, "predGEO_25b", overwrite = TRUE)
# system("d.mon x0");system("d.rast.leg predGEO_25b")

# add category names
system("r.support map=predGEO_25b raster=GEO_25")
# system("r.category predGEO_25b")
# system("d.mon x0");system("d.rast.leg predGEO_25b")

# use only predictions for the Botucatu Formation
# system("r.category predGEO_25b")
system("r.reclass --o in=predGEO_25b out=botucatu25 <<EOF
       1 4 = NULL
       2 = 2 Botucatu Formation
       end
       EOF")
system("r.support map=botucatu25 raster=GEO_25")
# system("d.mon x0");system("d.rast.leg botucatu25")

# update GEO_25 ================================================================

# rpart model
system("r.patch --o input=GEO_25,botucatu25,predGEO_25a out=predGEO_25")