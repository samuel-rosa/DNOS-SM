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
initGRASS(gisBase = "/usr/lib/grass64/", gisDbase = GRASSgisDbase, 
          location = "dnos-sm-rs",
          mapset = "predictions", pid = Sys.getpid(), override = TRUE)

cmd <- paste("r.what input=ELEV_10 east_north=", coordinates(cal_data)[1, 1], ",",
             coordinates(cal_data)[1, 2], sep = "")
a <- system(cmd, intern = TRUE)

a <- execGRASS(cmd = "r.what", parameters = list(input = "ELEV_10",
                                            east_north = coordinates(cal_data)[2, ]),
               intern = TRUE)
a
a <- as.numeric(unlist(strsplit(a, "|", fixed = TRUE)))[-3]

input <- c("ELEV_10", "ELEV_90")

query <-
  function (input, east_north, intern = TRUE) {
    in_names <- input
    if (length(input) > 1) {
      input <- paste(input, collapse = ",")
    }
    res <- execGRASS(cmd = "r.what", 
                     parameters = list(input = input, east_north = east_north),
                     intern = TRUE)
    res <- as.numeric(unlist(strsplit(res, "|", fixed = TRUE)))[-3]
    res <- matrix(res, ncol = 2 + length(in_names))
    colnames(res) <- c("x", "y", in_names)
    return (res)
  }

query(c("ELEV_10", "ELEV_90"), coordinates(cal_data)[1, ])
