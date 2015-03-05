# DESCRIPTION ##################################################################
# Source code used to develop the first chapter (article) of my PhD research 
# project. The reference for the article is as follows:
# Samuel-Rosa, A.; Heuvelink, G. B. M.; Vasques, G. M. & Anjos, L. H. C. Do more
# detailed environmental covariates deliver more accurate soil maps?. Geoderma,
# v. 243-244, p. 214-227, 2015. doi:10.1016/j.geoderma.2014.12.017
# SETTINGS #####################################################################
rm(list = ls())
gc()
options(device = x11, stringsAsFactors = FALSE)
library(rgdal)
require(spgrass6)
require(raster)
require(gstat)
library(geoR)
library(pedometrics)
require(car)
require(caret)
require(MASS)
library(lattice)
require(xtable)
library(latticeExtra)
library(grid)
library(gridExtra)
library(Hmisc)
library(plotKML)
require(stringr)
require(plyr)
require(pbapply)
library(mail)
load("sm-dnos-general.RData")
load("sm-dnos-phd-chap1.RData")
load("sm-dnos-phd-chap1-final-models.RData")
data(R_pal)
data <- cal_data@data
source("/home/alessandro/PROJECTS/pedometrics/pedometrics/pkg/pedometrics/R/stepVIF.R")
source("/home/alessandro/PROJECTS/pedometrics/pedometrics/cooking/isAliased.R")
source("/home/alessandro/PROJECTS/pedometrics/pedometrics/cooking/stepAliased.R")
source("/home/alessandro/PROJECTS/pedometrics/pedometrics/cooking/whichAliased.R")
source("/home/alessandro/PROJECTS/pedometrics/pedometrics/cooking/invBoxCox.R")
source("/home/alessandro/PROJECTS/pedometrics/pedometrics/pkg/pedometrics/R/buildMS.R")
source("/home/alessandro/PROJECTS/pedometrics/pedometrics/pkg/pedometrics/R/statsMS.R")
source("/home/alessandro/PROJECTS/pedometrics/pedometrics/pkg/pedometrics/R/plotMS.R")
source("/home/alessandro/PROJECTS/pedometrics/pedometrics/cooking/plotHD.R")
source("/home/alessandro/PROJECTS/pedometrics/pedometrics/cooking/plotESDA.R")
source("/home/alessandro/PROJECTS/pedometrics/pedometrics/cooking/looCV.R")
source("/home/alessandro/PROJECTS/pedometrics/pedometrics/cooking/krigeCV.R")
source("/home/alessandro/PROJECTS/pedometrics/pedometrics/cooking/cvKrigeCDF.R")
source("/home/alessandro/PROJECTS/pedometrics/pedometrics/cooking/cvStats.R")
source("/home/alessandro/PROJECTS/pedometrics/pedometrics/cooking/spredict.R")
source("/home/alessandro/PROJECTS/pedometrics/pedometrics/cooking/linesREML.R")
initGRASS(gisBase = "/usr/lib/grass64/", gisDbase = GRASSgisDbase,
          location = "dnos-sm-rs", mapset = "predictions",
          pid = Sys.getpid(), override = TRUE)
system("g.region rast=dnos.raster")

# SETUP COVARIATE DATABASE #####################################################
# database 1
soil1 <- c("SOIL_100b", "SOIL_100c",
           "SOIL_100d", "SOIL_100e", "SOIL_100f")
soil1 <- paste(soil1, collapse = " + ")
land1 <- c("LU1980a", "LU1980b")
land1 <- paste(land1, collapse = " + ")
geo1 <- c("GEO_50a", "GEO_50b", "GEO_50c")
geo1 <- paste(geo1, collapse = " + ")
sat1 <- c("BLUE_30", "GREEN_30", "RED_30", "NIR_30a",
          "NIR_30b", "MIR_30", "NDVI_30", "SAVI_30")
sat1 <- paste(sat1, collapse = " + ")
dem1 <- c("ELEV_90", "SLP_90_3", "SLP_90_7", "SLP_90_15", "SLP_90_31", 
          "SLP_90_63", "SLP_90_127", "SLP_90_255", "TPI_90_3", "TPI_90_7",
          "TPI_90_15", "TPI_90_31", "TPI_90_63", "TPI_90_127", "TPI_90_255",
          "NOR_90_3", "NOR_90_7", "NOR_90_15", "NOR_90_31", "NOR_90_63",
          "NOR_90_127", "NOR_90_255", "TWI_90_3", "TWI_90_7", "TWI_90_15",
          "TWI_90_31", "TWI_90_63", "TWI_90_127", "TWI_90_255", "SPI_90_3",
          "SPI_90_7", "SPI_90_15", "SPI_90_31", "SPI_90_63", "SPI_90_127",
          "SPI_90_255")
dem1 <- paste(dem1, collapse = " + ")
# database 2
soil2 <- c("SOIL_25a", "SOIL_25b", "SOIL_25c", "SOIL_25d", 
           "SOIL_25h", "SOIL_25i", "SOIL_25j")
soil2 <- paste(soil2, collapse = " + ")
land2 <- c("LU2009a", "LU2009b", "LU2009c", "LU2009d", "LUdiff")
land2 <- paste(land2, collapse = " + ")
geo2 <- c("GEO_25a", "GEO_25b", "GEO_25c", "GEO_25d")
geo2 <- paste(geo2, collapse = " + ")
sat2 <- c("BLUE_5", "GREEN_5", "RED_5", "EDGE_5", "NIR_5", "NDVI_5a",
          "NDVI_5b", "SAVI_5a", "SAVI_5b")
sat2 <- paste(sat2, collapse = " + ")
dem2 <- c("ELEV_10", "SLP_10_3", "SLP_10_7", "SLP_10_15", "SLP_10_31", 
          "SLP_10_63", "SLP_10_127", "SLP_10_255", "TPI_10_3", "TPI_10_7",
          "TPI_10_15", "TPI_10_31", "TPI_10_63", "TPI_10_127", "TPI_10_255",
          "NOR_10_3", "NOR_10_7", "NOR_10_15", "NOR_10_31", "NOR_10_63", 
          "NOR_10_127", "NOR_10_255", "TWI_10_3", "TWI_10_7", "TWI_10_15", 
          "TWI_10_31", "TWI_10_63", "TWI_10_127", "TWI_10_255", "SPI_10_3", 
          "SPI_10_7", "SPI_10_15", "SPI_10_31", "SPI_10_63", "SPI_10_127", 
          "SPI_10_255")
dem2 <- paste(dem2, collapse = " + ")
# SETUP CANDIDATE MODELS #######################################################
# combs
combs <- list()
combs$main <- expand.grid(c("soil1", "soil2"), c("land1", "land2"),
                          c("geo1", "geo2"), c("sat1", "sat2"),
                          c("dem1", "dem2"), stringsAsFactors = FALSE)

combs$main <- split(combs$main, seq(1, nrow(combs$main), 1))
combs$num <- expand.grid(c(1, 2), c(1, 2), c(1, 2), c(1, 2), c(1, 2))
colnames(combs$num) <- c("soil", "land", "geo", "sat", "dem")
combs$base <- list()
for (i in 1:length(combs$main[[1]])) {
  combs$base[[i]] <- combs$main[[1]][-i]
}
combs$fine <- list()
for (i in 1:length(combs$main[[32]])) {
  combs$fine[[i]] <- combs$main[[32]][-i]
}
# predictors
preds <- list()
preds$main <- lapply(combs$main, function(X){parse(text = X)})
preds$main <- lapply(preds$main, function(X){lapply(X, eval)})
preds$main <- lapply(preds$main, function(X) {
  res = paste(unlist(X), collapse = " + ")
  return(res)
  })
preds$base <- lapply(combs$base, function(X){parse(text = X)})
preds$base <- lapply(preds$base, function(X){lapply(X, eval)})
preds$base <- lapply(preds$base, function(X) {
  res = paste(unlist(X), collapse = " + ")
  return(res)
})
preds$fine <- lapply(combs$fine, function(X){parse(text = X)})
preds$fine <- lapply(preds$fine, function(X){lapply(X, eval)})
preds$fine <- lapply(preds$fine, function(X) {
  res = paste(unlist(X), collapse = " + ")
  return(res)
})
# formulas
forms <- list()
forms$soil_attrs <- c("clay_bc", "carbon_bc", "ecec_bc")
forms$main <- lapply(forms$soil_attrs, function(X){paste(X, " ~ ", preds$main)})
forms$clay <-
  lapply(unlist(forms$main[which(forms$soil_attrs == "clay_bc")]), as.formula)
forms$soc <- 
  lapply(unlist(forms$main[which(forms$soil_attrs == "carbon_bc")]), as.formula)
forms$ecec <- 
  lapply(unlist(forms$main[which(forms$soil_attrs == "ecec_bc")]), as.formula)
forms$base <- lapply(forms$soil_attrs, function(X){paste(X, " ~ ", preds$base)})
forms$fine <- lapply(forms$soil_attrs, function(X){paste(X, " ~ ", preds$fine)})
forms$clay_base <-
  lapply(unlist(forms$base[which(forms$soil_attrs == "clay_bc")]), as.formula)
forms$clay_fine <-
  lapply(unlist(forms$fine[which(forms$soil_attrs == "clay_bc")]), as.formula)
forms$soc_base <-
  lapply(unlist(forms$base[which(forms$soil_attrs == "carbon_bc")]), as.formula)
forms$soc_fine <-
  lapply(unlist(forms$fine[which(forms$soil_attrs == "carbon_bc")]), as.formula)
forms$ecec_base <-
  lapply(unlist(forms$base[which(forms$soil_attrs == "ecec_bc")]), as.formula)
forms$ecec_fine <-
  lapply(unlist(forms$fine[which(forms$soil_attrs == "ecec_bc")]), as.formula)

# LINEAR MODEL FITTING #########################################################
# CLAY -------------------------------------------------------------------------
# build model series using several strategies
clay_full     <- buildMS(forms$clay, data)
clay_vif      <- buildMS(forms$clay, data, vif = TRUE)
clay_both     <- buildMS(forms$clay, data, vif = TRUE, aic = TRUE)
clay_forward  <- buildMS(forms$clay, data, vif = TRUE, aic = TRUE,
                         aic.direction = "forward")
clay_backward <- buildMS(forms$clay, data, vif = TRUE, aic = TRUE, 
                         aic.direction = "backward")
# get statistics of model series
clay_full_stats     <- statsMS(clay_full, combs$num, "rmse")
clay_vif_stats      <- statsMS(clay_vif, combs$num, "rmse")
clay_both_stats     <- statsMS(clay_both, combs$num, "rmse")
clay_forward_stats  <- statsMS(clay_forward, combs$num, "rmse")
clay_backward_stats <- statsMS(clay_backward, combs$num, "rmse")
# plot and save all model series
grid <- c(2:6)
line <- "ADJ_r2"
ind  <- 2
color <- c("lightyellow", "palegreen")
a_plot <- plotMS(clay_full_stats, grid, line, ind, color = color, 
                 main = "full model")
b_plot <- plotMS(clay_vif_stats, grid, line, ind, color = color, 
                 main = "VIF selection")
c_plot <- plotMS(clay_forward_stats, grid, line, ind, color = color, 
                 main = "forward selection")
d_plot <- plotMS(clay_backward_stats, grid, line, ind, color = color, 
                 main = "backward selection")
e_plot <- plotMS(clay_both_stats, grid, line, ind, color = color,
                 main = "stepwise selection")
dev.off()
pdf(file = paste(update_covars_dir, "clay_model_series_all.pdf", sep = ""),
    width = 7, height = 15)
trellis.par.set(fontsize = list(text = 8, points = 6))
grid.arrange(a_plot, b_plot, c_plot, d_plot, e_plot, ncol = 1)
dev.off()
rm(grid, line, ind, color, a_plot, b_plot, c_plot, d_plot, e_plot)
# MODEL SERIES PLOT - stepwise variable selection
dev.off()
pdf(file = paste(update_covars_dir, "clay_models.pdf", sep = ""),
    width = 19/cm(1), height = 8/cm(1))
trellis.par.set(fontsize = list(text = 7, points = 5),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
plotMS(clay_both_stats, grid = c(2:6), line = "ADJ_r2", ind = 2, 
       color = c("lightyellow", "palegreen"),
       xlab = "CLAY model ranking", scales = list(cex = c(1, 1)))
dev.off()
# check the effect of the number of observations (200, 300)
data <- cal_data@data
data <- data[sample(c(1:350), size = 200), ]
tmp <- buildMS(forms$clay, data)
tmp <- statsMS(tmp, combs$num, "rmse")
grid <- c(2:6)
line <- "ADJ_r2"
ind  <- 2
color <- c("lightyellow", "palegreen")
plotMS(tmp, grid, line, ind, color = color)
# get base and best models
clay_sel <- list()
clay_sel$poor_lm <- clay_both[head(clay_both_stats, 1)$id][[1]]
clay_sel$base_lm <- clay_both[[1]]
clay_sel$fine_lm <- clay_both[[32]]
clay_sel$best_lm <- clay_both[tail(clay_both_stats, 1)$id][[1]]

clay_base_lm <- clay_both[[1]]
clay_best_lm <- rev(clay_both[tail(clay_both_stats, 1)$id])[[1]]

# CLAY - analysis of the residuals ---------------------------------------------
# BASE MODEL
data   <- cal_data
model  <- clay_base_lm
res    <- residuals(model)
lambda <- bc_lambda$clay
# residual plots
dev.off()
pdf(file = paste(update_covars_dir, "clay_base_lm_res.pdf", sep = ""),
    width = 7, height = 11)
par(mfrow = c(3, 2))
plot(model, which = c(1:6))
dev.off()
# exploratory spatial data analysis
dev.off()
pdf(file = paste(update_covars_dir, "clay_base_lm_esda.pdf", sep = ""))
plotESDA(res, lon = coordinates(data)[, 1], lat = coordinates(data)[, 2])
dev.off()

# BEST MODEL
data  <- cal_data
model <- clay_best_lm
res   <- residuals(model)
lambda <- bc_lambda$clay
# residual plots
dev.off()
pdf(file = paste(update_covars_dir, "clay_best_lm_res.pdf", sep = ""),
    width = 7, height = 11)
par(mfrow = c(3, 2))
plot(model, which = c(1:6))
dev.off()
# exploratory spatial data analysis
dev.off()
pdf(file = paste(update_covars_dir, "clay_best_lm_esda.pdf", sep = ""))
plotESDA(res, coordinates(data)[, 1], coordinates(data)[, 2],
         cutoff = 1000, width = 1000/10)
dev.off()

# CARBON -----------------------------------------------------------------------
data <- cal_data@data
# fit using several strategies
carbon_full     <- buildMS(forms$soc, data)
carbon_vif      <- buildMS(forms$soc, data, vif = TRUE)
carbon_both     <- buildMS(forms$soc, data, vif = TRUE, 
                                    aic = TRUE, aic.direction = "both")
carbon_forward  <- buildMS(forms$soc, data, vif = TRUE, 
                                    aic = TRUE, aic.direction = "forward")
carbon_backward <- buildMS(forms$soc, data, vif = TRUE, 
                                    aic = TRUE, aic.direction = "backward")
# get statistics of model series
carbon_full_stats     <- statsMS(carbon_full, combs$num, "rmse")
carbon_vif_stats      <- statsMS(carbon_vif, combs$num, "rmse")
carbon_both_stats     <- statsMS(carbon_both, combs$num, "rmse")
carbon_forward_stats  <- statsMS(carbon_forward, combs$num, "rmse")
carbon_backward_stats <- statsMS(carbon_backward, combs$num, "rmse")

# plot and save all model series
grid <- c(2:6)
line <- "ADJ_r2"
ind  <- 2
color <- c("lightyellow", "palegreen")
a_plot <- plotMS(carbon_full_stats, grid, line, ind, color = color,
                 main = "full model")
b_plot <- plotMS(carbon_vif_stats, grid, line, ind, color = color, main = "VIF selection")
c_plot <- plotMS(carbon_forward_stats, grid, line, ind, color = color, main = "forward selection")
d_plot <- plotMS(carbon_backward_stats, grid, line, ind, color = color, main = "backward selection")
e_plot <- plotMS(carbon_both_stats, grid, line, ind, color = color, main = "stepwise selection")
dev.off()
pdf(file = paste(update_covars_dir, "carbon_model_series_all.pdf", sep = ""),
    width = 7, height = 15)
trellis.par.set(fontsize = list(text = 8, points = 6))
grid.arrange(a_plot, b_plot, c_plot, d_plot, e_plot, ncol = 1)
dev.off()
rm(grid, line, ind, color, a_plot, b_plot, c_plot, d_plot, e_plot)
gc()
# MODEL SERIES PLOT - stepwise variable selection
dev.off()
pdf(file = paste(update_covars_dir, "carbon_models.pdf", sep = ""),
    width = 19/cm(1), height = 8/cm(1))
trellis.par.set(fontsize = list(text = 7, points = 5),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
plotMS(carbon_both_stats, grid = c(2:6), line = "ADJ_r2", ind = 2, 
       color = c("lightyellow", "palegreen"),
       xlab = "SOC model ranking", scales = list(cex = c(1, 1)))
dev.off()
# check the effect of the number of observations (200, 300)
data <- cal_data@data
data <- data[sample(c(1:350), size = 300), ]
tmp <- buildMS(forms$soc, data)
tmp <- statsMS(tmp, combs$num, "rmse")
grid <- c(2:6)
line <- "ADJ_r2"
ind  <- 2
color <- c("lightyellow", "palegreen")
plotMS(tmp, grid, line, ind, color = color)
# get base and best models
soc_sel <- list()
soc_sel$poor_lm <- carbon_both[head(carbon_both_stats, 1)$id][[1]]
soc_sel$base_lm <- carbon_both[[1]]
soc_sel$fine_lm <- carbon_both[[32]]
soc_sel$best_lm <- carbon_both[tail(carbon_both_stats, 1)$id][[1]]

carbon_base_lm <- carbon_both[[1]]
carbon_best_lm <- rev(carbon_both[tail(carbon_both_stats, 1)$id])[[1]]

# CARBON - analysis of the residuals -------------------------------------------

# BASE MODEL
data  <- cal_data
model <- carbon_base_lm
res   <- residuals(model)
lambda <- bc_lambda$carbon
# residual plots
dev.off()
pdf(file = paste(update_covars_dir, "carbon_base_lm_res.pdf", sep = ""),
    width = 7, height = 11)
par(mfrow = c(3, 2))
plot(model, which = c(1:6))
dev.off()
# exploratory spatial data analysis
dev.off()
pdf(file = paste(update_covars_dir, "carbon_base_lm_esda.pdf", sep = ""))
plotESDA(res, coordinates(data)[, 1], coordinates(data)[, 2],
         cutoff = 1000, width = 1000/10)
dev.off()

# BEST MODEL
data  <- cal_data
model <- carbon_best_lm
res   <- residuals(model)
lambda <- bc_lambda$carbon
# residual plots
dev.off()
pdf(file = paste(update_covars_dir, "carbon_best_lm_res.pdf", sep = ""),
    width = 7, height = 11)
par(mfrow = c(3, 2))
plot(model, which = c(1:6))
dev.off()
# exploratory spatial data analysis
dev.off()
pdf(file = paste(update_covars_dir, "carbon_best_lm_esda.pdf", sep = ""))
plotESDA(res, coordinates(data)[, 1], coordinates(data)[, 2],
         cutoff = 1000, width = 1000/10)
dev.off()

# ECEC -------------------------------------------------------------------------
data <- cal_data@data
# fit using several strategies
ecec_full     <- buildMS(forms$ecec, data)
ecec_vif      <- buildMS(forms$ecec, data, vif = TRUE)
ecec_both     <- buildMS(forms$ecec, data, vif = TRUE, aic = TRUE, 
                                  aic.direction = "both")
ecec_forward  <- buildMS(forms$ecec, data, vif = TRUE, aic = TRUE, 
                                  aic.direction = "forward")
ecec_backward <- buildMS(forms$ecec, data, vif = TRUE, aic = TRUE, 
                                  aic.direction = "backward")
# get statistics of model series
ecec_full_stats <- statsMS(ecec_full, combs$num, "rmse")
ecec_vif_stats <- statsMS(ecec_vif, combs$num, "rmse")
ecec_both_stats <- statsMS(ecec_both, combs$num, "rmse")
ecec_forward_stats <- statsMS(ecec_forward, combs$num, "rmse")
ecec_backward_stats <- statsMS(ecec_backward, combs$num, "rmse")

# plot and save all model series
grid <- c(2:6)
line <- "ADJ_r2"
ind  <- 2
color <- c("lightyellow", "palegreen")
a_plot <- plotMS(ecec_full_stats, grid, line, ind, color = color, main = "full model")
b_plot <- plotMS(ecec_vif_stats, grid, line, ind, color = color, main = "VIF selection")
c_plot <- plotMS(ecec_forward_stats, grid, line, ind, color = color, main = "forward selection")
d_plot <- plotMS(ecec_backward_stats, grid, line, ind, color = color, main = "backward selection")
e_plot <- plotMS(ecec_both_stats, grid, line, ind, color = color, main = "stepwise selection")
dev.off()
pdf(file = paste(update_covars_dir, "ecec_model_series_all.pdf", sep = ""),
    width = 7, height = 15)
trellis.par.set(fontsize = list(text = 8, points = 6))
grid.arrange(a_plot, b_plot, c_plot, d_plot, e_plot, ncol = 1)
dev.off()
rm(grid, line, ind, color, a_plot, b_plot, c_plot, d_plot, e_plot)
gc()
# MODEL SERIES PLOT - stepwise variable selection
dev.off()
pdf(file = paste(update_covars_dir, "ecec_models.pdf", sep = ""), 
    width = 19/cm(1), height = 8/cm(1))
trellis.par.set(fontsize = list(text = 7, points = 5),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
plotMS(ecec_both_stats, grid = c(2:6), line = "ADJ_r2", ind = 2, 
       color = c("lightyellow", "palegreen"),
       xlab = "ECEC model ranking", scales = list(cex = c(1, 1)))
dev.off()
# check the effect of the number of observations (200, 300)
data <- cal_data@data
data <- data[sample(c(1:350), size = 200), ]
tmp <- buildMS(forms$clay, data)
tmp <- statsMS(tmp, combs$num, "rmse")
grid <- c(2:6)
line <- "ADJ_r2"
ind  <- 2
color <- c("lightyellow", "palegreen")
plotMS(tmp, grid, line, ind, color = color)
# get base and best models
ecec_sel <- list()
ecec_sel$poor_lm <- ecec_both[head(ecec_both_stats, 1)$id][[1]]
ecec_sel$base_lm <- ecec_both[[1]]
ecec_sel$fine_lm <- ecec_both[[32]]
ecec_sel$best_lm <- ecec_both[tail(ecec_both_stats, 1)$id][[1]]

ecec_base_lm <- ecec_both[[1]]
ecec_best_lm <- rev(ecec_both[tail(ecec_both_stats, 1)$id])[[1]]

# ECEC - analysis of the residuals ---------------------------------------------

# BASE MODEL
data   <- cal_data
model  <- ecec_base_lm
res    <- residuals(model)
lambda <- bc_lambda$ecec
# residual plots
dev.off()
pdf(file = paste(update_covars_dir, "ecec_base_lm_res.pdf", sep = ""),
    width = 7, height = 11)
par(mfrow = c(3, 2))
plot(model, which = c(1:6))
dev.off()
# exploratory spatial data analysis
dev.off()
pdf(file = paste(update_covars_dir, "ecec_base_lm_esda.pdf", sep = ""))
plotESDA(res, coordinates(data)[, 1], coordinates(data)[, 2],
         cutoff = 1000, width = 1000/10)
dev.off()

# BEST MODEL
data  <- cal_data
model <- ecec_best_lm
res   <- residuals(model)
lambda <- bc_lambda$ecec
# residual plots
dev.off()
pdf(file = paste(update_covars_dir, "ecec_best_lm_res.pdf", sep = ""),
    width = 7, height = 11)
par(mfrow = c(3, 2))
plot(model, which = c(1:6))
dev.off()
# exploratory spatial data analysis
dev.off()
pdf(file = paste(update_covars_dir, "ecec_best_lm_esda.pdf", sep = ""))
plotESDA(res, coordinates(data)[, 1], coordinates(data)[, 2],
         cutoff = 1000, width = 1000/10)
dev.off()

# SENSITIVITY ANALYSIS #########################################################

# Effect of dropping one environmental covariate -------------------------------
deltaR2 <- 
  function (a, b, cols = c("r2", "adj_r2", "ADJ_r2"), 
            rows = c("soil", "land", "geo", "sat", "dem")) {
    a <- a[, c("r2", "adj_r2", "ADJ_r2")]
    b <- b[, c("r2", "adj_r2", "ADJ_r2")]
    ab <- list()
    for (i in 1:3){
      ab[[i]] <- as.numeric(a)[i] - as.numeric(b[, i])
    }
    ab <- as.data.frame(ab)
    colnames(ab) <- c("r2", "adj_r2", "ADJ_r2")
    rownames(ab) <- rows
    return (ab)
  }
drop <- list()
# clay
drop$clay$base_lm <- buildMS(forms$clay_base, data, vif = TRUE, aic = TRUE)
drop$clay$fine_lm <- buildMS(forms$clay_fine, data, vif = TRUE, aic = TRUE)
drop$clay$base_r2 <- statsMS(drop$clay$base_lm)
drop$clay$fine_r2 <- statsMS(drop$clay$fine_lm)
drop$clay$base_dr2 <- deltaR2(clay_both_stats[clay_both_stats$id == 1, ],
                              drop$clay$base_r2)
drop$clay$fine_dr2 <- deltaR2(clay_both_stats[clay_both_stats$id == 32, ],
                              drop$clay$fine_r2)
# soc
drop$soc$base_lm <- buildMS(forms$soc_base, data, vif = TRUE, aic = TRUE)
drop$soc$fine_lm <- buildMS(forms$soc_fine, data, vif = TRUE, aic = TRUE)
drop$soc$base_r2 <- statsMS(drop$soc$base_lm)
drop$soc$fine_r2 <- statsMS(drop$soc$fine_lm)
drop$soc$base_dr2 <- deltaR2(carbon_both_stats[carbon_both_stats$id == 1, ],
                             drop$soc$base_r2)
drop$soc$fine_dr2 <- deltaR2(carbon_both_stats[carbon_both_stats$id == 32, ],
                             drop$soc$fine_r2)
# ecec
drop$ecec$base_lm <- buildMS(forms$ecec_base, data, vif = TRUE, aic = TRUE)
drop$ecec$fine_lm <- buildMS(forms$ecec_fine, data, vif = TRUE, aic = TRUE)
drop$ecec$base_r2 <- statsMS(drop$ecec$base_lm)
drop$ecec$fine_r2 <- statsMS(drop$ecec$fine_lm)
drop$ecec$base_dr2 <- deltaR2(ecec_both_stats[ecec_both_stats$id == 1, ],
                              drop$ecec$base_r2)
drop$ecec$fine_dr2 <- deltaR2(ecec_both_stats[ecec_both_stats$id == 32, ],
                              drop$ecec$fine_r2)
# table
EnvCov <- data.frame(drop$clay$base_dr2$ADJ_r2, drop$clay$fine_dr2$ADJ_r2,
                     drop$soc$base_dr2$ADJ_r2, drop$soc$fine_dr2$ADJ_r2,
                     drop$ecec$base_dr2$ADJ_r2, drop$ecec$fine_dr2$ADJ_r2)
EnvCov <- round(EnvCov, 3)
rownames(EnvCov) <- c("\\texttt{soil}", "\\texttt{land}",
                      "\\texttt{geo}", "\\texttt{sat}", "\\texttt{dem}")
colnames(EnvCov) <- rep(c("less", "more"), 3)
long_cap <- "Variation of the adjusted R$^2$ ($\\Delta$R$^2_{adj}$) when dropping one environmental covariate (EnvCov) in the models built using only the less accurate or the more accurate version of all environmental covariates."
foot <- "* Environmental covariates (EnvCov): \\texttt{soil} - soil map, \\texttt{land} - land use map, \\texttt{geo} - geological map, \\texttt{sat} - satellite image, and \\texttt{dem} - digital elevation model. ** $\\Delta$R$^2_{adj} = R$^2_{adj}_{p=5} - R$^2_{adj}_{p=5-1}$."
file <-  paste(update_covars_dir, "drop-covars.tex", sep = "")
latex(EnvCov, file = file, label = "tab:drop", table.env = TRUE, 
      longtable = FALSE, cgroup = c("CLAY","SOC", "ECEC"), n.cgroup = c(2, 2, 2), 
      na.blank = TRUE, ctable = TRUE, caption = long_cap, where = NULL,
      size = "scriptsize", insert.bottom = foot, cgroupTexCmd = NULL,
      rgroupTexCmd = NULL)
rm(EnvCov, long_cap, foot, file)

# LINEAR MIXED MODEL FITTING ###################################################
data <- as.data.frame(cal_data)

# # CLAY - geodata ---------------------------------------------------------------
# clay_geodata <- list()
# data.col  <- which(colnames(data) == "clay")
# clay_geodata$clay <- as.geodata(data, coords.col = c(1, 2), data.col = data.col)
# # poor model
# covars <- attr(terms(clay_sel$poor_lm), "term.labels")
# covar.col <- which(match(colnames(data), covars) != "NA")
# clay_geodata$poor <- data[, covar.col]
# # base model
# covars <- attr(terms(clay_sel$base_lm), "term.labels")
# covar.col <- which(match(colnames(data), covars) != "NA")
# clay_geodata$base <- data[, covar.col]
# # fine model
# covars <- attr(terms(clay_sel$fine_lm), "term.labels")
# covar.col <- which(match(colnames(data), covars) != "NA")
# clay_geodata$fine <- data[, covar.col]
# # best model
# covars <- attr(terms(clay_sel$best_lm), "term.labels")
# covar.col <- which(match(colnames(data), covars) != "NA")
# clay_geodata$best <- data[, covar.col]

# CLAY - poor model ------------------------------------------------------------
# create geodata object
model     <- clay_sel$poor_lm
data      <- as.data.frame(cal_data)
covars    <- attr(terms(model), "term.labels")
covar.col <- which(match(colnames(data), covars) != "NA")
data.col  <- "clay"
data.col  <- which(colnames(data) == data.col)
clay_poor_geodata <- as.geodata(data, coords.col = c(1, 2), data.col = data.col,
                                covar.col = covar.col)
summary(clay_poor_geodata)
rm(model, data, covar.col, data.col)
# calculate empirical variogram
geodata         <- clay_poor_geodata
trend           <- formula(paste("~", paste(covars, collapse = " + ")))
lambda          <- bc_lambda$clay
breaks          <- seq(0, 6500, 100)
clay_poor_vario <- variog(geodata, trend = trend, lambda = lambda,
                          breaks = breaks)
plot(clay_poor_vario, col = "blue", pch = 20, type = "b")
# estimate model parameters using REML
ini.cov.pars  <- as.matrix(expand.grid(c(0.9, 1.0, 1.1), c(400, 500, 600)))
nugget        <- c(0.3, 0.4, 0.5)
clay_poor_lmm <- likfit(geodata, trend = trend, ini.cov.pars = ini.cov.pars,
                        nugget = nugget, lambda = lambda, lik.method = "REML")
summary(clay_poor_lmm)
rm(geodata, trend, lambda, nugget, covars, ini.cov.pars)
gc()

# CLAY - base model ------------------------------------------------------------
# create geodata object
model     <- clay_base_lm
data      <- as.data.frame(cal_data)
covars    <- attr(terms(model), "term.labels")
covar.col <- which(match(colnames(data), covars) != "NA")
data.col  <- "clay"
data.col  <- which(colnames(data) == data.col)
clay_base_geodata <- as.geodata(data, coords.col = c(1, 2), data.col = data.col,
                                covar.col = covar.col)
summary(clay_base_geodata)
rm(model, data, covar.col, data.col)
# calculate empirical variogram
geodata         <- clay_base_geodata
trend           <- formula(paste("~", paste(covars, collapse = " + ")))
lambda          <- bc_lambda$clay
breaks          <- seq(0, 6500, 100)
clay_base_vario <- variog(geodata, trend = trend, lambda = lambda,
                          breaks = breaks)
plot(clay_base_vario, col = "blue", pch = 20, type = "b")
# estimate model parameters using REML
ini.cov.pars  <- as.matrix(expand.grid(c(0.9, 1.0, 1.1), c(400, 500, 600)))
nugget        <- c(0.3, 0.4, 0.5)
clay_base_lmm <- likfit(geodata, trend = trend, ini.cov.pars = ini.cov.pars,
                        nugget = nugget, lambda = lambda, lik.method = "REML")
summary(clay_base_lmm)
rm(geodata, trend, lambda, nugget, covars, ini.cov.pars)
gc()

# CLAY - fine model ------------------------------------------------------------
# create geodata object
model     <- clay_sel$fine_lm
data      <- as.data.frame(cal_data)
covars    <- attr(terms(model), "term.labels")
covar.col <- which(match(colnames(data), covars) != "NA")
data.col  <- "clay"
data.col  <- which(colnames(data) == data.col)
clay_fine_geodata <- as.geodata(data, coords.col = c(1, 2), data.col = data.col,
                                covar.col = covar.col)
summary(clay_fine_geodata)
rm(model, data, covar.col, data.col)
# calculate empirical variogram
geodata         <- clay_fine_geodata
trend           <- formula(paste("~", paste(covars, collapse = " + ")))
lambda          <- bc_lambda$clay
breaks          <- seq(0, 6500, 100)
clay_fine_vario <- variog(geodata, trend = trend, lambda = lambda, 
                          breaks = breaks)
plot(clay_fine_vario, col = "blue", pch = 20, type = "b")
# estimate model parameters using REML
ini.cov.pars  <- as.matrix(expand.grid(c(1.0, 1.1, 1.2), c(300, 400, 500)))
nugget        <- c(0.1, 0.2, 0.3)
clay_fine_lmm <- likfit(geodata, trend = trend, ini.cov.pars = ini.cov.pars,
                        nugget = nugget, lambda = lambda, lik.method = "REML")
summary(clay_fine_lmm)
rm(geodata, trend, covars, lambda, nugget, ini.cov.pars)
gc()

# CLAY - best model ------------------------------------------------------------
# create geodata object
model     <- clay_best_lm
data      <- as.data.frame(cal_data)
covars    <- attr(terms(model), "term.labels")
covar.col <- which(match(colnames(data), covars) != "NA")
data.col  <- "clay"
data.col  <- which(colnames(data) == data.col)
clay_best_geodata <- as.geodata(data, coords.col = c(1, 2), data.col = data.col,
                                covar.col = covar.col)
summary(clay_best_geodata)
rm(model, data, covar.col, data.col)
# calculate empirical variogram
geodata         <- clay_best_geodata
trend           <- formula(paste("~", paste(covars, collapse = " + ")))
lambda          <- bc_lambda$clay
breaks          <- seq(0, 6500, 100)
clay_best_vario <- variog(geodata, trend = trend, lambda = lambda, 
                          breaks = breaks)
plot(clay_best_vario, col = "blue", pch = 20, type = "b")
# estimate model parameters using REML
ini.cov.pars  <- as.matrix(expand.grid(c(1.0, 1.1, 1.2), c(300, 400, 500)))
nugget        <- c(0.1, 0.2, 0.3)
clay_best_lmm <- likfit(geodata, trend = trend, ini.cov.pars = ini.cov.pars,
                        nugget = nugget, lambda = lambda, lik.method = "REML")
summary(clay_best_lmm)
rm(geodata, trend, covars, lambda, nugget, ini.cov.pars)
gc()

# CLAY - plot experimental variograms and fitted linear mixed models -----------
xlim <- c(0, 3000)
ylim <- max(#clay_poor_vario$v[clay_poor_vario$u <= max(xlim)],
            clay_base_vario$v[clay_base_vario$u <= max(xlim)],
            #clay_fine_vario$v[clay_fine_vario$u <= max(xlim)],
            clay_best_vario$v[clay_best_vario$u <= max(xlim)]) * 1.1
ylim <- c(0, round(ylim, 1))
#l1 <- linesREML(clay_poor_lmm, add = FALSE)
l2 <- linesREML(clay_base_lmm, add = FALSE)
#l3 <- linesREML(clay_fine_lmm, add = FALSE)
l4 <- linesREML(clay_best_lmm, add = FALSE)
# v1 <- xyplot(clay_poor_vario$v ~ clay_poor_vario$u, ylim = ylim, pch = 20,
#              col =  "black", scales = list(tick.number = 8), xlim = xlim,
#              panel = function(x, y, ...) {
#                panel.xyplot(x, y, ...)
#                panel.lines(x = l1$x, y = l1$y, lty = 2, col = "black")
#              })
clay_v2 <- xyplot(clay_base_vario$v ~ clay_base_vario$u, ylim = ylim, pch = 20,
             scales = list(tick.number = 8), xlim = xlim, col =  "black",
             panel = function(x, y, ...) {
               panel.xyplot(x, y, ...)
               panel.lines(x = l2$x, y = l2$y, lty = 2, col = "black")
             })
# v3 <- xyplot(clay_fine_vario$v ~ clay_fine_vario$u, ylim = ylim, pch = 20,
#              col =  "black", scales = list(tick.number = 8), xlim = xlim,
#              panel = function(x, y, ...) {
#                panel.xyplot(x, y, ...)
#                panel.lines(x = l3$x, y = l3$y, lty = 2, col = "black")
#              })
clay_v4 <- xyplot(clay_best_vario$v ~ clay_best_vario$u, ylim = ylim, pch = 20,
             scales = list(tick.number = 8), xlim = xlim, col =  "black", 
             panel = function(x, y, ...) {
               panel.xyplot(x, y, ...)
               panel.lines(x = l4$x, y = l4$y, lty = 2, col = "black")
             })
ylab <- expression(paste("Semivariance [(g kg"^"-1", ")"^"2", "]", sep = ""))
# v1 <- update(c(v1, v2, v3, v4), ylab = ylab, xlab = "Distance [m]", asp = 1,
#             scales = list(cex = c(0.9, 0.9)), layout = c(4, 1))
clay_v <- update(c(clay_v2, clay_v4), ylab = ylab, xlab = "Distance [m]", 
                 asp = 1, scales = list(cex = c(1, 1)), layout = c(2, 1))
# save plot
dev.off()
pdf(file = paste(update_covars_dir, "clay_lmm.pdf", sep = ""), 
    width = 9/cm(1), height = 6/cm(1))
trellis.par.set(fontsize = list(text = 7, points = 5),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
plot(clay_v)
dev.off()
rm(xlim, ylim, l1, l2, l3, l4, clay_v1, clay_v2, clay_v3, clay_v4, clay_v)
gc()

# CARBON - poor model ----------------------------------------------------------
# create geodata object
model     <- soc_sel$poor_lm
data      <- as.data.frame(cal_data)
covars    <- attr(terms(model), "term.labels")
covar.col <- which(match(colnames(data), covars) != "NA")
data.col  <- "carbon"
data.col  <- which(colnames(data) == data.col)
carbon_poor_geodata <- as.geodata(data, coords.col = c(1, 2), 
                                  data.col = data.col, covar.col = covar.col)
summary(carbon_poor_geodata)
rm(data, model, covar.col, data.col)
# calculate empirical variogram
geodata           <- carbon_poor_geodata
trend             <- formula(paste("~", paste(covars, collapse = " + ")))
lambda            <- bc_lambda$carbon
breaks            <- seq(0, 6500, 100)
carbon_poor_vario <- variog(geodata, trend = trend, lambda = lambda, 
                            breaks = breaks)
plot(carbon_poor_vario, col = "blue", pch = 20, type = "b")
# estimate model parameters using REML
ini.cov.pars     <- as.matrix(expand.grid(c(0.12, 0.14, 0.16), c(500, 600, 700)))
nugget           <- c(0.10, 0.12, 0.14)
carbon_poor_lmm  <- likfit(geodata, trend = trend, ini.cov.pars = ini.cov.pars,
                           nugget = nugget, lambda = lambda, lik.method = "REML")
summary(carbon_poor_lmm)
rm(covars, geodata, trend, breaks, lambda, nugget, ini.cov.pars)
gc()

# CARBON - base model ----------------------------------------------------------
# create geodata object
model     <- carbon_base_lm
data      <- as.data.frame(cal_data)
covars    <- attr(terms(model), "term.labels")
covar.col <- which(match(colnames(data), covars) != "NA")
data.col  <- "carbon"
data.col  <- which(colnames(data) == data.col)
carbon_base_geodata <- as.geodata(data, coords.col = c(1, 2), 
                                  data.col = data.col, covar.col = covar.col)
summary(carbon_base_geodata)
rm(data, model, covar.col, data.col)
# calculate empirical variogram
geodata           <- carbon_base_geodata
trend             <- formula(paste("~", paste(covars, collapse = " + ")))
lambda            <- bc_lambda$carbon
breaks            <- seq(0, 6500, 100)
carbon_base_vario <- variog(geodata, trend = trend, lambda = lambda, 
                            breaks = breaks)
plot(carbon_base_vario, col = "blue", pch = 20, type = "b")
# estimate model parameters using REML
ini.cov.pars     <- as.matrix(expand.grid(c(0.12, 0.14, 0.16), c(500, 600, 700)))
nugget           <- c(0.10, 0.12, 0.14)
carbon_base_lmm  <- likfit(geodata, trend = trend, ini.cov.pars = ini.cov.pars,
                           nugget = nugget, lambda = lambda, lik.method = "REML")
summary(carbon_base_lmm)
rm(covars, geodata, trend, breaks, lambda, nugget, ini.cov.pars, breaks)
gc()
# CARBON - fine model ----------------------------------------------------------
# create geodata object
model     <- soc_sel$fine_lm
data      <- as.data.frame(cal_data)
covars    <- attr(terms(model), "term.labels")
covar.col <- which(match(colnames(data), covars) != "NA")
data.col  <- "carbon"
data.col  <- which(colnames(data) == data.col)
carbon_fine_geodata <- as.geodata(data, coords.col = c(1, 2), 
                                  data.col = data.col, covar.col = covar.col)
summary(carbon_fine_geodata)
rm(data, model, covar.col, data.col)
# calculate empirical variogram
geodata           <- carbon_fine_geodata
trend             <- formula(paste("~", paste(covars, collapse = " + ")))
lambda            <- bc_lambda$carbon
breaks            <- seq(0, 6500, 100)
carbon_fine_vario <- variog(geodata, trend = trend, lambda = lambda, 
                            breaks = breaks)
plot(carbon_fine_vario, col = "blue", pch = 20, type = "b")
# estimate model parameters using REML
ini.cov.pars     <- as.matrix(expand.grid(c(0.12, 0.14, 0.16), c(500, 600, 700)))
nugget           <- c(0.10, 0.12, 0.14)
carbon_fine_lmm  <- likfit(geodata, trend = trend, ini.cov.pars = ini.cov.pars,
                           nugget = nugget, lambda = lambda, lik.method = "REML")
summary(carbon_fine_lmm)
rm(covars, geodata, trend, breaks, lambda, nugget, ini.cov.pars)
gc()
# CARBON - best model ----------------------------------------------------------
# create geodata object
model     <- carbon_best_lm
data      <- as.data.frame(cal_data)
covars    <- attr(terms(model), "term.labels")
covar.col <- which(match(colnames(data), covars) != "NA")
data.col  <- "carbon"
data.col  <- which(colnames(data) == data.col)
carbon_best_geodata <- as.geodata(data, coords.col = c(1, 2), 
                                  data.col = data.col, covar.col = covar.col)
summary(carbon_best_geodata)
rm(data, model, covar.col, data.col)
# calculate empirical variogram
geodata           <- carbon_best_geodata
trend             <- formula(paste("~", paste(covars, collapse = " + ")))
lambda            <- bc_lambda$carbon
breaks            <- seq(0, 6500, 100)
carbon_best_vario <- variog(geodata, trend = trend, lambda = lambda,
                            breaks = breaks)
plot(carbon_best_vario, col = "blue", pch = 20, type = "b")
# estimate model parameters using REML
ini.cov.pars    <- as.matrix(expand.grid(c(0.10, 0.12, 0.14), c(400, 500, 600)))
nugget          <- c(0.10, 0.12, 0.14)
carbon_best_lmm <- likfit(geodata, trend = trend, ini.cov.pars = ini.cov.pars,
                          nugget = nugget, lambda = lambda, lik.method = "REML")
summary(carbon_best_lmm)
# clean workspace
rm(covars, geodata, trend, breaks, lambda, nugget, ini.cov.pars)
gc()

# CARBON - plot experimental variograms and fitted linear mixed models ---------
xlim <- c(0, 3000)
ylim <- max(#carbon_poor_vario$v[carbon_poor_vario$u <= max(xlim)],
            carbon_base_vario$v[carbon_base_vario$u <= max(xlim)],
            #carbon_fine_vario$v[carbon_fine_vario$u <= max(xlim)],
            carbon_best_vario$v[carbon_best_vario$u <= max(xlim)]) * 1.1
ylim <- c(0, round(ylim, 2))
# l1 <- linesREML(carbon_poor_lmm, add = FALSE)
l2 <- linesREML(carbon_base_lmm, add = FALSE)
# l3 <- linesREML(carbon_fine_lmm, add = FALSE)
l4 <- linesREML(carbon_best_lmm, add = FALSE)
# v1 <- xyplot(carbon_poor_vario$v ~ carbon_poor_vario$u, ylim = ylim, pch = 20,
#              col =  "black", scales = list(tick.number = 8), xlim = xlim,
#              panel = function(x, y, ...) {
#                panel.xyplot(x, y, ...)
#                panel.lines(x = l1$x, y = l1$y, col = "black", lty = 2)
#              })
soc_v2 <- xyplot(carbon_base_vario$v ~ carbon_base_vario$u, ylim = ylim, pch = 20,
             col =  "black", scales = list(tick.number = 8), xlim = xlim,
             panel = function(x, y, ...) {
               panel.xyplot(x, y, ...)
               panel.lines(x = l2$x, y = l2$y, col = "black", lty = 2)
             })
# soc_v3 <- xyplot(carbon_fine_vario$v ~ carbon_fine_vario$u, ylim = ylim, pch = 20,
#              col =  "black", scales = list(tick.number = 8), xlim = xlim,
#              panel = function(x, y, ...) {
#                panel.xyplot(x, y, ...)
#                panel.lines(x = l3$x, y = l3$y, col = "black", lty = 2)
#              })
soc_v4 <- xyplot(carbon_best_vario$v ~ carbon_best_vario$u, ylim = ylim, pch = 20,
             col =  "black", scales = list(tick.number = 8), xlim = xlim,
             panel = function(x, y, ...) {
               panel.xyplot(x, y, ...)
               panel.lines(x = l4$x, y = l4$y, col = "black", lty = 2)
             })
ylab <- expression(paste("Semivariance [(g kg"^"-1", ")"^"2", "]", sep = ""))
# soc_v <- update(c(v1, v2, v3, v4), ylab = ylab, layout = c(4, 1), 
#             xlab = "Distance [m]", asp = 1, scales = list(cex = c(0.9, 0.9)))
soc_v <- update(c(soc_v2, soc_v4), ylab = ylab, layout = c(2, 1), 
            xlab = "Distance [m]", asp = 1, scales = list(cex = c(1, 1)))
# save plot
dev.off()
pdf(file = paste(update_covars_dir, "carbon_lmm.pdf", sep = ""), 
    width = 9/cm(1), height = 6/cm(1))
trellis.par.set(fontsize = list(text = 7, points = 5),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
plot(soc_v)
dev.off()
rm(xlim, ylim, l1, l2, soc_v1, soc_v2, soc_v3, soc_v4, soc_v)
gc()

# ECEC - poor model ------------------------------------------------------------
# create geodata object
model     <- ecec_sel$poor_lm
data      <- as.data.frame(cal_data)
covars    <- attr(terms(model), "term.labels")
covar.col <- which(match(colnames(data), covars) != "NA")
data.col  <- "ecec"
data.col  <- which(colnames(data) == data.col)
ecec_poor_geodata <- as.geodata(data, coords.col = c(1, 2), data.col = data.col,
                                covar.col = covar.col)
summary(ecec_poor_geodata)
rm(model, data, covar.col, data.col)
# calculate empirical variogram
geodata         <- ecec_poor_geodata
trend           <- formula(paste("~", paste(covars, collapse = " + ")))
lambda          <- bc_lambda$ecec
breaks          <- seq(0, 6500, 100)
ecec_poor_vario <- variog(geodata, trend = trend, breaks = breaks,
                          lambda = lambda)
plot(ecec_poor_vario, col = "blue", pch = 20, type = "b")
# estimate model parameters using REML
ini.cov.pars   <- as.matrix(expand.grid(c(0.30, 0.32, 0.34), c(500, 600, 700)))
nugget         <- c(0.20, 0.22, 0.24)
ecec_poor_lmm  <- likfit(geodata, trend = trend, ini.cov.pars = ini.cov.pars,
                         nugget = nugget, lambda = lambda, lik.method = "REML")
summary(ecec_poor_lmm)
# clean workspace
rm(covars, geodata, trend, breaks, lambda, nugget, ini.cov.pars)
gc()

# ECEC - base model ------------------------------------------------------------
# create geodata object
model     <- ecec_base_lm
data      <- as.data.frame(cal_data)
covars    <- attr(terms(model), "term.labels")
covar.col <- which(match(colnames(data), covars) != "NA")
data.col  <- "ecec"
data.col  <- which(colnames(data) == data.col)
ecec_base_geodata <- as.geodata(data, coords.col = c(1, 2), data.col = data.col,
                        covar.col = covar.col)
summary(ecec_base_geodata)
rm(model, data, covar.col, data.col)
# calculate empirical variogram
geodata         <- ecec_base_geodata
trend           <- formula(paste("~", paste(covars, collapse = " + ")))
lambda          <- bc_lambda$ecec
breaks          <- seq(0, 6500, 100)
ecec_base_vario <- variog(geodata, trend = trend, breaks = breaks,
                          lambda = lambda)
plot(ecec_base_vario, col = "blue", pch = 20, type = "b")
# estimate model parameters using REML
ini.cov.pars   <- as.matrix(expand.grid(c(0.30, 0.32, 0.34), c(500, 600, 700)))
nugget         <- c(0.20, 0.22, 0.24)
ecec_base_lmm  <- likfit(geodata, trend = trend, ini.cov.pars = ini.cov.pars,
                         nugget = nugget, lambda = lambda, lik.method = "REML")
summary(ecec_base_lmm)
# clean workspace
rm(covars, geodata, trend, breaks, lambda, nugget, ini.cov.pars)
gc()

# ECEC - fine model ------------------------------------------------------------
# create geodata object
model     <- ecec_sel$fine_lm
data      <- as.data.frame(cal_data)
covars    <- attr(terms(model), "term.labels")
covar.col <- which(match(colnames(data), covars) != "NA")
data.col  <- "ecec"
data.col  <- which(colnames(data) == data.col)
ecec_fine_geodata <- as.geodata(data, coords.col = c(1, 2), data.col = data.col,
                                covar.col = covar.col)
summary(ecec_fine_geodata)
rm(model, data, covar.col, data.col)
gc()
# calculate empirical variogram
geodata         <- ecec_fine_geodata
trend           <- formula(paste("~", paste(covars, collapse = " + ")))
lambda          <- bc_lambda$ecec
breaks          <- seq(0, 6500, 100)
ecec_fine_vario <- variog(geodata, trend = trend, lambda = lambda,
                          breaks = breaks)
plot(ecec_fine_vario, col = "blue", pch = 20, type = "b")
# estimate model parameters using REML
ini.cov.pars  <- as.matrix(expand.grid(c(0.20, 0.24, 0.28), c(500, 600, 700)))
nugget        <- c(0.16, 0.20, 0.24)
ecec_fine_lmm <- likfit(geodata, trend = trend, ini.cov.pars = ini.cov.pars,
                        nugget = nugget, lambda = lambda, lik.method = "REML")
summary(ecec_fine_lmm)
# clean workspace
rm(covars, geodata, trend, breaks, lambda, nugget, ini.cov.pars)
gc()

# ECEC - best model ------------------------------------------------------------
# create geodata object
model     <- ecec_best_lm
data      <- as.data.frame(cal_data)
covars    <- attr(terms(model), "term.labels")
covar.col <- which(match(colnames(data), covars) != "NA")
data.col  <- "ecec"
data.col  <- which(colnames(data) == data.col)
ecec_best_geodata <- as.geodata(data, coords.col = c(1, 2), data.col = data.col,
                        covar.col = covar.col)
summary(ecec_best_geodata)
rm(model, data, covar.col, data.col)
gc()
# calculate empirical variogram
geodata         <- ecec_best_geodata
trend           <- formula(paste("~", paste(covars, collapse = " + ")))
lambda          <- bc_lambda$ecec
breaks          <- seq(0, 6500, 100)
ecec_best_vario <- variog(geodata, trend = trend, lambda = lambda,
                          breaks = breaks)
plot(ecec_best_vario, col = "blue", pch = 20, type = "b")
# estimate model parameters using REML
ini.cov.pars  <- as.matrix(expand.grid(c(0.20, 0.24, 0.28), c(500, 600, 700)))
nugget        <- c(0.16, 0.20, 0.24)
ecec_best_lmm <- likfit(geodata, trend = trend, ini.cov.pars = ini.cov.pars,
                        nugget = nugget, lambda = lambda, lik.method = "REML")
summary(ecec_best_lmm)
# clean workspace
rm(covars, geodata, trend, breaks, lambda, nugget, ini.cov.pars)
gc()

# ECEC - plot experimental variograms and fitted linear mixed models -----------
xlim <- c(0, 3000)
ylim <- max(#ecec_poor_vario$v[ecec_poor_vario$u <= max(xlim)],
            ecec_base_vario$v[ecec_base_vario$u <= max(xlim)],
            #ecec_fine_vario$v[ecec_fine_vario$u <= max(xlim)],
            ecec_best_vario$v[ecec_best_vario$u <= max(xlim)]) * 1.1
ylim <- c(0, round(ylim, 2))
#l1 <- linesREML(ecec_poor_lmm, add = FALSE)
l2 <- linesREML(ecec_base_lmm, add = FALSE)
#l3 <- linesREML(ecec_fine_lmm, add = FALSE)
l4 <- linesREML(ecec_best_lmm, add = FALSE)
# v1 <- xyplot(ecec_poor_vario$v ~ ecec_poor_vario$u, ylim = ylim, pch = 20,
#              col =  "black", scales = list(tick.number = 8), xlim = xlim,
#              panel = function(x, y, ...) {
#                panel.xyplot(x, y, ...)
#                panel.lines(x = l1$x, y = l1$y, col = "black", lty = 2)
#              })
ecec_v2 <- xyplot(ecec_base_vario$v ~ ecec_base_vario$u, ylim = ylim, pch = 20,
             col =  "black", scales = list(tick.number = 8), xlim = xlim,
             panel = function(x, y, ...) {
               panel.xyplot(x, y, ...)
               panel.lines(x = l2$x, y = l2$y, col = "black", lty = 2)
             })
# v3 <- xyplot(ecec_fine_vario$v ~ ecec_fine_vario$u, ylim = ylim, pch = 20,
#              col =  "black", scales = list(tick.number = 8), xlim = xlim,
#              panel = function(x, y, ...) {
#                panel.xyplot(x, y, ...)
#                panel.lines(x = l3$x, y = l3$y, col = "black", lty = 2)
#              })
ecec_v4 <- xyplot(ecec_best_vario$v ~ ecec_best_vario$u, ylim = ylim, pch = 20,
             col =  "black", scales = list(tick.number = 8), xlim = xlim,
             panel = function(x, y, ...) {
               panel.xyplot(x, y, ...)
               panel.lines(x = l4$x, y = l4$y, col = "black", lty = 2)
             })
ylab <- expression(paste("Semivariance [(mmol kg"^"-1", ")"^"2", "]", sep = ""))
# v3 <- update(c(v1, v2, v3, v4), ylab = ylab, xlab = "Distance [m]", asp = 1,
#             layout = c(4, 1), scales = list(cex = c(0.9, 0.9)))
ecec_v <- update(c(ecec_v2, ecec_v4), ylab = ylab, xlab = "Distance [m]", asp = 1,
             layout = c(2, 1), scales = list(cex = c(1, 1)))
# save plot
dev.off()
pdf(file = paste(update_covars_dir, "ecec_lmm.pdf", sep = ""), 
    width = 9/cm(1), height = 6/cm(1))
trellis.par.set(fontsize = list(text = 7, points = 5),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
plot(ecec_v)
dev.off()
rm(xlim, ylim, l1, l2, l3, l4, ecec_v1, ecec_v2, ecec_v3, ecec_v4, ecec_v)
gc()

# LEAVE-ONE-OUT CROSS-VALIDATION ###############################################
# Start checking how many realizations are needed to stabilize the variance.
# Use the same number of realizations for all cases.
# CLAY
model           <- clay_base_lm
lambda          <- bc_lambda$clay
geodata         <- clay_base_geodata
clay_base_lm_cv <- looCV(model, geodata = geodata, lambda = lambda, 
                         simul.back = FALSE)
invBoxCox(mean = clay_base_lm_cv$pred, variance = c(clay_base_lm_cv$se.pred^2),
          lambda = lambda, profile = TRUE, simul.back = TRUE)
rm(model, lambda, geodata)
# SOC
model             <- carbon_base_lm
lambda            <- bc_lambda$carbon
geodata           <- carbon_base_geodata
carbon_base_lm_cv <- looCV(model, geodata = geodata, lambda = lambda, 
                         simul.back = FALSE)
invBoxCox(mean = carbon_base_lm_cv$pred, 
          variance = c(carbon_base_lm_cv$se.pred^2),
          lambda = lambda, profile = TRUE, simul.back = TRUE)
rm(model, lambda, geodata)
# ECEC
model           <- ecec_base_lm
lambda          <- bc_lambda$ecec
geodata         <- ecec_base_geodata
ecec_base_lm_cv <- looCV(model, geodata = geodata, lambda = lambda, 
                         simul.back = FALSE)
invBoxCox(mean = ecec_base_lm_cv$pred, variance = c(ecec_base_lm_cv$se.pred^2),
          lambda = lambda, profile = TRUE, simul.back = TRUE)
rm(model, lambda, geodata)

# CLAY - base linear model -----------------------------------------------------
model           <- clay_base_lm
lambda          <- bc_lambda$clay
geodata         <- clay_base_geodata
clay_base_lm_cv <- looCV(model, geodata = geodata, lambda = lambda, 
                         simul.back = TRUE, n.sim = 20000)
rm(model, lambda, geodata)

# CLAY - best linear model -----------------------------------------------------
model           <- clay_best_lm
lambda          <- bc_lambda$clay
geodata         <- clay_best_geodata
clay_best_lm_cv <- looCV(model, geodata = geodata, lambda = lambda, 
                         simul.back = TRUE, n.sim = 20000)
rm(model, lambda, geodata)

# CLAY - base linear mixed model -----------------------------------------------
geodata          <- clay_base_geodata
model            <- clay_base_lmm
clay_base_lmm_cv <- cvKrige(geodata, model = model, reestimate = TRUE,
                            output.reestimate = TRUE, n.sim = 20000)
rm(geodata, model)

# CLAY - best linear mixed model -----------------------------------------------
geodata          <- clay_best_geodata
model            <- clay_best_lmm
clay_best_lmm_cv <- cvKrige(geodata, model = model, reestimate = TRUE,
                            output.reestimate = TRUE, n.sim = 20000)
rm(geodata, model)

# CARBON - base linear model ---------------------------------------------------
model             <- carbon_base_lm
lambda            <- bc_lambda$carbon
geodata           <- carbon_base_geodata
carbon_base_lm_cv <- looCV(model, geodata = geodata, lambda = lambda, 
                           simul.back = TRUE, n.sim = 20000)
rm(model, lambda, geodata)

# CARBON - best linear model ---------------------------------------------------
model             <- carbon_best_lm
lambda            <- bc_lambda$carbon
geodata           <- carbon_best_geodata
carbon_best_lm_cv <- looCV(model, geodata = geodata, lambda = lambda, 
                           simul.back = TRUE, n.sim = 20000)
rm(model, lambda, geodata)

# CARBON - base linear mixed model ---------------------------------------------
geodata            <- carbon_base_geodata
model              <- carbon_base_lmm
carbon_base_lmm_cv <- cvKrige(geodata, model = model, reestimate = TRUE,
                              output.reestimate = TRUE, n.sim = 20000)
rm(geodata, model)

# CARBON - best linear mixed model ---------------------------------------------
geodata            <- carbon_best_geodata
model              <- carbon_best_lmm
carbon_best_lmm_cv <- cvKrige(geodata, model = model, reestimate = TRUE,
                              output.reestimate = TRUE, n.sim = 20000)
rm(geodata, model)

# ECEC - base linear model -----------------------------------------------------
model           <- ecec_base_lm
lambda          <- bc_lambda$ecec
geodata         <- ecec_base_geodata
ecec_base_lm_cv <- looCV(model, geodata = geodata, lambda = lambda, 
                         simul.back = TRUE, n.sim = 20000)
rm(model, lambda, geodata)

# ECEC - best linear model -----------------------------------------------------
model           <- ecec_best_lm
lambda          <- bc_lambda$ecec
geodata         <- ecec_best_geodata
ecec_best_lm_cv <- looCV(model, geodata = geodata, lambda = lambda, 
                         simul.back = TRUE, n.sim = 20000)
rm(model, lambda, geodata)

# ECEC - base linear mixed model -----------------------------------------------
geodata          <- ecec_base_geodata
model            <- ecec_base_lmm
ecec_base_lmm_cv <- cvKrige(geodata, model = model, reestimate = TRUE,
                            output.reestimate = TRUE, n.sim = 20000)
rm(geodata, model)
gc()

# ECEC - best linear mixed model -----------------------------------------------
geodata          <- ecec_best_geodata
model            <- ecec_best_lmm
ecec_best_lmm_cv <- cvKrige(geodata, model = model, reestimate = TRUE,
                            output.reestimate = TRUE, n.sim = 20000)
rm(geodata, model)
gc()

# SPATIAL PREDICTION - KRIGING #################################################
# Universal kriging, also known as kriging with external drift (trend)
# prepare base data
system("r.mask -o input=buffer_BASIN_10")
locations <- readRAST6("dnos.raster")

# CLAY - base linear mixed model -----------------------------------------------
# prepare base data
geodata <- clay_base_geodata
model   <- clay_base_lmm
file    <- "clay_base_lmm_krige.rda"
covars  <- colnames(geodata$covariate)
trend   <- formula(paste("~", paste(covars, collapse = " + ")))
covars  <- readRAST6(covars)
save(locations, geodata, model, file, covars, trend, spredict, file = "test-server.rda")
# make predictions
spredict(model = model, geodata = geodata, file = file, covars = covars, 
         simul.back = TRUE, n.sim = 20000, pred.loc = locations,
         n.tiles = 1000, email = "alessandrosamuelrosa@gmail.com")

# CLAY - best linear mixed model -----------------------------------------------
# prepare base data
geodata <- clay_best_geodata
model   <- clay_best_lmm
file    <- "clay_best_lmm_krige.rda"
covars  <- colnames(geodata$covariate)
trend   <- formula(paste("~", paste(covars, collapse = " + ")))
covars  <- readRAST6(covars)
# make predictions
spredict(model = model, geodata = geodata, file = file, covars = covars,
         simul.back = TRUE, n.sim = 20000, pred.loc = locations,
         n.tiles = 1000, email = "alessandrosamuelrosa@gmail.com")

# CLAY - kriging predictions ---------------------------------------------------
load("clay_base_lmm_krige.rda")
clay_base_lmm_krige <- data.frame(pred[1:4])
str(clay_base_lmm_krige)
rm(pred)
coordinates(clay_base_lmm_krige) <- ~ x.coord + y.coord
gridded(clay_base_lmm_krige) <- TRUE
proj4string(clay_base_lmm_krige) <- wgs1984utm22s
str(clay_base_lmm_krige)
image(clay_base_lmm_krige)
spplot(clay_base_lmm_krige, zcol = "krige.var", col.regions = bpy.colors(256))

load("clay_best_lmm_krige.rda")
clay_best_lmm_krige <- data.frame(pred[1:4])
str(clay_best_lmm_krige)
rm(pred)
coordinates(clay_best_lmm_krige) <- ~ x.coord + y.coord
gridded(clay_best_lmm_krige) <- TRUE
proj4string(clay_best_lmm_krige) <- wgs1984utm22s
str(clay_best_lmm_krige)
image(clay_best_lmm_krige)
spplot(clay_best_lmm_krige, zcol = "krige.var", col.regions = bpy.colors(256))

# prepare plots with predictions
breaks <- seq(min(cal_data$clay), max(cal_data$clay), by = 1)
max_z <- max(clay_base_lmm_krige$krige.pred, clay_best_lmm_krige$krige.pred)
if (max(breaks) < max_z) {
  breaks <- c(breaks, max_z)
}
psd.colors <- colorRampPalette(R_pal$tex_pal)
col <- psd.colors(length(breaks)-1)
clay_base_lmm_krige$krige.pred.cut <- cut(clay_base_lmm_krige$krige.pred, 
                                            breaks = breaks)
clay_best_lmm_krige$krige.pred.cut <- cut(clay_best_lmm_krige$krige.pred, 
                                            breaks = breaks)
p1 <- spplot(clay_base_lmm_krige, "krige.pred.cut", col.regions = col)
p2 <- spplot(clay_best_lmm_krige, "krige.pred.cut", col.regions = col)
p1$legend$right$args$key <- p1$legend$right$args$key[1:3]
p2$legend$right$args$key <- p2$legend$right$args$key[1:3]
# save plots with predictions
dev.off()
pdf(file = paste(update_covars_dir, "clay_base_lmm_krige.pdf", sep = ""),
    width = 6.3/cm(1), height = 6.3/cm(1))
trellis.par.set(fontsize = list(text = 7, points = 5),
                plot.line = list(lwd = 0.001),
                axis.line = list(lwd = 0.01),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
plot(p1)
dev.off()
pdf(file = paste(update_covars_dir, "clay_best_lmm_krige.pdf", sep = ""),
    width = 6.3/cm(1), height = 6.3/cm(1))
trellis.par.set(fontsize = list(text = 7, points = 5),
                plot.line = list(lwd = 0.001),
                axis.line = list(lwd = 0.01),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
plot(p2)
dev.off()
rm(p1, p2)
gc()
# prepare plots with prediction variance
# take the square root to improve plotting
clay_base_lmm_krige$krige.var <- sqrt(clay_base_lmm_krige$krige.var)
clay_best_lmm_krige$krige.var <- sqrt(clay_best_lmm_krige$krige.var)
max_z <- c(max(clay_base_lmm_krige$krige.var), 
           max(clay_best_lmm_krige$krige.var))
max_z <- max_z[which.min(max_z)]
breaks <- seq(0, max_z, by = 1)
max_z <- max(clay_base_lmm_krige$krige.var, clay_best_lmm_krige$krige.var)
if (max(breaks) < max_z) {
  breaks <- c(breaks, max_z)
}
col <- bpy.colors(length(breaks)-1)
clay_base_lmm_krige$krige.var.cut <- cut(clay_base_lmm_krige$krige.var, 
                                           breaks = breaks)
clay_best_lmm_krige$krige.var.cut <- cut(clay_best_lmm_krige$krige.var, 
                                           breaks = breaks)
p1 <- spplot(clay_base_lmm_krige, "krige.var.cut", col.regions = col)
p2 <- spplot(clay_best_lmm_krige, "krige.var.cut", col.regions = col)
p1$legend$right$args$key <- p1$legend$right$args$key[1:3]
p2$legend$right$args$key <- p2$legend$right$args$key[1:3]
# save plots with prediction variance
dev.off()
pdf(file = paste(update_covars_dir, "clay_base_lmm_krige_var.pdf", sep = ""),
    width = 6.3/cm(1), height = 6.3/cm(1))
trellis.par.set(fontsize = list(text = 7, points = 5),
                plot.line = list(lwd = 0.001),
                axis.line = list(lwd = 0.01),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
plot(p1)
dev.off()
pdf(file = paste(update_covars_dir, "clay_best_lmm_krige_var.pdf", sep = ""),
    width = 6.3/cm(1), height = 6.3/cm(1))
trellis.par.set(fontsize = list(text = 7, points = 5),
                plot.line = list(lwd = 0.001),
                axis.line = list(lwd = 0.01),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
plot(p2)
dev.off()
rm(p1, p2)
gc()

# CARBON - base linear mixed model ---------------------------------------------
# prepare base data
geodata <- carbon_base_geodata
model   <- carbon_base_lmm
file    <- "carbon_base_lmm_krige.rda"
covars  <- colnames(geodata$covariate)
trend   <- formula(paste("~", paste(covars, collapse = " + ")))
covars  <- readRAST6(covars)
# make predictions
spredict(model = model, geodata = geodata, file = file, covars = covars, 
         simul.back = TRUE, n.sim = 20000, pred.loc = locations,
         n.tiles = 1000, email = "alessandrosamuelrosa@gmail.com")

# CARBON - best linear mixed model ---------------------------------------------
# prepare base data
geodata <- carbon_best_geodata
model   <- carbon_best_lmm
file    <- "carbon_best_lmm_krige.rda"
covars  <- colnames(geodata$covariate)
trend   <- formula(paste("~", paste(covars, collapse = " + ")))
covars  <- readRAST6(covars)
# make predictions
spredict(model = model, geodata = geodata, file = file, covars = covars, 
         simul.back = TRUE, n.sim = 20000, pred.loc = locations,
         n.tiles = 1000, email = "alessandrosamuelrosa@gmail.com")

# CARBON - kriging predictions -------------------------------------------------
load("carbon_base_lmm_krige.rda")
carbon_base_lmm_krige <- data.frame(pred[1:4])
str(carbon_base_lmm_krige)
rm(pred)
coordinates(carbon_base_lmm_krige) <- ~ x.coord + y.coord
gridded(carbon_base_lmm_krige) <- TRUE
proj4string(carbon_base_lmm_krige) <- wgs1984utm22s
str(carbon_base_lmm_krige)
image(carbon_base_lmm_krige)

load("carbon_best_lmm_krige.rda")
carbon_best_lmm_krige <- data.frame(pred[1:4])
str(carbon_best_lmm_krige)
rm(pred)
coordinates(carbon_best_lmm_krige) <- ~ x.coord + y.coord
gridded(carbon_best_lmm_krige) <- TRUE
proj4string(carbon_best_lmm_krige) <- wgs1984utm22s
str(carbon_best_lmm_krige)
image(carbon_best_lmm_krige)
# prepare plots with predictions
breaks <- seq(min(cal_data$carbon), max(cal_data$carbon), by = 1)
max_z <- max(carbon_base_lmm_krige$krige.pred, carbon_best_lmm_krige$krige.pred)
if (max(breaks) < max_z) {
  breaks <- c(breaks, max_z)
}
soc.colors <- colorRampPalette(R_pal$soc_pal)
col <- soc.colors(length(breaks)-1)
carbon_base_lmm_krige$krige.pred.cut <- cut(carbon_base_lmm_krige$krige.pred, 
                                            breaks = breaks)
carbon_best_lmm_krige$krige.pred.cut <- cut(carbon_best_lmm_krige$krige.pred, 
                                            breaks = breaks)
p1 <- spplot(carbon_base_lmm_krige, "krige.pred.cut", col.regions = col)
p2 <- spplot(carbon_best_lmm_krige, "krige.pred.cut", col.regions = col)
p1$legend$right$args$key <- p1$legend$right$args$key[1:3]
p2$legend$right$args$key <- p2$legend$right$args$key[1:3]
# save plots with predictions
dev.off()
pdf(file = paste(update_covars_dir, "carbon_base_lmm_krige.pdf", sep = ""),
    width = 6.3/cm(1), height = 6.3/cm(1))
trellis.par.set(fontsize = list(text = 7, points = 5),
                plot.line = list(lwd = 0.001),
                axis.line = list(lwd = 0.01),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
plot(p1)
dev.off()
pdf(file = paste(update_covars_dir, "carbon_best_lmm_krige.pdf", sep = ""),
    width = 6.3/cm(1), height = 6.3/cm(1))
trellis.par.set(fontsize = list(text = 7, points = 5),
                plot.line = list(lwd = 0.001),
                axis.line = list(lwd = 0.01),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
plot(p2)
dev.off()
rm(p1, p2)
gc()
# prepare plots with prediction variance
# take the square root to improve plotting
carbon_base_lmm_krige$krige.var <- sqrt(carbon_base_lmm_krige$krige.var)
carbon_best_lmm_krige$krige.var <- sqrt(carbon_best_lmm_krige$krige.var)
max_z <- c(max(carbon_base_lmm_krige$krige.var), 
           max(carbon_best_lmm_krige$krige.var))
max_z <- max_z[which.min(max_z)]
breaks <- seq(0, max_z, by = 1)
max_z <- max(carbon_base_lmm_krige$krige.var, carbon_best_lmm_krige$krige.var)
if (max(breaks) < max_z) {
  breaks <- c(breaks, max_z)
}
col <- bpy.colors(length(breaks)-1)
carbon_base_lmm_krige$krige.var.cut <- cut(carbon_base_lmm_krige$krige.var, 
                                           breaks = breaks)
carbon_best_lmm_krige$krige.var.cut <- cut(carbon_best_lmm_krige$krige.var, 
                                           breaks = breaks)
p1 <- spplot(carbon_base_lmm_krige, "krige.var.cut", col.regions = col)
p2 <- spplot(carbon_best_lmm_krige, "krige.var.cut", col.regions = col)
p1$legend$right$args$key <- p1$legend$right$args$key[1:3]
p2$legend$right$args$key <- p2$legend$right$args$key[1:3]
# save plots with prediction variance
dev.off()
pdf(file = paste(update_covars_dir, "carbon_base_lmm_krige_var.pdf", sep = ""),
    width = 6.3/cm(1), height = 6.3/cm(1))
trellis.par.set(fontsize = list(text = 7, points = 5),
                plot.line = list(lwd = 0.001),
                axis.line = list(lwd = 0.01),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
plot(p1)
dev.off()
pdf(file = paste(update_covars_dir, "carbon_best_lmm_krige_var.pdf", sep = ""),
    width = 6.3/cm(1), height = 6.3/cm(1))
trellis.par.set(fontsize = list(text = 7, points = 5),
                plot.line = list(lwd = 0.001),
                axis.line = list(lwd = 0.01),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
plot(p2)
dev.off()
rm(p1, p2)
gc()

# extrapolation
a <- which(carbon_best_lmm_krige$krige.pred >= max(cal_data$carbon))
b <- coordinates(carbon_best_lmm_krige[a, ])
pts <- data.frame(b, carbon_best_lmm_krige$krige.pred[a])
coordinates(pts) <- ~ x.coord + y.coord
covar <- raster(readRAST6("TPI_10_15"))
image(covar, asp = 1)
points(pts, cex = 0.1, pch = 20)
pts$TPI_10_15 <- extract(x = covar, y = pts)
quantile(round(cal_data$TPI_10_15), prob = seq(0, 1, 0.1))
quantile(round(pts$TPI_10_15), prob = seq(0, 1, 0.1))

# ECEC - base linear mixed model -----------------------------------------------
# prepare base data
geodata <- ecec_base_geodata
model   <- ecec_base_lmm
file    <- "ecec_base_lmm_krige.rda"
covars  <- colnames(geodata$covariate)
trend   <- formula(paste("~", paste(covars, collapse = " + ")))
covars  <- readRAST6(covars)
# make predictions
spredict(model = model, geodata = geodata, file = file, covars = covars, 
         simul.back = TRUE, n.sim = 20000, pred.loc = locations,
         n.tiles = 1000, email = "alessandrosamuelrosa@gmail.com")

# ECEC - best linear mixed model -----------------------------------------------
# prepare base data
geodata <- ecec_best_geodata
model   <- ecec_best_lmm
file    <- "ecec_best_lmm_krige.rda"
covars  <- colnames(geodata$covariate)
trend   <- formula(paste("~", paste(covars, collapse = " + ")))
covars  <- readRAST6(covars)
# make predictions
spredict(model = model, geodata = geodata, file = file, covars = covars, 
         simul.back = TRUE, n.sim = 20000, pred.loc = locations,
         n.tiles = 1000, email = "alessandrosamuelrosa@gmail.com")

# ECEC - kriging predictions ---------------------------------------------------
load("ecec_base_lmm_krige.rda")
ecec_base_lmm_krige <- data.frame(pred[1:4])
str(ecec_base_lmm_krige)
rm(pred)
gc()
coordinates(ecec_base_lmm_krige) <- ~ x.coord + y.coord
gridded(ecec_base_lmm_krige) <- TRUE
proj4string(ecec_base_lmm_krige) <- wgs1984utm22s
str(ecec_base_lmm_krige)
image(ecec_base_lmm_krige)
spplot(ecec_base_lmm_krige, zcol = "krige.var", col.regions = bpy.colors(256))
summary(ecec_base_lmm_krige)

load("ecec_best_lmm_krige.rda")
ecec_best_lmm_krige <- data.frame(pred[1:4])
str(ecec_best_lmm_krige)
rm(pred)
gc()
coordinates(ecec_best_lmm_krige) <- ~ x.coord + y.coord
gridded(ecec_best_lmm_krige) <- TRUE
proj4string(ecec_best_lmm_krige) <- wgs1984utm22s
str(ecec_best_lmm_krige)
image(ecec_best_lmm_krige, zcol = "krige.pred")
spplot(ecec_best_lmm_krige, zcol = "krige.var", col.regions = bpy.colors(256))
summary(ecec_best_lmm_krige)
# prepare plots with predictions
breaks <- seq(min(cal_data$ecec), max(cal_data$ecec), by = 1)
max_z <- max(ecec_base_lmm_krige$krige.pred, ecec_best_lmm_krige$krige.pred)
if (max(breaks) < max_z) {
  breaks <- c(breaks, max_z)
}
ecec.colors <- colorRampPalette(R_pal$CEC_pal)
col <- ecec.colors(length(breaks)-1)
ecec_base_lmm_krige$krige.pred.cut <- cut(ecec_base_lmm_krige$krige.pred, 
                                          breaks = breaks)
ecec_best_lmm_krige$krige.pred.cut <- cut(ecec_best_lmm_krige$krige.pred, 
                                          breaks = breaks)
p1 <- spplot(ecec_base_lmm_krige, "krige.pred.cut", col.regions = col)
p2 <- spplot(ecec_best_lmm_krige, "krige.pred.cut", col.regions = col)
p1$legend$right$args$key <- p1$legend$right$args$key[1:3]
p2$legend$right$args$key <- p2$legend$right$args$key[1:3]
# save plots with predictions
dev.off()
pdf(file = paste(update_covars_dir, "ecec_base_lmm_krige.pdf", sep = ""),
    width = 6.3/cm(1), height = 6.3/cm(1))
trellis.par.set(fontsize = list(text = 7, points = 5),
                plot.line = list(lwd = 0.001),
                axis.line = list(lwd = 0.01),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
plot(p1)
dev.off()
pdf(file = paste(update_covars_dir, "ecec_best_lmm_krige.pdf", sep = ""),
    width = 6.3/cm(1), height = 6.3/cm(1))
trellis.par.set(fontsize = list(text = 7, points = 5),
                plot.line = list(lwd = 0.001),
                axis.line = list(lwd = 0.01),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
plot(p2)
dev.off()
rm(p1, p2)
gc()
# prepare plots with prediction variance
ecec_base_lmm_krige$krige.var <- sqrt(ecec_base_lmm_krige$krige.var)
ecec_best_lmm_krige$krige.var <- sqrt(ecec_best_lmm_krige$krige.var)
max_z <- c(max(ecec_base_lmm_krige$krige.var), 
           max(ecec_best_lmm_krige$krige.var))
max_z <- max_z[which.min(max_z)]
breaks <- seq(0, max_z, by = 1)
max_z <- max(ecec_base_lmm_krige$krige.var, ecec_best_lmm_krige$krige.var)
if (max(breaks) < max_z) {
  breaks <- c(breaks, max_z)
}
col <- bpy.colors(length(breaks)-1)
ecec_base_lmm_krige$krige.var.cut <- cut(ecec_base_lmm_krige$krige.var, 
                                         breaks = breaks)
ecec_best_lmm_krige$krige.var.cut <- cut(ecec_best_lmm_krige$krige.var, 
                                         breaks = breaks)
p1 <- spplot(ecec_base_lmm_krige, "krige.var.cut", col.regions = col)
p2 <- spplot(ecec_best_lmm_krige, "krige.var.cut", col.regions = col)
p1$legend$right$args$key <- p1$legend$right$args$key[1:3]
p2$legend$right$args$key <- p2$legend$right$args$key[1:3]
# save plots with prediction variance
dev.off()
pdf(file = paste(update_covars_dir, "ecec_base_lmm_krige_var.pdf", sep = ""),
    width = 6.3/cm(1), height = 6.3/cm(1))
trellis.par.set(fontsize = list(text = 7, points = 5),
                plot.line = list(lwd = 0.001),
                axis.line = list(lwd = 0.01),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
plot(p1)
dev.off()
pdf(file = paste(update_covars_dir, "ecec_best_lmm_krige_var.pdf", sep = ""),
    width = 6.3/cm(1), height = 6.3/cm(1))
trellis.par.set(fontsize = list(text = 7, points = 5),
                plot.line = list(lwd = 0.001),
                axis.line = list(lwd = 0.01),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
plot(p2)
dev.off()
rm(p1, p2)
gc()

# TABLES WITH RESULTS ##########################################################

# Tables - selected models and cross-validation statistics ---------------------
# selected models
covars <- list(CLAYbase = paste(colnames(clay_base_geodata$covariate), 
                                 collapse = "}, \texttt{"),
               CLAYbest = paste(colnames(clay_best_geodata$covariate), 
                                 collapse = "}, \texttt{"),
               SOCbase = paste(colnames(carbon_base_geodata$covariate), 
                                collapse = "}, \texttt{"),
               SOCbest = paste(colnames(carbon_best_geodata$covariate), 
                                collapse = "}, \texttt{"),
               ECECbase = paste(colnames(ecec_base_geodata$covariate), 
                                 collapse = "}, \texttt{"),
               ECECbest = paste(colnames(ecec_best_geodata$covariate), 
                                 collapse = "}, \texttt{"))
covars <- lapply(covars, 
                 function (X) {paste("\texttt{", X, "}", sep = "")})
covars <- data.frame(covars, stringsAsFactors = FALSE)
covars

# cross-validation statistics (ME, RMSE, SRMSE, R-square)
obs <- data.frame(cal_data$clay, cal_data$clay, cal_data$clay, cal_data$clay, 
                  cal_data$carbon, cal_data$carbon, cal_data$carbon,
                  cal_data$carbon, cal_data$ecec, cal_data$ecec, cal_data$ecec,
                  cal_data$ecec)
str(obs)               
pred <- data.frame(clay_base_lm_cv$back.transformed$mean,
                   clay_base_lmm_cv$back.transformed$predicted,
                   clay_best_lm_cv$back.transformed$mean,
                   clay_best_lmm_cv$back.transformed$predicted,
                   carbon_base_lm_cv$back.transformed$mean,
                   carbon_base_lmm_cv$back.transformed$predicted,
                   carbon_best_lm_cv$back.transformed$mean,
                   carbon_best_lmm_cv$back.transformed$predicted,
                   ecec_base_lm_cv$back.transformed$mean,
                   ecec_base_lmm_cv$back.transformed$predicted,
                   ecec_best_lm_cv$back.transformed$mean,
                   ecec_best_lmm_cv$back.transformed$predicted)
str(pred)
pev <- data.frame(clay_base_lm_cv$back.transformed$variance,
                  clay_base_lmm_cv$back.transformed$krige.var,
                  clay_best_lm_cv$back.transformed$variance,
                  clay_best_lmm_cv$back.transformed$krige.var,
                  carbon_base_lm_cv$back.transformed$variance,
                  carbon_base_lmm_cv$back.transformed$krige.var,
                  carbon_best_lm_cv$back.transformed$variance,
                  carbon_best_lmm_cv$back.transformed$krige.var,
                  ecec_base_lm_cv$back.transformed$variance,
                  ecec_base_lmm_cv$back.transformed$krige.var,
                  ecec_best_lm_cv$back.transformed$variance,
                  ecec_best_lmm_cv$back.transformed$krige.var)
str(pev)
models <- c("clay_base_lm", "clay_base_lmm", "clay_best_lm",
          "clay_best_lmm", "carbon_base_lm", "carbon_base_lmm",
          "carbon_best_lm", "carbon_best_lmm", "ecec_base_lm",
          "ecec_base_lmm", "ecec_best_lm", "ecec_best_lmm")
stats <- list()
for (i in 1:dim(obs)[2]) {
  stats[[i]] <- cvStats(obs[, i], pred[, i], pev[, i])
}
nam <- names(stats[[1]])
stats <- t(matrix(unlist(stats), ncol = 12))
colnames(stats) <- nam
stats <- data.frame(stats)
Structure <- rep(c("LM", "LMM"), 6)
Model <- data.frame(Structure, stats[, c("me", "mae", "rmse", "srmse", "r2")])
rgroup <- c("CLAY (g kg$^{-1}$)", "SOC (g kg$^{-1}$)", "ECEC (mmol kg$^{-1}$)")
colheads <- c("Type", "ME", "MAE", "RMSE", "SRMSE", "AVE")
rowname <- rep(c("Base", "", "Best", ""), 3)
long_cap <- "Statistics of the LOO-CV of \\textit{base} and \\textit{best} multiple linear regression models (LM) and linear mixed models (LMM)."
foot <- "Statistics: mean error (ME), mean absolute error (MAE), root-mean-squared error (RMSE), scaled root-mean-squared error (SRMSE, unitless), and amount of variance explained (AVE, percent)."
file <-  paste(update_covars_dir, "cv-stats.tex", sep = "")
digits <- c(0, 2, 1, 1, 2, 1)
latex(Model, ctable = TRUE, n.rgroup = c(4, 4, 4), rgroup = rgroup, 
      file = file, label = "tab:cv-stats", table.env = TRUE, 
      cgroupTexCmd = NULL, rgroupTexCmd = NULL, rowname = rowname,
      colheads = colheads, na.blank = TRUE, caption = long_cap, where = NULL,
      size = "scriptsize", insert.bottom = foot, cdec = digits)

\{1}{l}{Model}
# SAVE DATA ####################################################################
ls()
# selected models and geodata
save(clay_base_geodata, clay_best_geodata, clay_base_lmm, clay_base_vario,
     clay_base_lm, clay_best_lm, clay_best_lmm, clay_best_vario,
     carbon_base_geodata, carbon_best_geodata, carbon_base_lmm, 
     carbon_base_vario, carbon_base_lm, carbon_best_lm, carbon_best_lmm,
     carbon_best_vario,
     ecec_base_geodata, ecec_best_geodata, ecec_base_lmm, ecec_base_vario, 
     ecec_base_lm, ecec_best_lm, ecec_best_lmm, ecec_best_vario,
     clay_sel, soc_sel, ecec_sel,
     file = "sm-dnos-phd-chap1-final-models.RData")

# other data
save(soil1, soil2, land1, land2, geo1, geo2, dem1, dem2, sat1, sat2,
     combs, preds, forms, drop, deltaR2,
     # base models
     clay_full, clay_vif, clay_both, clay_forward, clay_backward, 
     carbon_full, carbon_vif, carbon_both, carbon_forward, carbon_backward, 
     ecec_full, ecec_vif, ecec_both, ecec_forward, ecec_backward, 
     # base models stats
     clay_full_stats, clay_vif_stats, clay_both_stats, clay_forward_stats,
     clay_backward_stats, carbon_full_stats, carbon_vif_stats, carbon_both_stats, 
     carbon_forward_stats, carbon_backward_stats, ecec_full_stats, 
     ecec_vif_stats, ecec_both_stats, ecec_forward_stats, ecec_backward_stats,
     # cross-validation
     clay_base_lm_cv,    clay_base_lm_cv_cdf,
     clay_base_lmm_cv,   clay_base_lmm_cv_cdf,
     clay_best_lm_cv,    clay_best_lm_cv_cdf,
     clay_best_lmm_cv,   clay_best_lmm_cv_cdf,
     carbon_base_lm_cv,  carbon_base_lm_cv_cdf,
     carbon_base_lmm_cv, carbon_base_lmm_cv_cdf,
     carbon_best_lm_cv,  carbon_best_lm_cv_cdf,
     carbon_best_lmm_cv, carbon_best_lmm_cv_cdf,
     ecec_base_lm_cv,    ecec_base_lm_cv_cdf,
     ecec_base_lmm_cv,   ecec_base_lmm_cv_cdf,
     ecec_best_lm_cv,    ecec_best_lm_cv_cdf,
     ecec_best_lmm_cv,   ecec_best_lmm_cv_cdf,
     file = "sm-dnos-phd-chap1.RData")
# End!