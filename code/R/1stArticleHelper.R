# USER DEFINED FUNCTIONS #######################################################

# Convert pdf figures to png
# This is to reduce the size of the document submitted to the review process.
# The published document is prepared with the high-quality pdf images.
pdf2png <- 
  function(file, density = 300, quality = 100) {
    cmd <- paste("convert -density ", density, " ", file, ".pdf -quality ",
                 quality, " ", file, ".png", sep = "")
    lapply(cmd, system)
  }

# Effect of dropping one environmental covariate
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

# Create an object of class 'geodata' from an object of class 'lm'
makeGeodata <-
  function (y, model, data) {
    covars <- attr(terms(model), "term.labels")
    covar_col <- which(match(colnames(data), covars) != "NA")
    data_col  <- which(colnames(data) == y)
    coords.col <- c(ncol(data) - 1, ncol(data))
    geod <- as.geodata(obj = data, coords.col = coords.col,
                       data.col = data_col, covar.col = covar_col)
    return (geod)
  }

# Fit an empirical variogram using an object of class 'lm'
fitVariog <-
  function (y, model, data, lambda, breaks) {
    covars <- attr(terms(model), "term.labels")
    geod <- makeGeodata(y = y, model = model, data = data)
    trend <- formula(paste("~", paste(covars, collapse = " + ")))
    vario <- variog(geodata = geod, trend = trend, lambda = lambda,
                    breaks = breaks)
    return (vario)
  }

# Fit a linear mixed model with REML using an object of class 'lm'
fitREML <-
  function (y, model, data, ini.cov.pars, nugget, lambda) {
    covars <- attr(terms(model), "term.labels")
    covar_col <- which(match(colnames(data), covars) != "NA")
    data_col  <- which(colnames(data) == y)
    coords.col <- c(ncol(data) - 1, ncol(data))
    geod <- as.geodata(data, coords.col = coords.col,
                       data.col = data_col, covar.col = covar_col)
    trend <- formula(paste("~", paste(covars, collapse = " + ")))
    res <- likfit(geodata = geod, trend = trend, ini.cov.pars = ini.cov.pars,
                  nugget = nugget, lambda = lambda, lik.method = "REML")
    return (res)
  }
# cross-validation of geostatistical models
krigeCV <- 
  function (model, geodata, back = TRUE, simul.back = TRUE, n.sim = 20000, 
            digits, ...) {
    # Cross-validation of geostatistical models
    #
    # Args:
    #   model:      Object of class likGRF and/or variomodel. The geostatistical
    #               model to be cross-validated.
    #   geodata:    An object of the class "geodata".
    #   back:       Logical for indicating if predicted values should be
    #               returned in both transformed and back-transformed scales.
    #               Defaults to back = TRUE.
    #   simul.back: Logical for indicating if back-transformation should 
    #               be done using Monte Carlo
    #               simulations of a Gaussian process.
    #   n.sim:      Integer to specify the number of
    #               Monte Carlo simulations to be used in the back-
    #               transformation.
    #   digits:     Integer indicating the number of 
    #               decimal places to be used to round cross-validation results.
    #   ...:        Further arguments passed to the function xvalid().
    #    
    # Returns:
    #   Cross-validation.
    #
    # check arguments ##########################################################
    if (missing(geodata)) {
      stop("<geodata> is a mandatory argument")
    }
    if (missing(model)) {
      stop("<model> is a mandatory argument")
    }
    if (any(class(model) != c("likGRF", "variomodel"))) {
      stop("<model> should be of class likGRF and/or variomodel")
    }
    if (class(geodata) != "geodata") {
      stop("<geodata> should be of class geodata")
    }
    if (!missing(back)) {
      if (class(back) != "logical") {
        stop("<back> should be a logical value")
      }
    }
    if (!missing(simul.back)) {
      if (class(simul.back) != "logical") {
        stop("<simul.back> should be a logical value")
      }
    }
    if (!missing(n.sim)) {
      if (round(n.sim) != n.sim || n.sim <= 1) {
        stop("<n.sim> should be a positive integer value larger than 1")
      }
    }
    if (!missing(digits)) {
      if (length(digits) > 1) {
        stop("<digits> should be a single integer value")
      }
      if (round(digits) != digits) {
        stop("<digits> should be an integer value")
      }
    }
    # back-transform cross-validation results ##################################
    if (back) {
      cv <- list(cv.results = NA, back.transformed = NA)
      geodata0 <- geodata
      geodata$data <- geoR::BCtransform(geodata$data, model$lambda)$data
      model$parameters.summary["lambda", 2] <- 1
      cv$cv.results <- geoR::xvalid(geodata, model = model,
                                    #reestimate = reestimate,
                                    variog.obj = variog.obj, 
                                    #output.reestimate = output.reestimate,
                                    #locations.xvalid = locations.xvalid,
                                    #data.xvalid = data.xvalid, 
                                    ...)
      a <- attributes(cv$cv.results) 
      a$names[1:6] <- c("obs", "pred", "pev", "error", "zscore", "prob")
      attributes(cv$cv.results) <- a
      # back-transformation using simulations ################################
      if (simul.back == TRUE) {
        message(paste("back-transforming cross-validation predictions using ",
                      n.sim, " Gaussian simulations", sep = ""))
        back_cv <- invBoxCox(lambda = model$lambda, 
                             mean = cv[[1]][["pred"]], 
                             variance = cv[[1]][["pev"]], 
                             n.sim = n.sim, simul.back = TRUE, 
                             profile = FALSE)
        # naive back-transformation ##########################################
      } 
      else {
        message(paste("back-transforming cross-validation predictions",
                      sep = ""))
        back_cv <- invBoxCox(lambda = model$lambda, 
                             mean = cv[[1]][["pred"]], 
                             variance = cv[[1]][["pev"]],
                             simul.back = FALSE, profile = FALSE)
      }
      # prepare output #######################################################
      obs <- geodata0$data
      pred <- back_cv$mean
      pev <- back_cv$variance
      if (!missing(digits)) {
        pred <- round(pred, digits)
        pev <- round(pev, digits)
      }
      error <- obs - pred
      zscore <- error / sqrt(pev)
      prob  <- pnorm(obs, mean = pred, sd = sqrt(pev))
      cv$back.transformed <- data.frame(obs, pred, pev, error, zscore, prob)
    }  
    else {
      cv <- geoR::xvalid(geodata, model = model,
                         #reestimate = reestimate,
                         variog.obj = variog.obj, 
                         #output.reestimate = output.reestimate,
                         #locations.xvalid = locations.xvalid,
                         #data.xvalid = data.xvalid, 
                         ...)
    }
    return (cv)
  }
# End!

# leave-one-out cross-validation of linear models (lm)
looCV <-
  function (model, back = FALSE, simul.back = TRUE, original, lambda,
            n.sim = 20000, digits) {
    # Leave-one-out cross-validation of linear models (lm)
    #
    # Args:
    #   model:      Object of class lm. The linear model to be cross-validated.
    #   back:       Logical for back-transforming cross-validation results.
    #               Defaults to back = FALSE. Arguments original and lambda are
    #               mandatory if back = TRUE.
    #   simul.back: Logical for back-transforming cross-validation results using
    #               Monte Carlo simulations. Defaults to simul.back = TRUE.
    #   original:   Vector object with original values (numeric or integer). Used
    #               when the linear model was fitted using Box-Cox transformed 
    #               values and cross-validation results should be back-transformed.
    #   lambda:     Lambda value used for the Box-Cox transformation.
    #   n.sim:      Number (integer) of simulations that should be used for
    #               back-transforming cross-validation results.
    #   digits:     Integer indicating the number of decimal places to be used to
    #               round cross-validation results. If back = TRUE, only back-
    #               transformed values are rounded. A vector with two integer 
    #               values can be passed to round both cross-validation results in
    #               transformed and back-transformed scales.
    #
    # Returns:
    #   Data frame containing cross-validation results: observed (obs) and 
    #   predicted (pred) values, plus the prediction error variance (pev). When
    #   cross-validation results are back-trasnformed, then a list with 
    #   cross-validation results in the transformed and back-trasnformed scales 
    #   is returned.
    #
    # check arguments ##########################################################
    if (missing(model)) {
      stop("<model> is a mandatory argument")
    }
    if (class(model) != "lm") {
      stop("<model> should be of class lm")
    }
    if (back) {
      if (missing(original)) {
        stop("<original> is a mandatory argument to perform back-transformation")
      }
      if (!any(class(original) == c("numeric", "integer"))) {
        stop("<original> should be of class numeric or integer")
      }
      if (missing(lambda)) {
        stop("<lambda> is a mandatory argument to perform back-transformation")
      }
      if (!any(class(lambda) == c("numeric", "integer"))) {
        stop("<lambda> should be of class numeric or integer")
      }
      if (missing(n.sim)) {
        stop("<n.sim> is a mandatory argument to perform back-transformation")
      }
      if (round(n.sim) != n.sim || n.sim <= 1) {
        stop("<n.sim> should be a positive integer value larger than 1")
      }
    }
    if (!missing(digits)) {
      if (length(digits) > 2) {
        stop("<digits> should be a vector with one or two integer values")
      }
      if (round(digits) != digits) {
        stop("<digits> should be an integer value")
      }
    }
    # prepare data #############################################################
    data <- model.frame(model)
    idx  <- seq(1, dim(data)[1], 1)
    data <- cbind(data, idx)
    obs  <- numeric()
    pred <- numeric()
    pev  <- numeric()
    # sd.mean <- numeric()
    # sigma   <- numeric()
    # run cross-validation
    for (i in 1:length(idx)) {
      cv.fit     <- lm(formula(model), subset = c(idx != 1), data = data)
      cv.pred    <- predict(cv.fit, data[i, ], se.fit = TRUE)
      obs[i]     <- data[i, 1]
      pred[i]    <- cv.pred$fit
      # standard deviation of the predicted mean value of Y
      # sd.mean[i] <- cv.pred$se.fit
      # standard deviation of a predicted values of an individual observation
      # also known as prediction error variance
      pev[i] <- (sqrt(1 + c(cv.pred$se.fit/cv.pred$residual.scale) ^ 2) *
                   cv.pred$residual.scale) ^ 2
      # residual sum of squares
      # sigma[i]   <- cv.pred$residual.scale
    }
    # res <- data.frame(obs, pred, sigma, sd.mean, sd.pred)
    res <- data.frame(obs, pred, pev)
    # round values if required
    if (!missing(digits)) {
      if (back == FALSE) {
        res$pred <- round(res$pred, digits)
        res$pev  <- round(res$pev, digits)
      }
      if (back == TRUE && length(digits) == 2) {
        res$pred <- round(res$pred, digits[1])
        res$pev  <- round(res$pev, digits[1])
      }
    }
    if (back) {
      # prepare data
      cv <- list(cv.results = NA, back.transformed = NA)
      cv$cv.results <- res
      # back-transformation
      inv <- invBoxCox(mean = res$pred, variance = pev, lambda = lambda, 
                       simul.back = simul.back, n.sim = n.sim)
      obs  <- original
      # round values if required
      if (!missing(digits)) {
        if (length(digits) == 2) {
          pred <- round(inv$mean, digits[2])
          pev  <- round(inv$variance, digits[2])
        }
        if (length(digits) == 1) {
          pred <- round(inv$mean, digits)
          pev  <- round(inv$variance, digits)
        }
      } else {
        pred <- inv$mean
        pev  <- inv$variance
      }
      inv  <- data.frame(obs, pred, pev)
      cv$back.transformed <- inv
    }
    if (back) {
      return (cv)
    } else {
      return (res)  
    }
  }
# End!
#  back-transformation of Box-Cox transformed values
invBoxCox <- 
  function (mean, variance, lambda, simul.back = TRUE, n.sim = 20000, 
            profile = FALSE, digits, ...) {
    # Back-transformation of Box-Cox transformed values
    #
    # Args:
    #   mean:       Numeric or integer value to be back-transformed.
    #   variance:   Variance of the distribution of the numeric or integer
    #               value to be back-transformed.
    #   lambda:     Lambda value used for the Box-Cox transformation.
    #   simul.back: Logical for back-transforming cross-validation results using
    #               Monte Carlo simulations. Defaults to simul.back = TRUE.
    #   n.sim:      Integer specifying the number monte Carlo simulations that
    #               should be used for back-transforming cross-validation
    #               results. If profile = TRUE, then a vector vith at least
    #               four integer values should be passed to the function.
    #               If not, default values are used (100 to 100100).
    #   profile:    Logical for ploting a profile showing the evolution of the
    #               variance with the number of simulations.
    #   digits:     Integer indicating the number of decimal places to be used 
    #               to round back-transformed values. Defaults to no rounding.
    #   ...:        Further arguments to the plotting function. Used only if
    #               profile = TRUE.
    #
    # Returns:
    #   A plot (profile = TRUE) or a data frame with back-transformed values:
    #   mean, variance, distribution, median, and uncertainty.
    #    
    # check arguments ##########################################################
    if (missing(mean)) {
      stop("<mean> is a mandatory argument")
    }
    if (!missing(mean)) {
      if (!any(class(mean) == c("numeric", "integer"))) {
        stop("<mean> should be of class numeric or integer")
      }
    }
    if (missing(variance)) {
      stop("<variance> is a mandatory argument")
    }
    if (!missing(variance)) {
      if (!any(class(variance) == c("numeric", "integer"))) {
        stop("<variance> should be of class numeric or integer")
      }
    }
    if (length(mean) != length(variance)) {
      stop("<mean> and <variance> must have the same length")
    }
    if (missing(lambda)) {
      stop("<lambda> a mandatory argument")
    }
    if (length(lambda) > 1) {
      stop("<lambda> should be a single numeric or integer value")
    }
    if (profile == FALSE && missing(n.sim) == FALSE) {
      if (round(n.sim) != n.sim || n.sim <= 1) {
        stop("<n.sim> should be a positive integer value larger than 1")
      }
    }
    # plot profile of variance against the number of simulations ###############
    if (profile) {
      if (length(n.sim)  <= 3) {
        message("<n.sim> must be a vector with at least four integer values")
        message("default range will be used (100 - 100100)")
        n.sim <- seq(100, 100000, 10000)
      }
      # make simulations and print progress bar
      s <- pblapply(X = n.sim, function (X) {
        geoR::backtransform.moments(lambda = lambda, mean = mean, 
                                    variance = variance, simul.back = TRUE, 
                                    n.sim = X)
      })
      idx <- seq(1, length(n.sim))
      var <- data.frame(sapply(X = idx, function (X) {s <- s[[X]]$variance}))
      var.max <- apply(var, 1, max)
      var <- sapply(X = c(1:length(mean)), function (X) {var[X, ] / var.max[X]})
      var <- as.data.frame(var)
      # plot
      plot(n.sim, var[, 1], ylim = c(0, 1), type = "l", ylab = "scaled variance",
           xlab = "number of simulations", xaxp = c(0, max(n.sim), 4), ...)
      sapply(X = c(2:length(mean)), function(X) {lines(n.sim, var[, X], col = X)})
    }
    # back-transformation ######################################################
    else {
      s <- geoR::backtransform.moments(lambda = lambda, mean = mean, 
                                       variance = variance, 
                                       simul.back = simul.back, n.sim = n.sim)
      if (!missing(digits)) {
        s$mean <- round(s$mean, digits)
        s$variance <- round(s$variance, digits)
      }
      return(s)
    }
  }
# Plot lines of LMM
linesREML <-
  function (x, max.dist, scaled = FALSE, ...)
  {
    my.l <- list()
    if(missing(max.dist)){
      my.l$max.dist <- x$max.dist
      if (is.null(my.l$max.dist)) 
        stop("argument max.dist needed for this object")
    }
    else
      my.l$max.dist <- max.dist
    if (any(x$cov.model == c("matern","powered.exponential",
                             "cauchy", "gencauchy", "gneiting.matern"))) 
      my.l$kappa <- x$kappa
    else kappa <- NULL
    if (is.vector(x$cov.pars)) 
      my.l$sill.total <- x$nugget + x$cov.pars[1]
    else my.l$sill.total <- x$nugget + sum(x$cov.pars[, 1])
    my.l$nugget <- x$nugget
    my.l$cov.pars <- x$cov.pars
    my.l$cov.model <- x$cov.model
    if (scaled){
      if(is.vector(x$cov.model))
        my.l$cov.pars[1] <-  my.l$cov.pars[1]/my.l$sill.total
      else my.l$cov.pars[,1] <-  my.l$cov.cov.pars[,1]/my.l$sill.total
      my.l$sill.total <- 1
    }
    gamma.f <- function(x, my.l){
      if(any(my.l$cov.model == c("linear", "power")))
        return(my.l$nugget + my.l$cov.pars[1] * (x ^ my.l$cov.pars[2]))
      else
        return(my.l$sill.total -
                 cov.spatial(x, cov.model = my.l$cov.model,
                             kappa = my.l$kappa,
                             cov.pars = my.l$cov.pars))
    }
    line <- curve(gamma.f(x, my.l = my.l), from = 0, to = my.l$max.dist, ...)
    return(line)
  }

