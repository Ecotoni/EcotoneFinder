################################### Slope calculations and curves ##########################################

################################# Curve function without plotting #########################################

#' Adaptation of the curve function (without plot).
#'
#' @param expr The name of a function, or a call or an expression written as a
#'   function of x which will evaluate to an object of the same length as x.
#' @param from the range over which the function will be plotted (start).
#' @param to the range over which the function will be plotted (end).
#' @param n integer; the number of x values at which to evaluate.
#' @param add logical; if TRUE add to an already existing plot; if NA start a
#'   new plot taking the defaults for the limits and log-scaling of the x-axis
#'   from the previous plot. Taken as FALSE (with a warning if a different value
#'   is supplied) if no graphics device is open.
#' @param type plot type: see plot.default.
#' @param xname character string giving the name to be used for the x axis.
#' @param xlab labels and graphical parameters.
#' @param ylab labels and graphical parameters.
#' @param log labels and graphical parameters. See ‘Details’ for the
#'   interpretation of the default for log.
#' @param xlim NULL or a numeric vector of length 2; if non-NULL it provides the
#'   defaults for c(from, to) and, unless add = TRUE, selects the x-limits of
#'   the plot – see plot.window.
#' @param ... Additional graphical arguments.
#'
#' @return A vector containing the y values of the gaussian along the gradient.
#'
#' @details Silently used by SyntheticData and SyntheticDataSeries. Equivalent
#'   to the curve function of the graphics package. See the details of the curve
#'   function in graphics package for more details.
#'
#' @export
#'
#' @examples
#' gaussian <- function(x) a*exp(-(((x-b)^2)/2*(c^2)))
#' a <- 60
#' b <- 250
#' c <- 0.4
#' Curve=curveNoPlot(gaussian, from = 1, to = 500, n = 500)
#' Curve$y
#'
#'
curveNoPlot <- function (expr, from = NULL, to = NULL, n = 101, add = FALSE,
                         type = "l", xname = "x", xlab = xname, ylab = NULL, log = NULL,
                         xlim = NULL, ...)
{
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(suppressWarnings(graphics::par(oldpar)))
  sexpr <- substitute(expr)
  if (is.name(sexpr)) {
    expr <- call(as.character(sexpr), as.name(xname))
  }
  else {
    if (!((is.call(sexpr) || is.expression(sexpr)) && xname %in%
          all.vars(sexpr)))
      stop(gettextf("'expr' must be a function, or a call or an expression containing '%s'",
                    xname), domain = NA)
    expr <- sexpr
  }
  if (grDevices::dev.cur() == 1L && !identical(add, FALSE)) {
    warning("'add' will be ignored as there is no existing plot")
    add <- FALSE
  }
  addF <- identical(add, FALSE)
  if (is.null(ylab))
    ylab <- deparse(expr)
  if (is.null(from) || is.null(to)) {
    xl <- if (!is.null(xlim))
      xlim
    else if (!addF) {
      pu <- graphics::par("usr")[1L:2L]
      if (graphics::par("xaxs") == "r")
        pu <- grDevices::extendrange(pu, f = -1/27)
      if (graphics::par("xlog"))
        10^pu
      else pu
    }
    else c(0, 1)
    if (is.null(from))
      from <- xl[1L]
    if (is.null(to))
      to <- xl[2L]
  }
  lg <- if (length(log))
    log
  else if (!addF && graphics::par("xlog"))
    "x"
  else ""
  if (length(lg) == 0)
    lg <- ""
  if (grepl("x", lg, fixed = TRUE)) {
    if (from <= 0 || to <= 0)
      stop("'from' and 'to' must be > 0 with log=\"x\"")
    x <- exp(seq.int(log(from), log(to), length.out = n))
  }
  else x <- seq.int(from, to, length.out = n)
  ll <- list(x = x)
  names(ll) <- xname
  y <- eval(expr, envir = ll, enclos = parent.frame())
  if (length(y) != length(x))
    stop("'expr' did not evaluate to an object of length 'n'")
  invisible(list(x = x, y = y))
}


################################# Slope calculations #####################################################

#'Method to calculate the derivative of irregular functions:
#'
#'@param ecotonefinder A list containing elements named in the same way than
#'  EcotoneFinder function outcomes, or a nammed list or numercic vector,
#'  if method is set to "list".
#'@param method The name of the method for which the slopes should be
#'  calculated. Correspond to the names of the list produced by the
#'  EcotoneFinder function, or can be set to "list" if another object
#'  is provided. See details.
#'@param dist When method = "list", if NULL (default), the derivatives are
#'  computed on an interval corresponding to the number of rows of the
#'  provided object. A numeric vector of the same length as the obhect can
#'  be provided instead.
#'@param window Must be an odd number. The interval to be used for slope
#'  calculation. The bigger the window, the more averaged the slope will be.
#'@param axis.number If "dca" is chosen, indicate the number of axis over which
#'  to calculate the slope (first axis, first and second axis,...)
#'@param groups If any clustering method is chosen, corresponds to the index of
#'  the cluster for which the slope shold be calculated. If "all", the slope
#'  will be calculated for all the clusters.
#'@param diversity If "diversity" is chosen in the method argument, define the
#'  diversity index for which to calculate the slope. "all" can be chosen.
#'
#'@details Slope calculations are done by moving window analysis. The width of
#'  the windows is defined by the window argument. For each window, the result
#'  of slope coefficient of a linear model (lm function of the stat package) is
#'  stored and used to draw the general slope along the gradient. The bigger the
#'  window, the more points will be used to compute the linear models, meaning
#'  the obtained slopes will be smoother. This also results in the addition of
#'  NAs at the ends of the gradient.
#'
#'  The first axis of DCA has been used as a beta-diversity index, and its
#'  derivative as a method to locate ecotones (see Brownstein et al., 2013). The
#'  Slope function provide the possibility of computing the slope of the other
#'  axis, to avoid the loss of information induced by the reduction of the
#'  dimentionality of the original data. Similarly, the slopes of the fuzzy
#'  clusters can be used to pinpoint the transitions between them. The value of
#'  the slopes can be an indicator of the relative sharpness of the transition
#'  area. Particularly, as the memberships of the fuzzy clusters range between 0
#'  and 1, these values can readily be compared between studies and datasets.
#'  These values vary depending on the window width and can be very sensible to
#'  noise in the original data. A reliable method to mathematically identify
#'  breaks is still needed and careful interpretation by the user is still
#'  required.
#'
#'  The function is primarily designed to work with the output produced by
#'  the Ecotonefinder function, but can accept other lists or numeric
#'  vectors, if method is set to "custom". The provided list may be named or
#'  not, but must consist of numeric dataframes or matrices. The derivative
#'  of all columns will be computed. By default, the derivatives will be
#'  computed on an interval c(1:nrow(x)), where `x` is the provided matrix
#'  (c(1:length(x)), if x is a numeric vector). If the default behavior is
#'  not appropriate, a numeric vector can be specified through the "dist"
#'  argument, and will be used to calibrate the windows.
#'
#'
#'@return A list of dataframes containing the slope values for the specified
#'  methods and the original data.
#'
#'@export
#'
#' @examples
#'  #### Artificial dataset:
#'  SyntheticTrial <- SyntheticData(SpeciesNum = 21, CommunityNum = 3,
#'                                  SpCo = NULL ,Length = 500,
#'                                  Parameters = list(a=rep(60, 3),
#'                                                    b=c(0,250,500),
#'                                                    c=rep(0.015,3)),
#'                                  pal = c("#008585", "#FBF2C4", "#C7522B"))
#'
#'  ## Analyses:
#'  SyntheticEcoFinder <- EcotoneFinder(data = SyntheticTrial[,-1],
#'                                      dist = SyntheticTrial$Distance,
#'                                      method = "all", groups = 3,
#'                                      standardize = "hellinger", diversity = "all")
#'
#'  ## Slope calculation:
#'  SyntheticSlope <- Slope(SyntheticEcoFinder, method = "all", axis.number = 2,
#'                          diversity = "all")
#'
#'  ## With another object type:
#'  # Numeric vector:
#'  x <- sort(rbeta(100, 30, 0.4))
#'  x_slope <- Slope(x, method = "custom", dist = NULL)
#'  x_slope <- Slope(x, method = "custom", dist = seq(1, 5, length.out = 100))
#'
#'  # list:
#'  l <- list("Beta" = data.frame("A" = sort(rbeta(100, 30, 0.4)),
#'                                "B" = sort(rbeta(100, 50, 0.5)),
#'                                "C" = sort(rbeta(100, 20, 0.1))),
#'            "NegBin" = data.frame("A" = sort(rnbinom(100, 10, 0.1)),
#'                                  "B" = sort(rnbinom(100, 50, 0.7)),
#'                                  "C" = sort(rnbinom(100, 70, 0.8))))
#'  l_slope <- Slope(l, method = "custom", dist = NULL)
#'
Slope <- function (ecotonefinder, method = c("dca", "fanny", "vegclust", "cmeans", "diversity", "custom", "all"), window = 3, axis.number = 1,
                   groups = as.numeric(ecotonefinder["groups"]), diversity = c("shannon", "richness", "expShannon", "pielou", "all"), dist = NULL)
{
  # Checks:
  if (is.null(ecotonefinder)) {
    stop("ecotonefinder must be provided")
  }
  if (any(is.na(names(ecotonefinder[method]))) && all(method != "all") && all(method != "custom")) {
    stop("ecotonefinder must contains the elements required in 'method'.")
  }
  list_output <- list()
  if (any(method == "diversity") || any(method == "all")) {
    list_diversity <- list()
  }
  if ((window %% 2) == 0) {
    stop("window needs to be odd-numbered")
  }
  list_output[["window"]] <- window
  endcut <- floor(window/2)
  mod <- list(axis_score ~ sub_dist)
  if (all(method != "custom")) {
    list_output[["distance"]] <- ecotonefinder$distance
    list_output[["data"]] <- ecotonefinder$data
    list_output[["groups"]] <- groups
    for (k in 1:(length(ecotonefinder$distance) - (window - 1))) {
      j <- k + (window - 1)
      sub_dist <- ecotonefinder$distance[k:j]
      if (any(method == "dca") || any(method == "all")) {
        axis_score <- ecotonefinder$dca$rproj[k:j, 1:axis.number]
        test_lm <- mapply(stats::lm, formula = mod)
        if (axis.number == "1") {
          list_output[["dca_slope"]][[endcut + k]] <- test_lm[,1]$coefficients[2]
        } else { list_output[["dca_slope"]][[endcut + k]] <- test_lm[,1]$coefficients[2,] }
      }
      if (any(method == "fanny") || any(method == "all")) {
        axis_score <- ecotonefinder$fanny$membership[k:j, 1:groups]
        test_lm <- mapply(stats::lm, formula = mod)
        if (groups == "1") {
          list_output[["fanny_slope"]][[endcut + k]] <- test_lm[,1]$coefficients[2]
        } else { list_output[["fanny_slope"]][[endcut + k]] <- test_lm[,1]$coefficients[2,] }
      }
      if (any(method == "vegclust") || any(method == "all")) {
        axis_score <- as.matrix(ecotonefinder$vegclust$memb[k:j, 1:groups])
        test_lm <- mapply(stats::lm, formula = mod)
        if (groups == "1") {
          list_output[["vegclust_slope"]][[endcut + k]] <- test_lm[,1]$coefficients[2]
        } else { list_output[["vegclust_slope"]][[endcut + k]] <- test_lm[,1]$coefficients[2,] }
      }
      if (any(method == "cmeans") || any(method == "all")) {
        axis_score <- as.matrix(ecotonefinder$cmeans$membership[k:j, 1:groups])
        test_lm <- mapply(stats::lm, formula = mod)
        if (groups == "1") {
          list_output[["cmeans_slope"]][[endcut + k]] <- test_lm[,1]$coefficients[2]
        } else { list_output[["cmeans_slope"]][[endcut + k]] <- test_lm[,1]$coefficients[2,] }
      }
      if (any(method == "diversity") || any(method == "all")) {
        if (any(diversity == "richness") || any(diversity == "all")) {
          axis_score <- ecotonefinder$diversity$SpeciesRichness[k:j]
          test_lm <- mapply(stats::lm, formula = mod)
          list_diversity[["SpeciesRichness_slope"]][[endcut + k]] <- test_lm[,1]$coefficients[2]
        }
        if (any(diversity == "shannon") || any(diversity == "all")) {
          axis_score <- ecotonefinder$diversity$Shannon[k:j]
          test_lm <- mapply(stats::lm, formula = mod)
          list_diversity[["Shannon_slope"]][[endcut + k]] <- test_lm[,1]$coefficients[2]
        }
        if (any(diversity == "expShannon") || any(diversity == "all")) {
          axis_score <- ecotonefinder$diversity$ExpShannon[k:j]
          test_lm <- mapply(stats::lm, formula = mod)
          list_diversity[["ExpShannon_slope"]][[endcut + k]] <- test_lm[,1]$coefficients[2]
        }
        if (any(diversity == "pielou") || any(diversity == "all")) {
          axis_score <- ecotonefinder$diversity$Pielou[k:j]
          test_lm <- mapply(stats::lm, formula = mod)
          list_diversity[["Pielou_slope"]][[endcut + k]] <- test_lm[,1]$coefficients[2]
        }
      }
    }
  }
  if (any(method == "custom")) {
    if (is.list(ecotonefinder) == FALSE && is.numeric(ecotonefinder) == FALSE) {
      stop("ecotonefinder must be either a list object or a numeric vector/matrix")
    }
    if (is.list(ecotonefinder) == TRUE) {
      loop <- length(ecotonefinder)
      list <- TRUE
    }
    if (is.numeric(ecotonefinder) == TRUE) {
      loop <- 1
      list <- FALSE
    }
    for (i in 1:loop) {
      dat <- NULL
      if (list == TRUE) { dat <- ecotonefinder[[i]] } else { dat <- ecotonefinder }
      if (is.null(dist) == FALSE && length(dist) != nrow(as.matrix(dat))) {
        stop("dist must be either NULL or the same lenght as the provided data")
      }
      if (is.null(dist)) { dist <- c(1:nrow(as.matrix(dat))) }
      for (k in 1:(length(dist) - (window - 1))) {
        j <- k + (window - 1)
        sub_dist <- dist[k:j]
        if (dim(as.matrix(dat))[2] == 1) { axis_score <- dat[k:j] } else { axis_score <- as.matrix(dat[k:j,]) }
        test_lm <- mapply(stats::lm, formula = mod)
        if (dim(as.matrix(dat))[2] == 1) {
          list_output[[paste(i)]][[endcut + k]] <- test_lm[,1]$coefficients[2]
        } else {
          list_output[[paste(i)]][[endcut + k]] <- test_lm[,1]$coefficients[2,]
        }
      }
      list_output[[paste(i)]] <- purrr::reduce(list_output[[paste(i)]], rbind)
    }
    if(is.null(names(ecotonefinder)) == FALSE) {
      names(list_output) <- c("window", paste("Slope_",names(ecotonefinder), sep=""))
    } else { names(list_output) <- c("window", paste("Slope_",c(1:loop), sep="")) }
  }
  if (any(method == "dca") || any(method == "all")) {
    if (axis.number == "1") {
      list_output[["dca_slope"]] <- list_output$dca_slope
    } else { list_output[["dca_slope"]] <- purrr::reduce(list_output[["dca_slope"]], rbind) }
  }
  if (any(method == "fanny") || any(method == "all")) {
    if (groups == "1") {
      list_output[["fanny_slope"]] <- list_output$fanny_slope
    } else { list_output[["fanny_slope"]] <- purrr::reduce(list_output[["fanny_slope"]], rbind) }
  }
  if (any(method == "vegclust") || any(method == "all")) {
    if (groups == "1") {
      list_output[["vegclust_slope"]] <- list_output$vegclust_slope
    } else { list_output[["vegclust_slope"]] <- purrr::reduce(list_output[["vegclust_slope"]], rbind) }
  }
  if (any(method == "cmeans") || any(method == "all")) {
    if (groups == "1") {
      list_output[["cmeans_slope"]] <- list_output$cmeans_slope
    } else { list_output[["cmeans_slope"]] <- purrr::reduce(list_output[["cmeans_slope"]], rbind) }
  }
  if (any(method == "diversity") || method == "all") {
    if (any(diversity == "richness") || any(diversity == "all")) { list_diversity[["SpeciesRichness_slope"]] <- list_diversity$SpeciesRichness_slope }
    if (any(diversity == "shannon") || any(diversity == "all")) { list_diversity[["Shannon_slope"]] <- list_diversity$Shannon_slope }
    if (any(diversity == "expShannon") || any(diversity == "all")) { list_diversity[["ExpShannon_slope"]] <- list_diversity$ExpShannon_slope }
    if (any(diversity == "pielou") || any(diversity == "all")) { list_diversity[["Pielou_slope"]] <- list_diversity$Pielou_slope }
    list_output[["diversity_slope"]] <- list_diversity
  }
  return(list_output)
}
