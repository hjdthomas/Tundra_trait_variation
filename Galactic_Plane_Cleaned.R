#The following is code for the manuscript
#"Global plant trait relationships extend to the climatic extremes of the tundra biome"
#Part 1: Trait space
#20th Jan 2020


#Load packages####
library(data.table)
library(plyr)
library(dplyr)
library(ggplot2)
library(ggbiplot)
library(gridExtra)
library(raster)
library(rgdal)
library(tidyr)
library(raster)


#Create custom bi-plot functions (for rotation)
#ggbiplot function####
ggbiplot2<-function (pcobj, choices = 1:2, scale = 1, pc.biplot = TRUE, 
                     obs.scale = 1 - scale, var.scale = scale, groups = NULL, 
                     ellipse = TRUE, ellipse.prob = 0.68, labels = NULL, labels.size = 3, 
                     alpha = 1, var.axes = TRUE, circle = FALSE, circle.prob = 0.69, 
                     varname.size = 3, varname.adjust = 1.5, varname.abbrev = FALSE, 
                     color = muted("red"), # <- add new arguments to the function
                     linetype = "solid",
                     alpha_arrow = 1) 
{
  library(ggplot2)
  library(plyr)
  library(scales)
  library(grid)
  stopifnot(length(choices) == 2)
  if (inherits(pcobj, "prcomp")) {
    nobs.factor <- sqrt(nrow(pcobj$x) - 1)
    d <- pcobj$sdev
    u <- sweep(pcobj$x, 2, 1/(d * nobs.factor), FUN = "*")
    v <- pcobj$rotation
  }
  else if (inherits(pcobj, "princomp")) {
    nobs.factor <- sqrt(pcobj$n.obs)
    d <- pcobj$sdev
    u <- sweep(pcobj$scores, 2, 1/(d * nobs.factor), FUN = "*")
    v <- pcobj$loadings
  }
  else if (inherits(pcobj, "PCA")) {
    nobs.factor <- sqrt(nrow(pcobj$call$X))
    d <- unlist(sqrt(pcobj$eig)[1])
    u <- sweep(pcobj$ind$coord, 2, 1/(d * nobs.factor), FUN = "*")
    v <- sweep(pcobj$var$coord, 2, sqrt(pcobj$eig[1:ncol(pcobj$var$coord), 
                                                  1]), FUN = "/")
  }
  else if (inherits(pcobj, "lda")) {
    nobs.factor <- sqrt(pcobj$N)
    d <- pcobj$svd
    u <- predict(pcobj)$x/nobs.factor
    v <- pcobj$scaling
    d.total <- sum(d^2)
  }
  else {
    stop("Expected a object of class prcomp, princomp, PCA, or lda")
  }
  choices <- pmin(choices, ncol(u))
  df.u <- as.data.frame(sweep(u[, choices], 2, d[choices]^obs.scale, 
                              FUN = "*"))
  v <- sweep(v, 2, d^var.scale, FUN = "*")
  df.v <- as.data.frame(v[, choices])
  names(df.u) <- c("xvar", "yvar")
  names(df.v) <- names(df.u)
  if (pc.biplot) {
    df.u <- df.u * nobs.factor
  }
  r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)
  v.scale <- rowSums(v^2)
  df.v <- r * df.v/sqrt(max(v.scale))
  if (obs.scale == 0) {
    u.axis.labs <- paste("standardized PC", choices, sep = "")
  }
  else {
    u.axis.labs <- paste("PC", choices, sep = "")
  }
  u.axis.labs <- paste(u.axis.labs, sprintf("(%0.1f%% explained var.)", 
                                            100 * pcobj$sdev[choices]^2/sum(pcobj$sdev^2)))
  if (!is.null(labels)) {
    df.u$labels <- labels
  }
  if (!is.null(groups)) {
    df.u$groups <- groups
  }
  if (varname.abbrev) {
    df.v$varname <- abbreviate(rownames(v))
  }
  else {
    df.v$varname <- rownames(v)
  }
  df.v$angle <- with(df.v, (180/pi) * atan(yvar/xvar))
  df.v$hjust = with(df.v, (1 - varname.adjust * sign(xvar))/2)
  g <- ggplot(data = df.u, aes(x = xvar, y = yvar)) + xlab(u.axis.labs[1]) + 
    ylab(u.axis.labs[2]) + coord_equal()
  if (var.axes) {
    if (circle) {
      theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, 
                                                length = 50))
      circle <- data.frame(xvar = r * cos(theta), yvar = r * 
                             sin(theta))
      g <- g + geom_path(data = circle, color = muted("white"), 
                         size = 1/2, alpha = 1/3)
    }
    g <- g + geom_segment(data = df.v, aes(x = 0, y = 0, 
                                           xend = xvar, yend = yvar*2), arrow = arrow(length = unit(1/2, 
                                                                                                    "picas")), color = color, linetype = linetype, alpha = alpha_arrow)
  }
  if (!is.null(df.u$labels)) {
    if (!is.null(df.u$groups)) {
      g <- g + geom_text(aes(label = labels, color = groups), 
                         size = labels.size)
    }
    else {
      g <- g + geom_text(aes(label = labels), size = labels.size)
    }
  }
  else {
    if (!is.null(df.u$groups)) {
      g <- g + geom_point(aes(color = groups), alpha = alpha)
    }
    else {
      g <- g + geom_point(alpha = alpha)
    }
  }
  if (!is.null(df.u$groups) && ellipse) {
    theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
    circle <- cbind(cos(theta), sin(theta))
    ell <- ddply(df.u, "groups", function(x) {
      if (nrow(x) <= 2) {
        return(NULL)
      }
      sigma <- var(cbind(x$xvar, x$yvar))
      mu <- c(mean(x$xvar), mean(x$yvar))
      ed <- sqrt(qchisq(ellipse.prob, df = 2))
      data.frame(sweep(circle %*% chol(sigma) * ed, 2, 
                       mu, FUN = "+"), groups = x$groups[1])
    })
    names(ell)[1:2] <- c("xvar", "yvar")
    g <- g + geom_path(data = ell, aes(color = groups, group = groups))
  }
  if (var.axes) {
    g <- g + geom_text(data = df.v, aes(label = varname, 
                                        x = xvar, y = yvar*2, angle = angle, hjust = hjust), 
                       color = "black", size = varname.size)
  }
  return(g)
}

#ggbipplot reverse y####
ggbiplot_reverse_y<-function (pcobj, choices = 1:2, scale = 1, pc.biplot = TRUE, 
                              obs.scale = 1 - scale, var.scale = scale, groups = NULL, 
                              ellipse = TRUE, ellipse.prob = 0.68, labels = NULL, labels.size = 3, 
                              alpha = 1, var.axes = TRUE, circle = FALSE, circle.prob = 0.69, 
                              varname.size = 3, varname.adjust = 1.5, varname.abbrev = FALSE, 
                              color = muted("red"), # <- add new arguments to the function
                              linetype = "solid",
                              alpha_arrow = 1) 
{
  library(ggplot2)
  library(plyr)
  library(scales)
  library(grid)
  stopifnot(length(choices) == 2)
  if (inherits(pcobj, "prcomp")) {
    nobs.factor <- sqrt(nrow(pcobj$x) - 1)
    d <- pcobj$sdev
    u <- sweep(pcobj$x, 2, 1/(d * nobs.factor), FUN = "*")
    v <- pcobj$rotation
  }
  else if (inherits(pcobj, "princomp")) {
    nobs.factor <- sqrt(pcobj$n.obs)
    d <- pcobj$sdev
    u <- sweep(pcobj$scores, 2, 1/(d * nobs.factor), FUN = "*")
    v <- pcobj$loadings
  }
  else if (inherits(pcobj, "PCA")) {
    nobs.factor <- sqrt(nrow(pcobj$call$X))
    d <- unlist(sqrt(pcobj$eig)[1])
    u <- sweep(pcobj$ind$coord, 2, 1/(d * nobs.factor), FUN = "*")
    v <- sweep(pcobj$var$coord, 2, sqrt(pcobj$eig[1:ncol(pcobj$var$coord), 
                                                  1]), FUN = "/")
  }
  else if (inherits(pcobj, "lda")) {
    nobs.factor <- sqrt(pcobj$N)
    d <- pcobj$svd
    u <- predict(pcobj)$x/nobs.factor
    v <- pcobj$scaling
    d.total <- sum(d^2)
  }
  else {
    stop("Expected a object of class prcomp, princomp, PCA, or lda")
  }
  choices <- pmin(choices, ncol(u))
  df.u <- as.data.frame(sweep(u[, choices], 2, d[choices]^obs.scale, 
                              FUN = "*"))
  v <- sweep(v, 2, d^var.scale, FUN = "*")
  df.v <- as.data.frame(v[, choices])
  names(df.u) <- c("xvar", "yvar")
  names(df.v) <- names(df.u)
  if (pc.biplot) {
    df.u <- df.u * nobs.factor
  }
  r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)
  v.scale <- rowSums(v^2)
  df.v <- r * df.v/sqrt(max(v.scale))
  if (obs.scale == 0) {
    u.axis.labs <- paste("standardized PC", choices, sep = "")
  }
  else {
    u.axis.labs <- paste("PC", choices, sep = "")
  }
  u.axis.labs <- paste(u.axis.labs, sprintf("(%0.1f%% explained var.)", 
                                            100 * pcobj$sdev[choices]^2/sum(pcobj$sdev^2)))
  if (!is.null(labels)) {
    df.u$labels <- labels
  }
  if (!is.null(groups)) {
    df.u$groups <- groups
  }
  if (varname.abbrev) {
    df.v$varname <- abbreviate(rownames(v))
  }
  else {
    df.v$varname <- rownames(v)
  }
  df.v$angle <- with(df.v, (225/pi) * atan(-yvar/xvar))
  df.v$hjust = with(df.v, (1 - varname.adjust * sign(xvar))/2)
  g <- ggplot(data = df.u, aes(x = xvar, y = yvar)) + xlab(u.axis.labs[1]) + 
    ylab(u.axis.labs[2]) + coord_equal()
  if (var.axes) {
    if (circle) {
      theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, 
                                                length = 50))
      circle <- data.frame(xvar = r * cos(theta), yvar = r * 
                             sin(theta))
      g <- g + geom_path(data = circle, color = muted("white"), 
                         size = 1/2, alpha = 1/3)
    }
    g <- g + geom_segment(data = df.v, aes(x = 0, y = 0, 
                                           xend = xvar, yend = yvar*2), arrow = arrow(length = unit(1/2, 
                                                                                                    "picas")), color = color, linetype = linetype, alpha = alpha_arrow)
  }
  if (!is.null(df.u$labels)) {
    if (!is.null(df.u$groups)) {
      g <- g + geom_text(aes(label = labels, color = groups), 
                         size = labels.size)
    }
    else {
      g <- g + geom_text(aes(label = labels), size = labels.size)
    }
  }
  else {
    if (!is.null(df.u$groups)) {
      g <- g + geom_point(aes(color = groups), alpha = alpha)
    }
    else {
      g <- g + geom_point(alpha = alpha)
    }
  }
  if (!is.null(df.u$groups) && ellipse) {
    theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
    circle <- cbind(cos(theta), sin(theta))
    ell <- ddply(df.u, "groups", function(x) {
      if (nrow(x) <= 2) {
        return(NULL)
      }
      sigma <- var(cbind(x$xvar, x$yvar))
      mu <- c(mean(x$xvar), mean(x$yvar))
      ed <- sqrt(qchisq(ellipse.prob, df = 2))
      data.frame(sweep(circle %*% chol(sigma) * ed, 2, 
                       mu, FUN = "+"), groups = x$groups[1])
    })
    names(ell)[1:2] <- c("xvar", "yvar")
    g <- g + geom_path(data = ell, aes(color = groups, group = groups))
  }
  if (var.axes) {
    g <- g + geom_text(data = df.v, aes(label = varname, 
                                        x = xvar, y = yvar*2, angle = angle, hjust = hjust), 
                       color = "black", size = varname.size)
  }
  return(g)
}

#ggbiplot function 2####
ggbiplot_reverse<-function (pcobj, choices = 1:2, scale = 1, pc.biplot = TRUE, 
                            obs.scale = 1 - scale, var.scale = scale, groups = NULL, 
                            ellipse = TRUE, ellipse.prob = 0.68, labels = NULL, labels.size = 3, 
                            alpha = 1, var.axes = TRUE, circle = FALSE, circle.prob = 0.69, 
                            varname.size = 3, varname.adjust = 1.5, varname.abbrev = FALSE, 
                            color = muted("red"), # <- add new arguments to the function
                            linetype = "solid",
                            alpha_arrow = 1) 
{
  library(ggplot2)
  library(plyr)
  library(scales)
  library(grid)
  stopifnot(length(choices) == 2)
  if (inherits(pcobj, "prcomp")) {
    nobs.factor <- sqrt(nrow(pcobj$x) - 1)
    d <- pcobj$sdev
    u <- sweep(pcobj$x, 2, 1/(d * nobs.factor), FUN = "*")
    v <- pcobj$rotation
  }
  else if (inherits(pcobj, "princomp")) {
    nobs.factor <- sqrt(pcobj$n.obs)
    d <- pcobj$sdev
    u <- sweep(pcobj$scores, 2, 1/(d * nobs.factor), FUN = "*")
    v <- pcobj$loadings
  }
  else if (inherits(pcobj, "PCA")) {
    nobs.factor <- sqrt(nrow(pcobj$call$X))
    d <- unlist(sqrt(pcobj$eig)[1])
    u <- sweep(pcobj$ind$coord, 2, 1/(d * nobs.factor), FUN = "*")
    v <- sweep(pcobj$var$coord, 2, sqrt(pcobj$eig[1:ncol(pcobj$var$coord), 
                                                  1]), FUN = "/")
  }
  else if (inherits(pcobj, "lda")) {
    nobs.factor <- sqrt(pcobj$N)
    d <- pcobj$svd
    u <- predict(pcobj)$x/nobs.factor
    v <- pcobj$scaling
    d.total <- sum(d^2)
  }
  else {
    stop("Expected a object of class prcomp, princomp, PCA, or lda")
  }
  choices <- pmin(choices, ncol(u))
  df.u <- as.data.frame(sweep(u[, choices], 2, d[choices]^obs.scale, 
                              FUN = "*"))
  v <- sweep(v, 2, d^var.scale, FUN = "*")
  df.v <- as.data.frame(v[, choices])
  names(df.u) <- c("xvar", "yvar")
  names(df.v) <- names(df.u)
  if (pc.biplot) {
    df.u <- df.u * nobs.factor
  }
  r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)
  v.scale <- rowSums(v^2)
  df.v <- r * df.v/sqrt(max(v.scale))
  if (obs.scale == 0) {
    u.axis.labs <- paste("standardized PC", choices, sep = "")
  }
  else {
    u.axis.labs <- paste("PC", choices, sep = "")
  }
  u.axis.labs <- paste(u.axis.labs, sprintf("(%0.1f%% explained var.)", 
                                            100 * pcobj$sdev[choices]^2/sum(pcobj$sdev^2)))
  if (!is.null(labels)) {
    df.u$labels <- labels
  }
  if (!is.null(groups)) {
    df.u$groups <- groups
  }
  if (varname.abbrev) {
    df.v$varname <- abbreviate(rownames(v))
  }
  else {
    df.v$varname <- rownames(v)
  }
  df.v$angle <- with(df.v, (270/pi) * atan(yvar/-xvar))
  df.v$hjust = with(df.v, (1 - varname.adjust * sign(-xvar))/2)
  g <- ggplot(data = df.u, aes(x = xvar, y = yvar)) + xlab(u.axis.labs[1]) + 
    ylab(u.axis.labs[2]) + coord_equal()
  if (var.axes) {
    if (circle) {
      theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, 
                                                length = 50))
      circle <- data.frame(xvar = r * cos(theta), yvar = r * 
                             sin(theta))
      g <- g + geom_path(data = circle, color = muted("white"), 
                         size = 1/2, alpha = 1/3)
    }
    g <- g + geom_segment(data = df.v, aes(x = 0, y = 0, 
                                           xend = xvar, yend = yvar*2), arrow = arrow(length = unit(1/2, 
                                                                                                    "picas")), color = color, linetype = linetype, alpha = alpha_arrow)
  }
  if (!is.null(df.u$labels)) {
    if (!is.null(df.u$groups)) {
      g <- g + geom_text(aes(label = labels, color = groups), 
                         size = labels.size)
    }
    else {
      g <- g + geom_text(aes(label = labels), size = labels.size)
    }
  }
  else {
    if (!is.null(df.u$groups)) {
      g <- g + geom_point(aes(color = groups), alpha = alpha)
    }
    else {
      g <- g + geom_point(alpha = alpha)
    }
  }
  if (!is.null(df.u$groups) && ellipse) {
    theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
    circle <- cbind(cos(theta), sin(theta))
    ell <- ddply(df.u, "groups", function(x) {
      if (nrow(x) <= 2) {
        return(NULL)
      }
      sigma <- var(cbind(x$xvar, x$yvar))
      mu <- c(mean(x$xvar), mean(x$yvar))
      ed <- sqrt(qchisq(ellipse.prob, df = 2))
      data.frame(sweep(circle %*% chol(sigma) * ed, 2, 
                       mu, FUN = "+"), groups = x$groups[1])
    })
    names(ell)[1:2] <- c("xvar", "yvar")
    g <- g + geom_path(data = ell, aes(color = groups, group = groups))
  }
  if (var.axes) {
    g <- g + geom_text(data = df.v, aes(label = varname, 
                                        x = xvar, y = yvar*2, angle = angle, hjust = hjust), 
                       color = "black", size = varname.size)
  }
  return(g)
}

#ggbipplot reverse xy####
ggbiplot_reverse_xy<-function (pcobj, choices = 1:2, scale = 1, pc.biplot = TRUE, 
                               obs.scale = 1 - scale, var.scale = scale, groups = NULL, 
                               ellipse = TRUE, ellipse.prob = 0.68, labels = NULL, labels.size = 3, 
                               alpha = 1, var.axes = TRUE, circle = FALSE, circle.prob = 0.69, 
                               varname.size = 3, varname.adjust = 1.5, varname.abbrev = FALSE, 
                               color = muted("red"), # <- add new arguments to the function
                               linetype = "solid",
                               alpha_arrow = 1) 
{
  library(ggplot2)
  library(plyr)
  library(scales)
  library(grid)
  stopifnot(length(choices) == 2)
  if (inherits(pcobj, "prcomp")) {
    nobs.factor <- sqrt(nrow(pcobj$x) - 1)
    d <- pcobj$sdev
    u <- sweep(pcobj$x, 2, 1/(d * nobs.factor), FUN = "*")
    v <- pcobj$rotation
  }
  else if (inherits(pcobj, "princomp")) {
    nobs.factor <- sqrt(pcobj$n.obs)
    d <- pcobj$sdev
    u <- sweep(pcobj$scores, 2, 1/(d * nobs.factor), FUN = "*")
    v <- pcobj$loadings
  }
  else if (inherits(pcobj, "PCA")) {
    nobs.factor <- sqrt(nrow(pcobj$call$X))
    d <- unlist(sqrt(pcobj$eig)[1])
    u <- sweep(pcobj$ind$coord, 2, 1/(d * nobs.factor), FUN = "*")
    v <- sweep(pcobj$var$coord, 2, sqrt(pcobj$eig[1:ncol(pcobj$var$coord), 
                                                  1]), FUN = "/")
  }
  else if (inherits(pcobj, "lda")) {
    nobs.factor <- sqrt(pcobj$N)
    d <- pcobj$svd
    u <- predict(pcobj)$x/nobs.factor
    v <- pcobj$scaling
    d.total <- sum(d^2)
  }
  else {
    stop("Expected a object of class prcomp, princomp, PCA, or lda")
  }
  choices <- pmin(choices, ncol(u))
  df.u <- as.data.frame(sweep(u[, choices], 2, d[choices]^obs.scale, 
                              FUN = "*"))
  v <- sweep(v, 2, d^var.scale, FUN = "*")
  df.v <- as.data.frame(v[, choices])
  names(df.u) <- c("xvar", "yvar")
  names(df.v) <- names(df.u)
  if (pc.biplot) {
    df.u <- df.u * nobs.factor
  }
  r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)
  v.scale <- rowSums(v^2)
  df.v <- r * df.v/sqrt(max(v.scale))
  if (obs.scale == 0) {
    u.axis.labs <- paste("standardized PC", choices, sep = "")
  }
  else {
    u.axis.labs <- paste("PC", choices, sep = "")
  }
  u.axis.labs <- paste(u.axis.labs, sprintf("(%0.1f%% explained var.)", 
                                            100 * pcobj$sdev[choices]^2/sum(pcobj$sdev^2)))
  if (!is.null(labels)) {
    df.u$labels <- labels
  }
  if (!is.null(groups)) {
    df.u$groups <- groups
  }
  if (varname.abbrev) {
    df.v$varname <- abbreviate(rownames(v))
  }
  else {
    df.v$varname <- rownames(v)
  }
  df.v$angle <- with(df.v, (225/pi) * atan(-yvar/-xvar))
  df.v$hjust = with(df.v, (1 - varname.adjust * sign(-xvar))/2)
  g <- ggplot(data = df.u, aes(x = xvar, y = yvar)) + xlab(u.axis.labs[1]) + 
    ylab(u.axis.labs[2]) + coord_equal()
  if (var.axes) {
    if (circle) {
      theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, 
                                                length = 50))
      circle <- data.frame(xvar = r * cos(theta), yvar = r * 
                             sin(theta))
      g <- g + geom_path(data = circle, color = muted("white"), 
                         size = 1/2, alpha = 1/3)
    }
    g <- g + geom_segment(data = df.v, aes(x = 0, y = 0, 
                                           xend = xvar, yend = yvar*2), arrow = arrow(length = unit(1/2, 
                                                                                                    "picas")), color = color, linetype = linetype, alpha = alpha_arrow)
  }
  if (!is.null(df.u$labels)) {
    if (!is.null(df.u$groups)) {
      g <- g + geom_text(aes(label = labels, color = groups), 
                         size = labels.size)
    }
    else {
      g <- g + geom_text(aes(label = labels), size = labels.size)
    }
  }
  else {
    if (!is.null(df.u$groups)) {
      g <- g + geom_point(aes(color = groups), alpha = alpha)
    }
    else {
      g <- g + geom_point(alpha = alpha)
    }
  }
  if (!is.null(df.u$groups) && ellipse) {
    theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
    circle <- cbind(cos(theta), sin(theta))
    ell <- ddply(df.u, "groups", function(x) {
      if (nrow(x) <= 2) {
        return(NULL)
      }
      sigma <- var(cbind(x$xvar, x$yvar))
      mu <- c(mean(x$xvar), mean(x$yvar))
      ed <- sqrt(qchisq(ellipse.prob, df = 2))
      data.frame(sweep(circle %*% chol(sigma) * ed, 2, 
                       mu, FUN = "+"), groups = x$groups[1])
    })
    names(ell)[1:2] <- c("xvar", "yvar")
    g <- g + geom_path(data = ell, aes(color = groups, group = groups))
  }
  if (var.axes) {
    g <- g + geom_text(data = df.v, aes(label = varname, 
                                        x = xvar, y = yvar*2, angle = angle, hjust = hjust), 
                       color = "black", size = varname.size)
  }
  return(g)
}

#biplot_flip####
ggbiplot_reverse_flip<-function (pcobj, choices = 1:2, scale = 1, pc.biplot = TRUE, 
                                 obs.scale = 1 - scale, var.scale = scale, groups = NULL, 
                                 ellipse = TRUE, ellipse.prob = 0.68, labels = NULL, labels.size = 3, 
                                 alpha = 1, var.axes = TRUE, circle = FALSE, circle.prob = 0.69, 
                                 varname.size = 3, varname.adjust = 1.5, varname.abbrev = FALSE, 
                                 color = muted("red"), # <- add new arguments to the function
                                 linetype = "solid",
                                 alpha_arrow = 1) 
{
  library(ggplot2)
  library(plyr)
  library(scales)
  library(grid)
  stopifnot(length(choices) == 2)
  if (inherits(pcobj, "prcomp")) {
    nobs.factor <- sqrt(nrow(pcobj$x) - 1)
    d <- pcobj$sdev
    u <- sweep(pcobj$x, 2, 1/(d * nobs.factor), FUN = "*")
    v <- pcobj$rotation
  }
  else if (inherits(pcobj, "princomp")) {
    nobs.factor <- sqrt(pcobj$n.obs)
    d <- pcobj$sdev
    u <- sweep(pcobj$scores, 2, 1/(d * nobs.factor), FUN = "*")
    v <- pcobj$loadings
  }
  else if (inherits(pcobj, "PCA")) {
    nobs.factor <- sqrt(nrow(pcobj$call$X))
    d <- unlist(sqrt(pcobj$eig)[1])
    u <- sweep(pcobj$ind$coord, 2, 1/(d * nobs.factor), FUN = "*")
    v <- sweep(pcobj$var$coord, 2, sqrt(pcobj$eig[1:ncol(pcobj$var$coord), 
                                                  1]), FUN = "/")
  }
  else if (inherits(pcobj, "lda")) {
    nobs.factor <- sqrt(pcobj$N)
    d <- pcobj$svd
    u <- predict(pcobj)$x/nobs.factor
    v <- pcobj$scaling
    d.total <- sum(d^2)
  }
  else {
    stop("Expected a object of class prcomp, princomp, PCA, or lda")
  }
  choices <- pmin(choices, ncol(u))
  df.u <- as.data.frame(sweep(u[, choices], 2, d[choices]^obs.scale, 
                              FUN = "*"))
  v <- sweep(v, 2, d^var.scale, FUN = "*")
  df.v <- as.data.frame(v[, choices])
  names(df.u) <- c("yvar", "xvar")
  names(df.v) <- names(df.u)
  if (pc.biplot) {
    df.u <- df.u * nobs.factor
  }
  r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)
  v.scale <- rowSums(v^2)
  df.v <- r * df.v/sqrt(max(v.scale))
  if (obs.scale == 0) {
    u.axis.labs <- paste("standardized PC", choices, sep = "")
  }
  else {
    u.axis.labs <- paste("PC", choices, sep = "")
  }
  u.axis.labs <- paste(u.axis.labs, sprintf("(%0.1f%% explained var.)", 
                                            100 * pcobj$sdev[choices]^2/sum(pcobj$sdev^2)))
  if (!is.null(labels)) {
    df.u$labels <- labels
  }
  if (!is.null(groups)) {
    df.u$groups <- groups
  }
  if (varname.abbrev) {
    df.v$varname <- abbreviate(rownames(v))
  }
  else {
    df.v$varname <- rownames(v)
  }
  df.v$angle <- with(df.v, (225/pi) * atan(-yvar/xvar))
  df.v$hjust = with(df.v, (1 - varname.adjust * sign(xvar))/2)
  g <- ggplot(data = df.u, aes(x = xvar, y = yvar)) + xlab(u.axis.labs[1]) + 
    ylab(u.axis.labs[2]) + coord_equal()
  if (var.axes) {
    if (circle) {
      theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, 
                                                length = 50))
      circle <- data.frame(xvar = r * cos(theta), yvar = r * 
                             sin(theta))
      g <- g + geom_path(data = circle, color = muted("white"), 
                         size = 1/2, alpha = 1/3)
    }
    g <- g + geom_segment(data = df.v, aes(x = 0, y = 0, 
                                           xend = xvar, yend = yvar*2), arrow = arrow(length = unit(1/2, 
                                                                                                    "picas")), color = color, linetype = linetype, alpha = alpha_arrow)
  }
  if (!is.null(df.u$labels)) {
    if (!is.null(df.u$groups)) {
      g <- g + geom_text(aes(label = labels, color = groups), 
                         size = labels.size)
    }
    else {
      g <- g + geom_text(aes(label = labels), size = labels.size)
    }
  }
  else {
    if (!is.null(df.u$groups)) {
      g <- g + geom_point(aes(color = groups), alpha = alpha)
    }
    else {
      g <- g + geom_point(alpha = alpha)
    }
  }
  if (!is.null(df.u$groups) && ellipse) {
    theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
    circle <- cbind(cos(theta), sin(theta))
    ell <- ddply(df.u, "groups", function(x) {
      if (nrow(x) <= 2) {
        return(NULL)
      }
      sigma <- var(cbind(x$xvar, x$yvar))
      mu <- c(mean(x$xvar), mean(x$yvar))
      ed <- sqrt(qchisq(ellipse.prob, df = 2))
      data.frame(sweep(circle %*% chol(sigma) * ed, 2, 
                       mu, FUN = "+"), groups = x$groups[1])
    })
    names(ell)[1:2] <- c("yvar", "xvar")
    g <- g + geom_path(data = ell, aes(color = groups, group = groups))
  }
  if (var.axes) {
    g <- g + geom_text(data = df.v, aes(label = varname, 
                                        x = xvar, y = yvar*2, angle = angle, hjust = hjust), 
                       color = "black", size = varname.size)
  }
  return(g)
}

#biplot_flip2####
ggbiplot_flip<-function (pcobj, choices = 1:2, scale = 1, pc.biplot = TRUE, 
                         obs.scale = 1 - scale, var.scale = scale, groups = NULL, 
                         ellipse = TRUE, ellipse.prob = 0.68, labels = NULL, labels.size = 3, 
                         alpha = 1, var.axes = TRUE, circle = FALSE, circle.prob = 0.69, 
                         varname.size = 3, varname.adjust = 1.5, varname.abbrev = FALSE, 
                         color = muted("red"), # <- add new arguments to the function
                         linetype = "solid",
                         alpha_arrow = 1) 
{
  library(ggplot2)
  library(plyr)
  library(scales)
  library(grid)
  stopifnot(length(choices) == 2)
  if (inherits(pcobj, "prcomp")) {
    nobs.factor <- sqrt(nrow(pcobj$x) - 1)
    d <- pcobj$sdev
    u <- sweep(pcobj$x, 2, 1/(d * nobs.factor), FUN = "*")
    v <- pcobj$rotation
  }
  else if (inherits(pcobj, "princomp")) {
    nobs.factor <- sqrt(pcobj$n.obs)
    d <- pcobj$sdev
    u <- sweep(pcobj$scores, 2, 1/(d * nobs.factor), FUN = "*")
    v <- pcobj$loadings
  }
  else if (inherits(pcobj, "PCA")) {
    nobs.factor <- sqrt(nrow(pcobj$call$X))
    d <- unlist(sqrt(pcobj$eig)[1])
    u <- sweep(pcobj$ind$coord, 2, 1/(d * nobs.factor), FUN = "*")
    v <- sweep(pcobj$var$coord, 2, sqrt(pcobj$eig[1:ncol(pcobj$var$coord), 
                                                  1]), FUN = "/")
  }
  else if (inherits(pcobj, "lda")) {
    nobs.factor <- sqrt(pcobj$N)
    d <- pcobj$svd
    u <- predict(pcobj)$x/nobs.factor
    v <- pcobj$scaling
    d.total <- sum(d^2)
  }
  else {
    stop("Expected a object of class prcomp, princomp, PCA, or lda")
  }
  choices <- pmin(choices, ncol(u))
  df.u <- as.data.frame(sweep(u[, choices], 2, d[choices]^obs.scale, 
                              FUN = "*"))
  v <- sweep(v, 2, d^var.scale, FUN = "*")
  df.v <- as.data.frame(v[, choices])
  names(df.u) <- c("yvar", "xvar")
  names(df.v) <- names(df.u)
  if (pc.biplot) {
    df.u <- df.u * nobs.factor
  }
  r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)
  v.scale <- rowSums(v^2)
  df.v <- r * df.v/sqrt(max(v.scale))
  if (obs.scale == 0) {
    u.axis.labs <- paste("standardized PC", choices, sep = "")
  }
  else {
    u.axis.labs <- paste("PC", choices, sep = "")
  }
  u.axis.labs <- paste(u.axis.labs, sprintf("(%0.1f%% explained var.)", 
                                            100 * pcobj$sdev[choices]^2/sum(pcobj$sdev^2)))
  if (!is.null(labels)) {
    df.u$labels <- labels
  }
  if (!is.null(groups)) {
    df.u$groups <- groups
  }
  if (varname.abbrev) {
    df.v$varname <- abbreviate(rownames(v))
  }
  else {
    df.v$varname <- rownames(v)
  }
  df.v$angle <- with(df.v, (225/pi) * atan(yvar/xvar))
  df.v$hjust = with(df.v, (1 - varname.adjust * sign(xvar))/2)
  g <- ggplot(data = df.u, aes(x = xvar, y = yvar)) + xlab(u.axis.labs[1]) + 
    ylab(u.axis.labs[2]) + coord_equal()
  if (var.axes) {
    if (circle) {
      theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, 
                                                length = 50))
      circle <- data.frame(xvar = r * cos(theta), yvar = r * 
                             sin(theta))
      g <- g + geom_path(data = circle, color = muted("white"), 
                         size = 1/2, alpha = 1/3)
    }
    g <- g + geom_segment(data = df.v, aes(x = 0, y = 0, 
                                           xend = xvar, yend = yvar*2), arrow = arrow(length = unit(1/2, 
                                                                                                    "picas")), color = color, linetype = linetype, alpha = alpha_arrow)
  }
  if (!is.null(df.u$labels)) {
    if (!is.null(df.u$groups)) {
      g <- g + geom_text(aes(label = labels, color = groups), 
                         size = labels.size)
    }
    else {
      g <- g + geom_text(aes(label = labels), size = labels.size)
    }
  }
  else {
    if (!is.null(df.u$groups)) {
      g <- g + geom_point(aes(color = groups), alpha = alpha)
    }
    else {
      g <- g + geom_point(alpha = alpha)
    }
  }
  if (!is.null(df.u$groups) && ellipse) {
    theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
    circle <- cbind(cos(theta), sin(theta))
    ell <- ddply(df.u, "groups", function(x) {
      if (nrow(x) <= 2) {
        return(NULL)
      }
      sigma <- var(cbind(x$xvar, x$yvar))
      mu <- c(mean(x$xvar), mean(x$yvar))
      ed <- sqrt(qchisq(ellipse.prob, df = 2))
      data.frame(sweep(circle %*% chol(sigma) * ed, 2, 
                       mu, FUN = "+"), groups = x$groups[1])
    })
    names(ell)[1:2] <- c("yvar", "xvar")
    g <- g + geom_path(data = ell, aes(color = groups, group = groups))
  }
  if (var.axes) {
    g <- g + geom_text(data = df.v, aes(label = varname, 
                                        x = xvar, y = yvar*2, angle = angle, hjust = hjust), 
                       color = "black", size = varname.size)
  }
  return(g)
}


#Calculate global 'galactic plane'#####

#Read in global data
#
t.fread<-fread("...",verbose=TRUE) #Data available from www.try-db.org 
nrow(t.fread)
head(t.fread)
se <- function(x) sd(x)/sqrt(length(x))
`%notin%` <- function(x,y) !(x %in% y)

#Read in global data
traits<-as.data.frame(unique(t.fread$TraitName))
load("...")  #Data available from github.com/TundraTraitTeam/TraitHub
try.ttt<-try.ttt.clean

(nrow(t.fread[t.fread$AccSpeciesName %in% try.ttt$AccSpeciesName,]) / nrow(t.fread))*100
try.ttt<-subset(try.ttt,TraitShort=="SeedMass"|TraitShort=="StemSpecificDensity"|TraitShort=="LeafN"|TraitShort=="SLA"|TraitShort=="LDMC"|TraitShort=="LeafArea"|TraitShort=="PlantHeight")
try.ttt[try.ttt$TraitShort=="StemSpecificDensity",2]<-"SSD"

#Traits of Interest:
#Seed dry mass
#Stem dry mass per stem fresh volume (stem specific density, SSD, wood density)
#Leaf nitrogen (N) content per leaf dry mass
#Leaf area per leaf dry mass (specific leaf area, SLA)
#Leaf area
#Leaf dry mass per leaf fresh mass (Leaf dry matter content, LDMC)
#Plant height

target<-c("Seed dry mass",
          "Leaf area",
          "Stem dry mass per stem fresh volume (stem specific density, SSD, wood density)",
          "Leaf nitrogen (N) content per leaf dry mass",
          "Leaf area per leaf dry mass (specific leaf area, SLA)",
          "Leaf dry mass per leaf fresh mass (Leaf dry matter content, LDMC)",
          "Plant height")
allspecies<-filter(t.fread, TraitName %in% target)

#Select species####

#1) TTT species list

load("...") #Data available from github.com/TundraTraitTeam/TraitHub
cover_species<-as.data.frame(unique(coverc.sub$Name))

names(cover_species)[1]<-"Name"

#2) Species from sites below 10oC July temp

ttt_only<-try.ttt %>%
  filter(Source == "TTT") %>%
  dplyr::select(Lon,Lat,AccSpeciesName) 

ttt_coords<-try.ttt %>%
  filter(Source == "TTT") %>%
  dplyr::select(Lon,Lat) %>%
  distinct()

ttt_only$latlon<-paste(ttt_only$Lat,ttt_only$Lon,sep="")

#Get temp data
July_temp <- raster("...") #Data available from chelsa-climate.org

#Extract cooridnates
ttt_climate <- raster::extract(July_temp, ttt_coords)

#Combine climate
ttt_coords$temp<-ttt_climate/10
ttt_coords$latlon<-paste(ttt_coords$Lat,ttt_coords$Lon,sep="")

#Combine with full object
ttt_only$Jul_temp<-ttt_coords$temp[match(ttt_only$latlon, ttt_coords$latlon)]

cover_species_ttt<-ttt_only %>%
  filter(Jul_temp <10) %>%
  dplyr::select(AccSpeciesName) %>%
  distinct(AccSpeciesName)

cover_species_ttt<-as.data.frame(cover_species_ttt)
names(cover_species_ttt)[1]<-"Name"

#3) ITEX Species list
load("...") #Data available from polardata.ca
cover_species_ITEX<-as.data.frame(unique(itex.vasccover.full$name))
names(cover_species_ITEX)[1]<-"Name"

#4) summit div species list
sdiv<-read.csv("...") #Data from Steinbauer M.J., Grytnes J.-A., Wipf S. (2018) 
#Accelerated increase in plant species richness on mountain summits is linked to warming. Nature v. 556, pages 231–234
cover_species_sdiv<-as.data.frame(unique(sdiv$Accepted_species_name))
names(cover_species_sdiv)[1]<-"Name"

#COMBINE species
species_list<-rbind(cover_species_ITEX, cover_species_ttt, cover_species, cover_species_sdiv)

species_list<-unique(species_list$Name)
try.ttt_yes<-filter(try.ttt, AccSpeciesName %in% species_list)
try.ttt_no<-filter(try.ttt, AccSpeciesName %notin% species_list)

this_species_list<-as.data.frame(unique(try.ttt_yes$AccSpeciesName))

try.ttt_all<-try.ttt
try.ttt<-try.ttt_yes

####Clean out extreme values####

#Correct a few final synonyms
try.ttt[try.ttt$AccSpeciesName=="Myosotis alpina",]$AccSpeciesName<-"Myosotis alpestris"
try.ttt[try.ttt$AccSpeciesName=="Alnus alnobetula",]$AccSpeciesName<-"Alnus viridis"

allspecies[allspecies$TraitName=="Leaf nitrogen (N) content per leaf dry mass" & allspecies$AccSpeciesName=="Astragalus armatus",]$StdValue<-allspecies[allspecies$TraitName=="Leaf nitrogen (N) content per leaf dry mass" & allspecies$AccSpeciesName=="Astragalus armatus",]$StdValue/10
try.ttt<-subset(try.ttt,DataContributor!="John Dickie"|StdValue<50)

# #Exclude trees
PH<-subset(try.ttt,TraitShort=="PlantHeight")
PHs<-PH %>%
  group_by(AccSpeciesName) %>%
  summarise(height = mean(StdValue))
PM_species<-unique(subset(PHs,height>=5)$AccSpeciesName)

#Note tree species
PM_species
#Add boreal trees
PM_all<-c(PM_species,c("Alnus viridis", "Picea glauca", "Picea abies", "Picea mariana", "Betula pendula", "Betula glandulosa"))

try.ttt.PM<-subset(try.ttt,AccSpeciesName%in%PM_all)
try.ttt<-subset(try.ttt,AccSpeciesName%notin%PM_all)

#Add metadata
metadata<-unique(as.data.frame(try.ttt.full[,c(1,13)]))

try.ttt$F_Group<-metadata$GrowthForm[match(try.ttt$AccSpeciesName,metadata$Name)]

missing<-subset(try.ttt,is.na(F_Group))
try.ttt<-subset(try.ttt,!is.na(F_Group))

F_Groups<-read.csv("...") #Data from github.com/hjdthomas/Tundra_functional_groups

metadata2<-F_Groups[,c(2,12)]
missing$F_Group<-metadata2$F_Group[match(missing$AccSpeciesName,metadata2$AccSpeciesName)]

try.ttt<-rbind(try.ttt,missing)

out<- try.ttt %>% 
  group_by(AccSpeciesName) %>%
  summarise(Obs = length(StdValue),
            Sites = length(unique(Lat)))

sum(out$Obs==1) /nrow(out)
mean(out$Obs)
mean(out$Sites)

#Cleaning: all plants####
allspecies[allspecies$TraitName=="Leaf dry mass per leaf fresh mass (Leaf dry matter content, LDMC)"&allspecies$StdValue>1&!is.na(allspecies$StdValue),]$StdValue<-1/allspecies[allspecies$TraitName=="Leaf dry mass per leaf fresh mass (Leaf dry matter content, LDMC)"&allspecies$StdValue>1&!is.na(allspecies$StdValue),]$StdValue

try.ttt$unique<-paste(try.ttt$AccSpeciesName,try.ttt$TraitShort,try.ttt$DataContributor,try.ttt$StdValue,sep="_")

try.ttt_del<-try.ttt %>% 
  group_by(TraitShort,Genus) %>% 
  filter(StdValue>mean(StdValue)+sd(StdValue)*4|StdValue<mean(StdValue)-sd(StdValue)*4)

try.ttt_del$unique<-paste(try.ttt_del$AccSpeciesName,try.ttt_del$TraitShort,try.ttt_del$DataContributor,try.ttt_del$StdValue,sep="_")
try.ttt<-subset(try.ttt,unique%notin%try.ttt_del$unique)


allspecies$Genus<-gsub("([A-Za-z]+).*", "\\1", allspecies$AccSpeciesName)

allspecies_del<-allspecies %>% 
  group_by(TraitName,Genus) %>% 
  filter(StdValue>mean(StdValue)+sd(StdValue)*4|StdValue<mean(StdValue)-sd(StdValue)*4)

allspecies<-subset(allspecies,ObservationID%notin%allspecies_del$ObservationID)

allspecies$F_Group<-"World"

# T_Table<-NULL
# S_Table<-NULL
# #Clean data
# 
# for (i in unique(allspecies$TraitName)){
#   k<-subset(allspecies,TraitName==i)
#   k<-subset(k,!is.na(StdValue))
#   #S_table<-NULL
#   #4sd
#   #for(m in unique(k$AccSpeciesName)){
#    # p<-subset(k,AccSpeciesName==m)
#   #change back to p
#     cleaned<-subset(k,StdValue<(mean(k$StdValue,)+(sd(k$StdValue,)*4)))
#     cleaned<-subset(cleaned,StdValue>(mean(k$StdValue,)-(sd(k$StdValue,)*4)))
#     S_Table<-rbind(S_Table,cleaned)
# }
#   #}
#   #T_table<-rbind(T_Table,S_Table)
# allspecies<-S_Table

#Convert SLA to LMA:

allspecies[allspecies$TraitName=="Leaf area per leaf dry mass (specific leaf area, SLA)",]$StdValue<-1/allspecies[allspecies$TraitName=="Leaf area per leaf dry mass (specific leaf area, SLA)",]$StdValue
try.ttt[try.ttt$TraitShort=="SLA",]$StdValue<-1/try.ttt[try.ttt$TraitShort=="SLA",]$StdValue

#Log10 transform data

#Keep try.ttt.part. for later on 

#All plants
allspecies$StdValue<-log10(allspecies$StdValue)
allspecies<-allspecies[which((allspecies)$StdValue!="-Inf"),]

#TTT
try.ttt$StdValue<-log10(try.ttt$StdValue)
try.ttt<-try.ttt[which((try.ttt)$StdValue!="-Inf"),]

try.ttt$latlon<-paste0(try.ttt$Lat,try.ttt$Lon)

length(unique(try.ttt$latlon))

####All plants####
trait.meanstable<-ddply(allspecies,.(TraitName,AccSpeciesName),summarise,
                        mean = mean(StdValue, na.rm=TRUE),
                        sd = sd(StdValue, na.rm=TRUE),
                        se = se(StdValue),
                        n= length(StdValue))

trait.meanstable[trait.meanstable$TraitName=="Seed dry mass",1]<-"SeedMass"
trait.meanstable[trait.meanstable$TraitName=="Leaf area",1]<-"LeafArea"
trait.meanstable[trait.meanstable$TraitName=="Stem dry mass per stem fresh volume (stem specific density, SSD, wood density)",1]<-"SSD"
trait.meanstable[trait.meanstable$TraitName=="Leaf nitrogen (N) content per leaf dry mass",1]<-"LeafN"
trait.meanstable[trait.meanstable$TraitName=="Leaf area per leaf dry mass (specific leaf area, SLA)",1]<-"SLA"
trait.meanstable[trait.meanstable$TraitName=="Leaf dry mass per leaf fresh mass (Leaf dry matter content, LDMC)",1]<-"LDMC"
trait.meanstable[trait.meanstable$TraitName=="Plant height",1]<-"PlantHeight"

trait.meanstable<- trait.meanstable[order(trait.meanstable$TraitName),]

#Create means
species.names <- as.data.frame(unique(allspecies$AccSpeciesName))
names(species.names)<-"Species"
species.names <- as.data.frame(species.names[order(species.names$Species),])
names(species.names)<-"Species"

testtable<-NULL
species.meanstable<-species.names
for (i in unique(trait.meanstable$TraitName)){
  testtable<-species.names
  names(testtable)<-"Species"
  testtable<-cbind(testtable,trait.meanstable[trait.meanstable$TraitName==i,]$mean[match(testtable$Species,trait.meanstable[trait.meanstable$TraitName==i,]$AccSpeciesName)])
  testtable<-as.data.frame(testtable[,c(2)])
  names(testtable)<-c(paste(i,"_mean",sep=""))
  species.meanstable<-cbind(species.meanstable,testtable)
}
species.meanstable$Group<-"World"

####Tundra####
trait.meanstable_t<-ddply(try.ttt,.(TraitShort,AccSpeciesName),summarise,
                          mean = mean(StdValue, na.rm=TRUE),
                          sd = sd(StdValue, na.rm=TRUE),
                          se = se(StdValue),
                          n= length(StdValue))

trait.meanstable_t<- trait.meanstable_t[order(trait.meanstable_t$TraitShort),]
trait.meanstable_t<-subset(trait.meanstable_t,TraitShort!="RSratio")

#Turn into table for bi-plots
species.names_t <- as.data.frame(unique(trait.meanstable_t$AccSpeciesName))
names(species.names_t)<-"Species"
species.names_t <- as.data.frame(species.names_t[order(species.names_t$Species),])
names(species.names_t)<-"Species"

testtable_t<-NULL
species.meanstable_t<-species.names_t
for (i in unique(trait.meanstable_t$TraitShort)){
  testtable_t<-species.names_t
  names(testtable_t)<-"Species"
  testtable_t<-cbind(testtable_t,trait.meanstable_t[trait.meanstable_t$TraitShort==i,]$mean[match(testtable_t$Species,trait.meanstable_t[trait.meanstable_t$TraitShort==i,]$AccSpeciesName)])
  testtable_t<-as.data.frame(testtable_t[,c(2)])
  names(testtable_t)<-c(paste(i,"_mean",sep=""))
  species.meanstable_t<-cbind(species.meanstable_t,testtable_t)
}
species.meanstable_t$Group<-"Tundra"

#Combine with TTT
all_no_t<-subset(species.meanstable, Species %notin% species.names_t$Species)
species.meanstable_wTTT<-rbind(all_no_t,species.meanstable_t)


#PCA all traits
PCA_table6<-species.meanstable_wTTT[complete.cases(species.meanstable_wTTT[,2:7]),]
PCA_tundra<-subset(PCA_table6,Group=="Tundra")

length(unique(PCA_world$Species))

PCA_tundra$F_Group<-try.ttt$F_Group[match(PCA_tundra$Species,try.ttt$AccSpeciesName)]
PCA_world<-subset(PCA_table6,Group=="World")
PCA_world$F_Group<-"World"
PCA_table6<-rbind(PCA_world,PCA_tundra)

pca.data6<-PCA_table6[,2:7]
pca.species6<-PCA_table6[,1]
pca.groups6<-PCA_table6[,9]

nrow(PCA_table6)

tundra.pca6 <- prcomp(pca.data6,
                      center = TRUE,
                      scale. = TRUE)

library(reshape2)
loadings_6<-as.data.frame(tundra.pca6$rotation[,c(1:2)])
loadings_6$trait<-row.names(loadings_6)
loadings_6<-melt(loadings_6, id.vars="trait", measure.vars=c("PC1","PC2"))
loadings_6$value<-as.numeric(as.character(loadings_6$value))
names(loadings_6)[2]<-"PC"

ggplot(loadings_6,aes(PC,value,fill=trait))+
  geom_bar(stat="identity", position = "dodge")+
  theme_bw()

g6 <- ggbiplot_reverse_xy(tundra.pca6, obs.scale = 1, var.scale = 1, ellipse = FALSE, 
                          circle = FALSE, alpha=0.8,group= pca.groups6,varname.size=0,varname.adjust = 1.5,color = "black")
g6 <- g6 + scale_colour_manual(values=c("blue","grey80"),name = '')+
  theme_bw()+
  scale_x_reverse( lim=c(5,-5))+
  scale_y_reverse( lim=c(5,-5))+
  theme(legend.title = element_text(size=14), legend.text = element_text(size=11), legend.position=c(0.1, 1))+
  theme(legend.background = element_rect(color = "black", fill = "transparent", size = 2, linetype="blank"))+
  theme(plot.title = element_text(size=15, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  axis.text=element_text(size=10), axis.title=element_text(size=11,vjust = -0.5),axis.text.x=element_text(vjust = -0.05),legend.key = element_blank(), plot.margin = unit(c(0.5,0.5,0.5,0.5), units = , "cm"))+
  theme(legend.direction = 'horizontal', 
        legend.position = 'top')+
  theme(aspect.ratio=1/1)+
  theme(axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"))+
  labs(x="PC1 - Plant Size\n(37.5% Variation)",y="PC2 - Resource Economics\n(29.1% Variation)")
g6$layers <- c(g6$layers, g6$layers[[1]])
print(g6)


##############################
#Calculate PCAs#
##############################

#Tundra only#####
PCA_tundra_only<-subset(PCA_table6,Group=="Tundra")

try.ttt_only<-subset(try.ttt,AccSpeciesName %in% PCA_tundra_only$Species)

CHELSA_temp<-raster("/Volumes/Teamshrub/Climate_Data/Chelsa/CHELSA_bio10_1.tif")
try.ttt_only$temp<-raster::extract(CHELSA_temp, try.ttt_only[,c(5,4)],df = TRUE)[,2]


try.ttt.notemp<-try.ttt_only[is.na(try.ttt_only$temp),]
try.ttt.temp<-try.ttt_only[!is.na(try.ttt_only$temp),]

Out<-NULL
for(i in unique(try.ttt.temp$AccSpeciesName)){
  sp_subset<-subset(try.ttt.temp,AccSpeciesName==i)
  table_subset<-subset(PCA_tundra_only,Species==i)
  table_subset$min_WQ<-mean(sp_subset$temp, na.rm = T)/10
  Out<-rbind(Out,table_subset)
}

PCA_tundra_only<-Out

PCA_tundra_only$min_WQ_cat<-ifelse(PCA_tundra_only$min_WQ<=(-1),1,
                                   ifelse(PCA_tundra_only$min_WQ<(1),2,3))

PCA_tundra_only[is.na(PCA_tundra_only$min_WQ_cat),]$min_WQ_cat<-"No_data"

pca.data_to<-PCA_tundra_only[,2:7]
pca.species_to<-PCA_tundra_only[,1]
pca.groups_to<-PCA_tundra_only$min_WQ_cat

tundra_only.pca <- prcomp(pca.data_to,
                          center = TRUE,
                          scale. = TRUE)

tundra_only.pca

loadings_to<-as.data.frame(tundra_only.pca$rotation[,c(1:2)])
loadings_to$trait<-row.names(loadings_to)
loadings_to<-melt(loadings_to, id.vars="trait", measure.vars=c("PC1","PC2"))

loadings_to$value<-as.numeric(as.character(loadings_to$value))
names(loadings_to)[2]<-"PC"

ggplot(loadings_to,aes(PC,value,fill=trait))+
  geom_bar(stat="identity", position = "dodge")+
  theme_bw()

#Combine
loadings_6$type<-"World"
loadings_to$type<-"Tundra"
#Swap axes
loadings_to_swap_1<-subset(loadings_to,PC=="PC1")
loadings_to_swap_2<-subset(loadings_to,PC=="PC2")
loadings_to_swap_1$PC<-"PC2"
loadings_to_swap_2$PC<-"PC1"
loadings_to<-rbind(loadings_to_swap_1, loadings_to_swap_2)

loadings_6$value <-loadings_6$value*-1

loadings_all<-rbind(loadings_6,loadings_to)
#Reorder factors so world first
loadings_all$type <- factor(loadings_all$type, levels=c("World", "Tundra"))
loadings_all$trait <- factor(loadings_all$trait, levels=c("Height", "Leaf Area", "Seed Mass", "Leaf N", "LMA", "LDMC"))


pdf(file="...", width = 7.5, height = 6)
ggplot(loadings_all, aes(PC,value,fill=trait, colour=trait, linetype=type, alpha=type))+
  geom_bar(stat="identity", position = "dodge")+
  theme_bw()+
  theme(plot.title = element_text(size=15, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black"), axis.text=element_text(size=10), axis.title=element_text(size=12,vjust = -0.5),legend.key = element_blank(), plot.margin = unit(c(1,1,1,1), units = , "cm"))+
  scale_alpha_manual(values = c(1,0.3), name = '', guide="none")+
  scale_fill_manual(values = c("#2C5197","#bd0016","grey30","#2E6444","orange","burlywood3"),
                    labels = c("Plant Height", "Leaf Area", "Seed Mass", "Leaf Nitrogen", "LMA", "LDMC"),
                      name = 'Trait')+
  scale_colour_manual(values = c("#2C5197","#bd0016","grey30","#2E6444","orange","burlywood3"),
                        guide="none")+
  scale_linetype(guide="none")+
  geom_hline(aes(yintercept=0))+
  scale_x_discrete(labels=c("PC1" = "Plant Size:\nGlobal PC1\nTundra PC2", "PC2" = "Resource Economics:\nGlobal PC2\nTundra PC1"))+
  labs(x="\nPCA axis",y="Trait Loading\n")
dev.off()

#labels = PCA_tundra_only$Species,
g_to <- ggbiplot_flip(tundra_only.pca, obs.scale = 1, var.scale = 1, ellipse = FALSE,
                      circle = FALSE,group= as.factor(pca.groups_to),alpha=0.8,varname.size=0,varname.adjust = 1.5,color = "black")
g_to <- g_to + scale_colour_manual(values=c("midnightblue","lightskyblue","darkorange1", "grey50"),
                                   labels = c("Cold", "Mid", "Warm", "No data"),
                                   name = '')


g_to <- g_to +
  theme_bw()+
  xlim(-5,5)+
  ylim(-5,5)+
  theme(legend.title = element_text(size=14), legend.text = element_text(size=11), legend.position=c(0.1, 1))+
  theme(legend.background = element_rect(color = "black", fill = "transparent", size = 2, linetype="blank"))+
  theme(plot.title = element_text(size=15, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black"), axis.text=element_text(size=10), axis.title=element_text(size=11,vjust = -0.5),legend.key = element_blank(), plot.margin = unit(c(0.5,0.5,0.5,0.5), units = , "cm"))+
  theme(legend.direction = 'horizontal', 
        legend.position = 'top')+
  theme(aspect.ratio=1/1)+
  theme(legend.title=element_text(size=11))+
  labs(y="Global PC2 - Resource Economics\n(37.6% Variation)",x="Global PC1 - Plant Size\n(26.9% Variation)")+
  theme(axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"))
g_to$layers <- c(g_to$layers, g_to$layers[[1]])
print(g_to)

pdf(file="...", width = 8, height = 4)
grid.arrange(g6,g_to,ncol=2)
dev.off()

g_to <- ggbiplot_flip(tundra_only.pca, obs.scale = 1, var.scale = 1, ellipse = FALSE,
                      circle = FALSE,alpha=0.8,varname.size=4.5,varname.adjust = 1.5,color = "black")

g_to <- g_to +
  geom_point(colour="blue") +
  theme_bw()+
  xlim(-5,5)+
  ylim(-5,5)+
  theme(legend.title = element_text(size=14), legend.text = element_text(size=11), legend.position=c(0.1, 1))+
  theme(legend.background = element_rect(color = "black", fill = "transparent", size = 2, linetype="blank"))+
  theme(plot.title = element_text(size=15, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black"), axis.text=element_text(size=10), axis.title=element_text(size=11,vjust = -0.5),axis.text.x=element_text(hjust = -0.05),legend.key = element_blank(), plot.margin = unit(c(1,1,1,1), units = , "cm"))+
  theme(legend.direction = 'horizontal', 
        legend.position = 'top')+
  theme(aspect.ratio=1/1)+
  theme(legend.title=element_text(size=11))+
  labs(y="Global PC1 - Plant Economics (37.6% Variation Explained)",x="Global PC2 - Plant Size (26.9% Variation Explained)")+
  theme(axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"))
g_to$layers <- c(g_to$layers, g_to$layers[[1]])
print(g_to)

##############


library("factoextra")
data(decathlon2)
# Extract active variables/individuals for PCA
decathlon2.active <- decathlon2[1:23, 1:10]
head(decathlon2.active[, 1:6])

library("FactoMineR")
res.pca <- PCA(decathlon2.active, graph = FALSE)
print(res.pca)

eigenvalues <- res.pca$eig
head(eigenvalues[, 1:2])

library("factoextra")
fviz_screeplot(res.pca, ncp=10)

head(res.pca$var$coord)

fviz_pca_var(res.pca)

head(res.pca$var$cos2)

fviz_pca_var(res.pca, col.var="cos2") +
  scale_color_gradient2(low="white", mid="blue", 
                        high="red", midpoint=0.5) + theme_minimal()

head(res.pca$var$contrib)

#World
(w1<-fviz_contrib(tundra.pca6, choice = "var", axes = 1)+
  ggtitle("a) Global PC axis 1 - Plant Size")+
  labs(y = "Contribution to PC axis (%)", x = "Trait")+
  ylim(0,35)+
  theme_classic())

(w2<-fviz_contrib(tundra.pca6, choice = "var", axes = 2)+
  ggtitle("b) Global PC axis 2- Resource Economics ")+
  labs(y = "Contribution to PC axis (%)", x = "Trait")+
  ylim(0,35)+
  theme_classic())

#Tundra
(t2<-fviz_contrib(tundra_only.pca, choice = "var", axes = 1)+
  ggtitle("d) Tundra PC axis 1 - Resource Economics \n(Global PC axis 2)")+
  labs(y = "Contribution to PC axis (%)", x = "Trait")+
  ylim(0,40)+
  theme_classic())
(t1<-fviz_contrib(tundra_only.pca, choice = "var", axes = 2)+
  ggtitle("c) Tundra PC axis 2- Plant Size \n(Global PC axis 1)")+
  labs(y = "Contribution to PC axis (%)", x = "Trait")+
  ylim(0,40)+
  theme_classic())

pdf(file="...", width = 8.5, height = 6)
grid.arrange(w1,w2,t1,t2)
dev.off()


#Draw distribution along axis####
# library(devtools)
#install_github("clauswilke/ggjoy")
library(ggjoy)

PCA_tundra_only$PC1<-(tundra_only.pca$x[,1])
PCA_tundra_only$PC2<-(tundra_only.pca$x[,2])

length(PCA_tundra_only[PCA_tundra_only$min_WQ_cat==1,])

library(ggplot2)
library(ggridges)

(b<-ggplot(PCA_tundra_only[PCA_tundra_only$min_WQ_cat!="No_data",], aes(x=PC1, y=as.character(min_WQ_cat), fill=factor(min_WQ_cat), height=..density..)) +
    #geom_joy(data=tundra_density, aes(x=PC1,y=min_WQ_cat,height=..density..),scale=8, adjust=1.5,fill="grey80")+
    geom_vline(xintercept=c(-5,0,5), linetype = "dashed", size = 0.4, colour = "black")+
    geom_joy(scale=1.7, bandwidth=1.75, alpha = 0.85)+
    xlim(-8,8)+
    scale_fill_manual(values=c("midnightb lue","lightskyblue","darkorange1"),name = '')+
    scale_y_discrete(expand = c(0.1, 0), 
                     labels=c("1" = "Tundra\n(Cold)", "2" = "Tundra\n(Mid)",
                              "3" = "Tundra\n(Warm)"))+
    theme_bw()+
    theme(plot.title = element_text(size=10, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black"), axis.text=element_text(size=10), axis.title=element_text(size=11,vjust = -0.5),axis.text.x=element_text(hjust = -0.05),legend.key = element_blank(), plot.margin = unit(c(1,1,1,1), units = , "cm"))+
    labs(x="\nSlower                  PC1 (Plant Economics)                  Faster", y="Temperature Class\n")+
    guides(fill=FALSE)+
    ggtitle("PC1 - Plant Economics"))

(a<-ggplot(PCA_tundra_only[PCA_tundra_only$min_WQ_cat!="No_data",], aes(x=PC2, y=as.character(min_WQ_cat), fill=factor(min_WQ_cat), height=..density..)) +
    #geom_joy(data=tundra_density, aes(x=PC2,y=min_WQ_cat,height=..density..),scale=8, adjust=1.5,fill="grey80")+
    geom_vline(xintercept=c(-5,0,5), linetype = "dashed", size = 0.4, colour = "black")+
    geom_joy(scale=1.7, bandwidth=1.75, alpha = 0.85)+
    xlim(-8,8)+
    scale_fill_manual(values=c("midnightb lue","lightskyblue","darkorange1"),name = '')+
    scale_y_discrete(expand = c(0.1, 0), 
                     labels=c("1" = "Tundra\n(Cold)", "2" = "Tundra\n(Mid)",
                              "3" = "Tundra\n(Warm)"))+
    theme_bw()+
    theme(plot.title = element_text(size=10, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black"), axis.text=element_text(size=10), axis.title=element_text(size=11,vjust = -0.5),axis.text.x=element_text(hjust = -0.05),legend.key = element_blank(), plot.margin = unit(c(1,1,1,1), units = , "cm"))+
    labs(x="\nSmaller              PC2 (Plant Size)              Larger", y="Temperature Class\n")+
    guides(fill=FALSE)+
    ggtitle("PC2 - Plant Size"))

pdf(file="...", width = 6, height = 7.5)
grid.arrange(b,a)
dev.off()

summary(aov(PC1 ~ min_WQ_cat, data = PCA_tundra_only))
summary(aov(PC2 ~ min_WQ_cat, data = PCA_tundra_only))


#Non-flipped
g_to <- ggbiplot_reverse(tundra_only.pca, obs.scale = 1, var.scale = 1, ellipse = FALSE,
                         circle = FALSE,group= pca.groups_to,alpha=0.9,varname.size=4.5,varname.adjust = 1.5,color = "black")
g_to <- g_to + scale_colour_gradient2(midpoint=-14,low="blue4", mid="blueviolet", high="darkorange2", na.value = "grey50", name = "Minimum annual \ntemperature (°C)\n")
g_to <- g_to +
  theme_bw()+
  ylim(-5,5)+
  scale_x_reverse(lim=c(5,-5))+
  theme(legend.title = element_text(size=14), legend.text = element_text(size=11), legend.position=c(0.1, 1))+
  theme(legend.background = element_rect(color = "black", fill = "transparent", size = 2, linetype="blank"))+
  theme(plot.title = element_text(size=15, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black"), axis.text=element_text(size=10), axis.title=element_text(size=12,vjust = -0.5),axis.text.x=element_text(hjust = -0.05),legend.key = element_blank(), plot.margin = unit(c(1,1,1,1), units = , "cm"))+
  theme(legend.direction = 'horizontal', 
        legend.position = 'top')+
  theme(aspect.ratio=1/1)+
  theme(legend.title=element_text(size=12))+
  labs(y="PC2 (21.8% Explained)",x="PC1 (39.2% Explained)")+
  theme(axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"))
g_to$layers <- c(g_to$layers, g_to$layers[[1]])
print(g_to)

pdf(file="...", width = 12, height = 6)
grid.arrange(g6,g_to,ncol=2)
dev.off()

#Whole plane only
pdf(file="...", width = 10, height = 6)
g6
dev.off()

#-----------------------------------------------------------------------------------
#Repeat with trees
#-----------------------------------------------------------------------------------


PCA_trees<-subset(PCA_table6, Species %in% PM_all)
PCA_trees$Group<-"Trees"
PCA_table6<-rbind(PCA_world,PCA_tundra, PCA_trees)

pca.data6<-PCA_table6[,2:7]
pca.species6<-PCA_table6[,1]
pca.groups6<-PCA_table6[,9]

nrow(PCA_table6)

tundra.pca6 <- prcomp(pca.data6,
                      center = TRUE,
                      scale. = TRUE)

library(reshape2)
loadings_6<-as.data.frame(tundra.pca6$rotation[,c(1:2)])
loadings_6$trait<-row.names(loadings_6)
loadings_6<-melt(loadings_6, id.vars="trait", measure.vars=c("PC1","PC2"))
loadings_6$value<-as.numeric(as.character(loadings_6$value))
names(loadings_6)[2]<-"PC"

ggplot(loadings_6,aes(PC,value,fill=trait))+
  geom_bar(stat="identity", position = "dodge")+
  theme_bw()

g6_trees <- ggbiplot_reverse_xy(tundra.pca6, obs.scale = 1, var.scale = 1, ellipse = FALSE, 
                          circle = FALSE, alpha=0.8,group= pca.groups6,varname.size=0,varname.adjust = 1.5,color = "black")
g6_trees <- g6_trees + scale_colour_manual(values=c("forestgreen","blue","grey80"),name = '')+
  theme_bw()+
  scale_x_reverse( lim=c(5,-5))+
  scale_y_reverse( lim=c(5,-5))+
  theme(legend.title = element_text(size=14), legend.text = element_text(size=11), legend.position=c(0.1, 1))+
  theme(legend.background = element_rect(color = "black", fill = "transparent", size = 2, linetype="blank"))+
  theme(plot.title = element_text(size=15, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  axis.text=element_text(size=10), axis.title=element_text(size=11,vjust = -0.5),axis.text.x=element_text(vjust = -0.05),legend.key = element_blank(), plot.margin = unit(c(0.5,0.5,0.5,0.5), units = , "cm"))+
  theme(legend.direction = 'horizontal', 
        legend.position = 'top')+
  theme(aspect.ratio=1/1)+
  theme(axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"))+
  labs(x="PC1 - Plant Size\n(37.5% Variation)",y="PC2 - Resource Economics\n(29.1% Variation)")
g6_trees$layers <- c(g6_trees$layers, g6_trees$layers[[1]])
print(g6_trees)

#Whole plane only
pdf(file="...", width = 8, height = 4)
g6_trees
dev.off()

#-----------------------------------------------------------------------------------
#Repeat without LDMC conversion
#-----------------------------------------------------------------------------------

#For median niche
#species.meanstable_wTTT2$Coord<-niches$mean_NT[match(species.meanstable_wTTT2$Species,niches$i)]
#species.meanstable_wTTT2$Coord<-as.numeric(as.character(species.meanstable_wTTT2$Coord))

names(species.meanstable_wTTT2)[7]<-"LMA"

#PCA all traits
PCA_table6<-species.meanstable_wTTT2[complete.cases(species.meanstable_wTTT2[,2:7]),]
PCA_tundra<-subset(PCA_table6,Group=="Tundra")
PCA_tundra$F_Group<-try.ttt$F_Group[match(PCA_tundra$Species,try.ttt$AccSpeciesName)]
PCA_world<-subset(PCA_table6,Group=="World")
PCA_world$F_Group<-"World"
PCA_table6<-rbind(PCA_world,PCA_tundra)


pca.data6<-PCA_table6[,2:7]
pca.species6<-PCA_table6[,1]
pca.groups6<-PCA_table6[,8]


tundra.pca6 <- prcomp(pca.data6,
                      center = TRUE,
                      scale. = TRUE)

g6 <- ggbiplot_reverse_xy(tundra.pca6, obs.scale = 1, var.scale = 1, ellipse = FALSE, 
                          circle = FALSE, alpha=0.8,group= pca.groups6,varname.size=4.5,varname.adjust = 1.5,color = "black")
g6 <- g6 + scale_colour_manual(values=c("blue","grey80"),name = '')+
  theme_bw()+
  scale_x_reverse( lim=c(5,-5))+
  scale_y_reverse( lim=c(5,-5))+
  theme(legend.title = element_text(size=14), legend.text = element_text(size=11), legend.position=c(0.1, 1))+
  theme(legend.background = element_rect(color = "black", fill = "transparent", size = 2, linetype="blank"))+
  theme(plot.title = element_text(size=15, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  axis.text=element_text(size=10), axis.title=element_text(size=12,vjust = -0.5),axis.text.x=element_text(hjust = -0.05),legend.key = element_blank(), plot.margin = unit(c(1,1,1,1), units = , "cm"))+
  theme(legend.direction = 'horizontal', 
        legend.position = 'top')+
  theme(aspect.ratio=1/1)+
  theme(axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"))+
  labs(x="PC1: Plant Size Spectrum (37.5% Explained)",y="PC2: Plant Economics Spectrum (29.0% Explained)")
g6$layers <- c(g6$layers, g6$layers[[1]])
print(g6)

##############################


#Tundra only#####
PCA_tundra_only_LDMC<-subset(PCA_table6,Group=="Tundra")

#Add class
PCA_tundra_only_LDMC$Class<-PCA_tundra_only$min_WQ_cat[match(PCA_tundra_only_LDMC$Species,PCA_tundra_only$Species)]


pca.data_to<-PCA_tundra_only_LDMC[,2:7]
pca.species_to<-PCA_tundra_only_LDMC[,1]
pca.groups_to<-PCA_tundra_only_LDMC[,10]

tundra_only.pca <- prcomp(pca.data_to,
                          center = TRUE,
                          scale. = TRUE)

summary(tundra_only.pca)
summary(tundra.pca6)

g_to <- ggbiplot_flip(tundra_only.pca, obs.scale = 1, var.scale = 1, ellipse = FALSE,
                      circle = FALSE,group= as.factor(pca.groups_to),alpha=0.8,varname.size=4.5,varname.adjust = 1.5,color = "black")
g_to <- g_to + scale_colour_manual(values=c("blue2","darkgoldenrod1","red2"),
                                   labels = c("Cold", "Mid", "Warm"),
                                   name = 'Temperature Class:')
g_to <- g_to +
  theme_bw()+
  xlim(-5,5)+
  ylim(-5,5)+
  theme(legend.title = element_text(size=14), legend.text = element_text(size=11), legend.position=c(0.1, 1))+
  theme(legend.background = element_rect(color = "black", fill = "transparent", size = 2, linetype="blank"))+
  theme(plot.title = element_text(size=15, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black"), axis.text=element_text(size=10), axis.title=element_text(size=11,vjust = -0.5),axis.text.x=element_text(hjust = -0.05),legend.key = element_blank(), plot.margin = unit(c(1,1,1,1), units = , "cm"))+
  theme(legend.direction = 'horizontal', 
        legend.position = 'top')+
  theme(aspect.ratio=1/1)+
  theme(legend.title=element_text(size=11))+
  labs(y="PC1 - Plant Economics (36.9% Variation Explained)",x="PC2 - Plant Size (27.0% Variation Explained)")+
  theme(axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"))
g_to$layers <- c(g_to$layers, g_to$layers[[1]])
print(g_to)

pdf(file="...", width = 12, height = 6)
grid.arrange(g6,g_to,ncol=2)
dev.off()