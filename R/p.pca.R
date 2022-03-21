##' Visualize predictors using principal coordinate analysis
##'
##' @title Visualize predictors using principal coordinate analysis
##' @param x A data matrix (row: samples, col: predictors).
##' @param y A vector of an environment in which the samples were collected.
##' @param method A string to specify the method of regression for calculating R-squared values.
##' "linear" (default), "quadratic" or "cubic" regression model can be specified.
##' @param lower.thr The lower threshold of R-squared value to be indicated in a PCA plot (default: 0).
##' @param n.pred The number of candidate predictors for PLORN model to be indicated in a PCA plot (default: ncol(x)).
##' @param size The size of symbols in a PCA plot (default: 1).
##' @importFrom ggplot2 ggplot aes geom_point xlab ylab scale_fill_manual scale_fill_gradient theme element_rect element_blank
##' @importFrom stats lm princomp
##' @return A PCA plot
##' @examples
##' data(Pinus)
##' train <- p.clean(Pinus$train)
##' target <- Pinus$target
##' p.pca(train, target)
##' @author Takahiko Koizumi
##' @export
p.pca <- function(x, y, method = "linear", lower.thr = 0, n.pred = ncol(x), size = 1){
  PC1 <- PC2 <- R2 <- selection <- NULL

  pc <- princomp(t(scale(x)), cor = TRUE, scores = TRUE)
  coords <- pc$scores
  sd <- summary(pc)$sdev
  pc1 <- (sd ^ 2 / sum(sd ^ 2))[1]
  pc2 <- (sd ^ 2 / sum(sd ^ 2))[2]

  degree <- switch(method,
                   "linear" = 1,
                   "quadratic" = 2,
                   "cubic" = 3,
                   stop("Select the methods linear, quadratic, or cubic")
  )

  if(lower.thr > 0 & n.pred < ncol(x)){
    stop("Don't specify <lower.thr> and <n.pred> at a time")
  }
  if(lower.thr < 0 | lower.thr > 1){
    stop("<lower.thr> should be within the range of 0-1")
  }
  if(n.pred < 0){
    stop("<n.pred> should not be a negative value")
  }else if(n.pred > ncol(x)){
    stop(paste("<n.pred> must not exceed", ncol(x), sep = " "))
  }

  ## calculate R2 values
  r2 <- rep(NA, ncol(x))
  rc <- rep(NA, ncol(x))
  for(i in 1:ncol(x)){
    rs <- summary(lm(y ~ poly(x[, i], degree = degree, raw = TRUE)))
    r2[i] <- rs$r.squared
    rc[i] <- rs$coefficients[2, 1]
  }

  h <- data.frame(
    R2 = r2,
    Coef = rc,
    selection = rep(NA, ncol(x))
  )
  row.names(h) <- colnames(x)
  coords <- cbind(coords[, 1:2], h)
  colnames(coords) <- c("PC1", "PC2", "R2", "Coef", "selection")
  coords <- coords[order(r2, decreasing = FALSE), ]
  if(lower.thr > 0){
    coords[coords$R2 >= lower.thr, "selection"] <- "selected"
    coords[coords$R2 < lower.thr, "selection"] <- "not"
  }else if(n.pred < ncol(x)){
    coords[, "selection"] <- "selected"
    coords[1:(nrow(coords) - n.pred), "selection"] <- "not"
  }else{
    coords <- coords
  }

  if(lower.thr > 0 | n.pred < ncol(x)){
    ggplot(coords, aes(x = PC1, y = PC2)) +
      geom_point(aes(fill = selection), colour = "white", shape = 21, alpha = 0.5, stroke = 0.1, size = size) +
      scale_fill_manual(values = c("gray", "red")) +
      xlab(paste("PC1 (", format(round(pc1 * 100, digits = 1), nsmall = 1), "%)", sep = "")) +
      ylab(paste("PC2 (", format(round(pc2 * 100, digits = 1), nsmall = 1), "%)", sep = "")) +
      theme(panel.background = element_rect(fill = "transparent", colour = "black"),
            panel.grid = element_blank(),
            strip.background = element_blank(),
            legend.position = "none")
  }else{
    coords[coords$Coef < 0, "R2"] <- -coords[coords$Coef < 0, "R2"]
    ggplot(coords, aes(x = PC1, y = PC2)) +
      geom_point(aes(fill = R2), colour = "white", shape = 21, alpha = 0.5, stroke = 0.1, size = size) +
      scale_fill_gradient(low = "gray", high = "red") +
      xlab(paste("PC1 (", format(round(pc1 * 100, digits = 1), nsmall = 1), "%)", sep = "")) +
      ylab(paste("PC2 (", format(round(pc2 * 100, digits = 1), nsmall = 1), "%)", sep = "")) +
      theme(panel.background = element_rect(fill = "transparent", colour = "black"),
            panel.grid = element_blank(),
            strip.background = element_blank(),
            legend.position = "none")
  }
}



