##' Sort and truncate predictors according to the strength of predictor-environment interaction
##'
##' @title Sort and truncate predictors according to the strength of predictor-environment interaction
##' @param x A data matrix (raw: samples, col: predictors).
##' @param y A vector of an environment in which the samples were collected.
##' @param method A string to specify the method of regression for calculating R-squared values.
##' "linear" (default), "quadratic" or "cubic" regression model can be specified.
##' @param n.pred The number of predictors to be included in PLORN model (default: ncol(x)).
##' @param trunc a threshold to be truncated (default: 1).
##' @importFrom stats lm
##' @return A data matrix (raw: samples, col: sorted predictors)
##' @examples
##' data(Pinus)
##' train <- p.clean(Pinus$train)
##' target <- Pinus$target
##' cor(target, train[, 1])
##'
##' train <- p.sort(train, target, trunc = 0.5)
##' cor(target, train[, 1])
##' @author Takahiko Koizumi
##' @export
p.sort <- function(x, y, method = "linear", n.pred = ncol(x), trunc = 1){
  degree <- switch(method,
                   "linear" = 1,
                   "quadratic" = 2,
                   "cubic" = 3,
                   stop("Select the <method> linear, quadratic, or cubic")
  )

  if(n.pred < ncol(x) & trunc < 1){
    stop("Don't specify <n.pred> and <trunc> at a time")
  }
  if(n.pred < 0){
    stop("<n.pred> should not be a negative value")
  }else if(n.pred > ncol(x)){
    stop(paste("<n.pred> must not exceed", ncol(x), sep = " "))
  }
  if(trunc < 0 | trunc > 1){
    stop("<trunc> should be within the range of 0-1")
  }

  ## calculate R2 values
  result <- rep(NA, ncol(x))
  for(i in 1:ncol(x)){
    result[i] <- summary(lm(y ~ poly(x[, i], degree = degree, raw = TRUE)))$r.squared
  }
  ## sort predictors in descending order of R2
  x <- x[, order(result, decreasing = TRUE)]

  ## extract predictors with higher R2
  if(n.pred <= ncol(x)){
    x <- x[, 1:n.pred]
  }else if(trunc < 1){
    x <- x[, 1:length(result[result >= trunc])]
  }else{
    x <- x
  }
}
