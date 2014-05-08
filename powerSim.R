require(MASS)
require(data.table)

#' @title This module performs a statistical power simulation for a non-binomial generalized linear model.
#' @author Jens Hooge
#' 
#' what is the interesting effect size
#' the estimator would relate to a 1% change of mutation load
#' so eg. exp(0.01)=1.01005, so an increase by 50% would be 1.01005^50=1.648708
#' i.e. comparing an indivdual with 0% AKT1mut and 100 reads for gene X to one with 50% AKTmut, the latter
#' has 164 reads of gene X
#' power as fx of effect size

#' Generates a feature matrix with random values.
#' Each column represents a feature (eg. gene) that is either non-mutated(x_ij = 0)
#' or mutated to some degree (eg.: 0 <= x_ij <= 100).
#'
#' @param nrow Number of rows
#' @param ncol Number of columns
#' @param dist Distribution to draw random samples from
#' 
#' @return X matrix{base}
generateFeatureMatrix <- function(nrow, ncol, dist=rnorm(n=nrow*ncol,mean=0,sd=1), seed=42) {
  set.seed(seed)
  X <- matrix(dist, nrow=nrow, ncol=ncol)
  return(X)
}

#' Generates random coefficients.
#'
#' @param dist Distribution to draw random samples from
#' 
#' @return beta numeric{base}
generateCoefficients <- function(dist, seed=42) {
  set.seed(seed)
  beta <- dist
  return(beta)
}

#' Generates the response vector
#' 
#' @param X The feature matrix [n x m]
#' @param beta The vector of coefficients
#' @param epsilon A scalar or vector defining additional noise
#' 
#' @return y matrix{base}
response <- function(X, beta, epsilon=0) {
  y <- ceiling((X %*% beta) + epsilon)
  return(y)
}

#' Takes a list of models and generates a data frame with the 
#' combined results of all models that have converged. Diverged
#' The resulting data frame includes a column with model indices,
#' the set of parameters, the model has been created with and 
#' a summary statistic of the model. The summary statistic consists
#' of the estimate, std error a z value and P(x > |z|).
#' 
#' @param models A list of summarized gbm.nb models
#' 
#' @return result data.frame{base}
postprocess <- function(models) {
  result <- data.frame()
  
  for(i in 1:length(models)) {
    model <- models[[i]]
    if(class(model)[2] != "error") {
      modelResult <- data.frame()
      
      ## remove rownames from coefficient data frame and 
      ## add these rownames as first column
      coefficients <- as.data.frame(model$coefficients)
      featureNames <- rownames(coefficients)
      coefficients <- cbind(featureNames, coefficients)
      rownames(coefficients) <- NULL
      
      ## replicate parameters, used to build this model, to add to result matrix
      modelParams <- rep(as.numeric(params[i, ]), nrow(coefficients))
      modelParams <- as.data.frame(matrix(modelParams, nrow=nrow(coefficients), byrow=TRUE))
      colnames(modelParams) <- colnames(params[i, ])
      
      ## replicate model number
      modelID <- rep(i, length(featureNames))
      
      # combine everythig in a resulting data frame
      modelResult <- cbind(modelID, modelParams)
      modelResult <- cbind(modelResult, coefficients)    
      result <- rbind(result, modelResult)
    }
  }
  return(result)
}

#' Computes the negative binomial generalized linear model
#' using variable effectSize. It uses helperfunctions \code{generateFeatureMatrix},
#' \code{generateCoefficients}, \code{response} to compute the respective
#' parts of the generalized linear model and returns the result data.frame of 
#' glm.nb
#'
#' @return out data.frame{base}
run <- function(effectSize, n, m, seed=42) {
  X <- generateFeatureMatrix(n, m, runif(n*m, min=0, max=100), seed=seed)
  beta <- generateCoefficients(rnorm(m, mean=0, sd=effectSize), seed=seed)
  y <- response(X, beta)
  
  out <- tryCatch({
    summary(glm.nb(y~X))
    },error=function(e){e})
  
  return(out)
}


###### One Example ######
# n <- 500
# m <- 2
# effectSize <- 0.005
# 
# X <- generateFeatureMatrix(n, m, runif(n*m, min=0, max=100))
# beta <- generateCoefficients(rnorm(m, mean=0, sd=effectSize))
# y <- response(X, beta)
# 
# out <- tryCatch({summary(glm.nb(as.numeric(y)~X))},error=function(e){e})
# out

###### expand.grid Example ######
seed=1337
params <- expand.grid(effectSize=seq(0.005, 0.030, 0.005), 
                      n=c(10, 20, 40, 80, 160, 320),
                      m=10)

models <- mapply(run, 
              effectSize=params$effectSize, 
              n=params$n, 
              m=params$m,
              seed=seed)

result <- postprocess(models)    
