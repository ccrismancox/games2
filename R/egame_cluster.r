egame12.Cluster <- function(model, cluster, cdjust=TRUE){
  bread <- vcov(model)
  cl <- model$call
  
  link <- model$link
  type <- model$type
  
  ## various sanity checks
  formulas <- model$formula
  ## make the model frame
  mf <- match(c("data", "subset", "na.action"), names(cl), 0L)
  mf <- cl[c(1L, mf)]
  mf$formula <- formulas
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  
  mfFull <- match(c("data", "subset"), names(cl), 0L)
  mfFull <- cl[c(1L, mfFull)]
  mfFull$formula <- formulas
  mfFull$drop.unused.levels <- TRUE
  mfFull$na.action <- "na.pass"
  mfFull[[1]] <- as.name("model.frame")
  mfFull <- eval(mfFull, parent.frame())
  
  clust.idx <- which( !(row.names(mfFull) %in% row.names(mf)))
  
  ## get the response and store it as factor (yf) and numeric (y)
  yf <- model.part(formulas, mf, lhs = 1, drop = TRUE)
  yf <- games:::makeResponse12(yf)
  y <- as.numeric(yf)
  
  
  
  cluster <- cluster[-clust.idx]
  G <- length(unique(cluster))
  
  ## makes a list of the 4 (or more, if variance formulas specified) matrices
  ## of regressors to be passed to estimation functions
  regr <- list()
  for (i in seq_len(length(formulas)[2]))
    regr[[i]] <- model.matrix(formulas, data = mf, rhs = i)
  rcols <- sapply(regr, ncol)
  
  b <- model$coefficients
  JAC <- games:::logLikGrad12(b, y, regr, link, type)
  meat <- Reduce(`+`,lapply(split.data.frame(JAC, f=cluster, drop=FALSE), crossprod))
  clusterVCOV <- bread %*% meat %*% bread
  if(cadjust){
    clusterVCOV <- clusterVCOV*(G/(G-1))
  }
  return(clusterVCOV)
  
  
  
}

