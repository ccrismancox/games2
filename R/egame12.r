##' @import utils
##' @import graphics
##' @import stats
##' @export
predict.egame12 <- function(object, newdata, type = c("outcome", "action"),
                            na.action = na.pass, ...)
{
  type <- match.arg(type)
  
  if (missing(newdata) || is.null(newdata)) {
    ## use original data if 'newdata' not supplied
    mf <- object$model
  } else {
    ## get rid of left-hand variables in the formula, since they're not
    ## needed for fitting
    formulas <- Formula(delete.response(terms(formula(object$formulas))))
    
    mf <- model.frame(formulas, data = newdata, na.action = na.action,
                      xlev = object$xlevels)
    
    ## check that variables are of the right classes
    Terms <- attr(object$model, "terms")
    if (!is.null(cl <- attr(Terms, "dataClasses")))
      .checkMFClasses(cl, mf)
  }
  
  regr <- list()
  for (i in seq_len(length(object$formulas)[2]))
    regr[[i]] <- model.matrix(object$formulas, data = mf, rhs = i)
  
  ## get action probabilities, as given by fitted model parameters
  ans <- makeProbs12(object$coefficients, regr = regr, link = object$link, type
                     = object$type)
  ans <- do.call(cbind, ans)
  
  if (type == "outcome") {
    ans <- data.frame(cbind(ans[, 1], ans[, 2] * ans[, 3],
                            ans[, 2] * ans[, 4]))
    names(ans) <- paste("Pr(", levels(object$y), ")", sep = "")
  } else {
    ans <- as.data.frame(ans)
    names(ans)[1:2] <- paste("Pr(", c("", "~"), levels(object$y)[1], ")",
                             sep = "")
    names(ans)[3:4] <- paste("Pr(", levels(object$y)[2:3], "|~",
                             levels(object$y)[1], ")", sep = "")
    names(ans) <- gsub("~~", "", names(ans))
  }
  
  return(ans)
}

sbi12 <- function(y, regr, link)
{
  names(regr) <- character(length(regr))
  names(regr)[1:4] <- c("X1", "X3", "X4", "Z")
  
  ## have to do this because binomial() issues warning if it's not directly
  ## passed a character string to its link argument
  if (link == "probit") {
    fam <- binomial(link = "probit")
  } else {
    fam <- binomial(link = "logit")
  }
  
  ## regression for player 2's choice
  Z2 <- regr$Z[y != 1, ]
  y2 <- as.numeric(y == 3)[y != 1]
  m2 <- suppressWarnings(glm.fit(Z2, y2, family = fam))
  p4 <- as.numeric(regr$Z %*% coef(m2))
  p4 <- if (link == "probit") pnorm(p4) else plogis(p4)
  
  ## regression for player 1's choice
  X1 <- cbind(-regr$X1, (1 - p4) * regr$X3, p4 * regr$X4)
  y1 <- as.numeric(y != 1)
  m1 <- suppressWarnings(glm.fit(X1, y1, family = fam))
  
  ## need to multiply by sqrt(2) because the standard glm assumes dispersion
  ## parameter 1, but agent error implies dispersion parameter sqrt(2)
  ans <- sqrt(2) * c(coef(m1), coef(m2))
  return(ans)
}

makeSDs12 <- function(b, regr, type)
{
  sds <- vector("list", 5)
  rcols <- sapply(regr, ncol)
  
  if (length(rcols) == 5) {  ## sdByPlayer == FALSE
    v <- exp(as.numeric(regr[[5]] %*% b))
    for (i in 1:5) sds[[i]] <- v
  } else {
    v1 <- exp(as.numeric(regr[[5]] %*% b[1:rcols[5]]))
    v2 <- exp(as.numeric(regr[[6]] %*% b[(rcols[5]+1):length(b)]))
    if (type == "agent") {
      sds[[1]] <- sds[[2]] <- v1
      sds[[3]] <- sds[[4]] <- v2
    } else {
      sds[[1]] <- sds[[2]] <- sds[[3]] <- v1
      sds[[4]] <- sds[[5]] <- v2
    }
  }
  
  return(sds)
}

makeProbs12 <- function(b, regr, link, type)
{
  utils <- makeUtils(b, regr, nutils = 4,
                     unames = c("u11", "u13", "u14", "u24"))
  
  ## length(utils$b) == 0 means no terms left for the variance components, so
  ## set these to 1
  if (length(utils$b) == 0) {
    sds <- as.list(rep(1, 5))
  } else {
    sds <- makeSDs12(utils$b, regr, type)
  }
  
  linkfcn <- switch(link,
                    logit = function(x, sd = 1) plogis(x, scale = sd),
                    probit = pnorm)
  
  if (type == "private") {
    sd4 <- sqrt(sds[[4]]^2 + sds[[5]]^2)
  } else {
    sd4 <- sqrt(sds[[3]]^2 + sds[[4]]^2)
  }
  p4 <- finiteProbs(linkfcn(utils$u24, sd = sd4))
  p3 <- 1 - p4
  
  if (type == "private") {
    sd2 <- sqrt(p3^2 * sds[[2]]^2 + p4^2 * sds[[3]]^2 + sds[[1]]^2)
  } else {
    sd2 <- sqrt(sds[[1]]^2 + sds[[2]]^2)
  }
  p2 <- p3 * utils$u13 + p4 * utils$u14 - utils$u11
  p2 <- finiteProbs(linkfcn(p2, sd = sd2))
  p1 <- 1 - p2
  
  return(list(p1 = p1, p2 = p2, p3 = p3, p4 = p4))
}

actionsToOutcomes12 <- function(probs, log.p = TRUE)
{
  probs <- do.call(cbind, probs)
  ans <- cbind(log(probs[, 1]),
               log(probs[, 2]) + log(probs[, 3]),
               log(probs[, 2]) + log(probs[, 4]))
  
  if (!log.p) ans <- exp(ans)
  return(ans)
}

logLik12Base <- function(b, y, regr, link, type, ...)
{
  names(regr) <- character(length(regr))
  names(regr)[1:4] <- c("X1", "X3", "X4", "Z")
  probs <- makeProbs12(b, regr, link, type)
  logProbs <- actionsToOutcomes12(probs, log.p = TRUE)
  ans <- logProbs[cbind(1:nrow(logProbs), y)]
  return(sum(ans))
}

logLik12 <- function(b, y, regr, link, type, logF, Cauchy, Firth, FirthExtra,...)
{
  names(regr) <- character(length(regr))
  names(regr)[1:4] <- c("X1", "X3", "X4", "Z")
  probs <- makeProbs12(b, regr, link, type)
  logProbs <- actionsToOutcomes12(probs, log.p = TRUE)
  ans <- logProbs[cbind(1:nrow(logProbs), y)]
  if(Cauchy){
    # ans <- sum(ans,  0.5* determinant(-numDeriv::hessian(func=logLik12Base, b, y=y, regr=regr, link=link, type=type))$m, na.rm=TRUE)
    # ans <- sum(ans,sum(dcauchy(b,scale=2.5, log = TRUE)))
    scale <- ifelse(str_detect(names(b), "Intercept"), 10,2.5)
    ans <- sum(ans,sum(-log((1+ (b/scale)^2))))
  }
  if(Firth){
    u <-  makeUtils(b, regr, nutils = 4,
                    unames = c("u11", "u13", "u14", "u24"))
    A <- -hessian12(b, y=y, regr=regr, link=link, type=type, FirthExtra=FirthExtra, p=probs, u=u)
    D <- determinant(A)
    if(D$sign<0){warning("The determinant of the Hessian is negative.  Consider switching to Cauchy.")}
    D <-  D$m/2
    ans <- sum(ans, D, na.rm=TRUE)
  }
  if(logF){ #b*m/2 - m* log(1+exp(b))
    m <- ifelse(str_detect(names(b), "Intercept"),1,1)
    ans <- sum(ans,sum(b*m/2 - m* log(1+exp(b))))
    
  }
  return(sum(ans))
}

logLikGrad12 <- function(b, y, regr, link, type, ...)
{
  names(regr) <- character(length(regr))
  names(regr)[1:4] <- c("X1", "X3", "X4", "Z")
  u <- makeUtils(b, regr, nutils = 4,
                 unames = c("u11", "u13", "u14", "u24"))
  p <- makeProbs12(b, regr, link, type)
  n <- nrow(regr$Z)
  rcols <- sapply(regr, ncol)
  
  if (link == "probit" && type == "private") {
    dp4db <- matrix(0L, nrow = n, ncol = sum(rcols[1:3]))
    dp4dg <- dnorm(u$u24 / sqrt(2)) * regr$Z / sqrt(2)
    dp4 <- cbind(dp4db, dp4dg)
    dp3 <- -dp4
    
    num2 <- p$p3 * u$u13 + p$p4 * u$u14 - u$u11
    denom2 <- sqrt(1 + p$p3^2 + p$p4^2)
    dn2 <- dnorm(num2 / denom2)
    dp2db1 <- dn2 * (-regr$X1 / denom2)
    dp2db3 <- dn2 * p$p3 * regr$X3 / denom2
    dp2db4 <- dn2 * p$p4 * regr$X4 / denom2
    dp2dg <- dn2 * ((u$u14-u$u13)*denom2 - (p$p4-p$p3)*num2/denom2)
    dp2dg <- (dp2dg * dp4dg) / (denom2^2)
    Dp2 <- cbind(dp2db1, dp2db3, dp2db4, dp2dg)
    Dp1 <- -Dp2
    Dp3 <- dp3
    Dp4 <- dp4
  } else if (type == "agent") {
    ## unfortunately these are not as consistent in the notation as the
    ## gradients in egame122, egame123, and ultimatum -- I intend to fix
    ## this for the sake of future maintainers' sanity, but the ugliness
    ## here doesn't cause any user-visible problems
    derivCDF <- switch(link,
                       logit = dlogis,
                       probit = dnorm)
    
    dp4 <- derivCDF(u$u24 / sqrt(2)) / sqrt(2)
    Dp4 <- cbind(matrix(0L, nrow = nrow(regr$X1), ncol = sum(rcols[1:3])),
                 dp4 * regr$Z)
    Dp3 <- -Dp4
    
    dp1 <- derivCDF((u$u11 - p$p3 * u$u13 - p$p4 * u$u14) / sqrt(2)) /
      sqrt(2)
    dbp1 <- dp1 * cbind(regr$X1, -p$p3 * regr$X3, -p$p4 * regr$X4)
    dgp1 <- dp1 * dp4 * (u$u13 - u$u14) * regr$Z
    Dp1 <- cbind(dbp1, dgp1)
    Dp2 <- -Dp1
  }
  
  dL1 <- Dp1 / p$p1
  dL3 <- Dp2 / p$p2 + Dp3 / p$p3
  dL4 <- Dp2 / p$p2 + Dp4 / p$p4
  
  ans <- matrix(NA, nrow = n, ncol = sum(rcols[1:4]))
  ans[y == 1, ] <- dL1[y == 1, ]
  ans[y == 2, ] <- dL3[y == 2, ]
  ans[y == 3, ] <- dL4[y == 3, ]
  
  return(ans)
}

makeResponse12 <- function(yf)
{
  if (length(dim(yf))) {  # response specified as dummies
    Y <- yf
    if (ncol(Y) > 2)
      warning("only first two columns of response will be used")
    
    Y <- Y[, 1:2]
    if (!all(unlist(yf) %in% c(0L, 1L)))
      stop("dummy responses must be dummy variables")
    
    y <- numeric(nrow(Y))
    y[Y[, 1] == 0] <- 1
    y[Y[, 1] == 1 & Y[, 2] == 0] <- 2
    y[Y[, 1] == 1 & Y[, 2] == 1] <- 3
    yf <- as.factor(y)
    levels(yf) <- c(paste("~", names(Y)[1], sep = ""),
                    paste(names(Y)[1], ",~", names(Y)[2], sep = ""),
                    paste(names(Y)[1], ",", names(Y)[2], sep = ""))
  } else {
    yf <- as.factor(yf)
    if (nlevels(yf) != 3) stop("dependent variable must have three values")
  }
  
  return(yf)
}

##' Strategic model with 2 players, 3 terminal nodes
##' 
##' Fits a strategic model with two players and three terminal nodes, as in the
##' game illustrated below in "Details".
##'
##' The model corresponds to the following extensive-form game, described in
##' Signorino (2003):
##' \preformatted{
##' .     1
##' .     /\
##' .    /  \
##' .   /    \ 2
##' .  u11   /\
##' .       /  \
##' .      /    \
##' .    u13    u14
##' .    0      u24}
##' 
##' If Player 1 chooses L, the game ends and Player 1 receives payoffs of u11.
##' (Player 2's utilities in this case cannot be identified in a statistical
##' model.)  If Player 1 chooses L, then Player 2 can choose L, resulting in
##' payoffs of u13 for Player 1 and 0 for Player 2, or R, with payoffs of u14
##' for 1 and u24 for 2.
##'
##' The four equations specified in the function's \code{formulas} argument
##' correspond to the regressors to be placed in u11, u13, u14, and u24
##' respectively.  If there is any regressor (including the constant) placed in
##' all of u11, u13, and u14, \code{egame12} will stop and issue an error
##' message, because the model is then unidentified (see Lewis and Schultz
##' 2003).  There are two equivalent ways to express the formulas passed to this
##' argument.  One is to use a list of four formulas, where the first contains
##' the response variable(s) (discussed below) on the left-hand side and the
##' other three are one-sided.  For instance, suppose:
##' \itemize{
##' \item u11 is a function of \code{x1}, \code{x2}, and a constant
##' \item u13 is set to 0
##' \item u14 is a function of \code{x3} and a constant
##' \item u24 is a function of \code{z} and a constant.}
##' The list notation would be \code{formulas = list(y ~ x1 + x2, ~ 0, ~ x3, ~
##' z)}.  The other method is to use the \code{\link{Formula}} syntax, with one
##' left-hand side and four right-hand sides (separated by vertical bars).  This
##' notation would be \code{formulas = y ~ x1 + x2 | 0 | x3 | z}.
##'
##' To fix a utility at 0, just use \code{0} as its equation, as in the example
##' just given.  To estimate only a constant for a particular utility, use
##' \code{1} as its equation.
##'
##' There are three equivalent ways to specify the outcome in \code{formulas}.
##' One is to use a numeric vector with three unique values, with their values
##' (from lowest to highest) corresponding with the terminal nodes of the game
##' tree illustrated above (from left to right).  The second is to use a factor,
##' with the levels (in order as given by \code{levels(y)}) corresponding to the
##' terminal nodes.  The final way is to use two indicator variables, with the
##' first standing for whether Player 1 moves L (0) or R (1), the second
##' standing for Player 2's choice if Player 1 moves R.  (The values of the
##' second when Player 1 moves L should be set to 0 or 1, \strong{not}
##' \code{NA}, in order to ensure that observations are not dropped from the
##' data when \code{na.action = na.omit}.)  The way to specify \code{formulas}
##' when using indicator variables is, for example, \code{y1 + y2 ~ x1 + x2 | 0
##' | x3 | z}.
##'
##' If \code{fixedUtils} or \code{sdformula} is specified, the estimated
##' parameters will include terms labeled \code{log(sigma)} (for probit links)
##' or \code{log(lambda)}.  These are the scale parameters of the stochastic
##' components of the players' utility.  If \code{sdByPlayer} is \code{FALSE},
##' then the variance of error terms (or the equation describing it, if
##' \code{sdformula} contains non-constant regressors) is assumed to be common
##' across all players.  If \code{sdByPlayer} is \code{TRUE}, then two variances
##' (or equations) are estimated: one for each player.  For more on the
##' interpretation of the scale parameters in these models and how it differs
##' between the agent error and private information models, see Signorino
##' (2003).
##'
##' The model is fit using \code{\link{maxLik}}, using the BFGS optimization
##' method by default (see \code{\link{maxBFGS}}).  Use the \code{method}
##' argument to specify an alternative from among those supplied by
##' \code{maxLik}.
##' @param formulas a list of four formulas, or a \code{\link{Formula}} object
##' with four right-hand sides.  See "Details" and the examples below.
##' @param data a data frame containing the variables in the model.
##' @param subset optional logical expression specifying which observations from
##' \code{data} to use in fitting.
##' @param na.action how to deal with \code{NA}s in \code{data}.  Defaults to
##' the \code{na.action} setting of \code{\link{options}}.  See
##' \code{\link{na.omit}}.
##' @param link whether to use a probit (default) or logit link structure,
##' @param type whether to use an agent-error ("agent", default) or
##' private-information ("private") stochastic structure.
##' @param startvals whether to calculate starting values for the optimization
##' using statistical backwards induction ("sbi", default), draw them from a
##' uniform distribution ("unif"), or to set them all to 0 ("zero")
##' @param fixedUtils numeric vector of values to fix for u11, u13, u14, and u24
##' respectively.  \code{NULL} (the default) indicates that these should be
##' estimated with regressors rather than fixed.
##' @param sdformula an optional list of formulas or a \code{\link{Formula}}
##' containing a regression equation for the scale parameter.  The formula(s)
##' should have nothing on the left-hand side; the right-hand side should have
##' one equation if \code{sdByPlayer} is \code{FALSE} and two equations if
##' \code{sdByPlayer} is \code{TRUE}.  See the examples below for how to specify
##' \code{sdformula}.
##' @param sdByPlayer logical: if scale parameters are being estimated (i.e.,
##' \code{sdformula} or \code{fixedUtils} is non-\code{NULL}), should a separate
##' one be estimated for each player?  This option is ignored unless
##' \code{fixedUtils} or \code{sdformula} is specified.
##' @param boot integer: number of bootstrap iterations to perform (if any).
##' @param bootreport logical: whether to print status bar when performing
##' bootstrap iterations.
##' @param profile output from running \code{\link{profile.game}} on a previous
##' fit of the model, used to generate starting values for refitting when an
##' earlier fit converged to a non-global maximum.
##' @param method character string specifying which optimization routine to use
##' (see \code{\link{maxLik}})
##' @param penalty type of penalty to use for penalized maximum likelihood estimation. 
##' efault is \code{"none"} (no penalty is used). Options include 
##' \item{\code{"Firth"}}{for Firth's (1993) Jeffreys prior penalty}}
##' \item{\code{"Cauchy"}}{for a Cauchy(0,2.5) penalty on non-intercepts and a Cauchy(0, 10) for intercepts}
##' \item{\code{"logF"}}{for a log-F(1,1) penalty}
##' See Crisman-Cox, Gasparyan, and Signorino (2022) for more details on the use of penalized estimation.
##' @param ... other arguments to pass to the fitting function (see
##' \code{\link{maxLik}}).
##' @return An object of class \code{c("game", "egame12")}. A
##' \code{game} object is a list containing: \describe{
##' \item{\code{coefficients}}{estimated parameters of the model.}
##' \item{\code{vcov}}{estimated variance-covariance matrix.  Cells referring to
##' a fixed parameter (e.g., a utility when \code{fixedUtils} is specified) will
##' contain \code{NA}s.}
##' \item{\code{log.likelihood}}{vector of individual log likelihoods (left
##' unsummed for use with non-nested model tests).}
##' \item{\code{call}}{the call used to produce the model.}
##' \item{\code{convergence}}{a list containing the optimization method used
##' (see argument \code{method}), the number of iterations to convergence, the
##' convergence code and message returned by \code{\link{maxLik}}, and an
##' indicator for whether the (analytic) gradient was used in fitting.}
##' \item{\code{formulas}}{the final \code{Formula} object passed to
##' \code{model.frame} (including anything specified for the scale parameters).}
##' \item{\code{link}}{the specified link function.}
##' \item{\code{type}}{the specified stochastic structure (i.e., agent error or
##' private information).}
##' \item{\code{model}}{the model frame containing all variables used in
##' fitting.}
##' \item{\code{xlevels}}{a record of the levels of any factor regressors.}
##' \item{\code{y}}{the dependent variable, represented as a factor.}
##' \item{\code{equations}}{names of each separate equation (e.g.,
##' "u1(sq)", "u1(cap)", etc.).}
##' \item{\code{fixed}}{logical vector specifying which parameter values, if
##' any, were fixed in the estimation procedure.}
##' \item{\code{boot.matrix}}{if \code{boot} was non-zero, a matrix of bootstrap
##' parameter estimates (otherwise \code{NULL}).}
##' \item{\code{localID}}{an indicator for whether the Hessian matrix is
##' negative definite, a sufficient condition for local identification of the
##' model parameters.}
##' }
##' The second class of the returned object, \code{egame12}, is for use in
##' generation of predicted probabilities.
##' @seealso \code{\link{summary.game}} and \code{\link{predProbs}} for
##' postestimation analysis; \code{\link{makeFormulas}} for formula
##' specification.
##' @export
##' @references Jeffrey B. Lewis and Kenneth A Schultz.  2003.  "Revealing
##' Preferences: Empirical Estimation of a Crisis Bargaining Game with
##' Incomplete Information."  \emph{Political Analysis} 11:345--367.
##'
##' Curtis S. Signorino.  2003.  "Structure and Uncertainty in Discrete Choice
##' Models."  \emph{Political Analysis} 11:316--344.
##' @author Casey Crisman-Cox (\email{ccrismancox@gmail.com}), Brenton Kenkel, and Curtis
##' S. Signorino
##' @example inst/examples/egame12.r
egame12 <- function(formulas, data, subset, na.action,
                    link = c("probit", "logit"),
                    type = c("agent", "private"),
                    startvals = c("sbi", "unif", "zero"),
                    fixedUtils = NULL,
                    sdformula = NULL,
                    sdByPlayer = FALSE,
                    boot = 0,
                    bootreport = TRUE,
                    profile,
                    method = "BFGS",
                    penalty= c("none", "Firth", "Cauchy", "logF"),
                    ...)
{
  cl <- match.call()
  
  link <- match.arg(link)
  type <- match.arg(type)
  penalty <- match.arg(penalty)
  Firth <- (penalty=="Firth")
  Cauchy <- (penalty=="Cauchy")
  logF <- (penalty=="logF")
  
  ## various sanity checks
  formulas <- checkFormulas(formulas)
  if (is.null(fixedUtils) && length(formulas)[2] != 4)
    stop("'formulas' should have four components on the right-hand side")
  
  if (!is.null(fixedUtils)) {
    if (length(fixedUtils) < 4)
      stop("fixedUtils must have 4 elements (u11, u13, u14, u24)")
    if (length(fixedUtils) > 4) {
      warning("only the first 4 elements of fixedUtils will be used")
      fixedUtils <- fixedUtils[1:4]
    }
    
    formulas <- update(formulas, . ~ 1 | 1 | 1 | 1)
    
    if (startvals == "sbi")
      startvals <- "zero"
    
    if (is.null(sdformula))
      sdformula <- if (sdByPlayer) Formula(~ 1 | 1) else Formula(~ 1)
  }
  
  if (!is.null(sdformula)) {
    sdformula <- checkFormulas(sdformula, argname = "sdformula")
    if (sdByPlayer && length(sdformula)[2] != 2)
      stop("'sdformula' should have two components (one for each player) on the right-hand side when sdByPlayer == TRUE")
    if (!sdByPlayer && length(sdformula)[2] != 1)
      stop("'sdformula' should have exactly one component on the right-hand side")
    formulas <- as.Formula(formula(formulas), formula(sdformula))
  }
  
  if (sdByPlayer && is.null(sdformula)) {
    warning("to estimate SDs by player, you must specify 'sdformula' or 'fixedUtils'")
    sdByPlayer <- FALSE
  }
  
  if (link == "logit" && type == "private") {
    warning("logit link cannot be used with private information model; changing to probit link")
    link <- "probit"
  }
  ###Firth Checks
  if(sum(Cauchy, Firth, logF) > 1){
    stop("Chose only one penalty Cauchy, logF, or Firth")
  }
  
  ## make the model frame
  mf <- match(c("data", "subset", "na.action"), names(cl), 0L)
  mf <- cl[c(1L, mf)]
  mf$formula <- formulas
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  
  ## get the response and store it as factor (yf) and numeric (y)
  yf <- model.part(formulas, mf, lhs = 1, drop = TRUE)
  yf <- makeResponse12(yf)
  y <- as.numeric(yf)
  
  ## makes a list of the 4 (or more, if variance formulas specified) matrices
  ## of regressors to be passed to estimation functions
  regr <- list()
  for (i in seq_len(length(formulas)[2]))
    regr[[i]] <- model.matrix(formulas, data = mf, rhs = i)
  rcols <- sapply(regr, ncol)
  
  ## makes starting values -- specify "unif" (a numeric vector of length two)
  ## to control the interval from which uniform values are drawn
  if (missing(profile) || is.null(profile)) {
    if(is.numeric(startvals)){
      sval <- startvals
    }else{
      startvals <- match.arg(startvals)
      if (startvals == "zero") {
        sval <- rep(0, sum(rcols))
      } else if (startvals == "unif") {
        if (!methods::hasArg(unif))
          unif <- c(-1, 1)
        sval <- runif(sum(rcols), unif[1], unif[2])
      } else if (startvals == "sbi") {
        sval <- sbi12(y, regr, link)
        sval <- c(sval, rep(0, sum(rcols) - length(sval)))
      }
    } 
  }else {
    sval <- svalsFromProfile(profile)
  }
  
  ## identification check
  varNames <- lapply(regr, colnames)
  idCheck <- (varNames[[1]] %in% varNames[[2]])
  idCheck <- idCheck & (varNames[[1]] %in% varNames[[3]])
  if (is.null(fixedUtils) && any(idCheck)) {
    stop("Identification problem: the following variables appear in all three of player 1's utility equations: ",
         paste(varNames[[1]][idCheck], collapse = ", "))
  }
  
  ## give names to starting values
  prefixes <- paste(c(rep("u1(", 3), "u2("), c(levels(yf), levels(yf)[3]),
                    ")", sep = "")
  sdterms <- if (!is.null(sdformula)) { if (sdByPlayer) 2L else 1L } else 0L
  utils <- if (is.null(fixedUtils)) 1:4 else numeric(0)
  varNames <- makeVarNames(varNames, prefixes, utils, link, sdterms)
  hasColon <- varNames$hasColon
  names(sval) <- varNames$varNames
  
  ## use the gradient iff no regressors in variance
  gr <- if (is.null(sdformula)) logLikGrad12 else NULL
  if(Cauchy |  Firth | logF){gr <- NULL}
  
  ## deal with fixed utilities
  fvec <- rep(FALSE, length(sval))
  names(fvec) <- names(sval)
  if (!is.null(fixedUtils)) {
    sval[1:4] <- fixedUtils
    fvec[1:4] <- TRUE
  }
  
  # if(Firth){
  names(regr) <- character(length(regr))
  names(regr)[1:4] <- c("X1", "X3", "X4", "Z")
  rcols <- sapply(regr, ncol)
  n <- nrow(regr$Z)
  
  FirthExtra <- list()
  FirthExtra$regr2 <- regr
  FirthExtra$regr2$Z <- cbind(matrix(0,
                                     nrow=n, ncol= sum(rcols[1:3])),
                              regr$Z)
  FirthExtra$regr2$X1 <- cbind(regr$X1,
                               matrix(0,
                                      nrow=n, ncol= sum(rcols[2:4])))
  FirthExtra$regr2$X3 <- cbind(matrix(0,
                                      nrow=n, ncol= sum(rcols[1])),
                               regr$X3,
                               matrix(0,
                                      nrow=n, ncol= sum(rcols[3:4])))
  FirthExtra$regr2$X4 <- cbind(matrix(0,
                                      nrow=n, ncol= sum(rcols[1:2])),
                               regr$X4,
                               matrix(0,
                                      nrow=n, ncol= sum(rcols[4])))
  # if(type=="private"){
  #   cx <- function(x){
  #     tcrossprod(x[1:length(sval)], x[(length(sval)+1):(2*length(sval))])+
  #       tcrossprod(x[(length(sval)+1):(2*length(sval))], x[1:length(sval)])
  #   }
  #   
  #   FirthExtra$x24x24 <- t(apply(FirthExtra$regr2$Z, 1, tcrossprod))
  #   FirthExtra$x11x11 <- t(apply(FirthExtra$regr2$X1, 1, tcrossprod))
  #   FirthExtra$x13x13 <- t(apply(FirthExtra$regr2$X3, 1, tcrossprod))
  #   FirthExtra$x14x14 <- t(apply(FirthExtra$regr2$X4, 1, tcrossprod))
  #   FirthExtra$x13x24 <- t(apply(cbind(FirthExtra$regr2$X3, FirthExtra$regr2$Z), 1, cx))
  #   FirthExtra$x14x24 <- t(apply(cbind(FirthExtra$regr2$X4, FirthExtra$regr2$Z), 1, cx))
  #   FirthExtra$x11x13 <- t(apply(cbind(FirthExtra$regr2$X1, FirthExtra$regr2$X3),1, cx))
  #   FirthExtra$x11x14 <- t(apply(cbind(FirthExtra$regr2$X1, FirthExtra$regr2$X4),1, cx))
  #   FirthExtra$x13x14 <- t(apply(cbind(FirthExtra$regr2$X3, FirthExtra$regr2$X4),1, cx))
  # }
  # }else{
  #   FirthExtra <- NULL
  # }
  # 
  
  results <- maxLik(logLik = logLik12, grad = gr, start = sval, fixed = fvec,
                    method = method, y = y, regr = regr, link = link, type =
                      type, logF=logF, Cauchy=Cauchy, Firth=Firth, FirthExtra=FirthExtra,...)
  
  ## some optimization routines in maxLik have different convergence codes for
  ## success, so we need to retrieve the proper code(s) for the supplied
  ## method and check against it/them
  cc <- convergenceCriterion(method)
  if (!(results$code %in% cc)) {
    warning("Model fitting did not converge\nCode:", results$code,
            "\nMessage: ", results$message)
  }
  
  ## check local identification
  lid <- checkLocalID(results$hessian, fvec)
  if(Cauchy | Firth | logF){
    probs <- makeProbs12(results$estimate, regr, link, type)
    u <-  makeUtils(results$estimate, regr, nutils = 4,
                    unames = c("u11", "u13", "u14", "u24"))
    H <- hessian12(results$estimate, y=y, regr=regr, link=link, type=type, 
                   FirthExtra=FirthExtra, p=probs, u=u)  
  }else{
    H <- results$hessian
  }
  
  if (!lid)
    warning("Hessian is not negative definite; coefficients may not be locally identified")
  
  if (boot > 0) {
    bootMatrix <-
      gameBoot(boot, report = bootreport, estimate = sval, y =
                 y, regr = regr, fn = logLik12, gr = gr, fixed = fvec,
               method = method, link = link, type = type, logF=logF,
               Cauchy = Cauchy, Firth = Firth, FirthExtra=FirthExtra,...)
  }
  
  ## create a 'game' object to store output
  ans <- list()
  ans$coefficients <- results$estimate
  ans$vcov <- getGameVcov(H, fvec)
  rownames(ans$vcov) <- colnames(ans$vcov) <- names(results$estimate)
  
  OPG <- FALSE
  if(anyNA(diag(ans$vcov)) || any(diag(ans$vcov)<0)){
    ans$vcov <- getGameVcov(-crossprod(logLikGrad12(b = results$estimate, y = y, regr = regr,link = link,type = type)), fvec)
    OPG <- TRUE
  }
  
  ans$log.likelihood <- logLik12(results$estimate, y = y, regr = regr, link =
                                   link, type = type, 
                                 logF=logF,
                                 Cauchy = Cauchy,
                                 Firth = Firth, FirthExtra=FirthExtra)
  # ##sanity check on the hessians##
  # ans$finalHessian <- results$hessian
  # 
  # probs <- makeProbs12(results$estimate, regr, link, type)
  # u <-  makeUtils(results$estimate, regr, nutils = 4,
  #                 unames = c("u11", "u13", "u14", "u24"))
  # H <- hessian12(results$estimate, y=y, regr=regr, link=link, type=type, FirthExtra=FirthExtra, p=probs, u=u)
  # ans$estHessian <- H
  # 
  # Grad <- numDeriv::grad(logLik12, results$estimate, y=y, regr=regr, link=link, type=type,
  #                        Cauchy=Cauchy, Firth=Firth, FirthExtra=FirthExtra,...)
  
  ans$call <- cl
  ans$convergence <- list(method = method, iter = nIter(results), code =
                            results$code, message = results$message, gradient =
                            !is.null(gr),
                          Penalty = c(logF,Cauchy, Firth),
                          # gradMaxLik = results$gradient,
                          # gradNumDeriv = Grad,
                          OPG=OPG)
  ans$formulas <- formulas
  ans$link <- link
  ans$type <- type
  ans$model <- mf
  ans$xlevels <- .getXlevels(attr(mf, "terms"), mf)
  ans$y <- yf
  ans$equations <- names(hasColon)
  attr(ans$equations, "hasColon") <- hasColon
  ans$fixed <- fvec
  if (boot > 0)
    ans$boot.matrix <- bootMatrix
  ans$localID <- lid
  
  class(ans) <- c("game", "egame12")
  
  return(ans)
}


hessian12 <- function(b, regr, y, type, link, FirthExtra, p, u){
  regr2 <- FirthExtra$regr2
  if(link=="probit" & type=="agent"){
    in2 <- (-u$u11 + u$u13*(p$p3) + u$u14*p$p4)
    d4 <- dnorm(u$u24/sqrt(2))
    d2 <- dnorm((in2/sqrt(2)))
    Dp4 <- regr2$Z  * d4/sqrt(2)
    
    DDB110 <- -regr2$X1 + regr2$X3*p$p3 + regr2$X4*p$p4 - ((u$u13)* Dp4) +  ((u$u14)* Dp4)
    # DDB11O <- t(apply(DDB110, 1, tcrossprod))
    # 
    # x24x24 <- t(apply(regr2$Z, 1, tcrossprod))
    # x24px24bx24 <- t(apply(regr2$Z, 1, function(x){x %*% t(x) %*%b%*%x}))
    # 
    # 
    # xlist <- list(x13x24=cbind(regr2$X3, regr2$Z),
    #               x14x24=cbind(regr2$X4, regr2$Z))
    # 
    # xlist <- lapply(xlist,
    #                 apply,
    #                 1,
    #                 function(x){
    #                   tcrossprod(x[1:length(b)], x[(length(b)+1):(2*length(b))])+
    #                     tcrossprod(x[(length(b)+1):(2*length(b))], x[1:length(b)])
    #                 }
    # )
    # xlist<-lapply(xlist,t)
    
    
    
    # DDp4Db.base <- ((-sqrt(2)/2) * p$p4 *d4 * x24px24bx24 - x24x24*(d4^2))/(2*p$p4^2) #OLD VECTOR FORMAT
    
    DDp4Db <- t(regr2$Z * (y==3)*drop(((-sqrt(2)/2) * p$p4 *d4 *(regr2$Z%*% b)- (d4^2))/(2*p$p4^2)  )  ) %*% regr2$Z
    
    ### D/Db of (Dp2/Db) / p2 is 3 parts####
    # DD1 <- -d2* in2  *  DDB11O/(2*sqrt(2) * p$p2)
    # DD2 <- ((d4* u$u13 *  x24px24bx24/(2*sqrt(2)) -
    #            d4* u$u14 * x24px24bx24/(2*sqrt(2)) -
    #            sqrt(2)/2 * xlist$x13x24 *d4+
    #            sqrt(2)/2 * xlist$x14x24 *d4)*
    #           dnorm(in2/sqrt(2)))/
    #   (sqrt(2)*p$p2)
    # DD3 <- -((DDB11O * dnorm(in2/sqrt(2))^2) / (2 * p$p2^2))
    # DDp2Db.base <- DD1+ DD2 +DD3
    
    DD1.y2 <-   t(DDB110* (y==2)* drop(-d2* in2/(2*sqrt(2) * p$p2))) %*% DDB110
    DD2.y2 <- ((t(regr2$Z * (y==2)* drop(d4* u$u13 *  (dnorm(in2/sqrt(2))/(sqrt(2)*p$p2)) *  (regr2$Z %*% b)/(2*sqrt(2)))) %*% regr2$Z - 
                  t(regr2$Z* (y==2)* drop(d4* u$u14 *  (dnorm(in2/sqrt(2))/(sqrt(2)*p$p2)) *  (regr2$Z %*% b)/(2*sqrt(2)))) %*% regr2$Z -
                  (t(regr2$X3 * (y==2)* drop(sqrt(2)/2 *  (dnorm(in2/sqrt(2))/(sqrt(2)*p$p2)) *d4)) %*% regr2$Z+
                     t(t(regr2$X3 * (y==2)* drop(sqrt(2)/2 *  (dnorm(in2/sqrt(2))/(sqrt(2)*p$p2)) *d4)) %*% regr2$Z))+
                  (t(regr2$X4 * (y==2)* drop(sqrt(2)/2  *  (dnorm(in2/sqrt(2))/(sqrt(2)*p$p2)) *d4)) %*% regr2$Z+ 
                     t(t(regr2$X4 * (y==2)* drop(sqrt(2)/2 *  (dnorm(in2/sqrt(2))/(sqrt(2)*p$p2)) *d4)) %*% regr2$Z))))
    DD3.y2 <- -t(DDB110* (y==2)* drop( (dnorm(in2/sqrt(2))^2) / (2 * p$p2^2))) %*% DDB110
    DDp2Db.y2 <- DD1.y2+ DD2.y2 +DD3.y2
    
    DD1.y3 <-   t(DDB110* (y==3)* drop(-d2* in2/(2*sqrt(2) * p$p2))) %*% DDB110
    DD2.y3 <- ((t(regr2$Z * (y==3)* drop(d4* u$u13 *  (dnorm(in2/sqrt(2))/(sqrt(2)*p$p2)) *  (regr2$Z %*% b)/(2*sqrt(2)))) %*% regr2$Z - 
                  t(regr2$Z * (y==3)* drop(d4* u$u14 *  (dnorm(in2/sqrt(2))/(sqrt(2)*p$p2)) *  (regr2$Z %*% b)/(2*sqrt(2)))) %*% regr2$Z -
                  (t(regr2$X3 * (y==3)* drop(sqrt(2)/2 *  (dnorm(in2/sqrt(2))/(sqrt(2)*p$p2)) *d4)) %*% regr2$Z+
                     t(t(regr2$X3 * (y==3)* drop(sqrt(2)/2 *  (dnorm(in2/sqrt(2))/(sqrt(2)*p$p2)) *d4)) %*% regr2$Z))+
                  (t(regr2$X4 * (y==3)*drop(sqrt(2)/2  *  (dnorm(in2/sqrt(2))/(sqrt(2)*p$p2)) *d4)) %*% regr2$Z+ 
                     t(t(regr2$X4 * (y==3)* drop(sqrt(2)/2 *  (dnorm(in2/sqrt(2))/(sqrt(2)*p$p2)) *d4)) %*% regr2$Z))))
    DD3.y3 <- -t(DDB110* (y==3)* drop( (dnorm(in2/sqrt(2))^2) / (2 * p$p2^2))) %*% DDB110
    DDp2Db.y3 <- DD1.y3+ DD2.y3 +DD3.y3
    
    
    ##2nd deriv of p3
    # DDp3Db.base <- ((1/(sqrt(2))) * p$p3 * d4 * x24px24bx24 -x24x24 * d4^2)/(2*p$p3^2)
    DDp3Db <- t(regr2$Z * (y==2)*drop( ((1/(sqrt(2))) * p$p3 * d4 * (regr2$Z%*% b) -  d4^2)/(2*p$p3^2) )  ) %*% regr2$Z
    
    ##2nd deriv of p1 is also 3 parts
    # DD1 <- d2 * in2  *  DDB11O/(2*sqrt(2) * p$p1)
    # DD2 <- -((d4* u$u13 *  x24px24bx24/(2*sqrt(2)) -
    #             d4* u$u14 *  x24px24bx24/(2*sqrt(2)) -
    #             sqrt(2)/2 * xlist$x13x24 *d4+
    #             sqrt(2)/2 * xlist$x14x24 *d4
    # ) * dnorm(in2/sqrt(2)))/
    #   (sqrt(2)*p$p1)
    # DD3 <- -(DDB11O * dnorm(in2/sqrt(2))^2) / (2 * p$p1^2)
    # DDp1Db.base <- DD1+ DD2 +DD3
    
    DD1 <-   t(DDB110* (y==1)* drop(d2* in2/(2*sqrt(2) * p$p1))) %*% DDB110
    DD2 <- -((t(regr2$Z * (y==1)* drop(d4* u$u13 *  (dnorm(in2/sqrt(2))/(sqrt(2)*p$p1)) *  (regr2$Z %*% b)/(2*sqrt(2)))) %*% regr2$Z - 
                t(regr2$Z * (y==1)* drop(d4* u$u14 *  (dnorm(in2/sqrt(2))/(sqrt(2)*p$p1)) *  (regr2$Z %*% b)/(2*sqrt(2)))) %*% regr2$Z -
                (t(regr2$X3 * (y==1)* drop(sqrt(2)/2 *  (dnorm(in2/sqrt(2))/(sqrt(2)*p$p1)) *d4)) %*% regr2$Z+
                   t(t(regr2$X3 * (y==1)*drop(sqrt(2)/2 *  (dnorm(in2/sqrt(2))/(sqrt(2)*p$p1)) *d4)) %*% regr2$Z))+
                (t(regr2$X4 * (y==1)* drop(sqrt(2)/2  *  (dnorm(in2/sqrt(2))/(sqrt(2)*p$p1)) *d4)) %*% regr2$Z+ 
                   t(t(regr2$X4 * (y==1)*drop(sqrt(2)/2 *  (dnorm(in2/sqrt(2))/(sqrt(2)*p$p1)) *d4)) %*% regr2$Z))))
    DD3 <- -t(DDB110* (y==1)* drop( (dnorm(in2/sqrt(2))^2) / (2 * p$p1^2))) %*% DDB110
    DDp1Db <- DD1+ DD2 +DD3
    
    # DDB.base <- (y==3)*(DDp4Db.base+DDp2Db.base)+(y==2)*(DDp3Db.base + DDp2Db.base)+(y==1)*(DDp1Db.base)
    DDB <- (DDp4Db+DDp2Db.y3) + (DDp3Db+DDp2Db.y2) +DDp1Db
    
  }else{
    if(type=="agent" & link=="logit"){
      d4 <- dlogis(u$u24/sqrt(2))
      In4 <- ((2 * d4 *p$p4 * exp(-u$u24/sqrt(2))) - d4)
      d2 <- dlogis((-u$u11 + u$u13*p$p3 + u$u14*p$p4)/sqrt(2))
      DDB110 <- -regr2$X1 + regr2$X3*p$p3 + regr2$X4*p$p4 - ((u$u13)* regr2$Z* d4/sqrt(2)) +  ((u$u14)* regr2$Z* d4/sqrt(2))
      DDB11O <- t(apply(DDB110, 1, tcrossprod))
      
      x24x24 <- t(apply(regr2$Z, 1, tcrossprod))
      x24px24bx24 <- t(apply(regr2$Z, 1, function(x){x %*% t(x) %*%b%*%x}))
      
      
      xlist <- list(x13x24=cbind(regr2$X3, regr2$Z),
                    x14x24=cbind(regr2$X4, regr2$Z))
      
      xlist <- lapply(xlist,
                      apply,
                      1,
                      function(x){
                        tcrossprod(x[1:length(b)], x[(length(b)+1):(2*length(b))])+
                          tcrossprod(x[(length(b)+1):(2*length(b))], x[1:length(b)])
                      }
      )
      xlist<-lapply(xlist,t)
      DDp4Db.base <- (In4 * x24x24)/ ( 2*p$p4) - (x24x24* d4^2)/ (2*p$p4^2)
      DDp4Db <- t(regr2$Z * (y==3)*drop( (In4)/( 2*p$p4)- (d4^2)/(2*p$p4^2)  )  ) %*% regr2$Z
      
      ###D2p2/Db is in 3 parts
      DD1 <- ((2 * d2 *p$p2 * exp(-(-u$u11 + u$u13*p$p3 + u$u14*p$p4)/sqrt(2)) - d2) * DDB11O)/(2*p$p2)
      DD2 <- (((-0.5 * In4 *(u$u13 * x24x24))+
                 (0.5 * In4 *u$u14 *  x24x24)-
                 sqrt(2) *xlist$x13x24 *d4/2+
                 sqrt(2) * xlist$x14x24 *d4/2)*
                d2)/
        (sqrt(2)*p$p2)
      DD3 <- DDB11O * (d2^2)/(2*(p$p2^2))
      DDp2Db.base <- DD1 + DD2- DD3
      
      DD1.y2 <- t(DDB110 * (y==2)* drop(((2 * d2 *p$p2 * exp(-(-u$u11 + u$u13*p$p3 + u$u14*p$p4)/sqrt(2)) - d2) )/(2*p$p2))) %*% DDB110
      DD2.y2 <-  t(regr2$Z * (y==2)*drop((-0.5 *(d2/(sqrt(2)*p$p2)) * In4 *(u$u13)))) %*% regr2$Z+
        t(regr2$Z * (y==2)*drop((0.5 *(d2/(sqrt(2)*p$p2)) * In4 *(u$u14)))) %*% regr2$Z-
        (t(regr2$X3 * (y==2)* drop(sqrt(2) *(d2/(sqrt(2)*p$p2)) * d4/2)) %*% regr2$Z+
           t(t(regr2$X3 * (y==2)* drop(sqrt(2) *(d2/(sqrt(2)*p$p2)) * d4/2)) %*% regr2$Z))+
        (t(regr2$X4 * (y==2)* drop(sqrt(2) *(d2/(sqrt(2)*p$p2)) * d4/2)) %*% regr2$Z+
           t(t(regr2$X4 * (y==2)* drop(sqrt(2) *(d2/(sqrt(2)*p$p2)) * d4/2)) %*% regr2$Z))
      
      
      DD3.y2 <-  t(DDB110 * (y==2)*drop((d2^2)/ (2*(p$p2^2)))) %*% DDB110
      DDp2Db.y2 <- DD1.y2 + DD2.y2- DD3.y2
      
      DD1.y3 <- t(DDB110 * (y==3)* drop(((2 * d2 *p$p2 * exp(-(-u$u11 + u$u13*p$p3 + u$u14*p$p4)/sqrt(2)) - d2) )/(2*p$p2))) %*% DDB110
      DD2.y3 <-  t(regr2$Z * (y==3)*drop((-0.5 *(d2/(sqrt(2)*p$p2)) * In4 *(u$u13)))) %*% regr2$Z+
        t(regr2$Z * (y==3)*drop((0.5 *(d2/(sqrt(2)*p$p2)) * In4 *(u$u14)))) %*% regr2$Z-
        (t(regr2$X3 * (y==3)* drop(sqrt(2) *(d2/(sqrt(2)*p$p2)) * d4/2)) %*% regr2$Z+
           t(t(regr2$X3 * (y==3)* drop(sqrt(2) *(d2/(sqrt(2)*p$p2)) * d4/2)) %*% regr2$Z))+
        (t(regr2$X4 * (y==3)* drop(sqrt(2) *(d2/(sqrt(2)*p$p2)) * d4/2)) %*% regr2$Z+
           t(t(regr2$X4 * (y==3)* drop(sqrt(2) *(d2/(sqrt(2)*p$p2)) * d4/2)) %*% regr2$Z))
      
      
      DD3.y3 <-  t(DDB110 * (y==3)*drop((d2^2)/ (2*(p$p2^2)))) %*% DDB110
      DDp2Db.y3 <- DD1.y3 + DD2.y3- DD3.y3
      
      ##p3
      # DDp3Db.base <- -(In4 * x24x24)/ ( 2*p$p3) - (x24x24 * d4^2)/ (2*p$p3^2)
      DDp3Db <- t(regr2$Z * (y==2)*drop( -(In4)/( 2*p$p3)- (d4^2)/(2*p$p3^2)  )  ) %*% regr2$Z
      
      ###D2p1/Db is in 3 parts
      # DD1 <- ((2 * d2 *p$p2 * exp(-(-u$u11 + u$u13*p$p3 + u$u14*p$p4)/sqrt(2)) - d2) * (DDB11O))/(2*p$p1)
      # DD2 <- (((-0.5 * In4 *u$u13 * x24x24)+
      #            (0.5 * In4 *u$u14 * x24x24)-
      #            sqrt(2) *xlist$x13x24  *d4/2+
      #            sqrt(2) * xlist$x14x24  *d4/2)*
      #           d2)/
      #   (sqrt(2)*p$p1)
      # DD3 <- ((DDB11O) * d2^2)/(2*(p$p1^2))
      # DDp1Db.base <-  -DD1 - DD2- DD3
      
      DD1 <- t(DDB110 * (y==1)* drop(((2 * d2 *p$p2 * exp(-(-u$u11 + u$u13*p$p3 + u$u14*p$p4)/sqrt(2)) - d2) )/(2*p$p1))) %*% DDB110
      DD2 <-  t(regr2$Z * (y==1)*drop((-0.5 *(d2/(sqrt(2)*p$p1)) * In4 *(u$u13)))) %*% regr2$Z+
        t(regr2$Z * (y==1)*drop((0.5 *(d2/(sqrt(2)*p$p1)) * In4 *(u$u14)))) %*% regr2$Z-
        (t(regr2$X3 * (y==1)* drop(sqrt(2) *(d2/(sqrt(2)*p$p1)) * d4/2)) %*% regr2$Z+
           t(t(regr2$X3 * (y==1)* drop(sqrt(2) *(d2/(sqrt(2)*p$p1)) * d4/2)) %*% regr2$Z))+
        (t(regr2$X4 * (y==1)* drop(sqrt(2) *(d2/(sqrt(2)*p$p1)) * d4/2)) %*% regr2$Z+
           t(t(regr2$X4 * (y==1)* drop(sqrt(2) *(d2/(sqrt(2)*p$p1)) * d4/2)) %*% regr2$Z))
      DD3 <-  t(DDB110 * (y==1)*drop((d2^2)/ (2*(p$p1^2)))) %*% DDB110
      DDp1Db <- -DD1 - DD2 - DD3
      
      
      # DDB.base <- (y==3)*(DDp4Db.base+DDp2Db.base)+(y==2)*(DDp3Db.base + DDp2Db.base)+(y==1)*(DDp1Db.base)
      DDB <- (DDp4Db+DDp2Db.y3) + (DDp3Db+DDp2Db.y2) +DDp1Db
      
    }else{
      
      
      in2 <- (-u$u11 + u$u13*(p$p3) + u$u14*p$p4)
      denom2 <- sqrt(1 + p$p3^2 + p$p4^2)
      d4 <- dnorm(u$u24/sqrt(2))
      d2 <- dnorm((in2/denom2))
      roots <- (-sqrt(2)*d4*p$p3*regr2$Z + sqrt(2)*d4*p$p4*regr2$Z)
      bigRoots <-  - (in2*roots/(2*denom2^3)) +
        (-d4*u$u13*regr2$Z/sqrt(2)+d4*u$u14*regr2$Z/sqrt(2))/denom2
      DU4 <- -d4*u$u13*regr2$Z/sqrt(2) + d4*u$u14*regr2$Z/sqrt(2)
      
      
      # roots.sq <- t(apply(roots, 1, tcrossprod))
      # Broots.sq <- t(apply(bigRoots, 1, tcrossprod))
      # x24px24bx24 <- t(apply(regr2$Z, 1, function(x){x %*% t(x) %*%b%*%x}))
      # 
      # 
      # xlist <- list(roots.x1=cbind(regr2$X1, roots),
      #               roots.x3=cbind(regr2$X3, roots),
      #               roots.x4=cbind(regr2$X4, roots),
      #               roots.d4=cbind(DU4, roots),
      #               roots.sq=cbind(roots, roots),
      #               Broots.x1=cbind(bigRoots, regr2$X1),
      #               Broots.x3=cbind(bigRoots, regr2$X3),
      #               Broots.x4=cbind(bigRoots, regr2$X4))
      # 
      # xlist <- lapply(xlist,
      #                 apply,
      #                 1,
      #                 function(x){
      #                   tcrossprod(x[1:length(b)], x[(length(b)+1):(2*length(b))])+
      #                     tcrossprod(x[(length(b)+1):(2*length(b))], x[1:length(b)])
      #                 }
      # )
      # xlist<-lapply(xlist,t)
      # xlist$roots.d4 <- xlist$roots.d4/2
      # 
      
      #2nd deriv of p4
      # DDp4Db <- ((-sqrt(2)/2) * p$p4 *d4 * x24px24bx24 - FirthExtra$x24x24*(d4^2))/(2*p$p4^2)
      DDp4Db <- t(regr2$Z * (y==3)*drop(((-sqrt(2)/2) * p$p4 *d4 * (regr2$Z %*% b) - (d4^2))/(2*p$p4^2))) %*% regr2$Z
      
      ##2nd deriv of p3 (same as before)
      # DDp3Db <-  ((1/(sqrt(2))) * p$p3 * d4 * x24px24bx24 -FirthExtra$x24x24 * d4^2)/(2*p$p3^2)
      DDp3Db <-  t(regr2$Z * (y==2)*drop(((1/(sqrt(2))) * p$p3 * d4 * (regr2$Z %*% b) - d4^2)/(2*p$p3^2))) %*% regr2$Z
      
      # ### D/Db of (Dp2/Db) / p2 is 10 parts####
      # ddp2db11 <- -(d2^2 * FirthExtra$x11x11)/(p$p2^2*denom2^2) - (d2*in2*FirthExtra$x11x11)/(p$p2 * denom2^3)
      # ddp2db12 <- (d2^2 * p$p3*FirthExtra$x11x13)/(p$p2^2*denom2^2) +
      #   (d2*p$p3*in2*FirthExtra$x11x13)/(p$p2 * denom2^3)
      # ddp2db13 <- (d2^2 * p$p4*FirthExtra$x11x14)/(p$p2^2*denom2^2) +
      #   (d2*p$p4*in2*FirthExtra$x11x14)/(p$p2 * denom2^3)
      # ddp2db14 <- d2*xlist$roots.x1 /(2*p$p2*denom2^3)+
      #   d2^2 *xlist$Broots.x1/(p$p2^2*denom2) +
      #   d2*in2*xlist$Broots.x1/(p$p2*denom2^2)
      # ddp2db22 <- -(d2^2*p$p3^2*FirthExtra$x13x13)/(p$p2^2*denom2^2) - (d2*p$p3^2*in2*FirthExtra$x13x13)/(p$p2*denom2^3)
      # ddp2db23 <- -(d2^2*p$p3*p$p4*FirthExtra$x13x14)/(p$p2^2*denom2^2) -
      #   (d2*p$p3*p$p4*in2*FirthExtra$x13x14)/(p$p2*denom2^3)
      # ddp2db24 <- -(d2*d4*FirthExtra$x13x24/(sqrt(2)*p$p2*denom2)) -
      #   d2*p$p3*xlist$roots.x3/(2*p$p2*denom2^3) -
      #   d2^2*p$p3 * xlist$Broots.x3/(p$p2^2*denom2) -
      #   d2*p$p3*in2* xlist$Broots.x3/(p$p2 * denom2^2)
      # ddp2db33 <- -(d2^2*p$p4^2*FirthExtra$x14x14)/(p$p2^2 * (denom2^2)) -
      #   (d2*p$p4^2 * in2* FirthExtra$x14x14)/(p$p2*denom2^3)
      # ddp2db34 <- (d2 * d4* FirthExtra$x14x24)/(sqrt(2)*p$p2*denom2)-
      #   (d2*p$p4*xlist$roots.x4)/(2*p$p2*denom2^3)-
      #   (d2^2*p$p4*xlist$Broots.x4)/(p$p2^2*denom2)-
      #   (d2*p$p4*in2*xlist$Broots.x4)/(p$p2*denom2^2)
      # ddp2db44 <- -(d2^2 * Broots.sq)/p$p2^2 -
      #   d2*in2*Broots.sq/(p$p2*denom2) +
      #   d2/p$p2 * ( (3*in2 * roots.sq)/(4 *denom2^5) -
      #                 xlist$roots.d4/(denom2^3)-
      #                 in2*(2*d4^2*FirthExtra$x24x24 + d4*p$p3*x24px24bx24/sqrt(2)-d4*p$p4*x24px24bx24/sqrt(2))/(2*denom2^3)+
      #                 (d4*u$u13*x24px24bx24/(2*sqrt(2)) - d4*u$u14*x24px24bx24/(2*sqrt(2)))/(denom2)
      #   )
      # DDp2Db <- ddp2db11 +
      #   ddp2db12 +
      #   ddp2db13 +
      #   ddp2db14 +
      #   ddp2db22 +
      #   ddp2db23 +
      #   ddp2db24 +
      #   ddp2db33 +
      #   ddp2db34 +
      #   ddp2db44
      
      ddp2db11 <- t(regr2$X1 * (y==2)*drop( -(d2^2)/(p$p2^2*denom2^2) - (d2*in2)/(p$p2 * denom2^3)))  %*% regr2$X1
      ddp2db12 <- t(regr2$X1 * (y==2)*drop( (d2^2 * p$p3)/(p$p2^2*denom2^2) +(d2*p$p3*in2)/(p$p2 * denom2^3))) %*% regr2$X3+
        t(t(regr2$X1 * (y==2)*drop( (d2^2 * p$p3)/(p$p2^2*denom2^2) +(d2*p$p3*in2)/(p$p2 * denom2^3))) %*% regr2$X3)
      ddp2db13 <- t(regr2$X1 * (y==2)*drop( (d2^2 * p$p4)/(p$p2^2*denom2^2) + (d2*p$p4*in2)/(p$p2 * denom2^3))) %*% regr2$X4+
        t(t(regr2$X1 * (y==2)*drop( (d2^2 * p$p4)/(p$p2^2*denom2^2) + (d2*p$p4*in2)/(p$p2 * denom2^3))) %*% regr2$X4)
      ddp2db14 <- t(regr2$X1 * (y==2)*drop(d2/(2*p$p2*denom2^3))) %*% roots + t(t(regr2$X1 * (y==2)*drop(d2/(2*p$p2*denom2^3))) %*% roots)+
        t(regr2$X1 * (y==2)*drop(d2^2/(p$p2^2*denom2))) %*% bigRoots + t( t(regr2$X1 * (y==2)*drop(d2^2/(p$p2^2*denom2))) %*% bigRoots)+
        t(regr2$X1 * (y==2)*drop(d2*in2/(p$p2*denom2^2))) %*% bigRoots + t(t(regr2$X1 * (y==2)*drop(d2*in2/(p$p2*denom2^2))) %*% bigRoots )
      ddp2db22 <- t(regr2$X3 * (y==2)*drop( -(d2^2*p$p3^2)/(p$p2^2*denom2^2) - (d2*p$p3^2*in2)/(p$p2*denom2^3))) %*% regr2$X3
      ddp2db23 <- t(regr2$X3 * (y==2)*drop(-(d2^2*p$p3*p$p4)/(p$p2^2*denom2^2) - (d2*p$p3*p$p4*in2)/(p$p2*denom2^3))) %*% regr2$X4 +
        t( t(regr2$X3 * (y==2)*drop(-(d2^2*p$p3*p$p4)/(p$p2^2*denom2^2) - (d2*p$p3*p$p4*in2)/(p$p2*denom2^3))) %*% regr2$X4 )
      ddp2db24 <- t(regr2$X3 * (y==2)*drop(-(d2*d4/(sqrt(2)*p$p2*denom2)))) %*% regr2$Z + t(t(regr2$X3 * (y==2)*drop(-(d2*d4/(sqrt(2)*p$p2*denom2)))) %*% regr2$Z) -
        (t(regr2$X3 * (y==2)*drop( d2*p$p3/(2*p$p2*denom2^3))) %*% roots +t( t(regr2$X3 * (y==2)*drop( d2*p$p3/(2*p$p2*denom2^3))) %*% roots) )-
        (t(regr2$X3 * (y==2)*drop( d2^2*p$p3/(p$p2^2*denom2) )) %*% bigRoots + t(t(regr2$X3 * (y==2)*drop( d2^2*p$p3/(p$p2^2*denom2) )) %*% bigRoots)) -
        (t(regr2$X3 * (y==2)*drop(d2*p$p3*in2/(p$p2 * denom2^2))) %*% bigRoots + t(t(regr2$X3 * (y==2)*drop(d2*p$p3*in2/(p$p2 * denom2^2))) %*% bigRoots))
      ddp2db33 <- t(regr2$X4 * (y==2)*drop( -(d2^2*p$p4^2)/(p$p2^2 * (denom2^2)) -  (d2*p$p4^2 * in2)/(p$p2*denom2^3))) %*% regr2$X4
      ddp2db34 <- (t(regr2$X4* (y==2)*drop( (d2 * d4)/(sqrt(2)*p$p2*denom2))) %*% regr2$Z + t(t(regr2$X4* (y==2)*drop( (d2 * d4)/(sqrt(2)*p$p2*denom2))) %*% regr2$Z) )-
        (t(regr2$X4 * (y==2)*drop((d2*p$p4)/(2*p$p2*denom2^3))) %*% roots + t( t(regr2$X4 * (y==2)*drop((d2*p$p4)/(2*p$p2*denom2^3))) %*% roots))-
        (t(regr2$X4 * (y==2)*drop((d2^2*p$p4)/(p$p2^2*denom2))) %*% bigRoots + t(t(regr2$X4 * (y==2)*drop((d2^2*p$p4)/(p$p2^2*denom2))) %*% bigRoots))-
        (t(regr2$X4 * (y==2)*drop((d2*p$p4*in2)/(p$p2*denom2^2))) %*% bigRoots + t(t(regr2$X4 * (y==2)*drop((d2*p$p4*in2)/(p$p2*denom2^2))) %*% bigRoots))
      ddp2db44 <- t(bigRoots * (y==2)*drop(-(d2^2)/p$p2^2 -d2*in2/(p$p2*denom2))) %*% bigRoots +
        ( t(roots * (y==2)*drop((3*in2*d2/p$p2)/(4 *denom2^5))) %*% roots -
            t(roots * (y==2)*drop( (d2/p$p2 )/(denom2^3))) %*% DU4 -
            t(regr2$Z * (y==2)*drop(d2/p$p2 *in2*(2*d4^2 + d4*p$p3*(regr2$Z%*%b)/sqrt(2)-d4*p$p4*(regr2$Z%*%b)/sqrt(2))/(2*denom2^3))) %*% regr2$Z +
            t(regr2$Z * (y==2)*drop(d2/p$p2 *(d4*u$u13*(regr2$Z %*% b)/(2*sqrt(2)) - d4*u$u14*(regr2$Z %*% b)/(2*sqrt(2)))/(denom2) )) %*% regr2$Z
        )
      DDp2Db.y2 <- ddp2db11 +
        ddp2db12 +
        ddp2db13 +
        ddp2db14 +
        ddp2db22 +
        ddp2db23 +
        ddp2db24 +
        ddp2db33 +
        ddp2db34 +
        ddp2db44
      
      
      
      ddp2db11 <- t(regr2$X1 * (y==3)*drop( -(d2^2)/(p$p2^2*denom2^2) - (d2*in2)/(p$p2 * denom2^3)))  %*% regr2$X1
      ddp2db12 <- t(regr2$X1 * (y==3)*drop( (d2^2 * p$p3)/(p$p2^2*denom2^2) +(d2*p$p3*in2)/(p$p2 * denom2^3))) %*% regr2$X3+
        t(t(regr2$X1 * (y==3)*drop( (d2^2 * p$p3)/(p$p2^2*denom2^2) +(d2*p$p3*in2)/(p$p2 * denom2^3))) %*% regr2$X3)
      ddp2db13 <- t(regr2$X1 * (y==3)*drop( (d2^2 * p$p4)/(p$p2^2*denom2^2) + (d2*p$p4*in2)/(p$p2 * denom2^3))) %*% regr2$X4+
        t(t(regr2$X1 * (y==3)*drop( (d2^2 * p$p4)/(p$p2^2*denom2^2) + (d2*p$p4*in2)/(p$p2 * denom2^3))) %*% regr2$X4)
      ddp2db14 <- t(regr2$X1 * (y==3)*drop(d2/(2*p$p2*denom2^3))) %*% roots + t(t(regr2$X1 * (y==3)*drop(d2/(2*p$p2*denom2^3))) %*% roots)+
        t(regr2$X1 * (y==3)*drop(d2^2/(p$p2^2*denom2))) %*% bigRoots + t( t(regr2$X1 * (y==3)*drop(d2^2/(p$p2^2*denom2))) %*% bigRoots)+
        t(regr2$X1 * (y==3)*drop(d2*in2/(p$p2*denom2^2))) %*% bigRoots + t(t(regr2$X1 * (y==3)*drop(d2*in2/(p$p2*denom2^2))) %*% bigRoots )
      ddp2db22 <- t(regr2$X3 * (y==3)*drop( -(d2^2*p$p3^2)/(p$p2^2*denom2^2) - (d2*p$p3^2*in2)/(p$p2*denom2^3))) %*% regr2$X3
      ddp2db23 <- t(regr2$X3 * (y==3)*drop(-(d2^2*p$p3*p$p4)/(p$p2^2*denom2^2) - (d2*p$p3*p$p4*in2)/(p$p2*denom2^3))) %*% regr2$X4 +
        t( t(regr2$X3 * (y==3)*drop(-(d2^2*p$p3*p$p4)/(p$p2^2*denom2^2) - (d2*p$p3*p$p4*in2)/(p$p2*denom2^3))) %*% regr2$X4 )
      ddp2db24 <- t(regr2$X3 * (y==3)*drop(-(d2*d4/(sqrt(2)*p$p2*denom2)))) %*% regr2$Z + t(t(regr2$X3 * (y==3)*drop(-(d2*d4/(sqrt(2)*p$p2*denom2)))) %*% regr2$Z) -
        (t(regr2$X3 * (y==3)*drop( d2*p$p3/(2*p$p2*denom2^3))) %*% roots +t( t(regr2$X3 * (y==3)*drop( d2*p$p3/(2*p$p2*denom2^3))) %*% roots) )-
        (t(regr2$X3 * (y==3)*drop( d2^2*p$p3/(p$p2^2*denom2) )) %*% bigRoots + t(t(regr2$X3 * (y==3)*drop( d2^2*p$p3/(p$p2^2*denom2) )) %*% bigRoots)) -
        (t(regr2$X3 * (y==3)*drop(d2*p$p3*in2/(p$p2 * denom2^2))) %*% bigRoots + t(t(regr2$X3 * (y==3)*drop(d2*p$p3*in2/(p$p2 * denom2^2))) %*% bigRoots))
      ddp2db33 <- t(regr2$X4 * (y==3)*drop( -(d2^2*p$p4^2)/(p$p2^2 * (denom2^2)) -  (d2*p$p4^2 * in2)/(p$p2*denom2^3))) %*% regr2$X4
      ddp2db34 <- (t(regr2$X4* (y==3)*drop( (d2 * d4)/(sqrt(2)*p$p2*denom2))) %*% regr2$Z + t(t(regr2$X4* (y==3)*drop( (d2 * d4)/(sqrt(2)*p$p2*denom2))) %*% regr2$Z) )-
        (t(regr2$X4 * (y==3)*drop((d2*p$p4)/(2*p$p2*denom2^3))) %*% roots + t( t(regr2$X4 * (y==3)*drop((d2*p$p4)/(2*p$p2*denom2^3))) %*% roots))-
        (t(regr2$X4 * (y==3)*drop((d2^2*p$p4)/(p$p2^2*denom2))) %*% bigRoots + t(t(regr2$X4 * (y==3)*drop((d2^2*p$p4)/(p$p2^2*denom2))) %*% bigRoots))-
        (t(regr2$X4 * (y==3)*drop((d2*p$p4*in2)/(p$p2*denom2^2))) %*% bigRoots + t(t(regr2$X4 * (y==3)*drop((d2*p$p4*in2)/(p$p2*denom2^2))) %*% bigRoots))
      ddp2db44 <- t(bigRoots * (y==3)*drop(-(d2^2)/p$p2^2 -d2*in2/(p$p2*denom2))) %*% bigRoots +
        ( t(roots * (y==3)*drop((3*in2*d2/p$p2)/(4 *denom2^5))) %*% roots -
            t(roots * (y==3)*drop( (d2/p$p2 )/(denom2^3))) %*% DU4 -
            t(regr2$Z * (y==3)*drop(d2/p$p2 *in2*(2*d4^2 + d4*p$p3*(regr2$Z%*%b)/sqrt(2)-d4*p$p4*(regr2$Z%*%b)/sqrt(2))/(2*denom2^3))) %*% regr2$Z +
            t(regr2$Z * (y==3)*drop(d2/p$p2 *(d4*u$u13*(regr2$Z %*% b)/(2*sqrt(2)) - d4*u$u14*(regr2$Z %*% b)/(2*sqrt(2)))/(denom2) )) %*% regr2$Z
        )
      DDp2Db.y3 <- ddp2db11 +
        ddp2db12 +
        ddp2db13 +
        ddp2db14 +
        ddp2db22 +
        ddp2db23 +
        ddp2db24 +
        ddp2db33 +
        ddp2db34 +
        ddp2db44
      
      
      
      
      
      
      
      
      ### D/Db of (Dp1/Db)/ p1 is also 10 parts####
      # ddp1db11 <- -(d2^2 * FirthExtra$x11x11)/(p$p1^2*denom2^2) + (d2*in2*FirthExtra$x11x11)/(p$p1 * denom2^3)
      # ddp1db12 <- (d2^2 * p$p3*FirthExtra$x11x13)/((p$p1^2) * (denom2^2)) -
      #   (d2*p$p3*in2*FirthExtra$x11x13)/(p$p1 * denom2^3) #poss
      # ddp1db13 <- (d2^2 * p$p4*FirthExtra$x11x14)/(p$p1^2*denom2^2) -
      #   (d2*p$p4*in2*FirthExtra$x11x14)/(p$p1 * denom2^3)
      # ddp1db14 <- -d2*xlist$roots.x1 /(2*p$p1*denom2^3) +
      #   d2^2 *xlist$Broots.x1/(p$p1^2*denom2) -
      #   d2*in2*xlist$Broots.x1/(p$p1*denom2^2)
      # ddp1db22 <- -(d2^2*p$p3^2*FirthExtra$x13x13)/(p$p1^2*denom2^2) +
      #   (d2*p$p3^2*in2*FirthExtra$x13x13)/(p$p1*denom2^3) #inter
      # ddp1db23 <- -(d2^2*p$p3*p$p4*FirthExtra$x13x14)/(p$p1^2*denom2^2) +
      #   (d2*p$p3*p$p4*in2*FirthExtra$x13x14)/(p$p1*denom2^3)
      # ddp1db24 <- (d2*d4*FirthExtra$x13x24/(sqrt(2)*p$p1*denom2)) +
      #   d2*p$p3*xlist$roots.x3/(2*p$p1*denom2^3) -
      #   d2^2*p$p3 * xlist$Broots.x3/(p$p1^2*denom2) +
      #   d2*p$p3*in2* xlist$Broots.x3/(p$p1 * denom2^2)
      # ddp1db33 <- -(d2^2*p$p4^2*FirthExtra$x14x14)/(p$p1^2 * (denom2^2)) +
      #   (d2*p$p4^2 * in2* FirthExtra$x14x14)/(p$p1*denom2^3)
      # ddp1db34 <- -(d2 * d4* FirthExtra$x14x24)/(sqrt(2)*p$p1*denom2)+
      #   (d2*p$p4*xlist$roots.x4)/(2*p$p1*denom2^3)-
      #   (d2^2*p$p4*xlist$Broots.x4)/(p$p1^2*denom2)+
      #   (d2*p$p4*in2*xlist$Broots.x4)/(p$p1*denom2^2)
      # ddp1db44 <- -(d2^2 * Broots.sq)/p$p1^2 +
      #   d2*in2*Broots.sq/(p$p1*denom2) -
      #   d2/p$p1 * ( (3*in2 * roots.sq)/(4 *denom2^5) -
      #                 xlist$roots.d4/(denom2^3)-
      #                 in2*(2*d4^2*FirthExtra$x24x24 + d4*p$p3*x24px24bx24/sqrt(2)-d4*p$p4*x24px24bx24/sqrt(2))/(2*denom2^3)+
      #                 (d4*u$u13*x24px24bx24/(2*sqrt(2)) - d4*u$u14*x24px24bx24/(2*sqrt(2)))/(denom2)
      #   )
      # DDp1Db <- ddp1db11 +
      #   ddp1db12 +
      #   ddp1db13 +
      #   ddp1db14 +
      #   ddp1db22 +
      #   ddp1db23 +
      #   ddp1db24 +
      #   ddp1db33 +
      #   ddp1db34 +
      #   ddp1db44
      ### D/Db of (Dp1/Db)/ p1 is also 10 parts####
      ddp1db11A <- t(regr2$X1 *(y==1)*drop(-(d2^2)/(p$p1^2*denom2^2) + (d2*in2)/(p$p1 * denom2^3))) %*% regr2$X1
      ddp1db12A <- t(regr2$X1 *(y==1)*drop((d2^2 * p$p3)/((p$p1^2) * (denom2^2)) - (d2*p$p3*in2)/(p$p1 * denom2^3)))  %*% regr2$X3 +
        t(t(regr2$X1 *(y==1)*drop((d2^2 * p$p3)/((p$p1^2) * (denom2^2)) - (d2*p$p3*in2)/(p$p1 * denom2^3)))  %*% regr2$X3)
      ddp1db13A <- t(regr2$X1 *(y==1)*drop((d2^2 * p$p4)/(p$p1^2*denom2^2) - (d2*p$p4*in2)/(p$p1 * denom2^3))) %*% regr2$X4+
        t(t(regr2$X1 *(y==1)*drop((d2^2 * p$p4)/(p$p1^2*denom2^2) - (d2*p$p4*in2)/(p$p1 * denom2^3))) %*% regr2$X4)
      ddp1db14A <- (t(regr2$X1 *(y==1)*drop(-d2/(2*p$p1*denom2^3))) %*% roots + t(t(regr2$X1 *(y==1)*drop(-d2/(2*p$p1*denom2^3))) %*% roots)) +
        (t(regr2$X1 *(y==1)*drop(d2^2/(p$p1^2*denom2) -d2*in2/(p$p1*denom2^2))) %*% bigRoots+
           t(t(regr2$X1 *(y==1)*drop(d2^2 /(p$p1^2*denom2) -d2*in2/(p$p1*denom2^2))) %*% bigRoots))
      ddp1db22A <- t(regr2$X3 *(y==1)*drop(-(d2^2*p$p3^2)/(p$p1^2*denom2^2) +(d2*p$p3^2*in2)/(p$p1*denom2^3))) %*% regr2$X3
      ddp1db23A <- (t(regr2$X3 *(y==1)*drop(-(d2^2*p$p3*p$p4)/(p$p1^2*denom2^2) + (d2*p$p3*p$p4*in2)/(p$p1*denom2^3))) %*% regr2$X4 +
                      t(t(regr2$X3 * (y==1)*drop(-(d2^2*p$p3*p$p4)/(p$p1^2*denom2^2) + (d2*p$p3*p$p4*in2)/(p$p1*denom2^3))) %*% regr2$X4))
      ddp1db24A <- (t(regr2$X3 * (y==1)*drop((d2*d4/(sqrt(2)*p$p1*denom2)))) %*% regr2$Z + t(t(regr2$X3 * (y==1)*drop((d2*d4/(sqrt(2)*p$p1*denom2)))) %*% regr2$Z )) +
        (t(regr2$X3 * (y==1)*drop(d2*p$p3/(2*p$p1*denom2^3))) %*% roots + t(t(regr2$X3 * (y==1)*drop(d2*p$p3/(2*p$p1*denom2^3))) %*% roots)) -
        (t(regr2$X3 * (y==1)*drop(d2^2*p$p3/(p$p1^2*denom2))) %*% bigRoots + t(t(regr2$X3 * (y==1)*drop(d2^2*p$p3/(p$p1^2*denom2))) %*% bigRoots)) +
        (t(regr2$X3 * (y==1)*drop(d2*p$p3*in2/(p$p1 * denom2^2)))%*% bigRoots + t(t(regr2$X3 * (y==1)*drop(d2*p$p3*in2/(p$p1 * denom2^2)))%*% bigRoots ))
      ddp1db33A <- t(regr2$X4 * (y==1)*drop(-(d2^2*p$p4^2)/(p$p1^2 * (denom2^2)) + (d2*p$p4^2 * in2)/(p$p1*denom2^3))) %*% regr2$X4
      ddp1db34A <- (t(regr2$X4* (y==1)*drop(-(d2 * d4)/(sqrt(2)*p$p1*denom2))) %*% regr2$Z+ t(t(regr2$X4* (y==1)*drop(-(d2 * d4)/(sqrt(2)*p$p1*denom2))) %*% regr2$Z))+
        (t(regr2$X4 * (y==1)*drop((d2*p$p4)/(2*p$p1*denom2^3))) %*% roots + t(t(regr2$X4 * (y==1)*drop((d2*p$p4)/(2*p$p1*denom2^3))) %*% roots))-
        (t(regr2$X4 * (y==1)*drop((d2^2*p$p4)/(p$p1^2*denom2))) %*% bigRoots + t(t(regr2$X4 * (y==1)*drop((d2^2*p$p4)/(p$p1^2*denom2))) %*% bigRoots))+
        (t(regr2$X4 * (y==1)*drop((d2*p$p4*in2)/(p$p1*denom2^2))) %*% bigRoots + t(t(regr2$X4 * (y==1)*drop((d2*p$p4*in2)/(p$p1*denom2^2))) %*% bigRoots))
      ddp1db44A <- t(bigRoots * (y==1)*drop(-(d2^2 )/p$p1^2)) %*% bigRoots+
        t(bigRoots* (y==1)*drop(d2*in2/(p$p1*denom2))) %*% bigRoots  +
        t(roots *(y==1)*drop( (-d2/p$p1 * 3*in2 )/(4 *denom2^5))) %*% roots +
        t(roots *(y==1)*drop( (d2/p$p1 )/(denom2^3))) %*% DU4 +
        t(regr2$Z * (y==1)*drop( d2/p$p1 * in2*(2*d4^2 + d4*p$p3*(regr2$Z%*%b)/sqrt(2)-d4*p$p4*(regr2$Z%*%b)/sqrt(2))/(2*denom2^3))) %*% regr2$Z+
        t(regr2$Z * (y==1)*drop(- d2/p$p1 * (d4*u$u13*(regr2$Z%*%b)/(2*sqrt(2)) - d4*u$u14*(regr2$Z%*%b)/(2*sqrt(2)))/(denom2))) %*% regr2$Z
      
      DDp1DbA <- ddp1db11A +
        ddp1db12A +
        ddp1db13A +
        ddp1db14A +
        ddp1db22A +
        ddp1db23A +
        ddp1db24A +
        ddp1db33A +
        ddp1db34A +
        ddp1db44A
      
      # DDB <- (y==3)*(DDp4Db+DDp2Db)+
      #   (y==2)*(DDp3Db+DDp2Db)+
      #   (y==1)*(DDp1Db)
      
      DDB <-(DDp4Db+DDp2Db.y3)+
        (DDp3Db+DDp2Db.y2)+
        (DDp1DbA)
      
    }
  }
  ans <- DDB
  # ans <- colSums(DDB)
  # ans <- matrix(ans, nrow=length(b))
  return(ans)
}

