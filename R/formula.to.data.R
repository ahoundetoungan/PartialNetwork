#' @importFrom Formula as.Formula
#' @importFrom formula.tools env
#' @importFrom stats model.frame
#' @importFrom stats terms
#' @importFrom stats update
#' @importFrom stats model.response
#' @importFrom stats model.matrix
#' @importFrom stats delete.response
formula.to.data <- function(formula, contextual, data, type = "outcome") {
  if (missing(data)) {
    data           <- env(formula)
  }
  if(is.null(data)) {
    data           <- env(formula)
  }
  
  ## Extract data from the formula
  formula          <- as.Formula(formula)
  
  if (type == "outcome") {
    stopifnot(length(formula)[1] == 1L, length(formula)[2] %in% 1:2)
  } else {
    if (!all(length(formula) == c(0, 1))) {
      stop("mlinks.formula is not defined as mlinks.formula = ~ X1 + X2 + ...")
    }
  }
  
  
  # try to handle dots in formula
  has_dot          <- function(formula) inherits(try(terms(formula), silent = TRUE), "try-error")
  if(has_dot(formula)) {
    f1             <- formula(formula, rhs = 1)
    f2             <- formula(formula, lhs = 0, rhs = 2)
    if(!has_dot(f1) & has_dot(f2)) {
      formula      <- Formula::as.Formula(f1, update(formula(formula, lhs = 0, rhs = 1), f2))
    }
  }
  
  ## call model.frame()
  mf               <- model.frame(formula, data = data)
  ## extract response, terms, model matrices
  y                <- model.response(mf, "numeric")
  mtXone           <- terms(formula, data = data, rhs = 1)
  Xone             <- model.matrix(mtXone, mf)
  Colnames.xone    <- colnames(Xone)
  intercept.one    <- "(Intercept)" %in% Colnames.xone
  
  mtX              <- NULL
  X                <- NULL
  if(length(formula)[2] >= 2L) {
    mtX            <- delete.response(terms(formula, data = data, rhs = 2))
    X              <- model.matrix(mtX, mf)
  }
  
  if (contextual) {
    if (!is.null(X)) {
      stop("contextual cannot be TRUE while contextual variable are declared after the pipe.")
    }
    X              <- Xone
    tmpx           <- as.character.default(formula(formula, lhs = 1, rhs = 1))
    formula        <- Formula::as.Formula(paste(c(tmpx[c(2, 1, 3)], "|", tmpx[3]), collapse = " "))
  } 
  
  Colnames.x       <- colnames(X)
  intercept        <- "(Intercept)" %in% Colnames.x
  
  if (intercept) {
    X              <- X[,-1, drop = FALSE]
  }  
  
  
  list("formula" = formula, 
       "Xone"    = Xone,
       "X"       = X, 
       "y"       = y)
}



formula.to.data.smm <- function(formula, data, fixed.effects) {
  if (missing(data)) {
    data           <- env(formula)
  }
  if(is.null(data)) {
    data           <- env(formula)
  }
  
  ftmp              <- gsub(" ", "", as.character(Reduce(paste, deparse(formula))), fixed = TRUE)
  pospipe           <- gregexpr("\\|", ftmp)[[1]]
  Gyobs             <- FALSE
  GXobs             <- FALSE
  if(length(pospipe) == 1){
    if(pospipe > 0) {
      Gyobs         <- TRUE
    }
  }
  if(length(pospipe) == 2){
    if((pospipe[2] - pospipe[1]) == 1) {
      GXobs         <- TRUE
      ftmp          <- paste0(substr(ftmp, 1, pospipe[1]), substring(ftmp, pospipe[2] + 1))
    } else {
      Gyobs         <- TRUE
      GXobs         <- TRUE
    }
  }
  formula          <- as.Formula(ftmp)
  stopifnot(length(formula)[1] == 1L, length(formula)[2] %in% 1:3)
  
  ## call model.frame()
  mf               <- model.frame(formula, data = data)
  ## extract response, terms, model matrices
  y                <- model.response(mf, "numeric")
  mtX              <- terms(formula, data = data, rhs = 1)
  X                <- model.matrix(mtX, mf)
  if(fixed.effects & ("(Intercept)" %in% colnames(X))){
    X              <- X[, -1, drop = FALSE]
  }
  GX               <- NULL
  Gy               <- NULL
  if(Gyobs){
    mtGy           <- delete.response(terms(formula, data = data, rhs = 2))
    Gy             <- model.matrix(mtGy, mf)
    if ("(Intercept)" %in% colnames(Gy)) {
      Gy           <- Gy[,-1]
    }
    if(!is.null(dim(Gy))) stop("Several variables declared as Gy")
    if(GXobs){
      mtGX         <- delete.response(terms(formula, data = data, rhs = 3))
      GX           <- model.matrix(mtGX, mf)
      if ("(Intercept)" %in% colnames(GX)) {
        GX         <- GX[,-1,drop = FALSE]
      }
    }
  } else{
    if(GXobs){
      mtGX         <- delete.response(terms(formula, data = data, rhs = 2))
      GX           <- model.matrix(mtGX, mf)
      if ("(Intercept)" %in% colnames(GX)) {
        GX         <- GX[,-1,drop = FALSE]
      }
      formula      <- as.Formula(paste0(substr(ftmp, 1, pospipe[1]), "|", substring(ftmp, pospipe[1] + 1)))
    }
  }
  
  list("formula" = formula, 
       "X"       = X,
       "y"       = y,
       "GX"      = GX,
       "Gy"      = Gy)
}


formula.to.data.bias <- function(formula, data) {
  formula          <- as.Formula(formula)

  ## call model.frame()
  mf               <- model.frame(formula, data = data)
  ## extract response, terms, model matrices
  y                <- model.response(mf, "numeric")
  mtX              <- terms(formula, data = data, rhs = 1)
  X                <- model.matrix(mtX, mf)
  cnames           <- colnames(X)
  
  mtZ              <- delete.response(terms(formula, data = data, rhs = 2))
  Z                <- model.matrix(mtZ, mf)

  list("S" = X, "Z" = Z)
}