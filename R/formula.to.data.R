#' @importFrom Formula as.Formula
#' @importFrom stats model.frame
formula.to.data <- function(formula, contextual, data) {
  
  ## Extract data from the formula
  cl               <- match.call()
  formula          <- as.Formula(formula)
  if (missing(data)) {
    data           <- environment(formula)
  }
  
  stopifnot(length(formula)[1] == 1L, length(formula)[2] %in% 1:2)
  
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
  mf               <- model.frame(formula)
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
    X              <- Xone
  } 
  
  Colnames.x       <- colnames(X)
  intercept        <- "(Intercept)" %in% Colnames.x
  
  if (intercept) {
    X              <- X[,-1]
  } 
  
  list("formula" = formula, 
       "Xone" = Xone,
       "X" = X, 
       "y" = y)
}