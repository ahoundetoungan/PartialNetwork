fcoefficients          <- function(coef, std) {
  cnames               <- names(coef)
  tval                 <- coef/std
  pval                 <- 2*(1 - pnorm(abs(tval)))
  
  pval_print           <- unlist(lapply(pval, function(u){
    ifelse(u < 2e-16, "<2e-16", format(u, digit = 3))
  }))
  
  refprob              <- c(0.001, 0.01, 0.05, 0.1)
  refstr               <- c("***",  "**", "*", ".", "")
  str                  <- unlist(lapply(pval, function(u) refstr[1 + sum(u > refprob)]))
  
  
  out_print            <- data.frame("X1" = round(coef, 6),
                                     "X2" = round(std, 6),
                                     "X3" = round(tval, 2),
                                     "X4" = pval_print,
                                     "X5" = str)
  
  out                  <- data.frame("X1" = coef,
                                     "X2" = std,
                                     "X3" = tval,
                                     "X4" = pval)
  
  rownames(out_print)  <- cnames
  colnames(out_print)  <- c("Estimate", "Robust SE", "t value", "Pr(>|t|)", "")
  rownames(out)        <- cnames
  colnames(out)        <- c("Estimate", "Robust SE", "t value", "Pr(>|t|)")
  
  list(out_print = out_print, out = out)
}
