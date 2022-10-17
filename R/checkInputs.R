# Purpose: Check the inputs of the functions (defensive programming)

.checkInput_p <- function(C, w, G, Sens) {

  if (missing(C)) stop("missing value for C is not allowed")
  if (missing(w)) stop("missing value for w is not allowed")
  if (missing(G)) stop("missing value for G is not allowed")
  if (missing(Sens)) stop("missing value for Sens is not allowed")
  if (!is.numeric(c(C, w, G, Sens))) stop("'C', 'w', 'G', & 'Sens' must be numeric")
  if (Sens >= 1 || Sens <= 0 ) stop("'Sens' must be between 0 and 1")
  if (C <= 0)  stop("'C' must be a positive number")
  if (w <= 0)  stop("'w' must be a positive number")
  if (G <= 0)  stop("'G' must be a positive number")
}

.checkInput_r <- function(r) {
  if (missing(r) || r<=0)  stop("Specify a positive value for r")
}

.checkInput_f <- function(f) {
  if (missing(f) || f > 1 || f <= 0) stop("Specify a value for f between 0 and 1")
}

.checkInput_prob_det<- function(prob_det){
  if (prob_det >= 1 || prob_det <= 0 ) stop("'prob_det' must be between 0 and 1")
}

.checkInput_dollar <- function(samp_dollar,lot_dollar) {
  if (!is.numeric(c(samp_dollar,lot_dollar))) stop("'samp_dollar', & 'lot_dollar' must be numeric")
  if (missing(samp_dollar) || samp_dollar < 0) stop("'samp_dollar' must be greater than or equal to 0")
  if (missing(lot_dollar) || (lot_dollar < 0)) stop("'lot_dollar' must be greater than or equal to 0")
}

