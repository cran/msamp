#' Calculate the probability of a single sample being contaminated
#'
#' The p() function calculates the probability of a single sample unit, with weight
#' w, and postulated contamination, C, being contaminated above a target level, G.
#'
#' @param C suspected lot contamination (CFU/g)
#' @param w weight of single sample unit (g)
#' @param G target value to detect (CFU/g)
#' @param Sens sensitivity of the analytical test (%)
#' @param D distribution of the bacteria in the product lot: "homogeneous", "heterogeneous", or "localized"
#' @param r for the heterogeneous case only, the degree of heterogeneity. r > 0
#' @param f for the localized case, r is further specified. 0 < f < 1
#'
#' @return A numeric value: the probability of a single sample unit being contaminated
#' above target level.
#'
#' @details Refer to vignette for details.
#'
#' @examples
#'
#' #A sample of 25 grams (w=25) is collected and analyzed using an analytical
#' #test with sensitivity of 90% (Sens=.9), to detect at least 5 CFU's/g (G=5).
#' #The suspected or postulated level of contamination in the lot is 4 CFU's/g (C=4)
#'
#'#homogeneous case
#'p(C=4,w=25,G=5,Sens=.9,D="homogeneous",r=NULL,f=NULL)
#'# 0.006117884
#'#heterogeneous case-- dispersion, r, is postulated as 2
#'p(C=4,w=25,G=5,Sens=.9,D="heterogeneous",r=2,f=NULL)
#'# 0.2576463
#'#localized case -- 30% of the lot is postulated to be contaminated
#'p(C=4,w=25,G=5,Sens=.9,D="localized",r=NULL,f=.3)
#'# 0.001835365
#'
#' @importFrom stats ppois
#' @importFrom stats pnbinom
#'
#' @export

p <- function(C, w, G, Sens, D=c("homogeneous","heterogeneous","localized"), r = NULL, f = NULL) {

match.arg(D)

.checkInput_p(C, w, G, Sens)

mu <- C * w
g  <- round(G * w,0)

if (D == 'homogeneous') {
    p1 <- Sens*stats::ppois(g, mu, lower.tail = FALSE, log.p = FALSE)
    return(p1)
} else if (D == 'heterogeneous') {
      .checkInput_r(r)
       p2 <- Sens*stats::pnbinom(g, mu=mu, size=r,lower.tail = FALSE, log.p = FALSE)
      return(p2)
} else if (D == 'localized') {
    .checkInput_f(f)
      p3 <- Sens*f*stats::ppois(g, mu, lower.tail = FALSE, log.p = FALSE)
    return(p3)
    }
}

#' Calculate the sample size necessary to detect contamination above
#' target level
#'
#' The n() function calculates the sample size,n, necessary to detect contamination above
#' a target level, G, in a product lot, where the probability of a single sample unit being above
#' the target level is calculated by the msamp function p(). The total cost, cost_tot,
#' associated with sample size is also output.
#'
#' @param C suspected lot contamination (CFU/g)
#' @param w weight of single sample unit (g)
#' @param G target value to detect (CFU/g)
#' @param Sens sensitivity of the analytical test (%)
#' @param D distribution of the bacteria in the product lot: "homogeneous", "heterogeneous", or "localized"
#' @param r for the heterogeneous case only, the degree of heterogeneity. r > 0
#' @param f for the localized case, r is further specified. 0 < f < 1
#' @param prob_det desired probability of detecting bacterial contamination above the
#' target level in the product lot. Set to 0.9 by default
#' @param samp_dollar cost per sample unit in $
#' @param lot_dollar  fixed cost (if any) of sampling the lot in $
#'
#' @return A list containing:
#'   \itemize{
#'     \item{\strong{n}: }{The sample size}
#'     \item{\strong{p}: }{The probability of a single sample unit being contaminated above the
#'                          target level}
#'     \item{\strong{cost_tot}: }{The total cost associated with sampling of
#'        \emph{n} samples}
#'   }
#'
#' @details Refer to vignette for details.
#'
#' @examples
#'
#' #A sample of 25 grams (w=25) is collected and analyzed using an analytical
#' #test with sensitivity of 90% (Sens=.9), to detect at least 5 CFU's/g (G=5).
#' #The suspected or postulated level of contamination in the lot is 4 CFU's/g (C=4).
#' #The desired probability of picking at least one sample unit contaminated above the target
#' #level is 0.9 (prob_det=0.9), the cost of a single sampling unit is $100 (samp_dollar=100),
#' #and the fixed cost for sampling the entire lot is $200 (lot_dollar=200).
#'
#'#homogeneous case
#'n(C=4,w=25,G=5,Sens=.9,D="homogeneous",r=NULL,f=NULL,prob_det=0.9,samp_dollar=100,lot_dollar=200)
#'# n=376, total cost=$37,722
#'#heterogeneous case
#'n(C=4,w=25,G=5,Sens=.9,D="heterogeneous",r=10,f=NULL,prob_det=0.9,samp_dollar=100,lot_dollar=200)
#'# n=12, total cost=$1,319
#'#localized case
#'n(C=4,w=25,G=5,Sens=.9,D="localized",r=NULL,f=.3,prob_det=0.9,samp_dollar=100,lot_dollar=200)
#'# n=1,254 , total cost=$125,541
#'
#' @importFrom stats ppois
#' @importFrom stats pnbinom
#'
#' @export


n <- function(C,w,G,Sens,D=c("homogeneous","heterogeneous","localized"),r=NULL,f=NULL,prob_det=0.9,samp_dollar,lot_dollar) {

  match.arg(D)

  .checkInput_p(C, w, G, Sens)
  .checkInput_prob_det(prob_det)
  .checkInput_dollar(samp_dollar,lot_dollar)
  if (D == 'heterogeneous') { .checkInput_r(r)}
  if (D == 'localized') { .checkInput_f(f)}

 p<- p(C,w,G,Sens,D,r,f)

 n<- ceiling(log(1 - prob_det) / log(1 - p))

 if( n > 100) warning('Large sample size (n > 100). Decreasing the weight, w, of
 the sample unit (g) or relaxing the probability of detection, prob_det results
 in smaller sample sizes')


 cost_tot<-  round(lot_dollar + (log(1 - prob_det) / log(1 - p)) * samp_dollar,0)

return(list(n=n,p=p,cost_tot=cost_tot))
}

#' Plots the relation between the probability of detection and the sample size, n
#'
#' The plotn() function examines the effect of increasing the probability of detection
#' on the sample size, n, where the probability of a single sample unit being contaminated
#' above the target limit is calculated from the msamp function p()
#'
#' @param C suspected lot contamination (CFU/g)
#' @param w weight of single sample unit (g)
#' @param G target value to detect (CFU/g)
#' @param Sens sensitivity of the analytical test (%)
#' @param D distribution of the bacteria in the product lot: "homogeneous", "heterogeneous", or "localized"
#' @param r for the heterogeneous case only, the degree of heterogeneity. r > 0
#' @param f for the localized case, r is further specified. 0 < f < 1
#'
#' @return A plot, of recordedplot class. The probability of detection is on the y-axis and
#' the  sample size n is on the x-axis. Overlaid at intersecting red dashed lines is the
#' sample size for probability of detection (prob_det) = 0.9.
#'
#' @details Refer to vignette for details.
#'
#' @importFrom stats ppois
#' @importFrom stats pnbinom
#' @importFrom graphics plot
#' @importFrom graphics mtext
#' @importFrom grDevices recordPlot
#' @importFrom graphics abline
#'
#' @export

plotn <- function(C,w,G,Sens,D=c("homogeneous","heterogeneous","localized"),r=NULL,f=NULL) {

  match.arg(D)

  .checkInput_p(C, w, G, Sens)
  if (D == 'heterogeneous') { .checkInput_r(r)}
  if (D == 'localized') { .checkInput_f(f)}

  prob<- p(C,w,G,Sens,D,r,f)
  n <- ceiling(log(1 - 0.9) / log(1 - prob))

  prob_seq<- as.vector(seq(from=0.1, to=1, by=0.01))
  n_seq<- as.vector(ceiling(log(1 - prob_seq) / log(1 - prob)))

  Plot<- graphics::plot(x=n_seq, y=prob_seq, type="l", ylim=c(0, 1), xlab="", ylab="")
 graphics::mtext(text="Desired Probability of Detecting Contamination (prob_det)",side=2,line=2.2, cex=0.8)
 graphics::mtext(text="Sample Size (n)", side=1, line=2.2, cex=0.8)
# add reference lines for prob_detect=0.9 and estimated n
  Plot + graphics::abline(h=0.9, col="red", lty=2, lwd=1)
  Plot + graphics::abline(v=n, col="red", lty=2, lwd=1)
rec= grDevices::recordPlot()
return(rec)

}

