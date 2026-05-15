#' @title Objective scores
#'
#' @description Extract the scores (optimal objective values) of the evaluated
#' DMUs from a stochastic DEA solution.
#' Note that these scores may not always be interpreted as efficiencies.
#'
#' @param x Object of class \code{dea_stoch} obtained with some of the stochastic DEA
#' \code{modelstoch_*} functions.
#' @param ... Other options (for compatibility reasons)
#'
#' @returns A matrix with the scores (optimal objective values).
#'
#' @author
#' \strong{Vicente Bolós} (\email{vicente.bolos@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' \strong{Rafael Benítez} (\email{rafael.suarez@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' \strong{Vicente Coll-Serrano} (\email{vicente.coll@@uv.es}).
#' \emph{Quantitative Methods for Measuring Culture (MC2). Applied Economics.}
#'
#' University of Valencia (Spain)
#'
#' @examples
#' \donttest{
#' # Example 1.
#' library(deaR)
#' data("Coll_Blasco_2006")
#' ni <- 2 # number of inputs
#' no <- 2 # number of outputs
#' data_example <- make_deadata(datadea = Coll_Blasco_2006,
#'                              ni = ni,
#'                              no = no)
#' nd <- length(data_example$dmunames) # number of DMUs
#' var_input <- matrix(1, nrow = ni, ncol = nd)
#' var_output <- matrix(1, nrow = no, ncol = nd)
#' data_stoch <- make_deadata_stoch(datadea = data_example,
#'                                  var_input = var_input,
#'                                  var_output = var_output)
#' # Evaluate the sixth DMU
#' Collstoch <- modelstoch_radial(data_stoch, dmu_eval = 6)
#' efficiencies(Collstoch)
#' }
#'
#' @export

efficiencies.dea_stoch <-
  function(x, ...) {
    deasol <- x
    if ("efficiency" %in% names(deasol$DMU[[1]])) {
      eff <- unlist(lapply(deasol$DMU, function(x)
        x$efficiency))
      if (length(deasol$DMU[[1]]$efficiency) > 1) {
        eff  <- do.call(rbind, lapply(deasol$DMU, function(x)
          x$efficiency))
        mean_eff <- unlist(lapply(deasol$DMU, function(x)
          x$mean_efficiency))
        eff <- cbind(eff, mean_eff)
      }
    } else if ("beta" %in% names(deasol$DMU[[1]])) {
      eff <- unlist(lapply(deasol$DMU, function(x)
        x$beta))
    } else if ("delta" %in% names(deasol$DMU[[1]])) {
      eff <- unlist(lapply(deasol$DMU, function(x)
        x$delta))
    } else if ("gamma" %in% names(deasol$DMU[[1]])) {
      eff <- unlist(lapply(deasol$DMU, function(x)
        x$gamma))
    } else if ("objval" %in% names(deasol$DMU[[1]])) {
      eff <- unlist(lapply(deasol$DMU, function(x)
        x$objval))
    } else {
      stop("No efficiency/beta/delta/objval parameters in this solution!")
    }

    return(round(eff, 6))

}
