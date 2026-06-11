#' @title Chance Constrained Radial Super-efficiency Models
#'
#' @description Solve chance constrained radial super-efficiency DEA models, based on the
#' Cooper et al. (2002) chance constrained radial efficiency models. Analogously to the deterministic case,
#' it removes the evaluated DMU from the set of reference DMUs \code{dmu_ref}
#' with respect to which it is evaluated.
#'
#' @usage
#' modelstoch_radial_supereff(datadea,
#'                            dmu_eval = NULL,
#'                            dmu_ref = NULL,
#'                            parallel = FALSE,
#'                            ...)
#'
#' @param datadea The data of class \code{deadata_stoch} with the expected values
#' of inputs and outputs.
#' @param dmu_eval A numeric vector containing which DMUs have to be evaluated.
#' If \code{NULL} (default), all DMUs are considered.
#' @param dmu_ref A numeric vector containing which DMUs are the evaluation reference set.
#' If \code{NULL} (default), all DMUs are considered.
#' @param parallel Logical, if \code{TRUE}, the DMUs are computed in parallel
#' (default \code{FALSE}).
#' @param ... Model parameters like \code{orientation} or \code{rts}, and other
#' parameters to be passed to the solver.
#'
#' @returns A list with the results for the evaluated DMUs and other parameters
#' for reproducibility.
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
#' @references
#' Cooper, W.W.; Deng, H.; Huang, Z.; Li, S.X. (2002). “Chance constrained programming
#' approaches to technical efficiencies and inefficiencies in stochastic data envelopment
#' analysis", Journal of the Operational Research Society, 53:12, 1347-1356.
#'
#' Land, K.C; Lovell, C.A.K.; Thore, S. (1993). "Chance-constrained data envelopment analysis",
#' Managerial and Decision Economics, Vol. 14, No. 6, 541-554.
#'
#' @note Radial super-efficiency chance constrained model under non constant
#' (vrs, nirs, ndrs, grs) returns to scale can be unfeasible for certain DMUs.
#'
#' @seealso \code{\link{modelstoch_radial}}
#'
#' @examples
#' \dontrun{
#' # Example 1. Replication of results in Land et al. (1993)
#' # with super-efficiencies and parallelized.
#'
#' library(deaR)
#' library(doParallel)
#' data("PFT1981")
#' # Selecting DMUs in Program Follow Through (PFT)
#' PFT <- PFT1981[1:49, ]
#' PFT <- make_deadata(PFT,
#'                     inputs = 2:6,
#'                     outputs = 7:9)
#' c <- 0.5
#' var_output <- matrix(c^2, nrow = 3, ncol = 49)
#' PFT_stoch <- make_deadata_stoch(datadea = PFT,
#'                                 var_output = var_output)
#' # Preparing parallelization
#' no_cores <- detectCores() - 1
#' cl <- makeCluster(no_cores)
#' registerDoParallel(no_cores)
#' # Evaluate all DMUs
#' res <- modelstoch_radial_supereff(PFT_stoch,
#'                                   parallel = TRUE)
#' stopImplicitCluster()
#' efficiencies(res)
#' }
#' @importFrom foreach foreach %dopar%
#'
#' @export

modelstoch_radial_supereff <-
  function(datadea,
           dmu_eval = NULL,
           dmu_ref = NULL,
           parallel = FALSE,
           ...) {

    # Checking whether datadea is of class "deadata_stoch" or not...
    if (!inherits(datadea, "deadata_stoch")) {
      stop("Data should be of class deadata_stoch. Run make_deadata_stoch function first!")
    }

  #optlist <- list(...)

  dmunames <- datadea$dmunames
  nd <- length(dmunames) # number of dmus

  if (is.null(dmu_eval)) {
    dmu_eval <- 1:nd
  } else if (!all(dmu_eval %in% (1:nd))) {
    stop("Invalid set of DMUs to be evaluated (dmu_eval).")
  }
  names(dmu_eval) <- dmunames[dmu_eval]
  nde <- length(dmu_eval)

  if (is.null(dmu_ref)) {
    dmu_ref <- 1:nd
  } else if (!all(dmu_ref %in% (1:nd))) {
    stop("Invalid set of reference DMUs (dmu_ref).")
  }
  names(dmu_ref) <- dmunames[dmu_ref]
  ndr <- length(dmu_ref)

  DMU <- vector(mode = "list", length = nde)

  # For loop ------
  if (parallel) {

    DMU <- foreach(i = seq_len(nde - 1)) %dopar% {

      ii <- dmu_eval[i]

      deasol <- modelstoch_radial(datadea,
                                  dmu_eval = ii,
                                  dmu_ref = dmu_ref[dmu_ref != ii],
                                  ...)

      DMU_i <- deasol$DMU[[1]]

      if ((ii %in% dmu_ref) && (!is.null(DMU_i$lambda))) {
        newlambda <- rep(0, ndr)
        # newlambda[dmu_ref == ii] <- 0
        newlambda[dmu_ref != ii] <- DMU_i$lambda
        names(newlambda) <- dmunames[dmu_ref]
        DMU_i$lambda <- newlambda
      }

      DMU_i
    }

    i <- nde
    ii <- dmu_eval[i]

    deasol <- modelstoch_radial(datadea,
                                dmu_eval = ii,
                                dmu_ref = dmu_ref[dmu_ref != ii],
                                ...)

    DMU[[i]] <- deasol$DMU[[1]]

    if ((ii %in% dmu_ref) && (!is.null(DMU[[i]]$lambda))) {
      newlambda <- rep(0, ndr)
      # newlambda[dmu_ref == ii] <- 0
      newlambda[dmu_ref != ii] <- DMU[[i]]$lambda
      names(newlambda) <- dmunames[dmu_ref]
      DMU[[i]]$lambda <- newlambda
    }

  } else {

    for (i in 1:nde) {

      ii <- dmu_eval[i]

      deasol <- modelstoch_radial(datadea,
                                  dmu_eval = ii,
                                  dmu_ref = dmu_ref[dmu_ref != ii],
                                  ...)

      DMU[[i]] <- deasol$DMU[[1]]

      if ((ii %in% dmu_ref) && (!is.null(DMU[[i]]$lambda))) {
        newlambda <- rep(0, ndr)
        # newlambda[dmu_ref == ii] <- 0
        newlambda[dmu_ref != ii] <- DMU[[i]]$lambda
        names(newlambda) <- dmunames[dmu_ref]
        DMU[[i]]$lambda <- newlambda
      }
    }
  }

  names(DMU) <- dmunames[dmu_eval]

  deaOutput <- deasol

  deaOutput$modelname <- "stoch_radial_supereff"
  deaOutput$DMU <- DMU
  deaOutput$dmu_eval <- dmu_eval
  deaOutput$dmu_ref <- dmu_ref

  return(deaOutput)

}
