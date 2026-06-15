#' @title Chance Constrained Additive E-model.
#'
#' @description It solves the chance constrained additive E-model based on the
#' deterministic additive model from Charnes et al. (1985), under constant and
#' non-constant returns to scale.
#'
#' Besides, the user can set weights for the input and/or output slacks. So, it
#' is also possible to solve chance constrained versions of weighted additive
#' models like Measure of Inefficiency Proportions (MIP) or Range Adjusted
#' Measure (RAM), see Cooper et al. (1999).
#'
#' We consider \eqn{\mathcal{D}=\left\{ \textrm{DMU}_1, \ldots ,\textrm{DMU}_n \right\} }
#' a set of \eqn{n} DMUs with \eqn{m} stochastic inputs and \eqn{s} stochastic outputs.
#' Matrices \eqn{\tilde{X}=(\tilde{x}_{ij})} and \eqn{\tilde{Y}=(\tilde{y}_{rj})}
#' are the input and output data matrices, respectively, where \eqn{\tilde{x}_{ij}}
#' and \eqn{\tilde{y}_{rj}} represent the \eqn{i}-th input and \eqn{r}-th output
#' of the \eqn{j}-th DMU. Moreover, we denote by \eqn{X=(x_{ij})} and \eqn{Y=(y_{rj})}
#' their expected values. In general, we denote vectors by bold-face letters and
#' they are considered as column vectors unless otherwise stated. The \eqn{0}-vector
#' is denoted by \eqn{\bm{0}} and the context determines its dimension.
#'
#' Given \eqn{0<\alpha <1}, the program for \eqn{\text{DMU}_o} with
#' constant returns to scale is given by
#' \deqn{\max \limits_{\bm{\lambda},\mathbf{s}^-,\mathbf{s}^+}\quad
#' \mathbf{w}^-\mathbf{s}^-+\mathbf{w}^+\mathbf{s}^+}
#' \deqn{\text{s.t.}\quad P\left\{ \left( \tilde{\mathbf{x}}_o-\tilde{X}
#' \bm{\lambda}-\mathbf{s}^-\right) _i\geq 0\right\}\geq 1-\alpha ,\qquad i=1,\ldots ,m,}
#' \deqn{P\left\{ \left( \tilde{Y}\bm{\lambda}-\tilde{\mathbf{y}}_o-\mathbf{s}^+\right) _r
#' \geq 0\right\}\geq 1-\alpha ,\qquad r=1,\ldots ,s,}
#' \deqn{\bm{\lambda}\geq \mathbf{0},\,\, \mathbf{s}^-\geq \mathbf{0},\,\,
#' \mathbf{s}^+\geq \mathbf{0},}
#' where \eqn{\bm{\lambda}=(\lambda_1,\ldots,\lambda_n)^\top}, \eqn{\tilde{\mathbf{x}}_o
#' =(\tilde{x}_{1o},\ldots,\tilde{x}_{mo})^\top} and \eqn{\tilde{\mathbf{y}}_o=
#' (\tilde{y}_{1o},\ldots,\tilde{y}_{so})^\top} are column vectors. Moreover,
#' \eqn{\mathbf{s}^-,\mathbf{s}^+} are column vectors with the slacks,
#' and \eqn{\mathbf{w}^-,\mathbf{w}^+} are positive row vectors with the weights
#' for the slacks.
#' Different returns to scale can be easily considered by adding the corresponding
#' constraints: \eqn{\mathbf{e}\bm{\lambda}=1} (VRS), \eqn{0\leq \mathbf{e}\bm{\lambda}
#' \leq 1} (NIRS), \eqn{\mathbf{e}\bm{\lambda}\geq 1} (NDRS) or \eqn{L\leq \mathbf{e}
#' \bm{\lambda}\leq U} (GRS), with \eqn{0\leq L\leq 1} and \eqn{U\geq 1}, where
#' \eqn{\mathbf{e}=(1,\ldots ,1)} is a row vector.
#'
#' The deterministic equivalent for a multivariate normal distribution of inputs/outputs
#' is given by
#' \deqn{\max \limits_{\bm{\lambda},\mathbf{s}^-,\mathbf{s}^+}\quad \mathbf{w}^-
#' \mathbf{s}^-+\mathbf{w}^+\mathbf{s}^+ }
#' \deqn{\text{s.t.}\quad \mathbf{x}_o-X\bm{\lambda}-\mathbf{s}^-+\Phi ^{-1}(\alpha)
#' \bm{\sigma} ^-(\bm{\lambda})\geq \mathbf{0},}
#' \deqn{Y\bm{\lambda}-\mathbf{y}_o-\mathbf{s}^++\Phi ^{-1}(\alpha)\bm{\sigma}^+
#' (\bm{\lambda})\geq \mathbf{0},}
#' \deqn{\bm{\lambda}\geq \mathbf{0},\,\, \mathbf{s}^-\geq \mathbf{0},\,\,
#' \mathbf{s}^+\geq \mathbf{0},}
#' where \eqn{\Phi } is the standard normal distribution, and
#' \deqn{\displaystyle \left( \sigma ^-_i\left( \bm{\lambda}\right)\right) ^2 =
#' \sum _{j,q=1}^n\lambda _j\lambda _q\mathrm{Cov}(\tilde{x}_{ij},\tilde{x}_{iq})
#' -2\sum _{j=1}^n\lambda _j\mathrm{Cov}(\tilde{x}_{ij},\tilde{x}_{io})}
#' \deqn{+\mathrm{Var}(\tilde{x}_{io}),\quad i=1,\ldots ,m,}
#' \deqn{\displaystyle \left( \sigma ^+_r\left( \bm{\lambda}\right)\right) ^2 =
#' \sum _{j,q=1}^n\lambda _j\lambda _q\mathrm{Cov}(\tilde{y}_{rj},\tilde{y}_{rq})
#' -2\sum _{j=1}^n\lambda _j\mathrm{Cov}(\tilde{y}_{rj},\tilde{y}_{ro})}
#' \deqn{+\mathrm{Var}(\tilde{y}_{ro}),\quad r=1,\ldots ,s.}
#'
#' @note A DMU is \eqn{\alpha}-stochastically efficient if and only if the optimal
#' objective value of the problem, (\code{objval}), is zero (or less than zero).
#'
#' @usage
#' modelstoch_additive(datadea,
#'                     alpha = 0.05,
#'                     dmu_eval = NULL,
#'                     dmu_ref = NULL,
#'                     orientation = NULL,
#'                     weight_slack_i = 1,
#'                     weight_slack_o = 1,
#'                     rts = c("crs", "vrs", "nirs", "ndrs", "grs"),
#'                     L = 1,
#'                     U = 1,
#'                     solver = c("alabama"),
#'                     X = NULL,
#'                     give_X = FALSE,
#'                     n_attempts_max = 5,
#'                     returnqp = FALSE,
#'                     parallel = FALSE,
#'                     ...)
#'
#' @param datadea The data of class \code{deadata_stoch}, including \code{n} DMUs,
#' and the expected values of \code{m} inputs and \code{s} outputs.
#' @param alpha A value for parameter alpha.
#' @param dmu_eval A numeric vector containing which DMUs have to be evaluated.
#' If \code{NULL} (default), all DMUs are considered.
#' @param dmu_ref A numeric vector containing which DMUs are the evaluation reference set.
#' If \code{NULL} (default), all DMUs are considered.
#' @param orientation This parameter is either \code{NULL} (default) or a string, equal to
#' "io" (input-oriented) or "oo" (output-oriented). It is used to modify the weight slacks.
#' If input-oriented, \code{weight_slack_o} are taken 0.
#' If output-oriented, \code{weight_slack_i} are taken 0.
#' @param weight_slack_i A value, vector of length \code{m}, or matrix \code{m} x
#' \code{ne} (where \code{ne} is the length of \code{dmu_eval})
#' with the weights of the input slacks. If 0, output-oriented.
#' @param weight_slack_o A value, vector of length \code{s}, or matrix \code{s} x
#' \code{ne} (where \code{ne} is the length of \code{dmu_eval})
#' with the weights of the output slacks. If 0, input-oriented.
#' @param rts A string, determining the type of returns to scale, equal to "crs" (constant),
#' "vrs" (variable), "nirs" (non-increasing), "ndrs" (non-decreasing) or "grs" (generalized).
#' @param L Lower bound for the generalized returns to scale (grs).
#' @param U Upper bound for the generalized returns to scale (grs).
#' @param solver Character string with the name of the solver used by function \code{solvecop}
#' from package \code{optiSolve}. Only alabama solver is currently operational, but more solvers
#' could be implemented in the future.
#' @param X Initial vector for the solver in the first attempt. The variables are
#' lambdas of reference DMUs, input slacks, output slacks, input sigmas and output sigmas. It can be a matrix
#' of \code{ne} rows (where \code{ne} is the length of \code{dmu_eval}), with different initial
#' vectors for each evaluated DMU. If it is a vector, the same initial vector is applied to all
#' evaluated DMUs.
#' @param give_X Logical. If it is \code{TRUE} and \code{X} is \code{NULL}, it makes
#' the initial vector to be the variables corresponding to each evaluated DMU.
#' If it is \code{FALSE} (default) and \code{X} is \code{NULL}, the initial
#' vector is given internally by the solver (usually randomly generated).
#' @param n_attempts_max A value with the maximum number of attempts if the solver
#' does not converge. Each attempt uses a different initial vector.
#' @param returnqp Logical. If it is \code{TRUE}, it returns the quadratic
#' problems (objective function and constraints).
#' @param parallel Logical, if \code{TRUE}, the DMUs are computed in parallel
#' (default \code{FALSE}).
#' @param ... Other parameters to be passed to the solver.
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
#' Charnes, A.; Cooper, W.W.; Golany, B.; Seiford, L.; Stuz, J. (1985) "Foundations
#' of Data Envelopment Analysis for Pareto-Koopmans Efficient Empirical Production
#' Functions", Journal of Econometrics, 30(1-2), 91-107.
#' \doi{10.1016/0304-4076(85)90133-2}
#'
#' Cooper, W.W.; Park, K.S.; Pastor, J.T. (1999). "RAM: A Range Adjusted Measure
#' of Inefficiencies for Use with Additive Models, and Relations to Other Models
#' and Measures in DEA". Journal of Productivity Analysis, 11, p. 5-42.
#' \doi{10.1023/A:1007701304281}
#'
#' @import optiSolve stats
#' @importFrom foreach foreach %dopar%
#'
#' @export

modelstoch_additive <-
  function(datadea,
           alpha = 0.05,
           dmu_eval = NULL,
           dmu_ref = NULL,
           orientation = NULL,
           weight_slack_i = 1,
           weight_slack_o = 1,
           rts = c("crs", "vrs", "nirs", "ndrs", "grs"),
           L = 1,
           U = 1,
           solver = c("alabama"),
           X = NULL,
           give_X = FALSE,
           n_attempts_max = 5,
           returnqp = FALSE,
           parallel = FALSE,
           ...) {

    # Checking whether datadea is of class "deadata_stoch" or not...
    if (!inherits(datadea, "deadata_stoch")) {
      stop("Data should be of class deadata_stoch. Run make_deadata_stoch function first!")
    }

    # Checking alpha
    if (!((alpha > 0) && (alpha < 1))) {
      stop("Parameter alpha must be between 0 and 1.")
    }
    qa <- qnorm(alpha)

    # Checking non-discretionary inputs/outputs
    if ((!is.null(datadea$nd_inputs)) || (!is.null(datadea$nd_outputs))) {
      warning("This model does not take into account the non-discretionary feature for inputs/outputs.")
    }

    # Checking undesirable inputs/outputs
    if (!is.null(datadea$ud_inputs) || !is.null(datadea$ud_outputs)) {
      warning("This model does not take into account the undesirable feature for inputs/outputs.")
    }

    # Checking rts
    rts <- tolower(rts)
    rts <- match.arg(rts)

    if (rts == "grs") {
      if (L > 1) {
        stop("L must be <= 1.")
      }
      if (U < 1) {
        stop("U must be >= 1.")
      }
    }

    # Checking solver
    solver <- tolower(solver)
    solver <- match.arg(solver)

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

    input <- datadea$input
    output <- datadea$output
    inputnames <- rownames(input)
    outputnames <- rownames(output)
    ni <- nrow(input) # number of  inputs
    no <- nrow(output) # number of outputs
    inputref <- matrix(input[, dmu_ref], nrow = ni)
    outputref <- matrix(output[, dmu_ref], nrow = no)

    covX <- array(0, dim = c(ni, nd, nd))
    for (i in 1:ni) {
      covX[i, , ] <- datadea$cov_II[i, i, , ]
    }
    covY <- array(0, dim = c(no, nd, nd))
    for (i in 1:no) {
      covY[i, , ] <- datadea$cov_OO[i, i, , ]
    }
    covXref <- array(covX[, dmu_ref, dmu_ref], dim = c(ni, ndr, ndr))
    covYref <- array(covY[, dmu_ref, dmu_ref], dim = c(no, ndr, ndr))

    nc_inputs <- datadea$nc_inputs
    nc_outputs <- datadea$nc_outputs
    nnci <- length(nc_inputs)
    nnco <- length(nc_outputs)

    # Checking weights
    if(is.null(weight_slack_i)){
      weight_slack_i <- 1
    }
    if(is.null(weight_slack_o)){
      weight_slack_o <- 1
    }

    if (is.matrix(weight_slack_i)) {
      if ((nrow(weight_slack_i) != ni) || (ncol(weight_slack_i) != nde)) {
        stop("Invalid weight input matrix (number of inputs x number of evaluated DMUs).")
      }
    } else if ((length(weight_slack_i) == 1) || (length(weight_slack_i) == ni)) {
      weight_slack_i <- matrix(weight_slack_i, nrow = ni, ncol = nde)
    } else {
      stop("Invalid weight input vector (number of inputs).")
    }
    if ((!is.null(orientation)) && (orientation == "oo")) {
      weight_slack_i <- matrix(0, nrow = ni, ncol = nde)
    }
    rownames(weight_slack_i) <- inputnames
    colnames(weight_slack_i) <- dmunames[dmu_eval]

    if (is.matrix(weight_slack_o)) {
      if ((nrow(weight_slack_o) != no) || (ncol(weight_slack_o) != nde)) {
        stop("Invalid weight output matrix (number of outputs x number of evaluated DMUs).")
      }
    } else if ((length(weight_slack_o) == 1) || (length(weight_slack_o) == no)) {
      weight_slack_o <- matrix(weight_slack_o, nrow = no, ncol = nde)
    } else {
      stop("Invalid weight output vector (number of outputs).")
    }
    if ((!is.null(orientation)) && (orientation == "io")) {
      weight_slack_o <- matrix(0, nrow = no, ncol = nde)
    }
    rownames(weight_slack_o) <- outputnames
    colnames(weight_slack_o) <- dmunames[dmu_eval]

    namevar <- c(paste("lambda", 1:ndr, sep = "_"),
                 paste("slack_I", 1:ni, sep = "_"),
                 paste("slack_O", 1:no, sep = "_"),
                 paste("sigma_I", 1:ni, sep = "_"),
                 paste("sigma_O", 1:no, sep = "_"))
    nvar <- ndr + 2 * (ni + no)

    # Checking initial vector X
    if (!is.null(X)) {
      if (is.numeric(X) && is.vector(X)) {
        if (length(X) != nvar)
          stop("Incorrect number of variables for initial vector X.")
        X <- matrix(rep(X, nde), nrow = nde, byrow = TRUE)
      } else if (is.numeric(X) && is.matrix(X)) {
        if (nrow(X) != nde || ncol(X) != nvar)
          stop("Incorrect matrix X for initial vectors.")
      } else {
        stop("Initial vector must be a numeric vector or matrix.")
      }
    }

    ###########################

    # Lower and upper bounds constraints
    lbcon1 <- lbcon(val = 0, id = namevar)
    ubcon1 <- NULL

    if (rts == "crs") {
      f.con.rs <- NULL
      f.dir.rs <- NULL
      f.rhs.rs <- NULL
    } else {
      f.con.rs <- cbind(matrix(1, nrow = 1, ncol = ndr), matrix(0, nrow = 1, ncol = 2 * (ni + no)))
      f.rhs.rs <- 1
      if (rts == "vrs") {
        f.dir.rs <- "=="
        ubcon1 <- ubcon(val = rep(1, ndr), id = namevar[1:ndr])
      } else if (rts == "nirs") {
        f.dir.rs <- "<="
        ubcon1 <- ubcon(val = rep(1, ndr), id = namevar[1:ndr])
      } else if (rts == "ndrs") {
        f.dir.rs <- ">="
      } else {
        f.con.rs <- rbind(f.con.rs, f.con.rs)
        f.dir.rs <- c(">=", "<=")
        f.rhs.rs <- c(L, U)
        ubcon1 <- ubcon(val = rep(U, ndr), id = namevar[1:ndr])
      }
    }

    # Linear constraints matrix
    f.con.1 <- cbind(inputref, diag(ni), matrix(0, nrow = ni, ncol = no),
                     diag(-qa, ni), matrix(0, nrow = ni, ncol = no))
    f.con.1[nc_inputs, ndr + nc_inputs] <- 0
    f.con.1[nc_inputs, ndr + ni + no + nc_inputs] <- 0
    f.con.2 <- cbind(outputref, matrix(0, nrow = no, ncol = ni), -diag(no),
                     matrix(0, nrow = no, ncol = ni), diag(qa, no))
    f.con.2[nc_outputs, ndr + ni + nc_outputs] <- 0
    f.con.2[nc_outputs, ndr + 2 * ni + no + nc_outputs] <- 0
    f.con.nc <- matrix(0, nrow = (nnci + nnco), ncol = (nvar))
    f.con.nc[, ndr + c(nc_inputs, ni + nc_outputs)] <- diag(nnci + nnco)
    f.con <- rbind(f.con.1, f.con.2, f.con.nc, f.con.rs)
    rownames(f.con) <- paste("lc", 1:nrow(f.con), sep = "") # to prevent names errors in lincon

    # Linear directions vector
    f.dir <- c(rep("==", ni + no + nnci + nnco), f.dir.rs)

    # Quadratic matrices
    Qin <- array(0, dim = c(ni, nvar, nvar))
    Qout <- array(0, dim = c(no, nvar, nvar))
    Qin[, 1:ndr, 1:ndr] <- covXref
    Qout[, 1:ndr, 1:ndr] <- covYref
    for (j in 1:ni) {
      Qin[j, ndr + ni + no + j, ndr + ni + no + j] <- -1
    }
    for (j in 1:no) {
      Qout[j, ndr + 2 * ni + no + j, ndr + 2 * ni + no + j] <- -1
    }

    # Linear part of quadratic constraints
    ain <- matrix(0, nrow = ni, ncol = nvar)
    aout <- matrix(0, nrow = no, ncol = nvar)

    qclist <- vector(mode = "list", length = 2 * (ni + no))

    # For loop ------
    if (parallel) {

      DMU <- foreach(i = 1:nde) %dopar% {

        ii <- dmu_eval[i]

        # Objective function coefficients
        f.obj <- linfun(a = c(rep(0, ndr), weight_slack_i[, i],
                              weight_slack_o[, i], rep(0, ni + no)),
                        id = namevar)

        # Right hand side vector
        f.rhs <- c(input[, ii], output[, ii], rep(0, nnci + nnco), f.rhs.rs)

        lincon <- lincon(A = f.con, dir = f.dir, val = f.rhs, id = namevar)

        covXo <- matrix(covX[, ii, dmu_ref], nrow = ni, ncol = ndr)
        covYo <- matrix(covY[, ii, dmu_ref], nrow = no, ncol = ndr)
        ain[, 1:ndr] <- -2 * covXo
        valin <- -covX[, ii, ii]
        aout[, 1:ndr] <- -2 * covYo
        valout <- -covY[, ii, ii]

        quadconina <- vector(mode = "list", length = ni)
        quadconinb <- quadconina
        quadconouta <- vector(mode = "list", length = no)
        quadconoutb <- quadconouta
        for (j in 1:ni) {
          if (j %in% nc_inputs) {
            quadconina[[j]] <- NULL
            quadconinb[[j]] <- NULL
          } else {
            quadconina[[j]] <- quadcon(Q = Qin[j, , ], a = ain[j, ], val = valin[j], id = namevar)
            quadconinb[[j]] <- quadcon(Q = -Qin[j, , ], a = -ain[j, ], val = -valin[j], id = namevar)
          }
        }
        for (j in 1:no) {
          if (j %in% nc_outputs) {
            quadconouta[[j]] <- NULL
            quadconoutb[[j]] <- NULL
          } else {
            quadconouta[[j]] <- quadcon(Q = Qout[j, , ], a = aout[j, ], val = valout[j], id = namevar)
            quadconoutb[[j]] <- quadcon(Q = -Qout[j, , ], a = -aout[j, ], val = -valout[j], id = namevar)
          }
        }

        mycop <- cop(f = f.obj, max = TRUE, lb = lbcon1, ub = ubcon1, lc = lincon)
        qclist <- c(quadconina, quadconinb, quadconouta, quadconoutb)
        auxkk <- which(sapply(qclist, is.null))
        if (length(auxkk) > 0) {
          qclist <- qclist[-auxkk] # remove NULLs from non-controllable
        }
        mycop$qc <- qclist

        if (returnqp) {

          mycop

        } else {

          n_attempts <- 1

          while (n_attempts <= n_attempts_max) {

            # Initial vector
            Xini <- NULL
            if (n_attempts == 1) {
              if (is.matrix(X)) {
                Xini <- X[i, ]
                names(Xini) <- namevar
              } else if (give_X && (ii %in% dmu_ref)) {
                Xini <- rep(0, nvar)
                Xini[which(dmu_ref == ii)] <- 1
                names(Xini) <- namevar
              }
            }

            res <- solvecop(op = mycop, solver = solver, quiet = TRUE, X = Xini, ...)

            if (res$status == "successful convergence") {
              n_attempts <- n_attempts_max
            }
            n_attempts <- n_attempts + 1

          }

          if (res$status == "successful convergence") {

            res <- res$x

            lambda <- res[1:ndr]
            names(lambda) <- dmunames[dmu_ref]
            slack_input <- res[(ndr + 1) : (ndr + ni)]
            names(slack_input) <- inputnames
            slack_output <- res[(ndr + ni + 1) : (ndr + ni + no)]
            names(slack_output) <- outputnames
            sigma_input <- res[(ndr + ni + no + 1):(ndr + ni + no + ni)]
            names(sigma_input) <- inputnames
            sigma_output <- res[(ndr + 2 * ni + no + 1):nvar]
            names(sigma_output) <- outputnames

            target_input <- as.vector(inputref %*% lambda)
            target_output <- as.vector(outputref %*% lambda)
            names(target_input) <- inputnames
            names(target_output) <- outputnames

            objval <- weight_slack_i[, i] %*% slack_input + weight_slack_o[, i] %*% slack_output

          } else {

            objval <- NA
            lambda <- NA
            sigma_input <- NA
            sigma_output <- NA
            slack_input <- NA
            slack_output <- NA
            target_input <- NA
            target_output <- NA

          }

          list(objval = objval,
               lambda = lambda,
               slack_input = slack_input, slack_output = slack_output,
               target_input = target_input, target_output = target_output,
               sigma_input = sigma_input, sigma_output = sigma_output)

        }
      }

    } else {

      DMU <- vector(mode = "list", length = nde)

      for (i in 1:nde) {

        ii <- dmu_eval[i]

        # Objective function coefficients
        f.obj <- linfun(a = c(rep(0, ndr), weight_slack_i[, i],
                              weight_slack_o[, i], rep(0, ni + no)),
                        id = namevar)

        # Right hand side vector
        f.rhs <- c(input[, ii], output[, ii], rep(0, nnci + nnco), f.rhs.rs)

        lincon <- lincon(A = f.con, dir = f.dir, val = f.rhs, id = namevar)

        covXo <- matrix(covX[, ii, dmu_ref], nrow = ni, ncol = ndr)
        covYo <- matrix(covY[, ii, dmu_ref], nrow = no, ncol = ndr)
        ain[, 1:ndr] <- -2 * covXo
        valin <- -covX[, ii, ii]
        aout[, 1:ndr] <- -2 * covYo
        valout <- -covY[, ii, ii]

        quadconina <- vector(mode = "list", length = ni)
        quadconinb <- quadconina
        quadconouta <- vector(mode = "list", length = no)
        quadconoutb <- quadconouta
        for (j in 1:ni) {
          if (j %in% nc_inputs) {
            quadconina[[j]] <- NULL
            quadconinb[[j]] <- NULL
          } else {
            quadconina[[j]] <- quadcon(Q = Qin[j, , ], a = ain[j, ], val = valin[j], id = namevar)
            quadconinb[[j]] <- quadcon(Q = -Qin[j, , ], a = -ain[j, ], val = -valin[j], id = namevar)
          }
        }
        for (j in 1:no) {
          if (j %in% nc_outputs) {
            quadconouta[[j]] <- NULL
            quadconoutb[[j]] <- NULL
          } else {
            quadconouta[[j]] <- quadcon(Q = Qout[j, , ], a = aout[j, ], val = valout[j], id = namevar)
            quadconoutb[[j]] <- quadcon(Q = -Qout[j, , ], a = -aout[j, ], val = -valout[j], id = namevar)
          }
        }

        mycop <- cop(f = f.obj, max = TRUE, lb = lbcon1, ub = ubcon1, lc = lincon)
        qclist <- c(quadconina, quadconinb, quadconouta, quadconoutb)
        auxkk <- which(sapply(qclist, is.null))
        if (length(auxkk) > 0) {
          qclist <- qclist[-auxkk] # remove NULLs from non-controllable
        }
        mycop$qc <- qclist

        if (returnqp) {

          DMU[[i]] <- mycop

        } else {

          n_attempts <- 1

          while (n_attempts <= n_attempts_max) {

            # Initial vector
            Xini <- NULL
            if (n_attempts == 1) {
              if (is.matrix(X)) {
                Xini <- X[i, ]
                names(Xini) <- namevar
              } else if (give_X && (ii %in% dmu_ref)) {
                Xini <- rep(0, nvar)
                Xini[which(dmu_ref == ii)] <- 1
                names(Xini) <- namevar
              }
            }

            res <- solvecop(op = mycop, solver = solver, quiet = TRUE, X = Xini, ...)

            if (res$status == "successful convergence") {
              n_attempts <- n_attempts_max
            }
            n_attempts <- n_attempts + 1

          }

          if (res$status == "successful convergence") {

            res <- res$x

            lambda <- res[1:ndr]
            names(lambda) <- dmunames[dmu_ref]
            slack_input <- res[(ndr + 1) : (ndr + ni)]
            names(slack_input) <- inputnames
            slack_output <- res[(ndr + ni + 1) : (ndr + ni + no)]
            names(slack_output) <- outputnames
            sigma_input <- res[(ndr + ni + no + 1):(ndr + ni + no + ni)]
            names(sigma_input) <- inputnames
            sigma_output <- res[(ndr + 2 * ni + no + 1):nvar]
            names(sigma_output) <- outputnames

            target_input <- as.vector(inputref %*% lambda)
            target_output <- as.vector(outputref %*% lambda)
            names(target_input) <- inputnames
            names(target_output) <- outputnames

            objval <- weight_slack_i[, i] %*% slack_input + weight_slack_o[, i] %*% slack_output

          } else {

            objval <- NA
            lambda <- NA
            sigma_input <- NA
            sigma_output <- NA
            slack_input <- NA
            slack_output <- NA
            target_input <- NA
            target_output <- NA

          }

          DMU[[i]] <- list(objval = objval,
                           lambda = lambda,
                           slack_input = slack_input, slack_output = slack_output,
                           target_input = target_input, target_output = target_output,
                           sigma_input = sigma_input, sigma_output = sigma_output)

        }
      }
    }

    names(DMU) <- dmunames[dmu_eval]

    # Checking if a DMU is in its own reference set (when rts = "grs")
    if (rts == "grs") {
      eps <- 1e-6
      for (i in 1:nde) {
        j <- which(dmu_ref == dmu_eval[i])
        if (length(j) == 1) {
          kk <- DMU[[i]]$lambda[j]
          kk2 <- sum(DMU[[i]]$lambda[-j])
          if ((kk > eps) && (kk2 > eps)) {
            warning(paste("Under generalized returns to scale,", dmunames[dmu_eval[i]],
                          "appears in its own reference set."))
          }
        }
      }
    }

    deaOutput <- list(modelname = "stoch_additive_e",
                      rts = rts,
                      L = L,
                      U = U,
                      DMU = DMU,
                      data = datadea,
                      alpha = alpha,
                      dmu_eval = dmu_eval,
                      dmu_ref = dmu_ref,
                      weight_slack_i = weight_slack_i,
                      weight_slack_o = weight_slack_o,
                      orientation = orientation)

    return(structure(deaOutput, class = "dea_stoch"))

  }
