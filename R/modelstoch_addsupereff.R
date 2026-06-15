#' @title Chance Constrained Super-efficiency Additive E-model.
#'
#' @description It solves the chance constrained super-efficiency additive
#' E-model based on the deterministic super-efficiency additive model from Du et
#' al. (2010), under constant and non-constant returns to scale.
#' Besides, the user can set weights for the input and/or output slacks.
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
#' \deqn{\min \limits_{\bm{\lambda},\mathbf{t}^-,\mathbf{t}^+}\quad
#' \mathbf{w}^-\mathbf{t}^-+\mathbf{w}^+\mathbf{t}^+}
#' \deqn{\text{s.t.}\quad P\left\{ \left( \tilde{\mathbf{x}}_o-\tilde{X}_{-o}
#' \bm{\lambda}+\mathbf{t}^-\right) _i\geq 0\right\}\geq 1-\alpha ,\qquad i=1,\ldots ,m,}
#' \deqn{P\left\{ \left( \tilde{Y}_{-o}\bm{\lambda}-\tilde{\mathbf{y}}_o+
#' \mathbf{t}^+\right) _r\geq 0\right\}\geq 1-\alpha ,\qquad r=1,\ldots ,s,}
#' \deqn{\bm{\lambda}\geq \mathbf{0},\,\, \mathbf{t}^-\geq \mathbf{0},\,\,
#' \mathbf{t}^+\geq \mathbf{0},}
#' where \eqn{\tilde{X}_{-o},\tilde{Y}_{-o}} are the input and output data matrices, respectively,
#' defined by \eqn{\mathcal{D}-\left\{ \textrm{DMU}_o\right\}}, \eqn{\bm{\lambda}
#' =(\lambda_1,\ldots ,\lambda_{o-1},\lambda_{o+1},\ldots ,\lambda_n)^{\top}},
#' \eqn{\tilde{\mathbf{x}}_o=(\tilde{x}_{1o},\ldots,\tilde{x}_{mo})^\top} and
#' \eqn{\tilde{\mathbf{y}}_o=(\tilde{y}_{1o},\ldots,\tilde{y}_{so})^\top} are
#' column vectors. Moreover, \eqn{\mathbf{s}^-,\mathbf{s}^+} are column vectors
#' with the slacks, and \eqn{\mathbf{w}^-,\mathbf{w}^+} are positive row vectors
#' with the weights for the slacks.
#' Different returns to scale can be easily considered by adding the corresponding
#' constraints: \eqn{\mathbf{e}\bm{\lambda}=1} (VRS), \eqn{0\leq \mathbf{e}\bm{\lambda}
#' \leq 1} (NIRS), \eqn{\mathbf{e}\bm{\lambda}\geq 1} (NDRS) or \eqn{L\leq \mathbf{e}
#' \bm{\lambda}\leq U} (GRS), with \eqn{0\leq L\leq 1} and \eqn{U\geq 1}, where
#' \eqn{\mathbf{e}=(1,\ldots ,1)} is a row vector.
#'
#' The deterministic equivalent for a multivariate normal distribution of inputs/outputs
#' is given by
#' \deqn{\min \limits_{\bm{\lambda},\mathbf{t}^-,\mathbf{t}^+}\quad
#' \mathbf{w}^-\mathbf{t}^-+\mathbf{w}^+\mathbf{t}^+}
#' \deqn{\text{s.t.}\quad \mathbf{x}_o-X_{-o}\bm{\lambda}+\mathbf{t}^-+\Phi ^{-1}
#' (\alpha)\bm{\sigma} ^-(\bm{\lambda})\geq \mathbf{0},}
#' \deqn{Y_{-o}\bm{\lambda}-\mathbf{y}_o+\mathbf{t}^++\Phi ^{-1}(\alpha)
#' \bm{\sigma}^+(\bm{\lambda})\geq \mathbf{0},}
#' \deqn{\bm{\lambda}\geq \mathbf{0},\,\, \mathbf{t}^-\geq \mathbf{0},\,\,
#' \mathbf{t}^+\geq \mathbf{0},}
#' where \eqn{\Phi } is the standard normal distribution, and
#' \deqn{\displaystyle \left( \sigma ^-_i\left( \bm{\lambda}\right)\right) ^2=\sum _
#' {\substack{j,q=1\\j,q\neq o}}^n\lambda _j\lambda _q\mathrm{Cov}(\tilde{x}_{ij},
#' \tilde{x}_{iq})-2\sum _{\substack{j=1\\j\neq o}}^n\lambda _j\mathrm{Cov}
#' (\tilde{x}_{ij},\tilde{x}_{io})}
#' \deqn{+\mathrm{Var}(\tilde{x}_{io}),\quad i=1,\ldots ,m,}
#' \deqn{\displaystyle \left( \sigma ^+_r\left( \bm{\lambda}\right)\right) ^2=\sum _
#' {\substack{j,q=1\\j,q\neq o}}^n\lambda _j\lambda _q\mathrm{Cov}
#' (\tilde{y}_{rj},\tilde{y}_{rq})-2\sum _{\substack{j=1\\j\neq o}}^n
#' \lambda _j\mathrm{Cov}(\tilde{y}_{rj},\tilde{y}_{ro})}
#' \deqn{+\mathrm{Var}(\tilde{y}_{ro}),\quad r=1,\ldots ,s.}
#'
#' @note A DMU is \eqn{\alpha}-stochastically efficient if and only if the optimal
#' objective value of the problem, (\code{objval}), is zero (or less than zero).
#'
#' @usage
#' modelstoch_addsupereff(datadea,
#'                        alpha = 0.05,
#'                        dmu_eval = NULL,
#'                        dmu_ref = NULL,
#'                        orientation = NULL,
#'                        weight_slack_i = NULL,
#'                        weight_slack_o = NULL,
#'                        rts = c("crs", "vrs", "nirs", "ndrs", "grs"),
#'                        L = 1,
#'                        U = 1,
#'                        solver = c("alabama"),
#'                        X = NULL,
#'                        n_attempts_max = 5,
#'                        compute_target = TRUE,
#'                        returnqp = FALSE,
#'                        parallel = FALSE,
#'                        ...)
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
#' with the weights of the input super-slacks (\code{t_input}).
#' If 0, output-oriented.
#' If \code{weight_slack_i} is the matrix of the inverses of inputs of DMUS in
#' \code{dmu_eval} (default), the model is unit invariant.
#' @param weight_slack_o A value, vector of length \code{s}, or matrix \code{s} x
#' \code{ne} (where \code{ne} is the length of \code{dmu_eval})
#' with the weights of the output super-slacks (\code{t_output}).
#' If 0, input-oriented.
#' If \code{weight_slack_o} is the matrix of the inverses of outputs of DMUS in
#' \code{dmu_eval} (default), the model is unit invariant.
#' @param rts A string, determining the type of returns to scale, equal to "crs" (constant),
#' "vrs" (variable), "nirs" (non-increasing), "ndrs" (non-decreasing) or "grs" (generalized).
#' @param L Lower bound for the generalized returns to scale (grs).
#' @param U Upper bound for the generalized returns to scale (grs).
#' @param solver Character string with the name of the solver used by function \code{solvecop}
#' from package \code{optiSolve}. Only alabama solver is currently operational, but more solvers
#' could be implemented in the future.
#' @param X Initial vector for the solver in the first attempt. The variables are
#' lambdas of reference DMUs, input t, output t, input sigmas and output sigmas. It can be a matrix
#' of \code{ne} rows (where \code{ne} is the length of \code{dmu_eval}), with different initial
#' vectors for each evaluated DMU. If it is a vector, the same initial vector is applied to all
#' evaluated DMUs.
#' @param n_attempts_max A value with the maximum number of attempts if the solver
#' does not converge. Each attempt uses a different initial vector.
#' @param compute_target Logical. If it is \code{TRUE}, it computes targets,
#' projections and slacks.
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
#' Du, J.; Liang, L.; Zhu, J. (2010). "A Slacks-based Measure of Super-efficiency
#' in Data Envelopment Analysis. A Comment", European Journal of Operational Research,
#' 204, 694-697. \doi{10.1016/j.ejor.2009.12.007}
#'
#' @import optiSolve stats
#' @importFrom foreach foreach %dopar%
#'
#' @export

modelstoch_addsupereff <-
  function(datadea,
           alpha = 0.05,
           dmu_eval = NULL,
           dmu_ref = NULL,
           orientation = NULL,
           weight_slack_i = NULL,
           weight_slack_o = NULL,
           rts = c("crs", "vrs", "nirs", "ndrs", "grs"),
           L = 1,
           U = 1,
           solver = c("alabama"),
           X = NULL,
           n_attempts_max = 5,
           compute_target = TRUE,
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
    if (is.null(weight_slack_i)) {
      weight_slack_i <- matrix(1 / input[, dmu_eval], nrow = ni) / (ni + no - nnci - nnco)
    } else {
      if (is.matrix(weight_slack_i)) {
        if ((nrow(weight_slack_i) != ni) || (ncol(weight_slack_i) != nde)) {
          stop("Invalid matrix of weights of the input slacks (number of inputs x number of evaluated DMUs).")
        }
      } else if ((length(weight_slack_i) == 1) || (length(weight_slack_i) == ni)) {
        weight_slack_i <- matrix(weight_slack_i, nrow = ni, ncol = nde)
      } else {
        stop("Invalid vector of weights of the input slacks.")
      }
    }
    weight_slack_i[nc_inputs, ] <- 0
    if ((!is.null(orientation)) && (orientation == "oo")) {
      weight_slack_i <- matrix(0, nrow = ni, ncol = nde)
    }
    rownames(weight_slack_i) <- inputnames
    colnames(weight_slack_i) <- dmunames[dmu_eval]

    if (is.null(weight_slack_o)) {
      weight_slack_o <- matrix(1 / output[, dmu_eval], nrow = no) / (ni + no - nnci - nnco)
    } else {
      if (is.matrix(weight_slack_o)) {
        if ((nrow(weight_slack_o) != no) || (ncol(weight_slack_o) != nde)) {
          stop("Invalid matrix of weights of the output slacks (number of outputs x number of evaluated DMUs).")
        }
      } else if ((length(weight_slack_o) == 1) || (length(weight_slack_o) == no)) {
        weight_slack_o <- matrix(weight_slack_o, nrow = no, ncol = nde)
      } else {
        stop("Invalid vector of weights of the output slacks.")
      }
    }
    weight_slack_o[nc_outputs, ] <- 0
    if ((!is.null(orientation)) && (orientation == "io")) {
      weight_slack_o <- matrix(0, nrow = no, ncol = nde)
    }
    rownames(weight_slack_o) <- outputnames
    colnames(weight_slack_o) <- dmunames[dmu_eval]

    namevar <- c(paste("lambda", 1:ndr, sep = "_"),
                 paste("t_I", 1:ni, sep = "_"),
                 paste("t_O", 1:no, sep = "_"),
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

    target_input <- NULL
    target_output <- NULL
    project_input <- NULL
    project_output <- NULL
    slack_input <- NULL
    slack_output <- NULL

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
    f.con.1 <- cbind(inputref, -diag(ni), matrix(0, nrow = ni, ncol = no),
                     diag(-qa, ni), matrix(0, nrow = ni, ncol = no))
    f.con.2 <- cbind(outputref, matrix(0, nrow = no, ncol = ni), diag(no),
                     matrix(0, nrow = no, ncol = ni), diag(qa, no))

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

        w0i <- which(weight_slack_i[, i] == 0)
        nw0i <- length(w0i)
        w0o <- which(weight_slack_o[, i] == 0)
        nw0o <- length(w0o)

        # Objective function coefficients
        f.obj <- linfun(a = c(rep(0, ndr), weight_slack_i[, i],
                              weight_slack_o[, i], rep(0, ni + no)),
                        id = namevar)

        # Linear constraints matrix
        f.con.se <- rep(0, ndr)
        f.con.se[dmu_ref == ii] <- 1
        f.con.se <- c(f.con.se, rep(0, 2 * (ni + no)))
        f.con.w0 <- matrix(0, nrow = (nw0i + nw0o), ncol = nvar)
        f.con.w0[, ndr + c(w0i, ni + w0o)] <- diag(nw0i + nw0o)
        f.con <- rbind(f.con.1, f.con.2, f.con.se, f.con.w0, f.con.rs)
        rownames(f.con) <- paste("lc", 1:nrow(f.con), sep = "") # to prevent names errors in lincon

        # Linear directions vector
        f.dir <- c(rep("<=", ni), rep(">=", no), rep("==", 1 + nw0i + nw0o), f.dir.rs)
        f.dir[nc_inputs] <- "=="
        f.dir[ni + nc_outputs] <- "=="

        # Linear right hand side vector
        f.rhs <- c(input[, ii], output[, ii], rep(0, 1 + nw0i + nw0o), f.rhs.rs)

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

        mycop <- cop(f = f.obj, max = FALSE, lb = lbcon1, ub = ubcon1, lc = lincon)
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
            if ((n_attempts == 1) && is.matrix(X)) {
              Xini <- X[i, ]
              names(Xini) <- namevar
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
            t_input <- res[(ndr + 1) : (ndr + ni)]
            names(t_input) <- inputnames
            t_output <- res[(ndr + ni + 1) : (ndr + ni + no)]
            names(t_output) <- outputnames
            sigma_input <- res[(ndr + ni + no + 1):(ndr + ni + no + ni)]
            names(sigma_input) <- inputnames
            sigma_output <- res[(ndr + 2 * ni + no + 1):nvar]
            names(sigma_output) <- outputnames

            delta_num <- 1 + sum(t_input / input[, ii]) / (ni - nnci)
            delta_den <- 1 - sum(t_output / output[, ii]) / (no - nnco)
            delta <- delta_num / delta_den

            if (compute_target) {
              target_input <- as.vector(inputref %*% lambda)
              names(target_input) <- inputnames
              target_output <- as.vector(outputref %*% lambda)
              names(target_output) <- outputnames

              project_input <- input[, ii] + t_input
              names(project_input) <- inputnames
              project_output <- output[, ii] - t_output
              names(project_output) <- outputnames

              slack_input <- project_input - target_input
              names(slack_input) <- inputnames
              slack_output <- target_output - project_output
              names(slack_output) <- outputnames
            }

            objval <- weight_slack_i[, i] %*% t_input + weight_slack_o[, i] %*% t_output

          } else {

            delta <- NA
            objval <- NA
            lambda <- NA
            t_input <- NA
            t_output <- NA
            sigma_input <- NA
            sigma_output <- NA
            if (compute_target) {
              target_input <- NA
              target_output <- NA
              project_input <- NA
              project_output <- NA
              slack_input <- NA
              slack_output <- NA
            }

          }

          list(delta = delta,
               objval = objval,
               lambda = lambda,
               t_input = t_input, t_output = t_output,
               slack_input = slack_input, slack_output = slack_output,
               project_input = project_input, project_output = project_output,
               target_input = target_input, target_output = target_output,
               sigma_input = sigma_input, sigma_output = sigma_output)

        }
      }

    } else {

      DMU <- vector(mode = "list", length = nde)

      for (i in 1:nde) {

        ii <- dmu_eval[i]

        w0i <- which(weight_slack_i[, i] == 0)
        nw0i <- length(w0i)
        w0o <- which(weight_slack_o[, i] == 0)
        nw0o <- length(w0o)

        # Objective function coefficients
        f.obj <- linfun(a = c(rep(0, ndr), weight_slack_i[, i],
                              weight_slack_o[, i], rep(0, ni + no)),
                        id = namevar)

        # Linear constraints matrix
        f.con.se <- rep(0, ndr)
        f.con.se[dmu_ref == ii] <- 1
        f.con.se <- c(f.con.se, rep(0, 2 * (ni + no)))
        f.con.w0 <- matrix(0, nrow = (nw0i + nw0o), ncol = nvar)
        f.con.w0[, ndr + c(w0i, ni + w0o)] <- diag(nw0i + nw0o)
        f.con <- rbind(f.con.1, f.con.2, f.con.se, f.con.w0, f.con.rs)
        rownames(f.con) <- paste("lc", 1:nrow(f.con), sep = "") # to prevent names errors in lincon

        # Linear directions vector
        f.dir <- c(rep("<=", ni), rep(">=", no), rep("==", 1 + nw0i + nw0o), f.dir.rs)
        f.dir[nc_inputs] <- "=="
        f.dir[ni + nc_outputs] <- "=="

        # Linear right hand side vector
        f.rhs <- c(input[, ii], output[, ii], rep(0, 1 + nw0i + nw0o), f.rhs.rs)

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

        mycop <- cop(f = f.obj, max = FALSE, lb = lbcon1, ub = ubcon1, lc = lincon)
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
            if ((n_attempts == 1) && is.matrix(X)) {
              Xini <- X[i, ]
              names(Xini) <- namevar
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
            t_input <- res[(ndr + 1) : (ndr + ni)]
            names(t_input) <- inputnames
            t_output <- res[(ndr + ni + 1) : (ndr + ni + no)]
            names(t_output) <- outputnames
            sigma_input <- res[(ndr + ni + no + 1):(ndr + ni + no + ni)]
            names(sigma_input) <- inputnames
            sigma_output <- res[(ndr + 2 * ni + no + 1):nvar]
            names(sigma_output) <- outputnames

            delta_num <- 1 + sum(t_input / input[, ii]) / (ni - nnci)
            delta_den <- 1 - sum(t_output / output[, ii]) / (no - nnco)
            delta <- delta_num / delta_den

            if (compute_target) {
              target_input <- as.vector(inputref %*% lambda)
              names(target_input) <- inputnames
              target_output <- as.vector(outputref %*% lambda)
              names(target_output) <- outputnames

              project_input <- input[, ii] + t_input
              names(project_input) <- inputnames
              project_output <- output[, ii] - t_output
              names(project_output) <- outputnames

              slack_input <- project_input - target_input
              names(slack_input) <- inputnames
              slack_output <- target_output - project_output
              names(slack_output) <- outputnames
            }

            objval <- weight_slack_i[, i] %*% t_input + weight_slack_o[, i] %*% t_output

          } else {

            delta <- NA
            objval <- NA
            lambda <- NA
            t_input <- NA
            t_output <- NA
            sigma_input <- NA
            sigma_output <- NA
            if (compute_target) {
              target_input <- NA
              target_output <- NA
              project_input <- NA
              project_output <- NA
              slack_input <- NA
              slack_output <- NA
            }

          }

          DMU[[i]] <- list(delta = delta,
                           objval = objval,
                           lambda = lambda,
                           t_input = t_input, t_output = t_output,
                           slack_input = slack_input, slack_output = slack_output,
                           project_input = project_input, project_output = project_output,
                           target_input = target_input, target_output = target_output,
                           sigma_input = sigma_input, sigma_output = sigma_output)

        }
      }
    }

    names(DMU) <- dmunames[dmu_eval]

    deaOutput <- list(modelname = "stoch_addsupereff",
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
                      orientation = NA)

    return(structure(deaOutput, class = "dea_stoch"))

  }
