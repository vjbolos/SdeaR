#' @title Chance Constrained Additive P-model.
#'
#' @description It solves the "max" version of the "almost 100% confidence"
#' chance constrained problem from Cooper et al. (1998), under constant and
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
#' Given \eqn{0<\alpha <1} and \eqn{\varepsilon} a positive non-Archimedean infinitesimal,
#' the program for \eqn{\text{DMU}_o} with constant returns to scale is given by
#' \deqn{\max \limits_{\bm{\lambda}}\quad P\left\{ \mathbf{w}^-\cdot \left(
#' \tilde{\mathbf{x}}_o-\tilde{X}\bm{\lambda}\right) +\mathbf{w}^+\cdot \left(
#' \tilde{Y}\bm{\lambda}-\tilde{\mathbf{y}}_o\right) > 0\right\} }
#' \deqn{\text{s.t.}\quad P\left\{ \left( \tilde{\mathbf{x}}_o-\tilde{X}
#' \bm{\lambda}\right) _i> 0\right\} \geq 1-\varepsilon,\quad i=1,\ldots ,m,}
#' \deqn{P\left\{ \left( \tilde{Y}\bm{\lambda}-\tilde{\mathbf{y}}_o\right) _r>
#' 0\right\} \geq 1-\varepsilon,\quad r=1,\ldots ,s,}
#' \deqn{\bm{\lambda}\geq \mathbf{0},}
#' where \eqn{\bm{\lambda}=(\lambda_1,\ldots,\lambda_n)^\top}, \eqn{\tilde{\mathbf{x}}_o
#' =(\tilde{x}_{1o},\ldots,\tilde{x}_{mo})^\top} and \eqn{\tilde{\mathbf{y}}_o=
#' (\tilde{y}_{1o},\ldots,\tilde{y}_{so})^\top} are column vectors. Moreover,
#' \eqn{\mathbf{w}^-,\mathbf{w}^+} are positive row vectors with the weights
#' for the slacks.
#' Different returns to scale can be easily considered by adding the corresponding
#' constraints: \eqn{\mathbf{e}\bm{\lambda}=1} (VRS), \eqn{0\leq \mathbf{e}\bm{\lambda}
#' \leq 1} (NIRS), \eqn{\mathbf{e}\bm{\lambda}\geq 1} (NDRS) or \eqn{L\leq \mathbf{e}
#' \bm{\lambda}\leq U} (GRS), with \eqn{0\leq L\leq 1} and \eqn{U\geq 1}, where
#' \eqn{\mathbf{e}=(1,\ldots ,1)} is a row vector.
#'
#' The deterministic equivalent for a multivariate normal distribution of inputs/outputs
#' is given by
#' \deqn{\max \limits_{\bm{\lambda}}\quad \mathbf{w}^-\cdot \left( \mathbf{x}_o-
#' X\bm{\lambda}\right)+\mathbf{w}^+\cdot \left( Y\bm{\lambda}-\mathbf{y}_o\right)
#' -\Phi ^{-1}(\alpha)\sigma(\bm{\lambda}) }
#' \deqn{\text{s.t.}\quad \mathbf{x}_o-X\bm{\lambda}+\Phi ^{-1}(\varepsilon)
#' \bm{\sigma} ^-(\bm{\lambda})\geq \mathbf{0},}
#' \deqn{Y\bm{\lambda}-\mathbf{y}_o+\Phi ^{-1}(\varepsilon)\bm{\sigma}^+
#' (\bm{\lambda})\geq \mathbf{0},}
#' \deqn{\bm{\lambda}\geq \mathbf{0},}
#' where \eqn{\Phi } is the standard normal distribution, and
#' \deqn{\left( \sigma \left( \bm{\lambda}\right)\right) ^2=\bm{\lambda}^\top
#' \left[ \displaystyle \sum _{i,k=1}^m\Delta ^{II}_{ik}+\sum _{r,p=1}^s
#' \Delta ^{OO}_{rp}-2\sum _{i=1}^m\sum _{r=1}^s \Delta ^{IO}_{ir}\right] \bm{\lambda} }
#' \deqn{\displaystyle +2\sum _{j=1}^n\lambda _j\left[ \sum _{i=1}^m\sum _{r=1}^s\left(
#' \mathrm{Cov}(\tilde{x}_{io},\tilde{y}_{rj})+\mathrm{Cov}(\tilde{x}_{ij},
#' \tilde{y}_{ro})\right) \right. }
#' \deqn{\displaystyle \left.-\sum _{i,k=1}^m\mathrm{Cov}(\tilde{x}_{io},\tilde{x}_{kj})
#' -\sum _{r,p=1}^s\mathrm{Cov}(\tilde{y}_{ro},\tilde{y}_{pj})\right] }
#' \deqn{\displaystyle +\sum _{i,k=1}^m\mathrm{Cov}(\tilde{x}_{io},\tilde{x}_{ko})
#' +\sum _{r,p=1}^s \mathrm{Cov}(\tilde{y}_{ro},\tilde{y}_{po})-2\sum _{i=1}^m
#' \sum _{r=1}^s \mathrm{Cov}(\tilde{x}_{io}, \tilde{y}_{ro}),}
#' \deqn{\displaystyle \left( \sigma ^-_i\left( \bm{\lambda}\right)\right) ^2=\bm{\lambda}^\top
#' \Delta ^{II}_{ii}\bm{\lambda}-2\sum _{j=1}^n\lambda _j\mathrm{Cov}(\tilde{x}_{ij},
#' \tilde{x}_{io})+\mathrm{Var}(\tilde{x}_{io}),\quad i=1,\ldots ,m,}
#' \deqn{\displaystyle \left( \sigma ^+_r\left( \bm{\lambda}\right)\right) ^2=\bm{\lambda}^\top
#' \Delta ^{OO}_{rr}\bm{\lambda}-2\sum _{j=1}^n\lambda _j\mathrm{Cov}(\tilde{y}_{rj},
#' \tilde{y}_{ro})+\mathrm{Var}(\tilde{y}_{ro}),\quad r=1,\ldots ,s.}
#' such that
#' \deqn{\left( \Delta ^{II}_{ik}\right) _{jq}=\mathrm{Cov}(\tilde{x}_{ij}, \tilde{x}_{kq}),}
#' \deqn{\left( \Delta ^{OO}_{rp}\right) _{jq}=\mathrm{Cov}(\tilde{y}_{rj}, \tilde{y}_{pq}),}
#' \deqn{\left( \Delta ^{IO}_{ir}\right) _{jq}=\left( \Delta ^{OI}_{ri}\right) _{qj}
#' =\mathrm{Cov}(\tilde{x}_{ij}, \tilde{y}_{rq}), }
#' with \eqn{1\leq j,q\leq n}.
#'
#' @note The model in this function is the "max" version of the model in Cooper
#' et al. (1998), in the sense that maximizes the sum of positive slacks like
#' the conventional additive model in Cooper et al. (1985). Hence, a DMU is
#' \eqn{\alpha}-stochastically efficient if and only if the optimal objective value of
#' the problem, (\code{objval}), is zero (or less than zero).
#'
#' @usage modelstoch_additive_p(datadea,
#'                alpha = 0.05,
#'                epsilon = 0.05,
#'                dmu_eval = NULL,
#'                dmu_ref = NULL,
#'                orientation = NULL,
#'                weight_slack_i = 1,
#'                weight_slack_o = 1,
#'                rts = c("crs", "vrs", "nirs", "ndrs", "grs"),
#'                L = 1,
#'                U = 1,
#'                solver = c("alabama", "cccp", "cccp2", "slsqp"),
#'                give_X = TRUE,
#'                n_attempts_max = 5,
#'                returnqp = FALSE,
#'                ...)
#'
#' @param datadea The data of class \code{deadata_stoch}, including \code{n} DMUs,
#' and the expected values of \code{m} inputs and \code{s} outputs.
#' @param alpha A value for parameter alpha.
#' @param epsilon A value for parameter epsilon.
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
#' from package \code{optiSolve}.
#' @param give_X Logical. If it is \code{TRUE}, it uses an initial vector (given by
#' the evaluated DMU) for the solver, except for "cccp". If it is \code{FALSE}, the initial vector is given
#' internally by the solver and it is usually randomly generated.
#' @param n_attempts_max A value with the maximum number of attempts if the solver
#' does not converge. Each attempt uses a different initial vector.
#' @param returnqp Logical. If it is \code{TRUE}, it returns the quadratic problems
#' (objective function and constraints).
#' @param ... Other parameters, like the initial vector \code{X}, to be passed to the solver.
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
#' Cooper, W.W.; Huang, Z.; Lelas, V.; Li, S.X.; Olesen, O.B. (1998) "Chance
#' Constrained Programming Formulations for Stochastic Characterizations of
#' Efficiency and Dominance in DEA", Journal of Productivity Analysis, 9, 53–79.
#' \doi{10.1023/A:1018320430249}
#'
#' Cooper, W.W.; Park, K.S.; Pastor, J.T. (1999). "RAM: A Range Adjusted Measure
#' of Inefficiencies for Use with Additive Models, and Relations to Other Models
#' and Measures in DEA". Journal of Productivity Analysis, 11, p. 5-42.
#' \doi{10.1023/A:1007701304281}
#'
#'
#'
#' @import optiSolve stats
#'
#' @export

modelstoch_additive_p <-
  function(datadea,
           alpha = 0.05,
           epsilon = 0.05,
           dmu_eval = NULL,
           dmu_ref = NULL,
           orientation = NULL,
           weight_slack_i = 1,
           weight_slack_o = 1,
           rts = c("crs", "vrs", "nirs", "ndrs", "grs"),
           L = 1,
           U = 1,
           solver = c("alabama", "cccp", "cccp2", "slsqp"),
           give_X = TRUE,
           n_attempts_max = 5,
           returnqp = FALSE,
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

    # Checking epsilon
    if (epsilon <= 0) {
      stop("Parameter epsilon must be greater than 0.")
    }
    qe <- qnorm(epsilon)

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

    cov_II <- datadea$cov_II
    cov_OO <- datadea$cov_OO
    cov_IO <- datadea$cov_IO
    covX <- array(0, dim = c(ni, nd, nd))
    for (i in 1:ni) {
      covX[i, , ] <- cov_II[i, i, , ]
    }
    covY <- array(0, dim = c(no, nd, nd))
    for (i in 1:no) {
      covY[i, , ] <- cov_OO[i, i, , ]
    }
    cov_IIref <- array(cov_II[, , dmu_ref, dmu_ref], dim = c(ni, ni, ndr, ndr))
    cov_OOref <- array(cov_OO[, , dmu_ref, dmu_ref], dim = c(no, no, ndr, ndr))
    cov_IOref <- array(cov_IO[, , dmu_ref, dmu_ref], dim = c(ni, no, ndr, ndr))
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
                 paste("sigma_I", 1:ni, sep = "_"),
                 paste("sigma_O", 1:no, sep = "_"),
                 "sigma")
    nvar <- ndr + ni + no + 1

    DMU <- vector(mode = "list", length = nde)
    names(DMU) <- dmunames[dmu_eval]

    ###########################

    # Lower and upper bounds constraints
    lbcon1 <- lbcon(val = 0, id = namevar)
    ubcon1 <- NULL

    if (rts == "crs") {
      f.con.rs <- NULL
      f.dir.rs <- NULL
      f.rhs.rs <- NULL
    } else {
      f.con.rs <- cbind(matrix(1, nrow = 1, ncol = ndr), matrix(0, nrow = 1, ncol = ni + no + 1))
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
    f.con.1 <- cbind(-inputref, diag(qe, ni), matrix(0, nrow = ni, ncol = no + 1))
    f.con.2 <- cbind(outputref, matrix(0, nrow = no, ncol = ni), diag(qe, no), matrix(0, nrow = no, ncol = 1))
    f.con.1[nc_inputs, ndr + nc_inputs] <- 0
    f.con.2[nc_outputs, ndr + ni + nc_outputs] <- 0
    f.con <- rbind(f.con.1, f.con.2, f.con.rs)
    rownames(f.con) <- paste("lc", 1:nrow(f.con), sep = "") # to prevent names errors in lincon

    # Linear directions vector
    f.dir <- c(rep(">=", ni + no), f.dir.rs)
    f.dir[c(nc_inputs, ni + nc_outputs)] <- "=="

    # Quadratic matrix for sigma constraint
    Q0 <- matrix(0, nrow = nvar, ncol = nvar)
    for (j1 in 1:ndr) {
      for (j2 in 1:ndr) {
        Q0[j1, j2] <- sum(cov_IIref[, , j1, j2]) + sum(cov_OOref[, , j1, j2]) - 2 * sum(cov_IOref[, , j1, j2])
      }
    }
    Q0[nvar, nvar] <- -1

    # Quadratic matrices
    Q1 <- array(0, dim = c(ni, nvar, nvar))
    Q2 <- array(0, dim = c(no, nvar, nvar))
    Q1[, 1:ndr, 1:ndr] <- covXref
    Q2[, 1:ndr, 1:ndr] <- covYref
    for (j in 1: ni) {
      Q1[j, ndr + j, ndr + j] <- -1
    }
    for (j in 1:no) {
      Q2[j, ndr + ni + j, ndr + ni + j] <- -1
    }

    # Linear part of quadratic constraints
    a1 <- matrix(0, nrow = ni, ncol = nvar)
    a2 <- matrix(0, nrow = no, ncol = nvar)

    qclist <- vector(mode = "list", length = 2 * (1 + ni + no))

    for (i in 1:nde) {

      ii <- dmu_eval[i]

      # Objective function coefficients
      f.obj <- linfun(a = c(colSums(outputref * weight_slack_o[, i]) -
                              colSums(inputref * weight_slack_i[, i]),
                            rep(0, ni + no), -qa),
                      d = input[, ii] %*% weight_slack_i[, i] -
                        output[, ii] %*% weight_slack_o[, i],
                      id = namevar)

      # Linear right hand side vector
      f.rhs <- c(-input[, ii], output[, ii], f.rhs.rs)

      # Linear constraints
      lincon1 <- lincon(A = f.con, dir = f.dir, val = f.rhs, id = namevar)

      # Quadratic sigma constraint
      a0 <- rep(0, nvar)
      for (j1 in 1:ndr) {
        a0[j1] <- 2 * (sum(cov_IOref[, , j1, ii]) + sum(cov_IOref[ , , ii, j1]) -
                         sum(cov_IIref[, , j1, ii]) - sum(cov_OOref[, , j1, ii]))
      }
      val0 <- 2 * sum(cov_IOref[, , ii, ii]) - sum(cov_IIref[, , ii, ii]) - sum(cov_OOref[, , ii, ii])
      quadcon0a <- quadcon(Q = Q0, a = a0, val = val0, id = namevar)
      quadcon0b <- quadcon(Q = -Q0, a = -a0, val = -val0, id = namevar)

      # Quadratic constraints
      a1[, 1:ndr] <- -2 * covX[, ii, dmu_ref]
      a2[, 1:ndr] <- -2 * covY[, ii, dmu_ref]
      val1 <- -covX[, ii, ii]
      val2 <- -covY[, ii, ii]
      quadcon1a <- vector(mode = "list", length = ni)
      quadcon1b <- quadcon1a
      quadcon2a <- vector(mode = "list", length = no)
      quadcon2b <- quadcon2a
      for (j in 1:ni) {
        quadcon1a[[j]] <- quadcon(Q = Q1[j, , ], a = a1[j, ], val = val1[j], id = namevar)
        quadcon1b[[j]] <- quadcon(Q = -Q1[j, , ], a = -a1[j, ], val = -val1[j], id = namevar)
      }
      for (j in 1:no) {
        quadcon2a[[j]] <- quadcon(Q = Q2[j, , ], a = a2[j, ], val = val2[j], id = namevar)
        quadcon2b[[j]] <- quadcon(Q = -Q2[j, , ], a = -a2[j, ], val = -val2[j], id = namevar)
      }

      mycop <- cop(f = f.obj, max = TRUE, lb = lbcon1, ub = ubcon1, lc = lincon1)
      qclist <- c(list(quadcon0a, quadcon0b), quadcon1a, quadcon1b, quadcon2a, quadcon2b)
      names(qclist) <- paste("qc", 1:length(qclist), sep = "")
      mycop$qc <- qclist

      if (returnqp) {

        DMU[[i]] <- mycop

      } else {

        n_attempts <- 1

        while (n_attempts <= n_attempts_max) {

          # Initial vector
          if ((n_attempts == 1) && give_X && (ii %in% dmu_ref)) {
            Xini <- rep(0, ndr + ni + no + 1)
            Xini[which(dmu_ref == ii)] <- 1
            names(Xini) <- namevar
          } else {
            Xini <- NULL
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
          sigma_input <- res[(ndr + 1):(ndr + ni)]
          names(sigma_input) <- inputnames
          sigma_output <- res[(ndr + ni + 1):(ndr + ni + no)]
          names(sigma_output) <- outputnames
          sigma <- res[nvar]

          target_input <- as.vector(inputref %*% lambda)
          target_output <- as.vector(outputref %*% lambda)
          names(target_input) <- inputnames
          names(target_output) <- outputnames

          slack_input <- input[, ii] - target_input
          slack_output <- target_output - output[, ii]
          names(slack_input) <- inputnames
          names(slack_output) <- outputnames

          objval <- weight_slack_i[, i] %*% slack_input + weight_slack_o[, i] %*% slack_output -
            qa * sigma

        } else {

          objval <- NA
          lambda <- NA
          sigma_input <- NA
          sigma_output <- NA
          sigma <- NA
          slack_input <- NA
          slack_output <- NA
          target_input <- NA
          target_output <- NA

        }

        DMU[[i]] <- list(objval = objval,
                         lambda = lambda,
                         slack_input = slack_input, slack_output = slack_output,
                         target_input = target_input, target_output = target_output,
                         sigma_input = sigma_input, sigma_output = sigma_output,
                         sigma = sigma)

      }

    }

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

    deaOutput <- list(modelname = "stoch_additive_p",
                      rts = rts,
                      L = L,
                      U = U,
                      DMU = DMU,
                      data = datadea,
                      alpha = alpha,
                      epsilon = epsilon,
                      dmu_eval = dmu_eval,
                      dmu_ref = dmu_ref,
                      weight_slack_i = weight_slack_i,
                      weight_slack_o = weight_slack_o,
                      orientation = orientation)

    return(structure(deaOutput, class = "dea_stoch"))

  }
