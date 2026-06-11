#' @title Chance Constrained Directional Models with Deterministic Directions
#'
#' @description It solves chance constrained directional models with deterministic directions,
#' under constant, variable, non-increasing, non-decreasing or generalized returns to
#' scale. Inputs and outputs must follow a multivariate normal distribution.
#' By default, models are solved in a two-stage process (slacks are maximized).
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
#' Given \eqn{0<\alpha <1}, the first stage program for \eqn{\text{DMU}_o} with
#' constant returns to scale is given by
#' \deqn{\max \limits_{\beta, \bm{\lambda}}\quad \beta}
#' \deqn{\text{s.t.} \quad P\left\{ \left( \tilde{\mathbf{x}}_o-\beta \mathbf{g}^--
#' \tilde{X}\bm{\lambda}\right) _i\geq 0\right\} \geq 1-\alpha,\quad i=1,\ldots ,m,}
#' \deqn{P\left\{ \left( \tilde{\mathbf{y}}_{o}+\beta \mathbf{g}^+-\tilde{Y}\bm{\lambda}
#' \right) _r\leq 0\right\} \geq 1-\alpha,\quad r=1,\ldots ,s,}
#' \deqn{\bm{\lambda}\geq \mathbf{0},}
#' where \eqn{\bm{\lambda}=(\lambda_1,\ldots,\lambda_n)^\top}, \eqn{\tilde{\mathbf{x}}_o
#' =(\tilde{x}_{1o},\ldots,\tilde{x}_{mo})^\top} and \eqn{\tilde{\mathbf{y}}_o=
#' (\tilde{y}_{1o},\ldots,\tilde{y}_{so})^\top} are column vectors, and
#' \eqn{\mathbf{g}=(-\mathbf{g}^-,\mathbf{g}^+)\neq \mathbf{0}} is a preassigned direction
#' (with \eqn{\mathbf{g}^-\in \mathbb{R}^m} and \eqn{\mathbf{g}^+\in\mathbb{R}^s}
#' non-negative column vectors).
#' Different returns to scale can be easily considered by adding the corresponding
#' constraints: \eqn{\mathbf{e}\bm{\lambda}=1} (VRS), \eqn{0\leq \mathbf{e}\bm{\lambda}
#' \leq 1} (NIRS), \eqn{\mathbf{e}\bm{\lambda}\geq 1} (NDRS) or \eqn{L\leq \mathbf{e}
#' \bm{\lambda}\leq U} (GRS), with \eqn{0\leq L\leq 1} and \eqn{U\geq 1}, where
#' \eqn{\mathbf{e}=(1,\ldots ,1)} is a row vector.
#'
#' The corresponding second stage program is given by
#' \deqn{\max \limits_{\bm{\lambda},\mathbf{s}^-,\mathbf{s}^+}\quad \mathbf{w}^-
#' \mathbf{s}^-+\mathbf{w}^+\mathbf{s}^+}
#' \deqn{\text{s.t.} \quad P\left\{ \left( \tilde{\mathbf{x}}_o-\beta ^*\mathbf{g}^--
#' \tilde{X}\bm{\lambda}-\mathbf{s}^-\right) _i\geq 0\right\} \geq 1-\alpha,\quad i=1,\ldots ,m,}
#' \deqn{P\left\{ \left( \tilde{\mathbf{y}}_{o}+\beta ^*\mathbf{g}^+-\tilde{Y}\bm{\lambda}
#' +\mathbf{s}^+\right) _r\leq 0\right\} \geq 1-\alpha,\quad r=1,\ldots ,s,}
#' \deqn{\bm{\lambda}\geq \mathbf{0},\,\, \mathbf{s}^-\geq \mathbf{0},\,\,
#' \mathbf{s}^+\geq \mathbf{0},}
#' where \eqn{\beta ^*} is the optimal objective function of the first stage
#' program, \eqn{\mathbf{s}^-,\mathbf{s}^+} are column vectors with the slacks,
#' and \eqn{\mathbf{w}^-,\mathbf{w}^+} are positive row vectors with the weights
#' for the slacks.
#'
#' The deterministic equivalents for a multivariate normal distribution of inputs/outputs
#' are given by
#' \deqn{\max \limits_{\beta, \bm{\lambda}} \quad \beta}
#' \deqn{\text{s.t.}\quad \beta \mathbf{g}^-+X\bm{\lambda}-\Phi ^{-1}
#' (\alpha)\bm{\sigma} ^-(\bm{\lambda}) \leq \mathbf{x}_o,}
#' \deqn{-\beta \mathbf{g}^++Y\bm{\lambda}+\Phi ^{-1}(\alpha)
#' \bm{\sigma} ^+(\bm{\lambda}) \geq \mathbf{y}_{o},}
#' \deqn{\bm{\lambda}\geq \mathbf{0},}
#' and for the second stage,
#' \deqn{\max \limits_{\bm{\lambda},\mathbf{s}^-,\mathbf{s}^+} \quad \mathbf{w}^-
#' \mathbf{s}^-+\mathbf{w}^+\mathbf{s}^+}
#' \deqn{\text{s.t.}\quad \beta ^*\mathbf{g}^-+X\bm{\lambda}+\mathbf{s}^--\Phi ^{-1}
#' (\alpha)\bm{\sigma} ^-(\bm{\lambda}) = \mathbf{x}_o,}
#' \deqn{-\beta ^*\mathbf{g}^++Y\bm{\lambda}-\mathbf{s}^++\Phi ^{-1}(\alpha)
#' \bm{\sigma} ^+(\bm{\lambda}) = \mathbf{y}_{o},}
#' \deqn{\bm{\lambda}\geq \mathbf{0},\,\, \mathbf{s}^-\geq \mathbf{0},\,\,
#' \mathbf{s}^+\geq \mathbf{0},}
#' where \eqn{\Phi } is the standard normal distribution, and
#' \deqn{\displaystyle \left( \sigma ^-_i\left( \bm{\lambda}\right)\right) ^2 =
#' \sum _{j,q=1}^n\lambda _j\lambda _q\mathrm{Cov}(\tilde{x}_{ij},\tilde{x}_{iq})-
#' 2\sum _{j=1}^n\lambda _j\mathrm{Cov}(\tilde{x}_{ij},\tilde{x}_{io})}
#' \deqn{+\mathrm{Var}(\tilde{x}_{io}),\quad i=1,\ldots ,m,}
#' \deqn{\displaystyle \left( \sigma ^+_r\left( \bm{\lambda}\right)\right) ^2 =
#' \sum _{j,q=1}^n\lambda _j\lambda _q\mathrm{Cov}(\tilde{y}_{rj},\tilde{y}_{rq})-
#' 2\sum _{j=1}^n\lambda _j\mathrm{Cov}(\tilde{y}_{rj},\tilde{y}_{ro})}
#' \deqn{+\mathrm{Var}(\tilde{y}_{ro}),\quad r=1,\ldots ,s.}
#'
#' @usage
#' modelstoch_dir_dd(datadea,
#'                   alpha = 0.05,
#'                   dmu_eval = NULL,
#'                   dmu_ref = NULL,
#'                   dir_input = NULL,
#'                   dir_output = NULL,
#'                   rts = c("crs", "vrs", "nirs", "ndrs", "grs"),
#'                   L = 1,
#'                   U = 1,
#'                   solver = c("alabama", "cccp", "cccp2", "slsqp"),
#'                   give_X = TRUE,
#'                   n_attempts_max = 5,
#'                   maxslack = FALSE,
#'                   weight_slack_i = 1,
#'                   weight_slack_o = 1,
#'                   compute_target = TRUE,
#'                   returnqp = FALSE,
#'                   silent_ud = FALSE,
#'                   parallel = FALSE,
#'                   ...)
#'
#' @param datadea The data of class \code{deadata_stoch}, including \code{n} DMUs,
#' and the expected values of \code{m} inputs and \code{s} outputs.
#' @param alpha A value for parameter alpha.
#' @param dmu_eval A numeric vector containing which DMUs have to be evaluated.
#' If \code{NULL} (default), all DMUs are considered.
#' @param dmu_ref A numeric vector containing which DMUs are the evaluation
#' reference set.
#' If \code{NULL} (default), all DMUs are considered.
#' @param dir_input A value, vector of length \code{m}, or matrix \code{m} x \code{ne}
#' (where \code{ne} is the length of \code{dmu_eval}) with the input directions.
#' If \code{dir_input} == input matrix (of DMUS in \code{dmu_eval}) and
#' \code{dir_output} == 0, it is equivalent to input oriented (\code{beta} = 1 -
#' \code{efficiency}). If \code{dir_input} is omitted, input matrix (of DMUS in
#' \code{dmu_eval}) is assigned.
#' @param dir_output A value, vector of length \code{s}, or matrix \code{s} x \code{ne}
#' (where \code{ne} is the length of \code{dmu_eval}) with the output directions.
#' If \code{dir_input} == 0 and \code{dir_output} == output matrix (of DMUS in
#' \code{dmu_eval}), it is equivalent to output oriented (\code{beta} = \code{efficiency} - 1).
#' If \code{dir_output} is omitted, output matrix (of DMUS in \code{dmu_eval}) is assigned.
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
#' @param maxslack Logical. If it is \code{TRUE}, it computes the max slack solution.
#' @param weight_slack_i A value, vector of length \code{m}, or matrix \code{m} x \code{ne}
#' (where \code{ne} is the length of \code{dmu_eval}) with the weights of the input slacks
#' for the max slack solution.
#' @param weight_slack_o A value, vector of length \code{s}, or matrix \code{s} x \code{ne}
#' (where \code{ne} is the length of \code{dmu_eval}) with the weights of the output
#' slacks for the max slack solution.
#' @param compute_target Logical. If it is \code{TRUE}, it computes targets of the
#' max slack solution.
#' @param returnqp Logical. If it is \code{TRUE}, it returns the quadratic problems
#' (objective function and constraints) of stage 1.
#' @param silent_ud Logical, to avoid warnings related with undesirable variables.
#' @param parallel Logical, if \code{TRUE}, the DMUs are computed in parallel
#' (default \code{FALSE}).
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
#' Bolós, V.J.; Benítez, R.; Coll-Serrano, V. (2024). “Chance constrained directional
#' models in stochastic data envelopment analysis", Operations Research Perspectives, 12, 100307..
#' \doi{10.1016/j.orp.2024.100307}
#'
#' @examples
#' \donttest{
#' # Example 1.
#'
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
#' Collstochdirdd <- modelstoch_dir_dd(data_stoch, dmu_eval = 6)
#' efficiencies(Collstochdirdd)
#'
#' }
#'
#' @import optiSolve stats
#' @importFrom foreach foreach %dopar%
#'
#' @export

modelstoch_dir_dd <-
  function(datadea,
           alpha = 0.05,
           dmu_eval = NULL,
           dmu_ref = NULL,
           dir_input = NULL,
           dir_output = NULL,
           rts = c("crs", "vrs", "nirs", "ndrs", "grs"),
           L = 1,
           U = 1,
           solver = c("alabama", "cccp", "cccp2", "slsqp"),
           give_X = TRUE,
           n_attempts_max = 5,
           maxslack = FALSE,
           weight_slack_i = 1,
           weight_slack_o = 1,
           compute_target = TRUE,
           returnqp = FALSE,
           silent_ud = FALSE,
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
    nc_inputs <- datadea$nc_inputs
    nc_outputs <- datadea$nc_outputs
    nd_inputs <- datadea$nd_inputs
    nd_outputs <- datadea$nd_outputs
    #ud_inputs <- datadea$ud_inputs
    #ud_outputs <- datadea$ud_outputs
    inputnames <- rownames(input)
    outputnames <- rownames(output)
    ni <- nrow(input)
    covX <- array(0, dim = c(ni, nd, nd))
    for (i in 1:ni) {
      covX[i, , ] <- datadea$cov_II[i, i, , ]
    }
    no <- nrow(output)
    covY <- array(0, dim = c(no, nd, nd))
    for (i in 1:no) {
      covY[i, , ] <- datadea$cov_OO[i, i, , ]
    }
    namevar1 <- c("beta",
                  paste("lambda", 1:ndr, sep = "_"),
                  paste("sigma_I", 1:ni, sep = "_"),
                  paste("sigma_O", 1:no, sep = "_"))
    namevar2 <- c(paste("lambda", 1:ndr, sep = "_"),
                  paste("slack_I", 1:ni, sep = "_"),
                  paste("slack_O", 1:no, sep = "_"),
                  paste("sigma_I", 1:ni, sep = "_"),
                  paste("sigma_O", 1:no, sep = "_"))

    if (is.null(dir_input)) {
      dir_input <- matrix(input[, dmu_eval], nrow = ni) # input of DMUs in dmu_eval
    } else {
      if (is.matrix(dir_input)) {
        if ((nrow(dir_input) != ni) || (ncol(dir_input) != nde)) {
          stop("Invalid input direction matrix (number of inputs x number of evaluated DMUs).")
        }
      } else if ((length(dir_input) == 1) || (length(dir_input) == ni)) {
        dir_input <- matrix(dir_input, nrow = ni, ncol = nde)
      } else {
        stop("Invalid input direction vector.")
      }
    }
    rownames(dir_input) <- inputnames
    colnames(dir_input) <- dmunames[dmu_eval]

    if (is.null(dir_output)) {
      dir_output <- matrix(output[, dmu_eval], nrow = no) # output of DMUs in dmu_eval
    } else {
      if (is.matrix(dir_output)) {
        if ((nrow(dir_output) != no) || (ncol(dir_output) != nde)) {
          stop("Invalid output direction matrix (number of outputs x number of evaluated DMUs).")
        }
      } else if ((length(dir_output) == 1) || (length(dir_output) == no)) {
        dir_output <- matrix(dir_output, nrow = no, ncol = nde)
      } else {
        stop("Invalid output direction vector.")
      }
    }
    rownames(dir_output) <- outputnames
    colnames(dir_output) <- dmunames[dmu_eval]

    inputref <- matrix(input[, dmu_ref], nrow = ni)
    outputref <- matrix(output[, dmu_ref], nrow = no)
    covXref <- array(covX[, dmu_ref, dmu_ref], dim = c(ni, ndr, ndr))
    covYref <- array(covY[, dmu_ref, dmu_ref], dim = c(no, ndr, ndr))

    ncd_inputs <- sort(unique(c(nc_inputs, nd_inputs)))
    ncd_outputs <- sort(unique(c(nc_outputs, nd_outputs)))
    nnci <- length(nc_inputs) # number of non-controllable inputs
    nnco <- length(nc_outputs) # number of non-controllable outputs

    target_input <- NULL
    target_output <- NULL
    orientation_param <- NULL

    ###########################

    # Objective function coefficients stage 1
    f.obj <- linfun(a = c(1, rep(0, ndr + ni + no)), id = namevar1)

    # Lower and upper bounds constraints
    lbcon1 <- lbcon(val = 0, id = namevar1)
    lbcon2 <- lbcon(val = 0, id = namevar2)
    ubcon12 <- NULL # 1st and 2nd stage

    if (rts == "crs") {
      f.con.rs <- NULL # Stage 1
      f.con2.rs <- NULL # Stage 2
      f.dir.rs <- NULL
      f.rhs.rs <- NULL
    } else {
      f.con.rs <- cbind(0, matrix(1, nrow = 1, ncol = ndr), matrix(0, nrow = 1, ncol = ni + no))
      f.con2.rs <- cbind(matrix(1, nrow = 1, ncol = ndr), matrix(0, nrow = 1, ncol = 2 * (ni + no)))
      f.rhs.rs <- 1
      if (rts == "vrs") {
        f.dir.rs <- "=="
        ubcon12 <- ubcon(val = rep(1, ndr), id = namevar2[1:ndr])
      } else if (rts == "nirs") {
        f.dir.rs <- "<="
        ubcon12 <- ubcon(val = rep(1, ndr), id = namevar2[1:ndr])
      } else if (rts == "ndrs") {
        f.dir.rs <- ">="
      } else {
        f.con.rs <- rbind(f.con.rs, f.con.rs)
        f.con2.rs <- rbind(f.con2.rs, f.con2.rs)
        f.dir.rs <- c(">=", "<=")
        f.rhs.rs <- c(L, U)
        ubcon12 <- ubcon(val = rep(U, ndr), id = namevar2[1:ndr])
      }
    }

    # Directions vector stage 1
    f.dir <- c(rep("<=", ni), rep(">=", no), f.dir.rs)
    f.dir[c(nc_inputs, ni + nc_outputs)] <- "=="

    if (maxslack && (!returnqp)) {

      # Checking weights
      if (is.matrix(weight_slack_i)) {
        if ((nrow(weight_slack_i) != ni) || (ncol(weight_slack_i) != nde)) {
          stop("Invalid weight input matrix (number of inputs x number of evaluated DMUs).")
        }
      } else if ((length(weight_slack_i) == 1) || (length(weight_slack_i) == ni)) {
        weight_slack_i <- matrix(weight_slack_i, nrow = ni, ncol = nde)
      } else {
        stop("Invalid weight input vector (number of inputs).")
      }
      rownames(weight_slack_i) <- inputnames
      colnames(weight_slack_i) <- dmunames[dmu_eval]
      weight_slack_i[nd_inputs, ] <- 0 # Non-discretionary io not taken into account for max slack solution

      if (is.matrix(weight_slack_o)) {
        if ((nrow(weight_slack_o) != no) || (ncol(weight_slack_o) != nde)) {
          stop("Invalid weight output matrix (number of outputs x number of evaluated DMUs).")
        }
      } else if ((length(weight_slack_o) == 1) || (length(weight_slack_o) == no)) {
        weight_slack_o <- matrix(weight_slack_o, nrow = no, ncol = nde)
      } else {
        stop("Invalid weight output vector (number of outputs).")
      }
      rownames(weight_slack_o) <- outputnames
      colnames(weight_slack_o) <- dmunames[dmu_eval]
      weight_slack_o[nd_outputs, ] <- 0 # Non-discretionary io not taken into account for max slack solution

      # Constraints matrix stage 2
      f.con2.1 <- cbind(inputref, diag(ni), matrix(0, nrow = ni, ncol = no),
                        diag(-qa, ni), matrix(0, nrow = ni, ncol = no))
      f.con2.1[nc_inputs, (ndr + 1) : (ndr + ni)] <- 0
      f.con2.1[nc_inputs, (ndr + ni + no + 1) : (ndr + 2 * ni + no)] <- 0

      f.con2.2 <- cbind(outputref, matrix(0, nrow = no, ncol = ni), -diag(no),
                        matrix(0, nrow = no, ncol = ni), diag(qa, no))
      f.con2.2[nc_outputs, (ndr + ni + 1) : (ndr + ni + no)] <- 0
      f.con2.2[nc_outputs, (ndr + 2 * ni + no + 1) : (ndr + 2 * (ni + no))] <- 0

      f.con2.nc <- matrix(0, nrow = (nnci + nnco), ncol = (ndr + 2 * (ni + no)))
      f.con2.nc[, ndr + c(nc_inputs, ni + nc_outputs)] <- diag(nnci + nnco)

      f.con2 <- rbind(f.con2.1, f.con2.2, f.con2.nc, f.con2.rs)

      # Directions vector stage 2
      f.dir2 <- c(rep("==", ni + no + nnci + nnco), f.dir.rs)

    }

    # Quadratic matrices
    Qin1 <- array(0, dim = c(ni, 1 + ndr + ni + no, 1 + ndr + ni + no))
    Qout1 <- array(0, dim = c(no, 1 + ndr + ni + no, 1 + ndr + ni + no))
    Qin2 <- array(0, dim = c(ni, ndr + 2 * (ni + no), ndr + 2 * (ni + no)))
    Qout2 <- array(0, dim = c(no, ndr + 2 * (ni + no), ndr + 2 * (ni + no)))
    Qin1[, 2:(ndr + 1), 2:(ndr + 1)] <- covXref
    Qout1[, 2:(ndr + 1), 2:(ndr + 1)] <- covYref
    Qin2[, 1:ndr, 1:ndr] <- covXref
    Qout2[, 1:ndr, 1:ndr] <- covYref
    for (j in 1: ni) {
      Qin1[j, 1 + ndr + j, 1 + ndr + j] <- -1
      Qin2[j, ndr + ni + no + j, ndr + ni + no + j] <- -1
    }
    for (j in 1:no) {
      Qout1[j, 1 + ndr + ni + j, 1 + ndr + ni + j] <- -1
      Qout2[j, ndr + 2 * ni + no + j, ndr + 2 * ni + no + j] <- -1
    }

    # Linear part of quadratic constraints
    ain1 <- matrix(0, nrow = ni, ncol = 1 + ndr + ni + no)
    aout1 <- matrix(0, nrow = no, ncol = 1 + ndr + ni + no)
    ain2 <- matrix(0, nrow = ni, ncol = ndr + 2 * (ni + no))
    aout2 <- matrix(0, nrow = no, ncol = ndr + 2 * (ni + no))

    qclist <- vector(mode = "list", length = 2 * (ni + no))

    # For loop ------
    if (parallel) {

      DMU <- foreach(i = 1:nde) %dopar% {

        ii <- dmu_eval[i]

        # Constraints matrix stage 1
        f.con.1 <- cbind(dir_input[, i], inputref, diag(-qa, ni),
                         matrix(0, nrow = ni, ncol = no))
        f.con.1[ncd_inputs, 1] <- 0
        f.con.1[nc_inputs, (ndr + 2) : (ndr + ni + 1)] <- 0
        f.con.2 <- cbind(-dir_output[, i], outputref,
                         matrix(0, nrow = no, ncol = ni), diag(qa, no))
        f.con.2[ncd_outputs, 1] <- 0
        f.con.2[nc_outputs, (ndr + 2) : (ndr + ni + 1)] <- 0
        f.con <- rbind(f.con.1, f.con.2, f.con.rs)

        # Right hand side vector stage 1
        f.rhs <- c(input[, ii], output[, ii], f.rhs.rs)

        rownames(f.con) <- paste("lc", 1:nrow(f.con), sep = "") # to prevent names errors in lincon
        lincon1 <- lincon(A = f.con, dir = f.dir, val = f.rhs, id = namevar1)

        # Linear part of quadratic constraints
        covXo <- matrix(covX[, ii, dmu_ref], nrow = ni, ncol = ndr)
        covYo <- matrix(covY[, ii, dmu_ref], nrow = no, ncol = ndr)
        ain1[, 2:(ndr + 1)] <- -2 * covXo
        aout1[, 2:(ndr + 1)] <- -2 * covYo

        # rhs quadratic constraints
        valin1 <- -covX[, ii, ii]
        valout1 <- -covY[, ii, ii]

        quadconin1a <- vector(mode = "list", length = ni)
        quadconin1b <- quadconin1a
        quadconout1a <- vector(mode = "list", length = no)
        quadconout1b <- quadconout1a
        for (j in 1:ni) {
          if (j %in% nc_inputs) {
            quadconin1a[[j]] <- NULL
            quadconin1b[[j]] <- NULL
          } else {
            quadconin1a[[j]] <- quadcon(Q = Qin1[j, , ], a = ain1[j, ], val = valin1[j], id = namevar1)
            quadconin1b[[j]] <- quadcon(Q = -Qin1[j, , ], a = -ain1[j, ], val = -valin1[j], id = namevar1)
          }
        }
        for (j in 1:no) {
          if (j %in% nc_outputs) {
            quadconout1a[[j]] <- NULL
            quadconout1b[[j]] <- NULL
          } else {
            quadconout1a[[j]] <- quadcon(Q = Qout1[j, , ], a = aout1[j, ], val = valout1[j], id = namevar1)
            quadconout1b[[j]] <- quadcon(Q = -Qout1[j, , ], a = -aout1[j, ], val = -valout1[j], id = namevar1)
          }
        }

        mycop <- cop(f = f.obj, max = TRUE, lb = lbcon1, ub = ubcon12, lc = lincon1)
        qclist <- c(quadconin1a, quadconin1b, quadconout1a, quadconout1b)
        names(qclist) <- paste("qc", 1:length(qclist), sep = "")
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
            if ((n_attempts == 1) && give_X && (ii %in% dmu_ref)) {
              Xini <- rep(0, ndr + ni + no + 1)
              Xini[which(dmu_ref == ii) + 1] <- 1
              names(Xini) <- namevar1
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
            beta <- res[1]

            if (maxslack) {

              # Objective function coefficients stage 2
              f.obj2 <- linfun(a = c(rep(0, ndr), weight_slack_i[, i],
                                     weight_slack_o[, i], rep(0, ni + no)),
                               id = namevar2)

              # Right hand side vector stage 2
              f.rhs2 <- c(input[, ii] - beta * dir_input[, i],
                          output[, ii] + beta * dir_output[, i],
                          rep(0, nnci + nnco), f.rhs.rs)
              f.rhs2[c(ncd_inputs, ni + ncd_outputs)] <- c(input[ncd_inputs, ii], output[ncd_outputs, ii])

              rownames(f.con2) <- paste("lc", 1:nrow(f.con2), sep = "") # to prevent names errors in lincon
              lincon2 <- lincon(A = f.con2, dir = f.dir2, val = f.rhs2, id = namevar2)

              ain2[, 1:ndr] <- -2 * covXo
              valin2 <- -covX[, ii, ii]
              aout2[, 1:ndr] <- -2 * covYo
              valout2 <- -covY[, ii, ii]

              quadconin2a <- vector(mode = "list", length = ni)
              quadconin2b <- quadconin2a
              quadconout2a <- vector(mode = "list", length = no)
              quadconout2b <- quadconout2a
              for (j in 1:ni) {
                if (j %in% nc_inputs) {
                  quadconin2a[[j]] <- NULL
                  quadconin2b[[j]] <- NULL
                } else {
                  quadconin2a[[j]] <- quadcon(Q = Qin2[j, , ], a = ain2[j, ], val = valin2[j], id = namevar2)
                  quadconin2b[[j]] <- quadcon(Q = -Qin2[j, , ], a = -ain2[j, ], val = -valin2[j], id = namevar2)
                }
              }
              for (j in 1:no) {
                if (j %in% nc_outputs) {
                  quadconout2a[[j]] <- NULL
                  quadconout2b[[j]] <- NULL
                } else {
                  quadconout2a[[j]] <- quadcon(Q = Qout2[j, , ], a = aout2[j, ], val = valout2[j], id = namevar2)
                  quadconout2b[[j]] <- quadcon(Q = -Qout2[j, , ], a = -aout2[j, ], val = -valout2[j], id = namevar2)
                }
              }

              mycop <- cop(f = f.obj2, max = TRUE, lb = lbcon2, ub = ubcon12, lc = lincon2)
              qclist <- c(quadconin2a, quadconin2b, quadconout2a, quadconout2b)
              auxkk <- which(sapply(qclist, is.null))
              if (length(auxkk) > 0) {
                qclist <- qclist[-auxkk] # remove NULLs from non-controllable
              }
              mycop$qc <- qclist

              n_attempts <- 1

              while (n_attempts <= n_attempts_max) {

                # Initial vector
                if ((n_attempts == 1) && give_X) {
                  lambda <- res[2 : (ndr + 1)]
                  sigma_input <- res[(ndr + 2):(ndr + ni + 1)]
                  sigma_output <- res[(ndr + ni + 2):(ndr + ni + no + 1)]
                  target_input <- as.vector(inputref %*% lambda)
                  target_output <- as.vector(outputref %*% lambda)
                  slack_input <- input[, ii] - beta * dir_input[, i] - target_input + qa * sigma_input
                  slack_input[nc_inputs] <- 0
                  slack_input[nd_inputs] <- input[nd_inputs, ii] - target_input[nd_inputs] +
                    qa * sigma_input[nd_inputs]
                  slack_output <- -output[, ii] - beta * dir_output[, i] + target_output + qa * sigma_output
                  slack_output[nc_outputs] <- 0
                  slack_output[nd_outputs] <- -output[nd_outputs, ii] + target_output[nd_outputs] +
                    qa * sigma_output[nd_outputs]
                  Xini <- c(lambda, slack_input, slack_output, sigma_input, sigma_output)
                  names(Xini) <- namevar2
                } else{
                  Xini <- NULL
                }

                res <- solvecop(mycop, solver = solver, quiet = TRUE, X = Xini, ...)

                if (res$status == "successful convergence") {
                  n_attempts <- n_attempts_max
                }
                n_attempts <- n_attempts + 1

              }

              if (res$status == "successful convergence") {
                res <- res$x
                lambda <- res[1 : ndr]
                names(lambda) <- dmunames[dmu_ref]

                slack_input <- res[(ndr + 1) : (ndr + ni)]
                names(slack_input) <- inputnames
                slack_output <- res[(ndr + ni + 1) : (ndr + ni + no)]
                names(slack_output) <- outputnames
                sigma_input <- res[(ndr + ni + no + 1):(ndr + ni + no + ni)]
                names(sigma_input) <- inputnames
                sigma_output <- res[(ndr + 2 * ni + no + 1):(ndr + 2 * (ni + no))]
                names(sigma_output) <- outputnames

                if (compute_target) {
                  target_input <- as.vector(inputref %*% lambda)
                  target_output <- as.vector(outputref %*% lambda)
                  names(target_input) <- inputnames
                  names(target_output) <- outputnames
                }

              } else {

                beta <- NA
                lambda <- NA
                slack_input <- NA
                slack_output <- NA
                sigma_input <- NA
                sigma_output <- NA
                if (compute_target) {
                  target_input <- NA
                  target_output <- NA
                }

              }

            } else {

              lambda <- res[2 : (ndr + 1)]
              names(lambda) <- dmunames[dmu_ref]
              sigma_input <- res[(ndr + 2):(ndr + ni + 1)]
              names(sigma_input) <- inputnames
              sigma_output <- res[(ndr + ni + 2):(ndr + ni + no + 1)]
              names(sigma_output) <- outputnames

              target_input <- as.vector(inputref %*% lambda)
              names(target_input) <- inputnames
              target_output <- as.vector(outputref %*% lambda)
              names(target_output) <- outputnames

              slack_input <- input[, ii] - beta * dir_input[, i] - target_input + qa * sigma_input
              slack_input[nc_inputs] <- 0
              slack_input[nd_inputs] <- input[nd_inputs, ii] - target_input[nd_inputs] +
                qa * sigma_input[nd_inputs]
              names(slack_input) <- inputnames
              slack_output <- -output[, ii] - beta * dir_output[, i] + target_output + qa * sigma_output
              slack_output[nc_outputs] <- 0
              slack_output[nd_outputs] <- -output[nd_outputs, ii] + target_output[nd_outputs] +
                qa * sigma_output[nd_outputs]
              names(slack_output) <- outputnames

            }

          } else {

            beta <- NA
            lambda <- NA
            slack_input <- NA
            slack_output <- NA
            sigma_input <- NA
            sigma_output <- NA
            if (compute_target) {
              target_input <- NA
              target_output <- NA
            }
          }

          list(beta = beta,
               lambda = lambda,
               slack_input = slack_input, slack_output = slack_output,
               sigma_input = sigma_input, sigma_output = sigma_output,
               target_input = target_input, target_output = target_output)

        }
      }

    } else {

      DMU <- vector(mode = "list", length = nde)

      for (i in 1:nde) {

        ii <- dmu_eval[i]

        # Constraints matrix stage 1
        f.con.1 <- cbind(dir_input[, i], inputref, diag(-qa, ni),
                         matrix(0, nrow = ni, ncol = no))
        f.con.1[ncd_inputs, 1] <- 0
        f.con.1[nc_inputs, (ndr + 2) : (ndr + ni + 1)] <- 0
        f.con.2 <- cbind(-dir_output[, i], outputref,
                         matrix(0, nrow = no, ncol = ni), diag(qa, no))
        f.con.2[ncd_outputs, 1] <- 0
        f.con.2[nc_outputs, (ndr + 2) : (ndr + ni + 1)] <- 0
        f.con <- rbind(f.con.1, f.con.2, f.con.rs)

        # Right hand side vector stage 1
        f.rhs <- c(input[, ii], output[, ii], f.rhs.rs)

        rownames(f.con) <- paste("lc", 1:nrow(f.con), sep = "") # to prevent names errors in lincon
        lincon1 <- lincon(A = f.con, dir = f.dir, val = f.rhs, id = namevar1)

        # Linear part of quadratic constraints
        covXo <- matrix(covX[, ii, dmu_ref], nrow = ni, ncol = ndr)
        covYo <- matrix(covY[, ii, dmu_ref], nrow = no, ncol = ndr)
        ain1[, 2:(ndr + 1)] <- -2 * covXo
        aout1[, 2:(ndr + 1)] <- -2 * covYo

        # rhs quadratic constraints
        valin1 <- -covX[, ii, ii]
        valout1 <- -covY[, ii, ii]

        quadconin1a <- vector(mode = "list", length = ni)
        quadconin1b <- quadconin1a
        quadconout1a <- vector(mode = "list", length = no)
        quadconout1b <- quadconout1a
        for (j in 1:ni) {
          if (j %in% nc_inputs) {
            quadconin1a[[j]] <- NULL
            quadconin1b[[j]] <- NULL
          } else {
            quadconin1a[[j]] <- quadcon(Q = Qin1[j, , ], a = ain1[j, ], val = valin1[j], id = namevar1)
            quadconin1b[[j]] <- quadcon(Q = -Qin1[j, , ], a = -ain1[j, ], val = -valin1[j], id = namevar1)
          }
        }
        for (j in 1:no) {
          if (j %in% nc_outputs) {
            quadconout1a[[j]] <- NULL
            quadconout1b[[j]] <- NULL
          } else {
            quadconout1a[[j]] <- quadcon(Q = Qout1[j, , ], a = aout1[j, ], val = valout1[j], id = namevar1)
            quadconout1b[[j]] <- quadcon(Q = -Qout1[j, , ], a = -aout1[j, ], val = -valout1[j], id = namevar1)
          }
        }

        mycop <- cop(f = f.obj, max = TRUE, lb = lbcon1, ub = ubcon12, lc = lincon1)
        qclist <- c(quadconin1a, quadconin1b, quadconout1a, quadconout1b)
        names(qclist) <- paste("qc", 1:length(qclist), sep = "")
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
            if ((n_attempts == 1) && give_X && (ii %in% dmu_ref)) {
              Xini <- rep(0, ndr + ni + no + 1)
              Xini[which(dmu_ref == ii) + 1] <- 1
              names(Xini) <- namevar1
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
            beta <- res[1]

            if (maxslack) {

              # Objective function coefficients stage 2
              f.obj2 <- linfun(a = c(rep(0, ndr), weight_slack_i[, i],
                                     weight_slack_o[, i], rep(0, ni + no)),
                               id = namevar2)

              # Right hand side vector stage 2
              f.rhs2 <- c(input[, ii] - beta * dir_input[, i],
                          output[, ii] + beta * dir_output[, i],
                          rep(0, nnci + nnco), f.rhs.rs)
              f.rhs2[c(ncd_inputs, ni + ncd_outputs)] <- c(input[ncd_inputs, ii], output[ncd_outputs, ii])

              rownames(f.con2) <- paste("lc", 1:nrow(f.con2), sep = "") # to prevent names errors in lincon
              lincon2 <- lincon(A = f.con2, dir = f.dir2, val = f.rhs2, id = namevar2)

              ain2[, 1:ndr] <- -2 * covXo
              valin2 <- -covX[, ii, ii]
              aout2[, 1:ndr] <- -2 * covYo
              valout2 <- -covY[, ii, ii]

              quadconin2a <- vector(mode = "list", length = ni)
              quadconin2b <- quadconin2a
              quadconout2a <- vector(mode = "list", length = no)
              quadconout2b <- quadconout2a
              for (j in 1:ni) {
                if (j %in% nc_inputs) {
                  quadconin2a[[j]] <- NULL
                  quadconin2b[[j]] <- NULL
                } else {
                  quadconin2a[[j]] <- quadcon(Q = Qin2[j, , ], a = ain2[j, ], val = valin2[j], id = namevar2)
                  quadconin2b[[j]] <- quadcon(Q = -Qin2[j, , ], a = -ain2[j, ], val = -valin2[j], id = namevar2)
                }
              }
              for (j in 1:no) {
                if (j %in% nc_outputs) {
                  quadconout2a[[j]] <- NULL
                  quadconout2b[[j]] <- NULL
                } else {
                  quadconout2a[[j]] <- quadcon(Q = Qout2[j, , ], a = aout2[j, ], val = valout2[j], id = namevar2)
                  quadconout2b[[j]] <- quadcon(Q = -Qout2[j, , ], a = -aout2[j, ], val = -valout2[j], id = namevar2)
                }
              }

              mycop <- cop(f = f.obj2, max = TRUE, lb = lbcon2, ub = ubcon12, lc = lincon2)
              qclist <- c(quadconin2a, quadconin2b, quadconout2a, quadconout2b)
              auxkk <- which(sapply(qclist, is.null))
              if (length(auxkk) > 0) {
                qclist <- qclist[-auxkk] # remove NULLs from non-controllable
              }
              mycop$qc <- qclist

              n_attempts <- 1

              while (n_attempts <= n_attempts_max) {

                # Initial vector
                if ((n_attempts == 1) && give_X) {
                  lambda <- res[2 : (ndr + 1)]
                  sigma_input <- res[(ndr + 2):(ndr + ni + 1)]
                  sigma_output <- res[(ndr + ni + 2):(ndr + ni + no + 1)]
                  target_input <- as.vector(inputref %*% lambda)
                  target_output <- as.vector(outputref %*% lambda)
                  slack_input <- input[, ii] - beta * dir_input[, i] - target_input + qa * sigma_input
                  slack_input[nc_inputs] <- 0
                  slack_input[nd_inputs] <- input[nd_inputs, ii] - target_input[nd_inputs] +
                    qa * sigma_input[nd_inputs]
                  slack_output <- -output[, ii] - beta * dir_output[, i] + target_output + qa * sigma_output
                  slack_output[nc_outputs] <- 0
                  slack_output[nd_outputs] <- -output[nd_outputs, ii] + target_output[nd_outputs] +
                    qa * sigma_output[nd_outputs]
                  Xini <- c(lambda, slack_input, slack_output, sigma_input, sigma_output)
                  names(Xini) <- namevar2
                } else{
                  Xini <- NULL
                }

                res <- solvecop(mycop, solver = solver, quiet = TRUE, X = Xini, ...)

                if (res$status == "successful convergence") {
                  n_attempts <- n_attempts_max
                }
                n_attempts <- n_attempts + 1

              }

              if (res$status == "successful convergence") {
                res <- res$x
                lambda <- res[1 : ndr]
                names(lambda) <- dmunames[dmu_ref]

                slack_input <- res[(ndr + 1) : (ndr + ni)]
                names(slack_input) <- inputnames
                slack_output <- res[(ndr + ni + 1) : (ndr + ni + no)]
                names(slack_output) <- outputnames
                sigma_input <- res[(ndr + ni + no + 1):(ndr + ni + no + ni)]
                names(sigma_input) <- inputnames
                sigma_output <- res[(ndr + 2 * ni + no + 1):(ndr + 2 * (ni + no))]
                names(sigma_output) <- outputnames

                if (compute_target) {
                  target_input <- as.vector(inputref %*% lambda)
                  target_output <- as.vector(outputref %*% lambda)
                  names(target_input) <- inputnames
                  names(target_output) <- outputnames
                }

              } else {

                beta <- NA
                lambda <- NA
                slack_input <- NA
                slack_output <- NA
                sigma_input <- NA
                sigma_output <- NA
                if (compute_target) {
                  target_input <- NA
                  target_output <- NA
                }

              }

            } else {

              lambda <- res[2 : (ndr + 1)]
              names(lambda) <- dmunames[dmu_ref]
              sigma_input <- res[(ndr + 2):(ndr + ni + 1)]
              names(sigma_input) <- inputnames
              sigma_output <- res[(ndr + ni + 2):(ndr + ni + no + 1)]
              names(sigma_output) <- outputnames

              target_input <- as.vector(inputref %*% lambda)
              names(target_input) <- inputnames
              target_output <- as.vector(outputref %*% lambda)
              names(target_output) <- outputnames

              slack_input <- input[, ii] - beta * dir_input[, i] - target_input + qa * sigma_input
              slack_input[nc_inputs] <- 0
              slack_input[nd_inputs] <- input[nd_inputs, ii] - target_input[nd_inputs] +
                qa * sigma_input[nd_inputs]
              names(slack_input) <- inputnames
              slack_output <- -output[, ii] - beta * dir_output[, i] + target_output + qa * sigma_output
              slack_output[nc_outputs] <- 0
              slack_output[nd_outputs] <- -output[nd_outputs, ii] + target_output[nd_outputs] +
                qa * sigma_output[nd_outputs]
              names(slack_output) <- outputnames

            }

          } else {

            beta <- NA
            lambda <- NA
            slack_input <- NA
            slack_output <- NA
            sigma_input <- NA
            sigma_output <- NA
            if (compute_target) {
              target_input <- NA
              target_output <- NA
            }
          }

          DMU[[i]] <- list(beta = beta,
                           lambda = lambda,
                           slack_input = slack_input, slack_output = slack_output,
                           sigma_input = sigma_input, sigma_output = sigma_output,
                           target_input = target_input, target_output = target_output)

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

    deaOutput <- list(modelname = "stoch_dir_dd",
                      dir_input = dir_input,
                      dir_output = dir_output,
                      rts = rts,
                      L = L,
                      U = U,
                      DMU = DMU,
                      data = datadea,
                      alpha = alpha,
                      dmu_eval = dmu_eval,
                      dmu_ref = dmu_ref,
                      maxslack = maxslack,
                      weight_slack_i = weight_slack_i,
                      weight_slack_o = weight_slack_o)

    return(structure(deaOutput, class = "dea_stoch"))

  }
