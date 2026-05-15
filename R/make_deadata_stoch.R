#' @title make_deadata_stoch
#'
#' @description This function creates, from a \code{deadata} object, a
#' \code{deadata_stoch} object by adding the corresponding covariance matrices.
#' These objects are prepared to be passed to a modelstoch_* function.
#'
#' We consider \eqn{\mathcal{D}=\left\{ \textrm{DMU}_1, \ldots ,\textrm{DMU}_n \right\} }
#' a set of \eqn{n} DMUs with \eqn{m} stochastic inputs and \eqn{s} stochastic outputs.
#' Matrices \eqn{\tilde{X}=(\tilde{x}_{ij})} and \eqn{\tilde{Y}=(\tilde{y}_{rj})}
#' are the input and output data matrices, respectively, where \eqn{\tilde{x}_{ij}}
#' and \eqn{\tilde{y}_{rj}} represent the \eqn{i}-th input and \eqn{r}-th output
#' of the \eqn{j}-th DMU. Moreover, we denote by \eqn{X=(x_{ij})} and \eqn{Y=(y_{rj})}
#' their expected values.
#' We suppose that \eqn{\tilde{x}_{ij}} and \eqn{\tilde{y}_{rj}} follow a multivariate
#' probability distribution with means \eqn{E(\tilde{x}_{ij})=x_{ij}}, \eqn{E(\tilde{y}_{rj})=y_{rj}}
#' and covariance matrix
#' \deqn{\Delta = \begin{pmatrix}
#' \Delta ^{II}_{11} & \cdots & \Delta ^{II}_{1m} & \Delta ^{IO}_{11} & \cdots & \Delta ^{IO}_{1s} \\
#' \vdots & & \vdots & \vdots & & \vdots \\
#' \Delta ^{II}_{m1} & \cdots & \Delta ^{II}_{mm} & \Delta ^{IO}_{m1} & \cdots & \Delta ^{IO}_{ms} \\ \\
#' \Delta ^{OI}_{11} & \cdots & \Delta ^{OI}_{1m} & \Delta ^{OO}_{11} & \cdots & \Delta ^{OO}_{1s} \\
#' \vdots & & \vdots & \vdots & & \vdots \\
#' \Delta ^{OI}_{s1} & \cdots & \Delta ^{OI}_{sm} & \Delta ^{OO}_{s1} & \cdots & \Delta ^{OO}_{ss} \\
#' \end{pmatrix}_{n(m+s)\times n(m+s)}}
#' where \eqn{\Delta ^{II}_{ik}}, \eqn{\Delta ^{OO}_{rp}}, \eqn{\Delta ^{IO}_{ir}}
#' and \eqn{\Delta ^{OI}_{ri}} are \eqn{n\times n} matrices, for \eqn{1\leq i,k\leq m}
#' and \eqn{1\leq r,p\leq s}, such that
#' \deqn{\left( \Delta ^{II}_{ik}\right) _{jq}=\mathrm{Cov}(\tilde{x}_{ij}, \tilde{x}_{kq}),}
#' \deqn{\left( \Delta ^{OO}_{rp}\right) _{jq}=\mathrm{Cov}(\tilde{y}_{rj}, \tilde{y}_{pq}),}
#' \deqn{\left( \Delta ^{IO}_{ir}\right) _{jq}=\left( \Delta ^{OI}_{ri}\right) _{qj}
#' =\mathrm{Cov}(\tilde{x}_{ij}, \tilde{y}_{rq}),}
#' with \eqn{1\leq j,q\leq n}.
#'
#' - If we have the covariances matrix in the general form above, it can be introduced
#' directly by parameter \code{cov_matrix}. Since this matrix is supposed to be
#' symmetric, only values above the diagonal are read, ignoring values below the diagonal.
#'
#' - Alternatively, we can introduce the covariances matrix using parameters \code{cov_II},
#' \code{cov_OO} and \code{cov_IO}, that are 4-dimensional arrays of size \eqn{m\times
#' m\times n\times n}, \eqn{s\times s\times n\times n} and \eqn{m\times s\times n\times n},
#' respectively, such that
#' \deqn{\texttt{cov\_II[i, k, , ]}=\Delta ^{II}_{ik},}
#' \deqn{\texttt{cov\_OO[r, p, , ]}=\Delta ^{OO}_{rp},}
#' \deqn{\texttt{cov\_IO[i, r, , ]}=\Delta ^{IO}_{ir},}
#' for \eqn{1\leq i,k\leq m} and \eqn{1\leq r,p\leq s}. Since matrices \eqn{\Delta ^{II}_{ii}} and \eqn{\Delta
#' ^{OO}_{rr}} are supposed to be symmetric, only values above the diagonal are
#' read, ignoring values below the diagonal. Moreover, since \eqn{\Delta ^{II}_{ki}}
#' is the transpose of \eqn{\Delta ^{II}_{ik}}, and \eqn{\Delta ^{OO}_{pr}} is
#' the transpose of \eqn{\Delta ^{OO}_{rp}}, only matrices \eqn{\Delta ^{II}_{ik}}
#' and \eqn{\Delta ^{OO}_{rp}} with \eqn{i\leq k} and \eqn{r\leq p} are necessary,
#' ignoring those with \eqn{i>k} and \eqn{r>p}.
#'
#' - If covariances between different inputs/outputs are zero, we can make use of
#' parameters \code{cov_input} and \code{cov_output}, that are arrays of size
#' \eqn{m\times n\times n} and \eqn{s\times n\times n}, respectively, such that
#' \deqn{\texttt{cov\_input[i, , ]}=\Delta ^{II}_{ii},}
#' \deqn{\texttt{cov\_output[r, , ]}=\Delta ^{OO}_{rr},}
#' for \eqn{1\leq i\leq m} and \eqn{1\leq r\leq s}. By symmetry of \eqn{\Delta
#' ^{II}_{ii}} and \eqn{\Delta ^{OO}_{rr}}, only values above the diagonal are
#' read, ignoring values below the diagonal.
#'
#' - Finally, if all the variables are independent then the covariances matrix is
#' diagonal. Hence, we might use parameters \code{var_input} and \code{var_output},
#' that are matrices of size \eqn{m\times n} and \eqn{s\times n}, respectively,
#' such that
#' \deqn{\texttt{var\_input[i, j]}=\mathrm{Var}\left( \tilde{x}_{ij}\right) ,}
#' \deqn{\texttt{var\_output[r, j]}=\mathrm{Var}\left( \tilde{y}_{rj}\right) ,}
#' for \eqn{1\leq i\leq m}, \eqn{1\leq r\leq s} and \eqn{1\leq j\leq n}.
#'
#' @usage make_deadata_stoch(datadea = NULL,
#'           var_input = NULL,
#'           var_output = NULL,
#'           cov_input = NULL,
#'           cov_output = NULL,
#'           cov_II = NULL,
#'           cov_OO = NULL,
#'           cov_IO = NULL,
#'           cov_matrix = NULL)
#'
#' @param datadea The \code{deadata} object.
#' @param var_input A matrix of size \code{m} x \code{n}, where \code{m} is the
#' number of inputs and \code{n} is the number of DMUs. It contains the variances
#' of each input of each DMU, such that \code{var_input[i, j]} is the variance
#' of the input \code{i} of DMU \code{j}.
#' Use this parameter if all covariances are 0.
#' @param var_output A matrix of size \code{s} x \code{n}, where \code{s} is the
#' number of outputs and \code{n} is the number of DMUs. It contains the variances
#' of each output of each DMU, such that \code{var_output[r, j]} is the variance
#' of the output \code{r} of DMU \code{j}.
#' Use this parameter if all covariances are 0.
#' @param cov_input An array of size \code{m} x \code{n} x \code{n} containing
#' the covariances of each input between DMUs, such that
#' \code{cov_input[i, j, k]} is the covariance between the input \code{i}
#' of DMU \code{j} and the input \code{i} of DMU \code{k}. The corresponding
#' variances are in the diagonal of each \code{n} x \code{n} submatrix. Since
#' these submatrices are supposed to be symmetric, only values above the diagonal
#' are read in order to reconstruct the symmetric submatrices, ignoring values
#' below the diagonal.
#' Use this parameter if covariances between different inputs are 0.
#' @param cov_output An array of size \code{s} x \code{n} x \code{n} containing
#' the covariances of each output between DMUs, such that
#' \code{cov_output[r, j, k]} is the covariance between the output \code{r}
#' of DMU \code{j} and the output \code{r} of DMU \code{k}. The corresponding
#' variances are in the diagonal of each \code{n} x \code{n} submatrix. Since
#' these submatrices are supposed to be symmetric, only values above the diagonal
#' are read in order to reconstruct the symmetric submatrices, ignoring values
#' below the diagonal.
#' Use this parameter if covariances between different outputs are 0.
#' @param cov_II An array of size \code{m} x \code{m} x \code{n} x \code{n}
#' containing the covariances between inputs and between DMUs, such that
#' \code{cov_II[i1, i2, j, k]} is the covariance between the input \code{i1}
#' of DMU \code{j} and the input \code{i2} of DMU \code{k}. For the input \code{i},
#' the corresponding variances are in the diagonal of the \code{n} x \code{n}
#' submatrices of the form \code{cov_II[i, i, , ]}. Since
#' these submatrices are supposed to be symmetric, only values above the diagonal
#' are read in order to reconstruct the symmetric submatrices, ignoring values
#' below the diagonal.
#' Moreover, since \code{cov_II[i2, i1, , ]} is the transpose of \code{cov_II[i1, i2, , ]},
#' only submatrices \code{cov_II[i1, i2, , ]} with \code{i1 <= i2} are necessary, ignoring
#' those with \code{i1 > i2}.
#' Use this parameter if covariances between inputs are nonzero.
#' @param cov_OO An array of size \code{s} x \code{s} x \code{n} x \code{n} containing
#' the covariances between outputs and between DMUs, such that
#' \code{cov_OO[r1, r2, j, k]} is the covariance between the output \code{r1}
#' of DMU \code{j} and the output \code{r2} of DMU \code{k}. For the output \code{r},
#' the corresponding variances are in the diagonal of the \code{n} x \code{n}
#' submatrices of the form \code{cov_OO[r, r, , ]}. Since
#' these submatrices are supposed to be symmetric, only values above the diagonal
#' are read in order to reconstruct the symmetric submatrices, ignoring values
#' below the diagonal.
#' Moreover, since \code{cov_OO[r2, r1, , ]} is the transpose of \code{cov_OO[r1, r2, , ]},
#' only submatrices \code{cov_OO[r1, r2, , ]} with \code{r1 <= r2} are necessary, ignoring
#' those with \code{r1 > r2}.
#' Use this parameter if covariances between outputs are nonzero.
#' @param cov_IO An array of size \code{m} x \code{s} x \code{n} x \code{n} containing
#' the covariances between inputs and outputs, and between DMUs, such that
#' \code{cov_IO[i, r, j, k]} is the covariance between the input \code{i}
#' of DMU \code{j} and the output \code{r} of DMU \code{k}.
#' Use this parameter if covariances between inputs and outputs are nonzero.
#' @param cov_matrix A matrix of size \code{n(m+s)} x \code{n(m+s)} following the
#' notation of Cooper et al. (1998). Since this matrix is supposed to be symmetric,
#' only values above the diagonal are read, ignoring values below the diagonal.
#'
#' @returns An object of class \code{deadata_stoch}.
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
#' Cooper, W.W.; Huang, Z.; Lelas, V.; Li, S.X.; Olesen, O.B. (1998). “Chance
#' Constrained Programming Formulations for Stochastic Characterizations of
#' Efficiency and Dominance in DEA", Journal of Productivity Analysis, 9, 53-79.
#'
#' Land, K.C; Lovell, C.A.K.; Thore, S. (1993). "Chance-constrained data envelopment analysis",
#' Managerial and Decision Economics, Vol. 14, No. 6, 541-554.
#'
#' @examples
#' # Example 1
#' library(deaR)
#' data("Coll_Blasco_2006")
#' ni <- 2 # number of inputs
#' no <- 2 # number of outputs
#' data_example <- make_deadata(datadea = Coll_Blasco_2006,
#'                              ni = ni,
#'                              no = no)
#' nd <- length(data_example$dmunames) # number of DMUs
#' # All variances are 1.
#' var_input <- matrix(1, nrow = ni, ncol = nd)
#' var_output <- matrix(1, nrow = no, ncol = nd)
#' data_stoch <- make_deadata_stoch(datadea = data_example,
#'                                  var_input = var_input,
#'                                  var_output = var_output)
#' # All covariances are 1.
#' cov_input <- array(1, dim = c(ni, nd, nd))
#' cov_output <- array(1, dim = c(no, nd, nd))
#' data_stoch2 <- make_deadata_stoch(datadea = data_example,
#'                                   cov_input = cov_input,
#'                                   cov_output = cov_output)
#'
#' # Example 2. Deterministic data with one stochastic input.
#' library(deaR)
#' dmunames <- c("A", "B", "C")
#' nd <- length(dmunames) # Number of DMUs
#' inputnames <- c("Input_1", "Input_2")
#' ni <- length(inputnames) # Number of Inputs
#' outputnames <- c("Output_1", "Output_2", "Output_3")
#' no <- length(outputnames) # Number of Outputs
#' X <- matrix(c(5, 14, 8, 15, 7, 12),
#'             nrow = ni, ncol = nd, dimnames = list(inputnames, dmunames))
#' Y <- matrix(c(9, 4, 16, 5, 7, 10, 4, 9, 13),
#'             nrow = no, ncol = nd, dimnames = list(outputnames, dmunames))
#' datadea <- make_deadata(inputs = X,
#'                         outputs = Y)
#' covX <- array(0, dim = c(2, 3, 3))
#' # The 2nd input is stochastic.
#' # Since the corresponding 3x3 covariances matrix is symmetric, only values
#' # above the diagonal are necessary.
#' covX[2, 1, ] <- c(1.4, 0.9, 0.6)
#' covX[2, 2, 2:3] <- c(1.5, 0.7)
#' covX[2, 3, 3] <- 1.2
#' # Alternatively (note that values below the diagonal are ignored).
#' covX[2, , ] <- matrix(c(1.4, 0.9, 0.6, 0, 1.5, 0.7, 0, 0, 1.2),
#'                       nrow = 3,
#'                       byrow = TRUE)
#' datadea_stoch <- make_deadata_stoch(datadea,
#'                                     cov_input = covX)
#'
#' # Example 3. Replication of Program Follow Through data in Land et al. (1993)
#' library(deaR)
#' data("PFT1981")
#' # Selecting DMUs in Program Follow Through (PFT)
#' PFT <- PFT1981[1:49, ]
#' PFT <- make_deadata(PFT,
#'                     inputs = 2:6,
#'                     outputs = 7:9)
#' c <- 0.5
#' var_output <- matrix(c^2, nrow = 3, ncol = 49)
#' PFT_stoch <- make_deadata_stoch(datadea = PFT, var_output = var_output)
#'
#' # Example 4. 5 random observations of 3 DMUs with 2 inputs and 2 outputs.
#' library(deaR)
#' # Generate random observations.
#' input1 <- data.frame(I1D1 = rnorm(5, mean = sample(5:10, 1)),
#'                      I1D2 = rnorm(5, mean = sample(5:10, 1)),
#'                      I1D3 = rnorm(5, mean = sample(5:10, 1)))
#' input2 <- data.frame(I2D1 = rnorm(5, mean = sample(5:10, 1)),
#'                      I2D2 = rnorm(5, mean = sample(5:10, 1)),
#'                      I2D3 = rnorm(5, mean = sample(5:10, 1)))
#' output1 <- data.frame(O1D1 = rnorm(5, mean = sample(5:10, 1)),
#'                       O1D2 = rnorm(5, mean = sample(5:10, 1)),
#'                       O1D3 = rnorm(5, mean = sample(5:10, 1)))
#' output2 <- data.frame(O2D1 = rnorm(5, mean = sample(5:10, 1)),
#'                       O2D2 = rnorm(5, mean = sample(5:10, 1)),
#'                       O2D3 = rnorm(5, mean = sample(5:10, 1)))
#' # Generate deadata with means of observations.
#' inputs <- matrix(mapply(mean, cbind(input1, input2)),
#'                  nrow = 2,
#'                  ncol = 3,
#'                  byrow = TRUE,
#'                  dimnames = list(c("I1", "I2"), c("D1", "D2", "D3")))
#' outputs <- matrix(mapply(mean, cbind(output1, output2)),
#'                   nrow = 2,
#'                   ncol = 3,
#'                   byrow = TRUE,
#'                   dimnames = list(c("O1", "O2"), c("D1", "D2", "D3")))
#' datadea <- make_deadata(inputs = inputs,
#'                         outputs = outputs)
#' # Generate covariances matrix cov_matrix.
#' cov_matrix <- cov(cbind(input1, input2, output1, output2))
#' # Generate deadata_stoch
#' datadea_stoch <- make_deadata_stoch(datadea,
#'                                     cov_matrix = cov_matrix)
#'
#' @export

make_deadata_stoch <- function(datadea = NULL,
                      var_input = NULL,
                      var_output = NULL,
                      cov_input = NULL,
                      cov_output = NULL,
                      cov_II = NULL,
                      cov_OO = NULL,
                      cov_IO = NULL,
                      cov_matrix = NULL) {

  if (!inherits(datadea, "deadata")) {
    stop("Data should be of class deadata. Run make_deadata function first!")
  }

  dmunames <- datadea$dmunames
  nd <- length(dmunames) # number of dmus
  inputnames <- rownames(datadea$input)
  outputnames <- rownames(datadea$output)
  ni <- nrow(datadea$input)
  no <- nrow(datadea$output)

  if (is.null(cov_matrix)) {

    if (is.null(cov_II)) {
      cov_II <- array(0, dim = c(ni, ni, nd, nd))
      if (is.null(cov_input)) {
        if (!is.null(var_input)) {
          if (length(dim(var_input)) != 2) {
            stop("var_input must be a matrix of size [number of inputs] x [number of DMUs].")
          }
          if (any(dim(var_input) != c(ni, nd))) {
            stop("var_input must be a matrix of size [number of inputs] x [number of DMUs].")
          }
          for (i in 1:ni) {
            for (j in 1:nd) {
              cov_II[i, i, j, j] <- var_input[i, j]
            }
          }
        }
      } else {
        if (length(dim(cov_input)) != 3) {
          stop("cov_input must be an array of size [number of inputs] x [number of DMUs] x [number of DMUs].")
        }
        if (any(dim(cov_input) != c(ni, nd, nd))) {
          stop("cov_input must be an array of size [number of inputs] x [number of DMUs] x [number of DMUs].")
        }
        # Reconstruct symmetrically covariances submatrices
        for (i in 1:ni) {
          if (!isSymmetric(cov_input[i, , ])) {
            for (j in 2:nd) {
              for (k in 1:(j - 1)) {
                cov_input[i, j, k] <- cov_input[i, k, j]
              }
            }
          }
          cov_II[i, i, , ] <- cov_input[i, , ]
        }
      }
    } else {
      if (length(dim(cov_II)) != 4) {
        stop("cov_II must be an array of size [number of inputs] x [number of inputs] x
           [number of DMUs] x [number of DMUs].")
      }
      if (any(dim(cov_II) != c(ni, ni, nd, nd))) {
        stop("cov_II must be an array of size [number of inputs] x [number of inputs] x
           [number of DMUs] x [number of DMUs].")
      }
      # Reconstruct symmetrically covariances submatrices
      for (i in 1:ni) {
        if (!isSymmetric(cov_II[i, i, , ])) {
          for (j in 2:nd) {
            for (k in 1:(j - 1)) {
              cov_II[i, i, j, k] <- cov_II[i, i, k, j]
            }
          }
        }
      }
      # Reconstruct transpose matrices
      if (ni > 1) {
        for (i in 2:ni) {
          for (j in 1:(i - 1)) {
            cov_II[i, j, , ] <- t(cov_II[j, i, , ])
          }
        }
      }
    }

    if (is.null(cov_OO)) {
      cov_OO <- array(0, dim = c(no, no, nd, nd))
      if (is.null(cov_output)) {
        if (!is.null(var_output)) {
          if (length(dim(var_output)) != 2) {
            stop("var_output must be a matrix of size [number of outputs] x [number of DMUs].")
          }
          if (any(dim(var_output) != c(no, nd))) {
            stop("var_output must be a matrix of size [number of outputs] x [number of DMUs].")
          }
          for (i in 1:no) {
            for (j in 1:nd) {
              cov_OO[i, i, j, j] <- var_output[i, j]
            }
          }
        }
      } else {
        if (length(dim(cov_output)) != 3) {
          stop("cov_output must be an array of size [number of outputs] x [number of DMUs] x [number of DMUs].")
        }
        if (any(dim(cov_output) != c(no, nd, nd))) {
          stop("cov_output must be an array of size [number of outputs] x [number of DMUs] x [number of DMUs].")
        }
        # Reconstruct symmetrically covariances submatrices
        for (i in 1:no) {
          if (!isSymmetric(cov_output[i, , ])) {
            for (j in 2:nd) {
              for (k in 1:(j - 1)) {
                cov_output[i, j, k] <- cov_output[i, k, j]
              }
            }
          }
          cov_OO[i, i, , ] <- cov_output[i, , ]
        }
      }
    } else {
      if (length(dim(cov_OO)) != 4) {
        stop("cov_OO must be an array of size [number of outputs] x [number of outputs] x
           [number of DMUs] x [number of DMUs].")
      }
      if (any(dim(cov_OO) != c(no, no, nd, nd))) {
        stop("cov_OO must be an array of size [number of outputs] x [number of outputs] x
           [number of DMUs] x [number of DMUs].")
      }
      # Reconstruct symmetrically covariances submatrices
      for (i in 1:no) {
        if (!isSymmetric(cov_OO[i, i, , ])) {
          for (j in 2:nd) {
            for (k in 1:(j - 1)) {
              cov_OO[i, i, j, k] <- cov_OO[i, i, k, j]
            }
          }
        }
      }
      # Reconstruct transpose matrices
      if (no > 1) {
        for (i in 2:no) {
          for (j in 1:(i - 1)) {
            cov_OO[i, j, , ] <- t(cov_OO[j, i, , ])
          }
        }
      }
    }

    if (is.null(cov_IO)) {
      cov_IO <- array(0, dim = c(ni, no, nd, nd))
    } else {
      if (length(dim(cov_IO)) != 4) {
        stop("cov_IO must be an array of size [number of inputs] x [number of outputs] x
           [number of DMUs] x [number of DMUs].")
      }
      if (any(dim(cov_IO) != c(ni, no, nd, nd))) {
        stop("cov_IO must be an array of size [number of inputs] x [number of outputs] x
           [number of DMUs] x [number of DMUs].")
      }
    }

  } else {

    if (length(dim(cov_matrix)) != 2) {
      stop("cov_matrix must be a square matrix of order (number of DMUs) * (number
      of inputs + number of outputs).")
    }
    ndio <- nd * (ni + no)
    if (any(dim(cov_matrix) != c(ndio, ndio))) {
      stop("cov_matrix must be a square matrix of order (number of DMUs) * (number
      of inputs + number of outputs).")
    }

    cov_II <- array(0, dim = c(ni, ni, nd, nd))
    cov_OO <- array(0, dim = c(no, no, nd, nd))
    cov_IO <- array(0, dim = c(ni, no, nd, nd))

    if (!isSymmetric(cov_matrix)) {
      cov_matrix[lower.tri(cov_matrix)] <- cov_matrix(upper.tri(cov_matrix))
    }

    for (i1 in 1:ni) {
      cov_II[i1, i1, , ] <- cov_matrix[(nd * (i1 - 1) + 1):(nd * i1),
                                       (nd * (i1 - 1) + 1):(nd * i1)]
      if (i1 < ni) {
        for (i2 in (i1 + 1):ni) {
          cov_II[i1, i2, , ] <- cov_matrix[(nd * (i1 - 1) + 1):(nd * i1),
                                           (nd * (i2 - 1) + 1):(nd * i2)]
          cov_II[i2, i1, , ] <- t(cov_II[i1, i2, , ])
        }
      }
    }

    for (i1 in 1:no) {
      cov_OO[i1, i1, , ] <- cov_matrix[(nd * (ni + i1 - 1) + 1):(nd * (ni + i1)),
                                       (nd * (ni + i1 - 1) + 1):(nd * (ni + i1))]
      if (i1 < no) {
        for (i2 in (i1 + 1):no) {
          cov_OO[i1, i2, , ] <- cov_matrix[(nd * (ni + i1 - 1) + 1):(nd * (ni + i1)),
                                           (nd * (ni + i2 - 1) + 1):(nd * (ni + i2))]
          cov_OO[i2, i1, , ] <- t(cov_OO[i1, i2, , ])

        }
      }
    }

    for (i1 in 1:ni) {
      for (i2 in 1:no) {
        cov_IO[i1, i2, , ] <- cov_matrix[(nd * (i1 - 1) + 1):(nd * i1),
                                         (nd * (ni + i2 - 1) + 1):(nd * (ni + i2))]
      }
    }

  }

  dimnames(cov_II) <- list(inputnames, inputnames, dmunames, dmunames)
  dimnames(cov_OO) <- list(outputnames, outputnames, dmunames, dmunames)
  dimnames(cov_IO) <- list(inputnames, outputnames, dmunames, dmunames)
  datadea$cov_II <- cov_II
  datadea$cov_OO <- cov_OO
  datadea$cov_IO <- cov_IO

  return(structure(datadea, class = c("deadata", "deadata_stoch")))

}
