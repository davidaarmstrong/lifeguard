#' Create matrix of lags
#'
#' Creates a matrix of lagged values of a variable.
#' @param x A vector of values to be lagged.
#' @param n An integer vector giving the lags to be created.  Including 0 will inculde the time t value of `x`.
#' @param rep0 If `TRUE`, replaces `NA` values with 0.  Default is `FALSE`.
#' @param ... Additional arguments, currently not implemented.
#' @return A matrix with the lagged values of `x`.  The first column is the time t value of `x` (if `0` is included in `n`), and subsequent columns are the lagged values.
#' @importFrom dplyr lag
#' @examples
#' lag_n(1:10, 0:3)
#'
#' @export
lag_n <- function(x, n=1, rep0 = FALSE, ...){
  lags <- sapply(n, \(t){
    out <- dplyr::lag(x, t)
    if(rep0)out <- ifelse(is.na(out), 0, out)
    out
  })
  if(!is.matrix(lags))lags <- matrix(lags, ncol=1)
  colnames(lags) <- paste0(substitute(x), "_L", n)
  lags
}

#' Create matrix of differences
#'
#' Creates a matrix of differenced values of a variable.
#' @param x A vector of values to be lagged.
#' @param n An integer vector of at least two values giving the differences to be created. For two values `n[1]` and `n[2]`, the function creates (`x` at time `n[2]`) - (`x` at time `n[1]`).
#' @param ... Additional arguments, currently not implemented.
#' @return A matrix with the differenced values of `x`.
#' @examples
#' diff_n(1:10, 0:3)
#'
#' @export
diff_n <- function(x, n=0:1, ...){
  if(length(n) < 2){stop("Differences require at least two lags\n")}
  lags <- lag_n(x, n)
  diffs <- lapply(1:(ncol(lags) - 1), \(i){
    lags[,i] - lags[,(i+1)]
  })
  diffs <- do.call(cbind, diffs)
  colnames(diffs) <- paste0(substitute(x), "_D", n[1:(length(n)-1)])
  diffs
}

#' CDF of Hinkley Distribution
#'
#' Calculates the CDF of the Hinkley distribution.  This is particularly useful for the ratio of normal variates.
#' @param r A numeric vector of ratios for which to calculate the CDF.
#' @param mu_a The mean of the numerator normal distribution.
#' @param mu_b The mean of the denominator normal distribution.
#' @param sigma_a The standard deviation of the numerator normal distribution.
#' @param sigma_b The standard deviation of the denominator normal distribution.
#' @param sigma_ab The covariance between the numerator and denominator normal distributions.
#' @param ... Additional arguments, currently not implemented.
#' @importFrom stats pnorm
#' @importFrom mvtnorm pmvnorm
#' @return A numeric vector of the CDF values for the given ratios.
#' @examples
#' x <- rnorm(100)
#' y <- x + rnorm(100, 0, .5)
#' v <- var(cbind(x,y))
#' phinkley(mean(x)/mean(y), mean(x), mean(y), v[1,1], v[2,2], v[1,2])
#' @export
phinkley <- function(r, mu_a, mu_b, sigma_a, sigma_b, sigma_ab, ...) {
  stopifnot(sigma_a > 0, sigma_b > 0, abs(sigma_ab) < sigma_a * sigma_b)
  r <- as.numeric(r)
  Phi  <- pnorm

  Phi2 <- function(x, y, rho) {
    mvtnorm::pmvnorm(lower = c(-Inf, -Inf),
                     upper = c(x, y),
                     mean  = c(0, 0),
                     sigma = matrix(c(1, rho, rho, 1), 2, 2))[1]
  }

  vapply(
    r,
    FUN = function(ri) {
      sigma_r <- sqrt(sigma_a^2 - 2 * ri * sigma_ab + (ri^2) * sigma_b^2)
      m_u     <- (mu_a - ri * mu_b) / sigma_r
      m_v     <-  mu_b / sigma_b
      rho_r   <- (sigma_ab - ri * sigma_b^2) / (sigma_r * sigma_b)
      Phi(-m_u) + Phi(-m_v) - 2 * Phi2(-m_u, -m_v, rho_r)
    },
    FUN.VALUE = numeric(1)
  )
}

#' Quantile function for the Hinkley distribution
#'
#' Estimates the quantiles of the Hinkley distribution.
#' @param p A numeric vector of quantile values.
#' @param mu_a The mean of the numerator normal distribution.
#' @param mu_b The mean of the denominator normal distribution.
#' @param sigma_a The standard deviation of the numerator normal distribution.
#' @param sigma_b The standard deviation of the denominator normal distribution.
#' @param sigma_ab The covariance between the numerator and denominator normal distributions.
#' @param tol The tolerance for the root-finding algorithm. Default is `1e-8`.
#' @param maxiter The maximum number of iterations for the root-finding algorithm. Default is `60`.
#' @param ... Additional arguments, currently not implemented.
#' @importFrom stats uniroot
#' @returns A numeric vector of quantiles corresponding to the input probabilities.
#' @examples
#' x <- rnorm(100)
#' y <- x + rnorm(100, 0, .5)
#' v <- var(cbind(x,y))
#' qhinkley(.975, 0, mean(y), v[1,1], v[2,2], v[1,2])
#' @export
qhinkley <- function(p, mu_a, mu_b, sigma_a, sigma_b, sigma_ab, tol = 1e-8, maxiter = 60, ...) {
  stopifnot(all(is.finite(p)), all(p >= 0 & p <= 1))
  q_one <- function(pi) {
    if (pi == 0) return(-Inf)
    if (pi == 1) return(Inf)

    # Find a bracket [L, U] with F(L) <= p <= F(U)
    L <- -1; U <- 1
    FL <- phinkley(L, mu_a, mu_b, sigma_a, sigma_b, sigma_ab)
    FU <- phinkley(U, mu_a, mu_b, sigma_a, sigma_b, sigma_ab)

    it <- 0
    while (FL > pi && it < maxiter) {
      U <- L; FU <- FL
      L <- 2 * L
      FL <- phinkley(L, mu_a, mu_b, sigma_a, sigma_b, sigma_ab)
      it <- it + 1
    }
    while (FU < pi && it < 2*maxiter) {
      L <- U; FL <- FU
      U <- 2 * U
      FU <- phinkley(U, mu_a, mu_b, sigma_a, sigma_b, sigma_ab)
      it <- it + 1
    }
    if (!(FL <= pi && pi <= FU)) {
      stop("Failed to bracket the quantile; try larger maxiter.", call. = FALSE)
    }

    uniroot(
      f = function(r) phinkley(r, mu_a, mu_b, sigma_a, sigma_b, sigma_ab) - pi,
      lower = L, upper = U, tol = tol
    )$root
  }
  vapply(p, q_one, numeric(1))
}

#' Confidence Intervals for Long-run Multipliers
#'
#' Calculates confidence intervals for long-run multipliers from a time-series analysis.
#' @details The function allows the use of either autoregressive distributed lag (ADL) or error correction models (ECM).  The confidence intervals are made from both the
#' t-distribution (using the delta method for standard errors) and with the Hinkley distribution.
#' Assuming `b` is the vector of model coefficients, the LRM for the ADL is calculated with
#' `sum(b[x_coefs[[i]])/(1-sum(b[y_coefs]))` for each of the `i` elements of `x_coefs`.  For the ECM, the LRM is calculated with `-b[x_coefs[[i]]/b[y_coefs]`.
#' @param obj A linear regression model in ADL form.
#' @param y_coefs A vector of index values for the lag-y coefficients from the model.
#' @param x_coefs For the ADL, a list of vectors of index values for the level and lag-x coefficients from the model. For the ECM, it should be a list, but each
#' element of the list should have a single value giving the coefficient on the lagged value of each x of interest.  Each
#' element of the list is assumed to be a different variable whose LRM is to be calculated.  If the list is
#' named, the output uses those names.
#' @param qb_dist The distribution for the quasi-Bayesian simulation either `"normal"` or `"t"` (default is `"t"`).
#' @param level The confidence level for the intervals.  Default is `0.95`.
#' @param t_min Minimum adjusted degrees of freedom in the t-distribution simulation.  Default is 2.
#' @param modtype The type of model, either `"adl"` (default) or `"ecm"`.
#' @param ... Additional arguments, currently not implemented.
#' @importFrom stats coef vcov terms qt quantile
#' @importFrom mvtnorm rmvnorm rmvt
#' @examples
#' x <- y <- c(0, rep(NA, 99))
#' for(t in 2:100){
#'   x[t] <- .5*x[(t-1)] + rnorm(1,0,1)
#' }
#'
#' for(t in 2:100){
#' y[t] <- .7*y[(t-1)] + .5*x[t] + .25*x[(t-1)] + rnorm(1,0,1)
#' }
#'
#' d <- data.frame(time=1:100, x=x, y=y)
#' adl <- lm(y ~ lag_n(y,1) + lag_n(x, c(0,1)), data=d)
#' lrm_ci(adl, y_coefs = 2, x_coefs = list(x=3:4))
#' @returns A data frame with the following variables:
#' * `vbl` - The name of the variable for which the long-run multiplier is calculated
#' * `est` - The estimate of the long-run multiplier
#' * `delta_se` - The approximate standard error of the estimate calculated by the delta method
#' * `lwr_h`, `upr_h` - The confidence bounds calculated from the Hinkley distribution.
#' * `lwr_qbn`, `upr_qbn` - The confidence bounds calculated from a quasi-Bayesian simulation using a multivariate normal distribution.
#' * `lwr_qbt`, `upr_qbt` - The confidence bounds calculated from a quasi-Bayesian simulation using a multivariate t distribution.
#' * `lwr_t`, `upr_t` - The confidence bounds calculated from the t-distribution using the delta method SE
#' * `lwr_f1`, `upr_f1`, `lwr_f2`, `upr_f2` - The confidence bounds calculated from Fieller's theorem.  If the confidence set is unbounded, both lower and both upper bounds are `NA`.  If the confidence interval is bounded or half-line, the `lwr_f2` = `upr_f2` = `NA`.
#' @export
lrm_ci <- function(obj, y_coefs, x_coefs, qb_dist = c("t", "normal"), level=.95, t_min = 2, modtype = c("adl", "ecm"), ...){
  dist <- match.arg(qb_dist, several.ok = TRUE)
  mt <- match.arg(modtype)
  delta_ratio_sd <- function(mu, Sigma) {
    # mu: vector c(mu_a, mu_b)
    # Sigma: 2x2 covariance matrix
    grad <- c(1 / mu[2], -mu[1] / (mu[2]^2))
    sqrt(t(grad) %*% Sigma %*% grad)
  }

  ll <- (1-level)/2
  ul <- 1-ll
  b <- coef(obj)
  v <- vcov(obj)
  if(is.null(names(x_coefs))){
    names(x_coefs) <- paste0("x_", seq_along(x_coefs))
  }
  if(mt == "adl"){
    if(attr(terms(obj), "intercept") == 1){
      b[1] <- 1
      v[1,] <- v[,1] <- 0
    }else{
      y_coefs <- y_coefs + 1
      x_coefs <- lapply(x_coefs, \(x)x+1)
      b <- c(1, b)
      v <- cbind(0, v)
      v <- rbind(0, v)
    }
    zv <- rep(0, length(b))
    xL <- lapply(seq_along(x_coefs), \(i){
      zv[x_coefs[[i]]] <- 1
      zv
    })
    yL <- zv
    yL[1] <- 1
    yL[y_coefs] <- -1
    yL <- matrix(yL, nrow=1)
    xL <- do.call(rbind, xL)
    L <- rbind(yL, xL)
    ests <- b %*% t(L)
    vcv <- L %*% v %*% t(L)
  }else{
    zv <- rep(0, length(b))
    xL <- lapply(seq_along(x_coefs), \(i){
      zv[x_coefs[[i]]] <- -1
      zv
    })
    yL <- zv
    yL[y_coefs] <- 1
    yL <- matrix(yL, nrow=1)
    xL <- do.call(rbind, xL)
    L <- rbind(yL, xL)
    ests <- b %*% t(L)
    vcv <- L %*% v %*% t(L)
  }
  delta_sds <- sapply(2:length(ests), \(i){
    delta_ratio_sd(c(ests[i], ests[1]), vcv[c(i, 1), c(i, 1)])
  })
  rho <- sum(b[y_coefs])
  inflation_factor = (1+rho)/(1-rho)
  if("normal" %in% dist){
    qb_ci_norm <- t(sapply(2:length(ests), \(i){
      B <- rmvnorm(10000, mean = c(ests[i], ests[1]), sigma = vcv[c(i, 1), c(i, 1)])
      good <- which(is.finite(B[,1]/B[,2]))
      if(length(good) > 5000){
        quantile(B[good[1:5000],1]/B[good[1:5000],2], probs = c(ll, ul))
      }else{
        matrix(c(NA, NA), ncol=2)
      }
    }))
    colnames(qb_ci_norm) <- c("lwr_qbn", "upr_qbn")
  }else{
    qb_ci_norm <- NULL
  }
  if("t" %in% dist){
    qb_ci_t <- t(sapply(2:length(ests), \(i){
      B <- rmvt(10000, sigma = vcv[c(i, 1), c(i, 1)], df = pmax(t_min,obj$df.residual/inflation_factor))
      B[,1] <- B[,1] + ests[i]
      B[,2] <- B[,2] + ests[1]
      good <- which(is.finite(B[,1]/B[,2]))
      if(length(good) > 5000){
        quantile(B[good[1:5000],1]/B[good[1:5000],2], probs = c(ll, ul))
      }else{
        matrix(c(NA, NA), ncol=2)
      }
    }))
    colnames(qb_ci_t) <- c("lwr_qbt", "upr_qbt")
  }else{
    qb_ci_norm <- NULL
  }
  cis <- t(sapply(2:length(ests), \(i){
    qhinkley(p=c(ll, ul), mu_a=ests[i], mu_b=ests[1], sigma_a=sqrt(vcv[i,i]), sigma_b=sqrt(vcv[1,1]), sigma_ab=vcv[i,1])
  }))
  colnames(cis) <- c("lwr_h", "upr_h")
  fcis <- t(sapply(2:length(ests), \(i){
    fi <- fieller_ci(b1=ests[i], b2=ests[1], V11=vcv[i,i], V22=vcv[1,1], V12=vcv[i,1], df=pmax(t_min,obj$df.residual/inflation_factor), alpha = 1-level)
    if(fi$type == "bounded"){
      out <- c(fi$ci, NA, NA)
    }else if(fi$type == "half-line"){
      out <- c(fi$ci, NA, NA)
    }else if(fi$type == "complement"){
      out <- c(fi$ci[[1]], fi$ci[[2]])
    }else if(fi$type == "all-real"){
      out <- c(-Inf, Inf, NA, NA)
    }else{
      out <- c(NA, NA, NA, NA)
    }
    out
  }))
  colnames(fcis) <- c("lwr_f1", "upr_f1", "lwr_f2", "upr_f2")
  out <- cbind(data.frame(vbl = names(x_coefs), est= ests[-1]/ests[1], delta_se = delta_sds), cis)
  out$lwr_t <- out$est - qt(ul, df = pmax(t_min,obj$df.residual/inflation_factor)) * delta_sds
  out$upr_t <- out$est + qt(ul, df = pmax(t_min,obj$df.residual/inflation_factor)) * delta_sds
  out <- cbind(out, qb_ci_norm, qb_ci_t, fcis)
  return(out)
}

#' Calculate Feiler Confidence Intervals
#'
#' @param b1 Estimate of the numerator.
#' @param b2 Estimate of the denominator.
#' @param V11 Variance of the numerator estimate.
#' @param V22 Variance of the denominator estimate.
#' @param V12 Covariance between the numerator and denominator estimates.
#' @param df Degrees of freedom for the t-distribution, potentially adjusted for persistence.
#' @param alpha Desired type I error rate (default is 0.05 for a 95% CI).
#' @param tol Numerical tolerance for determining special cases (default is 1e-12).
#' @param ... Additional arguments, currently not implemented.
#'
#' @details The degrees of freedom calculation could account for persistence in the data.  For an AR(1) process where $rho$ is the coefficient on y at time (t-1), the Inflation Factor is (1+rho)/(1-rho).  The adjusted degrees of freedom would be residual_df/Inflation Factor.
#' care should be taken to ensure that the adjusted degrees of freedom is not too small.  For example, setting a floor of 3 or more would work.
#'
#' @return A list with the following components:
#' * `estimate`: The point estimate of the ratio (b1/b2).
#' * `ci`: The confidence interval(s) for the ratio. This can be a vector of length 2 for a bounded interval, or a list of two vectors for a complement interval.
#' * `type`: A character string indicating the type of confidence interval: "bounded", "complement", "half-line", "all-real", or "empty".
#' @export
fieller_ci <- function(b1, b2, V11, V22, V12, df, alpha = 0.05, tol = 1e-12, ...) {
  stopifnot(is.finite(b1), is.finite(b2), is.finite(V11), is.finite(V22), is.finite(V12), df > 0)
  ccrit <- qt(1 - alpha/2, df)

  A <- b2^2 - ccrit^2 * V22
  B <- 2 * (b1 * b2 - ccrit^2 * V12)
  C <- b1^2 - ccrit^2 * V11
  D <- B^2 - 4 * A * C

  est <- b1 / b2

  # Numerical guards
  if (abs(A) < tol) A <- 0
  if (D < -tol) D <- -abs(D) # just to be explicit

  out <- list(estimate = est, A = A, B = B, C = C, disc = D, df = df, alpha = alpha)

  if (D < 0) {
    out$type <- "all-real"     # confidence set is the entire real line
    out$ci   <- c(-Inf, Inf)
    return(out)
  }

  # Compute roots
  sqrtD <- sqrt(max(D, 0))
  r1 <- (-B - sqrtD) / (2 * A)
  r2 <- (-B + sqrtD) / (2 * A)
  lo <- min(r1, r2); hi <- max(r1, r2)

  if (A > 0) {
    out$type <- "bounded"
    out$ci   <- c(lo, hi)      # (lo, hi)
  } else if (A < 0) {
    out$type <- "complement"   # (-Inf, lo] U [hi, Inf)
    out$ci <- list(c(-Inf, lo), c(hi,  Inf))
  } else {                     # A == 0: linear inequality → half-line
    # When A==0, inequality reduces to B*theta + C <= 0
    if (abs(B) < tol) {
      out$type <- if (C <= 0) "all-real" else "empty"
      out$ci   <- if (C <= 0) c(-Inf, Inf) else c(NA, NA)
    } else {
      thr <- -C / B
      if (B > 0) { out$type <- "half-line"; out$ci <- c(-Inf, thr) }
      else       { out$type <- "half-line"; out$ci <- c(thr,  Inf) }
    }
  }
  out
}


