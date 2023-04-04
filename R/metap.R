# copy from metap package: https://CRAN.R-project.org/package=metap
#' @importFrom stats na.fail
#' @importFrom stats pchisq
#' @importFrom stats pnorm
#' @importFrom stats qnorm
sumz <- function(p, weights = NULL, data = NULL, subset = NULL,
   na.action = na.fail, log.p = FALSE, log.input = FALSE)  {
   if(is.null(data)) data <- sys.frame(sys.parent())
   mf <- match.call()
   mf$data <- NULL
   mf$subset <- NULL
   mf$na.action <- NULL
   mf[[1]] <- as.name("data.frame")
   mf <- eval(mf, data)
   if(!is.null(subset)) mf <- mf[subset,]
   mf <- na.action(mf)
   p <- as.numeric(mf$p)
   weights <- mf$weights
   noweights <- is.null(weights)
   if(noweights) weights <- rep(1, length(p))
   if(length(p) != length(weights)) warning("Length of p and weights differ")
   if(log.input) {
      keep <- p < 0
   } else {
      keep <- (p > 0) & (p < 1)
   }
   invalid <- sum(1L * keep) < 2
   if(invalid) {
      warning("Must have at least two valid p values")
      res <- list(z = NA_real_, p = NA_real_,
         validp = p[keep], weights = weights)
   } else {
      if(sum(1L * keep) != length(p)) {
         warning("Some studies omitted")
         omitw <- weights[!keep]
         if((sum(1L * omitw) > 0) & !noweights)
            warning("Weights omitted too")
      }
      zp <- (qnorm(p[keep], lower.tail = FALSE, log.p = log.input) %*%
         weights[keep]) / sqrt(sum(weights[keep]^2))
      res <- list(z = zp, p = pnorm(zp, lower.tail = FALSE,
            log.p = log.p),
         validp = p[keep], weights = weights)
   }
   class(res) <- c("sumz", "metap")
   res
}

# copy from metap package: https://CRAN.R-project.org/package=metap
print.sumz <- function(x, ...) {
   cat("sumz = ", x$z, "p = ", x$p, "\n")
   invisible(x)
}

# copy from metap package: https://CRAN.R-project.org/package=metap
sumlog <- function(p, log.p = FALSE) {
   keep <- (p > 0) & (p <= 1)
   invalid <- sum(1L * keep) < 2
   if(invalid) {
      warning("Must have at least two valid p values")
      res <- list(chisq = NA_real_, df = NA_integer_,
         p = NA_real_, validp = p[keep])
   } else {
      lnp <- log(p[keep])
      chisq <- (-2) * sum(lnp)
      df <- 2 * length(lnp)
      if(length(lnp) != length(p)) {
         warning("Some studies omitted")
      }
      res <- list(chisq = chisq, df = df,
         p = pchisq(chisq, df, lower.tail = FALSE,
            log.p = log.p), validp = p[keep])
    }
   class(res) <- c("sumlog", "metap")
   res
}

# copy from metap package: https://CRAN.R-project.org/package=metap
print.sumlog <- function(x, ...) {
   cat("chisq = ", x$chisq, " with df = ", x$df, " p = ", x$p, "\n")
   invisible(x)
}


