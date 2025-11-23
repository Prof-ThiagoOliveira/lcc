#######################################################################
# Internal helpers for longitudinal components (shared by LCC/LPC/LA) #
#######################################################################

##' @title Internal helper to precompute longitudinal quantities
##' @description Precomputes quantities (polynomial bases, tGt, deltas, etc.)
##'   that are shared across LCC, LPC and LA for a given model and time grid.
##' @keywords internal
.precompute_longitudinal <- function(model, tk, q_f, q_r, basis = NULL) {
  # Polynomial bases (reuse if provided and compatible)
  if (!is.null(basis) &&
      identical(basis$tk, tk) &&
      isTRUE(all.equal(basis$q_f, q_f)) &&
      isTRUE(all.equal(basis$q_r, q_r)) &&
      !is.null(basis$Tk_f) && !is.null(basis$Tk_r)) {
    Tk_r <- basis$Tk_r
    Tk_f <- basis$Tk_f
  } else {
    Tk_r <- outer(tk, 0:q_r, `^`)
    Tk_f <- outer(tk, 0:q_f, `^`)
  }
  
  # Random-effect variance part: tGt = diag(Tk_r %*% G %*% t(Tk_r))
  G  <- getVarCov(model)
  AG <- Tk_r %*% G
  tGt <- rowSums(AG * Tk_r)
  
  # Residual variance
  sig2_epsilon <- model$sigma^2
  
  # Variance structure and delta parameters
  varcomp   <- summary(model)
  varStruct <- varcomp$modelStruct$varStruct
  deltas    <- getDelta(model = model)
  
  list(
    tk            = tk,
    q_f           = q_f,
    q_r           = q_r,
    Tk_r          = Tk_r,
    Tk_f          = Tk_f,
    tGt           = tGt,
    sig2_epsilon  = sig2_epsilon,
    varStruct     = varStruct,
    delta         = deltas$delta,
    deltal        = deltas$deltal,
    g             = deltas$g
  )
}

##' @title Internal helper to compute LCC list for one diffbeta
##' @keywords internal
.compute_LCC <- function(pre, diffbeta) {
  Tk_f         <- pre$Tk_f
  tGt          <- pre$tGt
  sig2_epsilon <- pre$sig2_epsilon
  varStruct    <- pre$varStruct
  delta        <- pre$delta
  deltal       <- pre$deltal
  g            <- pre$g
  tk           <- pre$tk
  
  # S^2 = (X_f * diffbeta)^2
  p_len    <- length(diffbeta)
  Tk_f_sub <- Tk_f[, seq_len(p_len), drop = FALSE]
  S2       <- (Tk_f_sub %*% diffbeta)^2
  
  # Determine variance structure type
  varType <- if (is.null(varStruct)) "NULL" else class(varStruct)[1L]
  
  computeRho <- function(gd, gdl) {
    den <- tGt + 0.5 * (sig2_epsilon * (gd + gdl) + S2)
    as.numeric(tGt / den)
  }
  
  if (varType == "varIdent") {
    gd  <- g(delta)
    gdl <- g(deltal)
    rho <- lapply(gdl, function(gdlVal) computeRho(gd, gdlVal))
    
  } else if (varType == "varExp") {
    form     <- attr(varStruct, "formula")
    form_chr <- gsub(" ", "", paste(deparse(form), collapse = ""))
    
    gd <- g(delta, tk)
    
    if (form_chr == "~time") {
      gdl <- g(deltal, tk)
      rho <- list(computeRho(gd, gdl), NA)
      
    } else if (form_chr == "~time|method") {
      gdl_list <- lapply(deltal, function(dl) g(dl, tk))
      rho <- lapply(gdl_list, function(gdlVal) computeRho(gd, gdlVal))
      
    } else {
      stop("Method not implemented for the specified varStruct formula.",
           call. = FALSE)
    }
    
  } else if (varType == "NULL") {
    # No variance structure: gd = gdl = 1
    rho <- list(
      as.numeric(tGt / (tGt + sig2_epsilon + 0.5 * S2)),
      NA
    )
    
  } else {
    stop("Unsupported variance structure type: ", varType, call. = FALSE)
  }
  
  rho
}

##' @title Internal helper to compute LPC list (does not depend on diffbeta)
##' @keywords internal
.compute_LPC <- function(pre) {
  tGt          <- pre$tGt
  sig2_epsilon <- pre$sig2_epsilon
  varStruct    <- pre$varStruct
  delta        <- pre$delta
  deltal       <- pre$deltal
  g            <- pre$g
  tk           <- pre$tk
  
  varType <- if (is.null(varStruct)) "NULL" else class(varStruct)[1L]
  
  calculateRhoPearson <- function(gd, gdl) {
    den1 <- tGt + sig2_epsilon * gd
    den2 <- tGt + sig2_epsilon * gdl
    as.numeric(tGt / sqrt(den1 * den2))
  }
  
  if (varType == "varIdent") {
    gd  <- g(delta)
    gdl <- g(deltal)
    rho.pearson <- lapply(gdl, function(gdli) calculateRhoPearson(gd, gdli))
    
  } else if (varType == "varExp") {
    form     <- attr(varStruct, "formula")
    form_chr <- gsub(" ", "", paste(deparse(form), collapse = ""))
    
    if (form_chr == "~time") {
      gd  <- g(delta, tk)
      gdl <- g(deltal, tk)
      rho.pearson <- list(calculateRhoPearson(gd, gdl), NA)
      
    } else if (form_chr == "~time|method") {
      gd   <- g(delta, tk)
      gdlL <- lapply(deltal, function(d) g(d, tk))
      rho.pearson <- lapply(gdlL, function(gdli) calculateRhoPearson(gd, gdli))
      
    } else {
      print("Method not implemented yet")
      rho.pearson <- list()
    }
    
  } else if (varType == "NULL") {
    rho.pearson <- list(
      as.numeric(tGt / (tGt + sig2_epsilon)),
      NA
    )
    
  } else {
    stop("Unsupported variance structure type: ", varType, call. = FALSE)
  }
  
  rho.pearson
}


##' @title Internal helper to compute LA list for one diffbeta
##' @keywords internal
.compute_LA <- function(pre, diffbeta) {
  Tk_f         <- pre$Tk_f
  tGt          <- pre$tGt
  sig2_epsilon <- pre$sig2_epsilon
  varStruct    <- pre$varStruct
  delta        <- pre$delta
  deltal       <- pre$deltal
  g            <- pre$g
  tk           <- pre$tk
  
  # S^2 = (X_f * diffbeta)^2
  p_len    <- length(diffbeta)
  Tk_f_sub <- Tk_f[, seq_len(p_len), drop = FALSE]
  S2       <- (Tk_f_sub %*% diffbeta)^2
  sqrt_S2  <- sqrt(S2)
  
  varType <- if (is.null(varStruct)) "NULL" else class(varStruct)[1L]
  
  computeLA <- function(gd, gdl) {
    den_gd  <- tGt + sig2_epsilon * gd
    den_gdl <- tGt + sig2_epsilon * gdl
    
    v <- sqrt(den_gd / den_gdl)
    u <- sqrt_S2 / (den_gd * den_gdl)^(1/4)
    
    as.numeric(2 / (v + 1 / v + u^2))
  }
  
  if (varType == "varIdent") {
    gd  <- g(delta)
    gdl <- g(deltal)
    LA  <- lapply(gdl, function(gdlVal) computeLA(gd, gdlVal))
    
  } else if (varType == "varExp") {
    form     <- attr(varStruct, "formula")
    form_chr <- gsub(" ", "", paste(deparse(form), collapse = ""))
    
    gd <- g(delta, tk)
    
    if (form_chr == "~time") {
      gdl <- g(deltal, tk)
      LA  <- list(computeLA(gd, gdl), NA)
      
    } else if (form_chr == "~time|method") {
      gdl_list <- lapply(deltal, function(dl) g(dl, tk))
      LA <- lapply(gdl_list, function(gdlVal) computeLA(gd, gdlVal))
      
    } else {
      stop("Method not implemented for the specified varStruct formula.",
           call. = FALSE)
    }
    
  } else if (varType == "NULL") {
    LA <- list(computeLA(1, 1), NA)
    
  } else {
    stop("Unsupported variance structure type: ", varType, call. = FALSE)
  }
  
  LA
}



#######################################################################
# lccWrapper / lpcWrapper / laWrapper                                 #
#######################################################################

##' @keywords internal
lccWrapper <- function(model, q_f, tk, diffbeta, n.delta) {
  G   <- getVarCov(model)
  q_r <- nrow(G) - 1L
  
  pre <- .precompute_longitudinal(model, tk, q_f = q_f, q_r = q_r)
  rho_list <- .compute_LCC(pre, diffbeta = diffbeta)
  
  # Original selection logic: if only one element or second is NA -> use first
  if (length(rho_list) == 1L || sum(is.na(rho_list[[2L]])) != 0) {
    rho_list[[1L]]
  } else {
    rho_list[[n.delta]]
  }
}

##' @keywords internal
lpcWrapper <- function(model, q_f, tk, n.delta) {
  G   <- getVarCov(model)
  q_r <- nrow(G) - 1L
  
  pre <- .precompute_longitudinal(model, tk, q_f = q_f, q_r = q_r)
  rho_pearson_list <- .compute_LPC(pre)
  
  rho_pearson_list[[n.delta]]
}

##' @keywords internal
laWrapper <- function(model, q_f, tk, diffbeta, n.delta) {
  G   <- getVarCov(model)
  q_r <- nrow(G) - 1L
  
  pre <- .precompute_longitudinal(model, tk, q_f = q_f, q_r = q_r)
  LA_list <- .compute_LA(pre, diffbeta = diffbeta)
  
  LA_list[[n.delta]]
}

##' @title Internal Function to Estimate the Longitudinal Concordance Correlation
##' @description Thin wrapper around .precompute_longitudinal() and .compute_LCC().
##' @keywords internal
lccBuilder <- function(G, diffbeta, tk, q_r, q_f, g,
                       sig2_epsilon, delta, deltal, model) {
  # We ignore G, g, sig2_epsilon, delta, deltal here and
  # recompute everything from 'model', 'tk', 'q_f', 'q_r'
  pre <- .precompute_longitudinal(
    model = model,
    tk    = tk,
    q_f   = q_f,
    q_r   = q_r
  )
  
  .compute_LCC(pre, diffbeta = diffbeta)
}


##' @title Internal Function to Estimate the Longitudinal Pearson Correlation
##' @description Thin wrapper around .precompute_longitudinal() and .compute_LPC().
##' @keywords internal
lpcBuilder <- function(G, tk, q_r, q_f, g,
                       sig2_epsilon, delta, deltal, model) {
  # Same idea: use model/tk/q_f/q_r via the helper
  pre <- .precompute_longitudinal(
    model = model,
    tk    = tk,
    q_f   = q_f,
    q_r   = q_r
  )
  
  .compute_LPC(pre)
}


##' @title Internal Function to Estimate the Longitudinal Accuracy (Bias Corrector)
##' @description Thin wrapper around .precompute_longitudinal() and .compute_LA().
##' @keywords internal
laBuilder <- function(G, diffbeta, tk, q_r, q_f, g,
                      sig2_epsilon, delta, deltal, model) {
  pre <- .precompute_longitudinal(
    model = model,
    tk    = tk,
    q_f   = q_f,
    q_r   = q_r
  )
  
  .compute_LA(pre, diffbeta = diffbeta)
}


#######################################################################
