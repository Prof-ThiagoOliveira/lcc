#######################################################################
#                                                                     #
# Package: lcc                                                        #
#                                                                     #
# File: summary_outputs.R                                             #
# Contains: lccSummary and fittedBuilder functions                    #
#                                                                     #
# Written by Thiago de Paula Oliveira                                 #
# copyright (c) 2017-18, Thiago P. Oliveira                           #
#                                                                     #
# First version: 11/10/2017                                           #
# Last update: 22/11/2025                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

##' @title Internal Function to Summarize Fitted and Sampled Values for \code{lcc} Objects
##'
##' @description Internally called function for summarizing fitted and sampled values, and the 
##'   concordance correlation coefficient between them for \code{lcc} objects.
##'
##' @usage NULL
##'
##' @details Returns a summary of fitted and sampled values and their concordance correlation.
##'
##' @author Thiago de Paula Oliveira, \email{thiago.paula.oliveira@@alumni.usp.br}
##'
##' @importFrom stats predict
##'
##' @keywords internal
lccSummary <- function(model, q_f, diffbeta, tk,
                       tk.plot, tk.plot2, rho, ENV.LCC,
                       rho.pearson, ENV.LPC, Cb, ENV.Cb,
                       ldb, ci, components) {
  
  # -------------------------------------------------------------------
  # Common pieces
  # -------------------------------------------------------------------
  method_levels <- levels(model$data$method)
  
  # Comparison labels: string for ldb == 1, list for ldb > 1
  if (ldb == 1L) {
    comp <- paste0(method_levels[2L], " vs. ", method_levels[1L])
  } else {
    comp <- vector("list", ldb)
    for (i in seq_len(ldb)) {
      comp[[i]] <- paste0(method_levels[i + 1L], " vs. ", method_levels[1L])
    }
  }
  
  # Goodness of fit: use CCC() function (not CCC_lin result)
  GF <- CCC(stats::predict(model), model$data$resp)
  
  # -------------------------------------------------------------------
  # No components: only LCC
  # -------------------------------------------------------------------
  if (!components) {
    # Longitudinal CCC by time (from data)
    CCC_vals <- CCC_lin(
      dataset = model$data,
      resp    = "resp",
      subject = "subject",
      method  = "method",
      time    = "time"
    )
    
    if (!ci) {
      # -------------------------------
      # components == FALSE, ci == FALSE
      # -------------------------------
      if (ldb == 1L) {
        LCC.data <- data.frame(
          "Time" = tk.plot,
          "LCC"  = rho
        )
        
        CCC.data <- data.frame(
          "Time" = tk.plot2,
          "CCC"  = CCC_vals
        )
        colnames(CCC.data) <- c("Time", "CCC")
        
        plot.data <- list(
          "fitted"  = LCC.data,
          "sampled" = CCC.data,
          "gof"     = GF,
          "comp"    = comp
        )
      } else {
        # ldb > 1: LCC per method, one CCC curve
        LCC.data <- vector("list", ldb)
        for (i in seq_len(ldb)) {
          LCC.data[[i]] <- data.frame(
            "Time" = tk.plot,
            "LCC"  = rho[[i]]
          )
        }
        
        CCC.data <- data.frame(
          "Time" = tk.plot2,
          "CCC"  = CCC_vals
        )
        colnames(CCC.data) <- c("Time", "CCC")
        
        plot.data <- list(
          "fitted"  = LCC.data,
          "sampled" = CCC.data,
          "gof"     = GF,
          "comp"    = comp
        )
      }
      
    } else {
      # -------------------------------
      # components == FALSE, ci == TRUE
      # -------------------------------
      if (ldb == 1L) {
        LCC.data <- data.frame(
          "Time"  = tk.plot,
          "LCC"   = rho,
          "Lower" = ENV.LCC[1, ],
          "Upper" = ENV.LCC[2, ]
        )
        
        CCC.data <- data.frame(
          "Time" = tk.plot2,
          "CCC"  = CCC_vals
        )
        colnames(CCC.data) <- c("Time", "CCC")
        
        plot.data <- list(
          "fitted"  = LCC.data,
          "sampled" = CCC.data,
          "gof"     = GF,
          "comp"    = comp
        )
      } else {
        LCC.data <- vector("list", ldb)
        for (i in seq_len(ldb)) {
          LCC.data[[i]] <- data.frame(
            "Time"  = tk.plot,
            "LCC"   = rho[[i]],
            "Lower" = ENV.LCC[[i]][1, ],
            "Upper" = ENV.LCC[[i]][2, ]
          )
        }
        
        CCC.data <- data.frame(
          "Time" = tk.plot2,
          "CCC"  = CCC_vals
        )
        colnames(CCC.data) <- c("Time", "CCC")
        
        plot.data <- list(
          "fitted"  = LCC.data,
          "sampled" = CCC.data,
          "gof"     = GF,
          "comp"    = comp
        )
      }
    }
    
  } else {
    # -----------------------------------------------------------------
    # components == TRUE: LCC, LPC, LA
    # -----------------------------------------------------------------
    CCC_vals <- CCC_lin(
      dataset = model$data,
      resp    = "resp",
      subject = "subject",
      method  = "method",
      time    = "time"
    )
    Pearson_vals <- Pearson(
      dataset = model$data,
      resp    = "resp",
      subject = "subject",
      method  = "method",
      time    = "time"
    )
    
    if (!ci) {
      # -------------------------------
      # components == TRUE, ci == FALSE
      # -------------------------------
      if (ldb == 1L) {
        LA <- CCC_vals[[1L]] / Pearson_vals[[1L]]
        
        LCC.data <- data.frame(
          "Time" = tk.plot,
          "LCC"  = rho,
          "LPC"  = rho.pearson,
          "LA"   = Cb
        )
        
        CCC.data <- data.frame(
          "Time"    = tk.plot2,
          "CCC"     = CCC_vals,
          "Pearson" = Pearson_vals,
          "Cb"      = LA
        )
        colnames(CCC.data) <- c("Time", "CCC", "Pearson", "Cb")
        
        plot.data <- list(
          "fitted"  = LCC.data,
          "sampled" = CCC.data,
          "gof"     = GF,
          "comp"    = comp
        )
      } else {
        LCC.data <- vector("list", ldb)
        CCC.data <- vector("list", ldb)
        LA       <- vector("list", ldb)
        
        for (i in seq_len(ldb)) {
          LA[[i]] <- CCC_vals[[i]] / Pearson_vals[[i]]
          
          LCC.data[[i]] <- data.frame(
            "Time" = tk.plot,
            "LCC"  = rho[[i]],
            "LPC"  = rho.pearson[[i]],
            "LA"   = Cb[[i]]
          )
          
          CCC.data[[i]] <- data.frame(
            "Time"    = tk.plot2,
            "CCC"     = CCC_vals[[i]],
            "Pearson" = Pearson_vals[[i]],
            "Cb"      = LA[[i]]
          )
          colnames(CCC.data[[i]]) <- c("Time", "CCC", "Pearson", "Cb")
        }
        
        plot.data <- list(
          "fitted"  = LCC.data,
          "sampled" = CCC.data,
          "gof"     = GF,
          "comp"    = comp
        )
      }
      
    } else {
      # -------------------------------
      # components == TRUE, ci == TRUE
      # -------------------------------
      if (ldb == 1L) {
        LA <- CCC_vals[[1L]] / Pearson_vals[[1L]]
        
        LCC.data <- data.frame(
          "Time"  = tk.plot,
          "LCC"   = rho,
          "Lower" = ENV.LCC[1, ],
          "Upper" = ENV.LCC[2, ]
        )
        
        LPC.data <- data.frame(
          "Time"  = tk.plot,
          "LPC"   = rho.pearson,
          "Lower" = ENV.LPC[1, ],
          "Upper" = ENV.LPC[2, ]
        )
        
        LA.data <- data.frame(
          "Time"  = tk.plot,
          "LA"    = Cb,
          "Lower" = ENV.Cb[1, ],
          "Upper" = ENV.Cb[2, ]
        )
        
        CCC.data <- data.frame(
          "Time"    = tk.plot2,
          "CCC"     = CCC_vals,
          "Pearson" = Pearson_vals,
          "Cb"      = LA
        )
        colnames(CCC.data) <- c("Time", "CCC", "Pearson", "Cb")
        
        fit <- list(
          "LCC" = LCC.data,
          "LPC" = LPC.data,
          "LA"  = LA.data
        )
        
        plot.data <- list(
          "fitted"  = fit,
          "sampled" = CCC.data,
          "gof"     = GF,
          "comp"    = comp
        )
        
      } else {
        LCC.data <- vector("list", ldb)
        LPC.data <- vector("list", ldb)
        LA.data  <- vector("list", ldb)
        CCC.data <- vector("list", ldb)
        LA       <- vector("list", ldb)
        
        for (i in seq_len(ldb)) {
          LA[[i]] <- CCC_vals[[i]] / Pearson_vals[[i]]
          
          LCC.data[[i]] <- data.frame(
            "Time"  = tk.plot,
            "LCC"   = rho[[i]],
            "Lower" = ENV.LCC[[i]][1, ],
            "Upper" = ENV.LCC[[i]][2, ]
          )
          
          LPC.data[[i]] <- data.frame(
            "Time"  = tk.plot,
            "LPC"   = rho.pearson[[i]],
            "Lower" = ENV.LPC[[i]][1, ],
            "Upper" = ENV.LPC[[i]][2, ]
          )
          
          LA.data[[i]] <- data.frame(
            "Time"  = tk.plot,
            "LA"    = Cb[[i]],
            "Lower" = ENV.Cb[[i]][1, ],
            "Upper" = ENV.Cb[[i]][2, ]
          )
          
          CCC.data[[i]] <- data.frame(
            "Time"    = tk.plot2,
            "CCC"     = CCC_vals[[i]],
            "Pearson" = Pearson_vals[[i]],
            "LA"      = LA[[i]]
          )
          colnames(CCC.data[[i]]) <- c("Time", "CCC", "Pearson", "Cb")
        }
        
        fit <- list(
          "LCC" = LCC.data,
          "LPC" = LPC.data,
          "LA"  = LA.data
        )
        
        plot.data <- list(
          "fitted"  = fit,
          "sampled" = CCC.data,
          "gof"     = GF,
          "comp"    = comp
        )
      }
    }
  }
  
  invisible(plot.data)
}

#######################################################################
#                                                                     #
# Package: lcc                                                        #
#                                                                     #
# Section: fittedBuilder                                              #
#                                                                     #
# Written by Thiago de Paula Oliveira                                 #
# copyright (c) 2017-18, Thiago P. Oliveira                           #
#                                                                     #
# First version: 11/10/2017                                           #
# Last update: 22/11/2025                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

##' @title Internal Function to Build Fitted Values for
##'   \code{lcc} Objects
##'
##' @description This is an internally called function used to build
##'   fitted values.
##'
##' @usage NULL
##'
##' @author Thiago de Paula Oliveira, \email{thiago.paula.oliveira@@alumni.usp.br}
##'
##' @keywords internal
fittedBuilder <- function(object, type) {
  # Map type -> index and column name
  .form <- switch(type,
                  "lcc" = 1L,
                  "lpc" = 2L,
                  "la"  = 3L)
  .name <- switch(type,
                  "lcc" = "LCC",
                  "lpc" = "LPC",
                  "la"  = "LA")
  
  comp   <- object$Summary.lcc$comp
  fitted <- object$Summary.lcc$fitted
  
  #------------------------------------------------------------------
  # Single comparison (ldb == 1): comp is a character vector
  #------------------------------------------------------------------
  if (inherits(comp, "character")) {
    # fitted is either:
    #  - data.frame (no components or ci==FALSE)
    #  - list(LCC = df, LPC = df, LA = df) when components=TRUE & ci=TRUE
    df <- if (inherits(fitted, "data.frame")) {
      fitted
    } else {
      fitted[[.form]]
    }
    
    ret <- data.frame(
      Methods = comp,
      Time    = df[, "Time"],
      LCC     = df[, .name]
    )
    
  } else {
    #----------------------------------------------------------------
    # Multiple comparisons (ldb > 1): comp is a list of labels
    #----------------------------------------------------------------
    # 'fitted' is either:
    #  - list of data.frames (no components or ci==FALSE)
    #  - list(LCC = list(df_i), LPC = list(df_i), LA = list(df_i))
    #    when components=TRUE & ci=TRUE
    if (is.null(fitted$LCC)) {
      # No named components: list of data.frames for this type
      df_list <- fitted
    } else {
      # Named components: pick the list corresponding to this type
      df_list <- fitted[[.form]]
    }
    
    ncomp <- length(comp)
    fit   <- vector("list", ncomp)
    rn    <- vector("list", ncomp)
    met   <- vector("list", ncomp)
    
    for (i in seq_len(ncomp)) {
      df_i <- df_list[[i]]
      
      fit[[i]] <- data.frame(LCC = df_i[, .name])
      rn[[i]]  <- data.frame(Time = df_i[, "Time"])
      met[[i]] <- data.frame(
        Methods = rep(comp[[i]], nrow(df_i))
      )
    }
    
    ret <- data.frame(
      do.call(rbind.data.frame, met),
      do.call(rbind.data.frame, rn),
      do.call(rbind.data.frame, fit)
    )
  }
  
  # Keep the final column name convention identical to original
  colnames(ret)[colnames(ret) == "LCC"] <- paste0("fitted.", .name)
  
  ret
}

