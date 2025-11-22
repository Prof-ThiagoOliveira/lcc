#######################################################################
#                                                                     #
# Package: lcc                                                        #
#                                                                     #
# File: fittedBuilder.R                                               #
# Contains: fittedBuilder function                                    #
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

