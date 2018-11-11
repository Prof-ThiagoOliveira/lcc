#######################################################################
#                                                                     #
# Package: lcc                                                        #
#                                                                     #
# File: dataBuilder.R                                                #
# Contains: dataBuilder function                                     #
#                                                                     #
# Written by Thiago de Paula Oliveira                                 #
# copyright (c) 2017-18, Thiago P. Oliveira                           #
#                                                                     #
# First version: 11/10/2017                                           #
# Last update: 18/06/2018                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################
##' @title Internal function to prepare the dataset for \code{lcc} objects
##'
##' @description This is an internally called function used to prepare the dataset for \code{lcc} objects
##'
##' @usage NULL
##'
##' @author Thiago de Paula Oliveira, \email{thiago.paula.oliveira@@usp.br}
##'
##' @keywords internal
dataBuilder <- function(dataset, resp, subject, method, time){
  Data <- data.frame(dataset)
  Data <- try(rename.vars(Data, from = c(resp, subject, method, time),
                      to = c("y", "ind", "FacA", "time"),
                      info = FALSE), TRUE)
  Data$ind <- as.factor(Data$ind)
  Data$FacA <- as.factor(Data$FacA)
  Data$time <- as.numeric(Data$time)
  Data$y <- as.numeric(Data$y)
  return(Data)
}
