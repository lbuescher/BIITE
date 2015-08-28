##############
# HYPOTHESIS #
##############

#' Generate a random hypothesis
#'
#' Returns a random hypothesis, used to initiate the MH chain. Based on
#' @param m Number of molecules in ELISpot population


# 1. generate random hypotheses
generate_random_hypothesis <- function(m, upper_lim=1){
  # return as string
  return(paste(runif(m, 0, upper_lim), collapse=","))
}