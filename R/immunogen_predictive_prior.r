####################
# PREDICTIVE PRIOR #
####################

# One DF with three cols, (Molec, Shape1 and Shape2) which contain the shape parameters 
# alpha resp beta for each molecule in the order of molecs

# 1. get_prior_hyp :: Evaluate the function at a given hypothesis
# 2. get_pred_prior :: gets the prediction prior shape params (based on NetMHCIIPan)


########
# CODE #
########

# 1. Evaluate the function at a given hypothesis
#' P(H)
#'
#' Computes the log of the prior of a given hypothesis. Returns 0 when using a uniform prior.
#' @param molecs Vector with names of the molecules in the cohort.
#' @param h Hypothesis we want the prior of
#' @param shape.df Dataframe with shape parameters for the prior Beta distributions. If NULL, a uniform, non-informative prior is used.
get_prior_hyp <- function(molecs, h, shape.df=NULL){
  if ( is.null(shape.df) ){
    # uniform prior on [0,1]^m
    return(0)
  }
  # h is a string; molecs and h are in same order, but shape.df might not be
  else {
    h_vec <- as.double(unlist(strsplit(h, ",")))
    log_prior <- 0
    for ( i in 1:length(h_vec) ){
      log_prior <- log_prior + 
        log(dbeta(h_vec[i], shape.df[shape.df$Molec==molecs[i],]$Shape1, 
                  shape.df[shape.df$Molec==molecs[i],]$Shape2))
    }
    return (log_prior)
  }
}

# 2. gets the prediction prior shape params (based on NetMHCIIPan)
# not used?
sd_fun <- function(x, sd, m){
  b <- x/m*(1-m) + 2 - 1/m
  x*b - sd^2*(x+b)^2*(x+b+1)
}


#' Beta shape parameters 
#'
#' Given a desired mode and standard deviation, compute the shape parameters for the Beta distribution
#' @param M desired mode
#' @param SD desired standard deviation

get_shape_params <- function(M, SD){
  # M is the mode 
  # SD is the standard deviation
  alpha <- uniroot(sd_fun, sd=SD, m=M, interval=c(1,10^5))$root
  beta <- alpha*(1-M)/M -1/M +2
  return(c(alpha,beta))
}

#' Beta shape parameters for (one peptide):(all HLA) combinations
#'
#' For a given peptide, and for all HLAs, compute the shape parameters of the prior distribution for that pHLA combination.
#' @param netmhc.df Dataframe with booleans for each pHLA combination: do we think it is immunogenic or not? This can come from predicitve data (hence the name).
#' @param pep Peptide currently being processed
#' @param mode_F mode for the Beta distribution for a pHLA combination that has a negative prior (pHLA considered to be non-immunogenic)
#' @param sd_F standard deviation for a pHLA combination that has a negative prior (pHLA considered to be non-immunogenic)
#' @param mode_T mode for the Beta distribution for a pHLA combination that has a positive prior (pHLA considered to be immunogenic)
#' @param sd_T standard deviation for a pHLA combination that has a positive prior (pHLA considered to be immunogenic)

get_shape_df <- function(netmhc.df, pep, mode_F=0.3, sd_F=0.1, mode_T=0.7, sd_T=0.2){
  shape.df <- data.frame(Molec=colnames(netmhc.df)[-1])
  shapes_T <- get_shape_params(mode_T, sd_T)
  shapes_F <- get_shape_params(mode_F, sd_F)
  shape.df$Shape1 <- shapes_T[1]
  shape.df$Shape2 <- shapes_T[2]
  shape.df[which(netmhc.df[netmhc.df[,1]==pep,-1]==F),]$Shape1 <- shapes_F[1]
  shape.df[which(netmhc.df[netmhc.df[,1]==pep,-1]==F),]$Shape2 <- shapes_F[2]
  return(shape.df)
}

#' Plotting a beta distribution with a given mode and standard deviation
#'
#' @param M mode
#' @param SD standard deviation

plot_beta <- function(M, SD){
  s <- get_shape_params(M,SD)
  x <- seq(0,1,0.001)
  y <- dbeta(x,s[1],s[2])
  print(ggplot(data.frame(cbind(x=x,  y=y))) + geom_line(aes(x=x,y=y)) +
          ggtitle(paste("M = ", M, " , SD = ", SD, sep="")))
}