#######################
# METROPOLIS-HASTINGS #
#######################

## sample through the population
# 1. get_proposal_unif :: find a hypothesis we might jump to
# 2. get_q_h1_given_h2 :: get Q(h1|h2)
# 3. one_time_step :: exactly what it says :: with evidence_1
# 4. mh_chain :: initialisation + running of mh algorithm :: with evidence_1
# 5. get_acc_rate_chain :: Get the acceptance rate of a chain; used to tune the radius

########
# CODE #
########


# 1. get_proposal :: find a hypothesis we might jump to
#' Choosing a proposal state
#'
#' Given the current state of the MH chain, pick a proposal state. Also returns Q(proposal state|current state)
#' @param h current State of the MH chain
#' @param unif.prop If TRUE, proposal is chosen from [0,1]^m. If false, only looks in the overlap of [0,1]^m with a cube around h with side length 2*radius
#' @param radius If unif.prop is FALSE, we look for a proposal state in the overlap of [0.1]^m and [h-radius, h+radius]. 

get_proposal <- function(h, unif.prop=T, radius){
  # output: vector
  # first el: the proposal hyp
  if ( unif.prop ){
    # just take a random point in unit hypercube
    return(c(paste(runif(length(unlist(strsplit(h, ",")))), collapse=","),0))
  }
  else {  # non-independent sampling
    h <- as.double(unlist(strsplit(h, ",")))
    lower <- pmax(0, h - radius)
    upper <- pmin(1, h + radius)
    g <- lower + unlist(lapply(rep(1, length(h)), runif))*(upper-lower)
    return(c(paste(g, collapse=","), -sum(log(upper-lower))))
  }
}

# 2. get_q_h1_given_h2 :: get Q(h1|h2)
#' Evaluate Q function

get_q_h1_given_h2 <- function(h1, h2, unif.prop=T, radius){
  if ( unif.prop ){
    return(0)
  }
  else{
    h2 <- as.double(unlist(strsplit(h2, ",")))
    -sum(log( pmin(1, h2 + radius) -  pmax(0, h2 - radius) ))
  }
}

# 3. one_time_step :: exactly what it says :: with evidence_1
#' One update of MH chain
#'
#' Performs one time step in the MH chain: 
#' @param h Hypothesis which is the current state of the chain 
#' @param loglik_h Log likelihood of h
#' @param eli.dat ELISpot dataset
#' @param pep Peptide currently being processed
#' @param mol.names Names of the HLA molecules in the ELISpot dataset
#' @param unif.prop If TRUE, choose proposal state from [0,1]^m
#' @param radius If unif.prop is FALSE, we look for a proposal state in a hypercube around h restricted by radius
#' @param p.df Dataframe which describes the parameters of the prior distributions for each HLA

one_time_step <- function(h, loglik_h, eli.dat, pep, mol.names, unif.prop, 
                          radius, p.df){
  ## get a proposal for next hypothesis
  new_vec <- get_proposal(h, unif.prop=unif.prop, radius)
  prop_h <- new_vec[1]
  q_new_given_old <- as.double(new_vec[2])
  # this is in log
  loglik_prop_h <- get_evidence_1(eli.dat, prop_h, mol.names, pep) +
                 get_prior_hyp(molecs, prop_h, shape.df=p.df)

  q_old_given_new <- get_q_h1_given_h2(h, prop_h, unif.prop=unif.prop, radius)
  # this is in log
  acc_rate <- min(1, exp(- loglik_h + loglik_prop_h)*exp(q_old_given_new - q_new_given_old) )
  if ( runif(1) < acc_rate ) {
    return(c(prop_h, loglik_prop_h))
  }
  else {
    return(c(h, loglik_h))
  }
}

# 4. mh_chain :: initialisation + running of mh algorithm :: with evidence_1
#' Metropolis-Hasting for a single peptide
#'
#' Implementation of the Metropolis-Hastings algorithm to find posterior distributions of pHLA immunogenicity for a single peptide
#' @param eli.dat ELISpot dataset
#' @param mol.names Names of the HLA molecules in the ELISpot dataset
#' @param init_h Initial state of the MH chain, can be produced randomly by generate_random_hypothesis
#' @param max_steps The required length of the chain. The default value, 5000, is too low but this was done for debugging reasons.
#' @param pep Peptide currently being processed
#' @param unif.prop If TRUE, choose proposal states from [0,1]^m
#' @param radius If unif.prop is FALSE, we look for a proposal state in a hypercube around h restricted by radius. Should be trained such that the acceptance rate of the chain is around 50%, see example.r.
#' @param p.df Dataframe which describes the parameters of the prior distributions for each HLA


mh_chain <- function(eli.dat, mol.names, init_h, max_steps=5000, pep, unif.prop=T, 
                     radius=0.1, p.df=NULL){
  mh_out <- data.frame(matrix(nrow=max_steps, ncol=2))
  colnames(mh_out) <- c("Hyp", "LogLik")
  mh_out [1,] <- c(init_h, get_evidence_1(eli.dat, init_h, mol.names, pep) +
                     get_prior_hyp(molecs, init_h, shape.df=p.df))
  t <- 2
  while ( ! t > max_steps){
    out_t <- one_time_step(h=mh_out[t-1,1], 
                           loglik_h=as.double(as.character(mh_out[t-1,2])), 
                           eli.dat, pep, mol.names, unif.prop, radius, p.df)  
    mh_out[t,] <- c(paste(out_t[1], collapse=","), out_t[2])
    t <- t + 1
  }
  mh_out$LogLik <- as.double(as.character(mh_out$LogLik))
  mh_out
}

# 5. Get the acceptance rate of a chain; used to tune the radius
#' Acceptance Rate of an MH chain
#'
#' Computes the acceptance rate of a chain by comparing each state to the next one.
#' @param chain Which chain do you want to compute the acceptance rate of?
get_acc_rate_chain <- function(chain){
  100*sum(!chain[1:(dim(chain)[1]-1),]$Hyp == chain[2:dim(chain)[1],]$Hyp)/dim(chain)[1]
}
