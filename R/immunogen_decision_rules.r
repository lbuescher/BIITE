##################
# Decision rules #
##################


# 1. get_mode :: wrapper for modeest, really
# 2a. get_modes_chain :: modes for posterior distributions of each pHLA combo
# 2b. get_meanss_chain :: same but means
# 3. get_info_df_single_chain :: Modes and means for a single peptide in a single dataframe
# 4. get_KL_vs_unif :: KL divergence of a sample wrt sample from a distribution
# 5. get_KL_vs_prior  :: KL divergence of a sample wrt sample from prior beta distribution
#    not used 
# 6. get_tpr_fpr :: Get TPR/FPR at different cut-offs
# 7. get_auc :: From TPR/FPR dataframe, get the AUC
# 8. 

library(modeest)
library(entropy)

# 1. wrapper for modeest, really
#' Estimate sample mode
#'
#' This function is a wrapper for the mlv function from the modeest package
#' @param sample Sample to get an estimated mode of

get_mode <- function(sample){
  modeest::mlv(sample, method = "hrm", bw=0.5)
}

# 2. for a peptide, get all modes/means (one for each HLA)
# 2a. get mode
#' Estimate posterior modes in the MH chain
#'
#' Returns a dataframe with the posterior mode of the pHLA immunogenicity, for each HLA.
#' @param chain Metropolis-Hastings chain, as obtained by mh_chain
#' @param molecs vector containing all the HLA-II molecules in your ELISpot subjects

get_modes_chain <- function(chain, molecs){
  out.df <- data.frame(matrix(nrow=length(molecs),ncol=2))
  colnames(out.df) <- c("Molec", "Mode")
  out.df$Molec <- molecs
  for ( i in 1:dim(out.df)[1] ){
    out.df$Mode[i] <- (get_mode(as.numeric(chain[,colnames(chain)==molecs[i]]))$M)
  }
  return(out.df)
}

# 2b. get means of a single chain
#' Estimate posterior modes in the MH chain
#'
#' Returns a dataframe with the posterior mean of the pHLA immunogenicity, for each HLA.
#' @param chain Metropolis-Hastings chain, as obtained by mh_chain
#' @param molecs vector containing all the HLA-II molecules in your ELISpot subjects

get_means_chain <- function(chain, molecs){
  out.df <- data.frame(matrix(nrow=length(molecs),ncol=2))
  colnames(out.df) <- c("Molec", "Mean")
  out.df$Molec <- molecs
  for ( i in 1:dim(out.df)[1] ){
    out.df$Mean[i] <- mean(as.numeric(chain[,colnames(chain)==molecs[i]]))
  }
  return(out.df)
}


# 3. Modes and means for a single peptide in a single dataframe
get_info_df_single_chain <- function(chain, molecs){
  # get the modes
  modes.df <- get_modes_chain(chain, molecs)
  # get the means
  modes.df <- merge(modes.df, get_means_chain(chain, molecs), by="Molec")
  # difference between mode and mean
  modes.df$Diff_mode_mean <- abs(modes.df$Mode - modes.df$Mean)
  return(modes.df)
}

# 4. KL divergence of a sample wrt sample from a distribution
#' Kullback-Leibler divergence of a sample (with respect to the uniform distribution)
#'
#' Returns the KL divergence of a sample wrt the uniform distribution
#' @param sample Sample to compute KL divergence of (in a vector)

get_KL_vs_unif <- function(sample){
  x1 <- density(sample, from=0, to=1)$y
  x2 <- rep(1/length(x1), length(x1))
  entropy::KL.plugin(x1, x2)
}

# 5. KL divergence of a sample wrt sample from prior beta distribution
# not used 
get_KL_vs_prior <- function(sample, M , SD){
  # get the shape parameters needed in order to get a beta with mean M and sd SD
  s <- get_shape_params(M, SD)
  # get density of the sample
  d <- density(sample, from=0, to=1)
  # get y values of this density
  x1 <- d$y
  # compute the density of a beta at the same points
  x2 <- dbeta(d$x,s[1],s[2])
  # get KL
  KL.plugin(x1[!x2==0],x2[!x2==0])
}

# 6. Get TPR/FPR at different cut-offs
get_tpr_fpr <- function(df, against="tg", dont_order=F){
  # test against would be either ba, tg or pred
  # dont_order = T if df is already in the order we want it to be
  if (dont_order==F) { df <- df[order(df$Value, decreasing=F),] }
  df <- df[!is.na(df[,against]),]
  TP <- c()
  TN <- c()
  FP <- c()
  FN <- c()
  for ( i in 1:((dim(df)[1]-1)) ){
    # count FP, TP, FN, TN if we put the threshold on Mode[i] (assumed you're working with modes)
    TP[i] <- sum(df[(i+1):dim(df)[1],against]==1)
    FP[i] <- sum(df[(i+1):dim(df)[1],against]==0)
    TN[i] <- sum(df[1:i,against]==0)
    FN[i] <- sum(df[1:i,against]==1)
  }
  TP <- c(TP,0)
  FP <- c(FP,0)
  TN <- c(TN, sum(df[,against]==0))
  FN <- c(FN, sum(df[,against]==1))
  
  TPR <- c(1,TP/(TP + FN))
  FPR <- c(1,FP/(TN + FP))
  
  TPR <- rev(TPR)
  FPR <- rev(FPR)
  
  return(data.frame(TPR, FPR))
}

# 7. From TPR/FPR dataframe, get the AUC
get_auc <- function(r.df){
  r2.df <-  r.df[!duplicated(r.df[,-1]),]
  dx <- r2.df$FPR[2:dim(r2.df)[1]] - r2.df$FPR[1:(dim(r2.df)[1]-1)]
  dy <- r2.df$TPR[2:dim(r2.df)[1]] - r2.df$TPR[1:(dim(r2.df)[1]-1)]
  auc <- sum(dx*dy)/2 + sum(dx*r2.df$TPR[2:dim(r2.df)[1]])
  return(auc)
}

# 8. Create dataframe with posterior modes, means, medians
#' Overview of posterior modes, means, medians for all pHLA combinations
#'
#' Returns a dataframe with the marginal posterior mode, mean and median of the pHLA immunogenicity, for each HLA.
#' Mode, mean and median are all on a separate line, so each pHLA combination occurs three times
#' Also gives the KL divergence of the marginal posterior with respect to the uniform distribution
#' @param peps_for_analysis Which peptides do you want to include
#' @param chainDir directory that contains the MH chains 
#' @param molecs A vector listing the molecules that are present in the ELISpot subjects

get_overview_df <- function(peps_for_analysis, chainDir, molecs){
  res <- data.frame(matrix(nrow=length(peps_for_analysis)*length(molecs)*3, ncol=5))
  colnames(res) <- c("Pep", "Molec", "What", "Value", "DKL")
  r <- 1
  for ( pep in peps_for_analysis ){
    try( ch <- read.csv(paste(chainDir, pep, "_full_chain.txt", sep=""), skip=1, header=F, sep=",") )
    if (exists("ch")){
      # load molecs for the peptide
      molecs_pep <- read.csv(paste0(chainDir, "molecs_", pep, ".txt"))
      colnames(ch) <- c(molecs_pep, "LL")
      for ( molec in molecs_pep ){
        res$Molec[r] <- molec
        res$Pep[r] <- pep
        res$What[r] <- "Mean"
        res$Value[r] <-  mean(ch[,molec])
        res$DKL[r] <- get_KL_vs_unif(ch[,molec])
        r <- r + 1
        
        res$Molec[r] <- molec
        res$Pep[r] <- pep
        res$What[r] <- "Mode"
        res$Value[r] <- get_mode(ch[,molec])$M
        res$DKL[r] <- get_KL_vs_unif(ch[,molec])
        r <- r + 1
        
        res$Molec[r] <- molec
        res$Pep[r] <- pep
        res$What[r] <- "Median"
        res$Value[r] <- median(ch[,molec])
        res$DKL[r] <- get_KL_vs_unif(ch[,molec])
        r <- r + 1
      }            
      rm(ch)
      rm(molecs_pep)
    }
  }
  res <- res[!is.na(res$Pep),]
  return(res)
}

# 9. Find 'explanatory' pHLA combinations (done for each peptide separately)
#' Find pHLA combinations that explain positive ELISpot results in the cohort.
#'
#' For each peptide, we rank the HLAs according to their posterior mode. 
#' For each subject, we calculate which of their HLAs is the highest ranked one. 
#' In return, for each HLA, we count in how many subjects this HLA is present and in how many subjects it is the highest ranked HLA.
#' This calculation is then repeated, taking into account only subjets with positive ELISpot for the peptide under consideration.
#' This allows the user to see which HLA 'explains' the positive ELISpot results.
#' @param peps_for_analysis Which peptides do you want to include
#' @param res The results table as generated by get_overview_df
#' @param molecs A vector listing the molecules that are present in the ELISpot subjects
#' @param elidat The original patient and ELISpot data
#' @param what The summary value used to rank the HLA. Defaults to "Mode", but can also be "Mean" or "Median"

get_hla_ranking <- function(peps_for_analysis, res, molecs, eli.dat, what="Mode"){
  for (pep in peps_for_analysis){
    eli.copy <- eli.dat[!is.na(eli.dat[,pep]) & as.character(eli.dat[,pep]) %in% c("0", "1", "TRUE", "FALSE"),]
    r1 <- res[res$What==what & res$Pep==pep,]
    r1 <- r1[order(r1$Value, decreasing=T),]
    r1$n <- unlist(lapply(as.character(r1$Molec), function(x){sum(eli.copy[,x]>0)}))
    subs <- eli.copy[,colnames(eli.copy) %in% molecs]
    eli.copy$strongest <- unlist(lapply(1:dim(eli.copy)[1],function(x){
      molec_person <- colnames(subs)[(subs[x,]>0)]
      as.character(r1[r1$Molec %in% molec_person,]$Molec[1])
    }))
    r1$strongest <- unlist(lapply(r1$Molec, function(x){
      sum(eli.copy[eli.copy[,colnames(eli.copy)==as.character(x)]>0,]$strongest==as.character(x))
    }))
    r1$npos <- unlist(lapply(r1$Molec, function(x){
      sum(eli.copy[eli.copy[,colnames(eli.copy)==as.character(x)]>0,pep]==T)
   }))
    r1$posstrongest <- unlist(lapply(r1$Molec, function(x){
      sum(eli.copy[eli.copy[,colnames(eli.copy)==as.character(x)]>0 & eli.copy[,pep]==T,]$strongest==as.character(x))
    }))
    r1 <- r1[,!colnames(r1)=="What"]
    colnames(r1) <- c("Pep",	"Molec",	"Posterior Mode",	"DKL",	"# Carriers",	"# Carriers with highest rank	#",
    "# Carriers with positive ELISpot",	"# Carriers with positive ELISpot and highest rank")

    write.csv(r1, paste0("overview_table_", pep, ".csv"))
  }
}
