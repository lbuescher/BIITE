library(BIITE)

################
# User-defined #
################

eliFileName <- "./Data/human_elispot.txt"
# location of ELISpot data

pep_names <- paste("pep", 1:17, sep="")
# name of your peptides

max_steps <- 1000#00
# how long do you want your MH chains to be?

outDir <- "/media/lboelen/MyPassport/MHCII/CODE/"
# what is the output directory?
# make sure it ends in "/" - or "\\" for windows

peps_for_analysis <- pep_names
# change this if you don't want to run all peptides

use_prior <- FALSE
# change to TRUE if you want a prior

predFileName <- NULL
# predFileName <- "./Data/netmhciipan_pred.txt"
# change if you want to use predictions

print_loglik_evol <- TRUE
# change to no if you don't want a plot of the evolution of the likelihood with each time step

molecs <- NULL
# change this if the code below in 3. doesn't grep your names for the molecules; or change it there!

########
# CODE #
########


# 1. import predictions if using a prior
if ( use_prior == T){
  pred_data <- read.table(predFileName, header=T)
}  

# 2. import elispot data
eli.dat <- read.table(eliFileName, sep="\t", header=T)   

# 3. extract names of the molecules in the population :: molecs
molecs <- colnames(eli.dat)[(grepl("DRB1_", colnames(eli.dat)) | grepl("DQB1_", colnames(eli.dat))) & !grepl("al", colnames(eli.dat))]
write.table(molecs, paste(outDir, "molecs.txt"), row.names=F, col.names=F, sep="\t")


# 4. for loop to analyze each peptide
for ( pep in peps_for_analysis ){
  cat(pep)
  cat("\n")
  # get initial hypotheses
  init <- unlist(lapply(rep(length(molecs), 1), generate_random_hypothesis))
  # get the dataframe with shape parameters if you are using a prior
  if ( use_prior == T){
    p.df <- get_shape_df(pred_data, pep, mode_F=0.001, sd_F=0.15, mode_T=0.35, sd_T=0.2)
  }
  else { p.df <- NULL }
  # We need to find out which radius to use
  # wrapper function so we can find the right radius to get an acceptance rate of about 50%
  get_acc_rate_wrap <- function(rad){
    get_acc_rate_chain(mh_chain(eli.dat,molecs,init,7500,
                                pep,unif.prop=F,radius=rad, p.df=p.df))-52
  }
  RAD <- uniroot(get_acc_rate_wrap, lower=0.0001, upper=0.5)
  cat("\t   got radius")
  cat("\n")
  # Now move on to the actual chain
  mh_out <- mh_chain(eli.dat, molecs, init_h=init, max_steps, pep, 
                     unif.prop=F, radius=RAD$root, p.df=p.df)
  mh_out$LogLik <- as.double(as.character(mh_out$LogLik))
  cat("\t   got chain")
  cat("\n")
  
  # saving chain
  write.csv(mh_out, paste(outDir, pep, "_full_chain.txt", sep=""), row.names=F, quote=F)
  
  # LogLik evol
  if ( print_loglik_evol==T ){
    mh_out$time <- 1:max_steps
    pic <- ggplot(mh_out) + geom_line(aes(x=time, y=LogLik)) + ggtitle(paste(pep, "chain", RAD))
    png(paste(outDir, pep, "_LogLik_chain.png", sep=""), width=3000)
    print(pic)
    dev.off()
    rm(pic)
  }
  # read in the chain in nicer form
  mh_out <- read.csv(paste(outDir, pep, "_full_chain.txt", sep=""), skip=1)
  colnames(mh_out) <- c(molecs, "LL")
  plot_posteriors(mh_out, nCol=3, fileName=paste(outDir, pep, "_posteriors.png"))
}

# 5. Output dataframe with mode, mean, median and DKL-vs-uniform for each pHLA combo
res <- get_overview_df(peps_for_analysis, chainDir=outDir, molecs)
write.table(res, paste(outDir, "results_table.txt"), col.names=T, row.names=F, sep="\t")
