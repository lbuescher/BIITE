\name{BIITE-package}
\alias{BIITE-package}
\alias{BIITE}
\docType{package}
\title{
  Bayesian Immunogenicity Inference Tool for ELISpot
}
\description{
  A Bayesian method to infer peptide:HLA-II immunogenicity from ELISpot data
}
\details{
  \tabular{ll}{
    Package: \tab BIITE\cr
    Type: \tab Package\cr
    Version: \tab 0.0.0.9000\cr
    Date: \tab 2015-09-01\cr
    License: \tab GPL (>= 2)\cr
    Depends: \tab ggplot2\cr \tab reshape2\cr \tab modeest\cr \tab entropy\cr
  }
}

\author{
  Lies Boelen <l.boelen@imperial.ac.uk> and Becca Asquith <b.asquith@imperial.ac.uk>
    
    Maintainer: Lies Boelen <l.boelen@imperial.ac.uk>
}
\references{
  To be confirmed.
}

\keyword{ package }
\examples{
  library(BIITE)
  library(ggplot2)
  data("simul_data")
  eli.dat <- simul_data
  pep_names <- paste0("pep_", 1:17)
  max_steps <- 1000
  outDir <- "./"
  peps_for_analysis <- pep_names[1:2]
  use_prior <- FALSE
  print_loglik_evol <- TRUE
  molecs <- colnames(eli.dat)[(grepl("DRB1_", colnames(eli.dat)) | grepl("DQB1_", colnames(eli.dat))) & !grepl("al", colnames(eli.dat))]
  write.table(molecs, paste0(outDir, "molecs.txt"), row.names=F, col.names=F, sep="\t")
  # for loop to analyze each peptide
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
    # Normally, we need to find out which radius to use; fixed now. 
    # See example.r for how this is normally handled
    RAD <- list(root=0.01)
    mh_out <- mh_chain(eli.dat, molecs, init_h=init, max_steps, pep,
                       unif.prop=F, radius=RAD$root, p.df=p.df)
    mh_out$LogLik <- as.double(as.character(mh_out$LogLik))
    # saving chain
    write.csv(mh_out, paste0(outDir, pep, "_full_chain.txt", sep=""), row.names=F, quote=F)
    # LogLik evol
    if ( print_loglik_evol==T ){
      mh_out$time <- 1:max_steps
      pic <- ggplot(mh_out) + geom_line(aes(x=time, y=LogLik)) + ggtitle(paste(pep, "chain", RAD))
      png(paste0(outDir, pep, "_LogLik_chain.png", sep=""), width=3000)
      print(pic)
      dev.off()
      rm(pic)
    }
    # read in the chain in nicer form
    mh_out <- read.csv(paste0(outDir, pep, "_full_chain.txt", sep=""), skip=1)
    colnames(mh_out) <- c(molecs, "LL")
    plot_posteriors(mh_out, nCol=3, fileName=paste0(outDir, pep, "_posteriors.png"))
  }
  # Output dataframe with mode, mean, median and DKL-vs-uniform for each pHLA combo
  res <- get_overview_df(peps_for_analysis, chainDir=outDir, molecs)
  write.table(res, paste0(outDir, "results_table.txt"), col.names=T, row.names=F, sep="\t")
}
