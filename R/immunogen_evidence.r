############
# EVIDENCE #
############

#' Evidence function
#'
#' Computes P(D|H)
#' @param eli.dat the dataframe containing the outcome of the ELISpot experiments
#' @param H hypothesis, a string of values between 0 and 1 separated by commas. There should be length(mol.names) values
#' @param mol.names Names of the molecules in your ELISpot subjects
#' @param pep Name of the peptide currently being processed

get_evidence_1 <- function(eli.dat, H, mol.names, pep){
  # H = hypothesis = vector of values between 0 and 1 as a string
  # output: log(evidence)
  H <- as.double(unlist(strsplit(H, ","))) # need to get individual values from string
  L_neg <- rep(0, times=dim(eli.dat)[1]) # vector: each entry is a person
  for ( i in 1:length(mol.names) ){
    L_neg <- L_neg + log(1-H[i])*eli.dat[,colnames(eli.dat)==mol.names[i]]
    ## this eli.dat value could also be 2 in case of homozygosity
  }
  ev <- sum(L_neg[eli.dat[,colnames(eli.dat)==pep]==F]) +
    sum (log(1-exp(L_neg[eli.dat[,colnames(eli.dat)==pep]==T])))
  return(ev)
}



