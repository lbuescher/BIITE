###################
# POSTERIOR PLOTS #
###################

library(ggplot2)
library(reshape2)

#' Get a plot of the marginal posteriors for each pHLA combination (all HLAs, one peptide p)
#'
#' Implementation of the Metropolis-Hastings algorithm to find posterior distributions of pHLA immunogenicity for a single peptide
#' @param chain The MH chain in question
#' @param nCol The number of columns in the output figure
#' @param fileName Filename for the figure

plot_posteriors <- function(chain, nCol=3, fileName){
  M <- reshape::melt(chain[,!colnames(chain)=="LL"])

  k <- ggplot(M)
  k <- k + geom_density(aes(x=value), size=1, show_guide=F)
  k <- k + stat_density(aes(x=value,), size=0, geom="line", position="identity")
  k <- k + facet_wrap(~variable, ncol=nCol)
  k <- k + theme_bw()
  k <- k + theme(panel.grid.major = element_line(colour = "black", size = 0.2),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black")
  )
  k <- k + theme(legend.justification = c(1,0), legend.position = c(1,0), legend.title = element_text(size=15),
                 legend.text = element_text(size=15))
  k <- k + scale_x_continuous(name="Immunogenicity", expand=c(0,0))
  k <- k + scale_y_continuous(name="Probability Density", expand=c(0,0))
  k <- k + theme(axis.text = element_text(size=17, colour="black"), axis.title.x = element_text(size=20), 
                 axis.title.y = element_text(vjust=0.8, size=20))
  k <- k + theme(strip.text = element_text(size=20),
                 legend.key = element_blank())
  k
  
  png(fileName,500,1000)
  print(k)
  dev.off()
}