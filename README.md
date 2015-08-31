# BIITE

This repository contains the code for the BIITE algorithm, which infers peptide:HLA immunogenicity from ELISpot data. 

To install, you can use devtools

```
#if you don't have devtools, install it
install.packages("devtools")
library(devtools)
devtools::install_github("liesb/BIITE")
```

Example data is available in Data; an example script is available in Example/example.r. This example script uses the example data. Before you run the example script, make sure to use 
```
setwd("some/directory")
```
to make sure the output goes where you want it to go.

To run it with your own data, you can adapt the example script; I've tried to indicate what should be (un)commented and changed by the user. 

We have submitted the BIITE paper for publication in a peer-reviewed journal. Please email me (lies.boelen@gmail.com) for a reference (I will add one as soon as the paper has been accepted).

I'm happy to provide support, just drop me a line!
