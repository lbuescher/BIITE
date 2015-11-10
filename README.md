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

To run it with your own data, you can adapt the example script; I've tried to indicate what should be (un)commented and changed by the user. Input data should at least contain:
- for each HLA: a column with the copy number per patient
- for each peptide: a Boolean column indicating whether the ELISpot for this peptide came up positive (TRUE) or negative, for each patient. 

A small example of input data:

|Patient|	DQB1_02 |DQB1_03	|DQB1_04	|DQB1_05	|DQB1_06	|pep_1	|pep_2	|pep_3|
|-------|---------|---------|---------|---------|---------|-------|-------|-----|
|Pat_1	|1	|0	|0	|0	|1	|TRUE	|TRUE	|TRUE|
|Pat_2	|1	|1	|0	|0	|0	|TRUE	|TRUE	|TRUE|
|Pat_3	|0	|2	|0	|0	|0	|TRUE	|FALSE	|FALSE|
|Pat_4	|0	|0	|0	|2	|0	|TRUE	|TRUE	|TRUE|
|Pat_5	|0	|0	|2	|0	|0	|TRUE	|TRUE	|FALSE|
|Pat_6	|1	|0	|0	|0	|1	|FALSE	|FALSE	|TRUE|



We have submitted the BIITE paper for publication in a peer-reviewed journal. Please email me (lies.boelen@gmail.com) for a reference (I will add one as soon as the paper has been accepted).

I'm happy to provide support, just drop me a line!
