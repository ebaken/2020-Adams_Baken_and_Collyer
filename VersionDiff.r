# Script to compare two Rmd files
 #Note, requires installation of 'latexdiff' in your tex components & a version of Perl

# NOTE: need to install latexdiffr from: https://github.com/hughjonesd/latexdiffr

library(latexdiffr)
latexdiff("2020-Evol-MultiRR_Vsn6.Rmd", "2020-Evol-MultiRR_Vsn6_AK.Rmd", clean = TRUE)

