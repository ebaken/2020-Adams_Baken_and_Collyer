# Script to compare two Rmd files
 #Note, requires installation of 'latexdiff' in your tex components & a version of Perl

# NOTE: need to install latexdiffr from: https://github.com/hughjonesd/latexdiffr

library(latexdiffr)
latexdiff("Manuscript-5.Rmd", "Manuscript-4-EKB.Rmd", clean = TRUE)

