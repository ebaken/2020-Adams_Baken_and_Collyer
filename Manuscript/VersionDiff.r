# Script to compare two Rmd files
 #Note, requires installation of 'latexdiff' in your tex components & a version of Perl

# NOTE: need to install latexdiffr from: https://github.com/hughjonesd/latexdiffr

library(latexdiffr)
latexdiff("Manuscript-6-DCA.Rmd", "Manuscript-6-EKB.Rmd", clean = TRUE, output = "V6-V6EB.diff.pdf")

