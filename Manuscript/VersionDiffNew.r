# Script to compare two Rmd files

# step by step, as 'latexdiffr' no longer working (breaks in the system2(latexdiff call)

#1: generate tex from Rmd
rmarkdown::render(input =  "2020-Evol-MultiRR_Submitted.Rmd")
rmarkdown::render(input = "2020-Evol-MultiRR_Rvsn.Rmd")

#2: latexdiff of tex files
shell(cmd="latexdiff.exe 2020-Evol-MultiRR_Submitted.tex 2020-Evol-MultiRR_Rvsn.tex > diff.tex")

#3: PDF of diff.tex
tinytex::latexmk("diff.tex") 












###########
#Note, requires installation of 'latexdiff' in your tex components & a version of Perl

# NOTE: need to install latexdiffr from: https://github.com/hughjonesd/latexdiffr

##OLD WAY: not working now
#library(latexdiffr)
#latexdiff("2020-Evol-MultiRR_Submitted.Rmd", "2020-Evol-MultiRR_Rvsn.Rmd", clean = TRUE)

