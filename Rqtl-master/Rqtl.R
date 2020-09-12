//Set the working directory
setwd("~/Documents/Hufford/Sarah")

//Get package needed for analysis
install.packages("qtl")
library(qtl)

//Read in the data
DMz18_205.BLAST <- read.cross("csv", "/Users/User/Documents/Hufford/Sarah","test.csv")