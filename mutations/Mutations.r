library(httr)
library(jsonlite)
library(xml2)
library(stringr)
library(seqinr)

source("Mutations.functions.r")
source("Run.mutation.r")
# Run 

Run.mutation("COSMIC-haematopoietic-point-mutations.csv")