#!/usr/bin/env Rscript
# Script to convert XGMix msp putput file to bed like format
# this output bed file is intended to be used for running
# Tracts: 
'Usage:
   global-ancestry-xgmix.R [-m <msp> -i <individual> -o <basenameoutput>]
   global-ancestry-xgmix.R (-h | --help)
Options:
   -m input msp file
   -i Indivisual the individual (sample) to extract
   -o output file basename,  -<individual>.tsv will be appened
 ]' -> doc

suppressPackageStartupMessages(library(docopt))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))

opts <- docopt(doc)

# functions ---------------------------------------------------------------


load_population_code <- function(msp_file) {
  
  pop_code_raw <- read_lines(msp_file, n_max = 1)
  
  pop_code_raw <- 
    pop_code_raw %>% 
    str_replace("#Subpopulation order/codes: ", "") %>% 
    str_split(pattern = "\\t") %>% 
    .[[1]]
  
  tibble(
    assignment = str_extract(pop_code_raw, pattern = "^\\w*"),
    ancestry = str_extract(pop_code_raw, pattern = "\\d$")
  ) %>% 
    mutate(ancestry = as.numeric(ancestry))
  
}


load_msp_as_bed <- function(msp_file, individual) {
  
  msp <- 
    read_tsv(msp_file, skip = 1) %>% 
    rename(chm = `#chm`, n_snps = `n snps`)
  
  # select individual
  
  msp <- msp %>% 
    select(chm, spos, epos, sgpos, egpos, n_snps, contains(individual))
  
  # the number of columns should be 8
  # if there are less than 8 coulumns it means the individual was not found
  
  if (ncol(msp) != 8) stop("individual not present in msp input")
  

  # pivot to long format
  msp <- 
    msp %>% 
    pivot_longer(
      cols = -c(chm, spos, epos, sgpos, egpos, n_snps),
      names_to = "sample_haplo",
      values_to = "ancestry"
    ) %>% 
    select(-n_snps) %>% 
    separate(sample_haplo, into = c("Individual", "Haplotype"), sep = "\\.")
  
  msp <- 
    msp %>% 
    inner_join(load_population_code(msp_file)) %>% 
    select(-ancestry) %>% 
    select(chm, spos, epos, assignment, sgpos, egpos, Individual, Haplotype)
  
  msp %>% 
    arrange(Individual, spos)
  
} 


# run ---------------------------------------------------------------------


msp_file <- opts$m
individual <- opts$i
out_file <- opts$o
out_file <- paste0(out_file, "-", individual, ".tsv")


load_msp_as_bed(msp_file, individual) %>% 
  write_tsv(out_file)