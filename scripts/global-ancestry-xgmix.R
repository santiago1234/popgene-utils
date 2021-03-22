#!/usr/bin/env Rscript
# Script to compute global ancestry for each individual
# based on local ancestry assignments (msp file).
# This script takes as input the msp output file generated
# by XGMix: https://github.com/AI-sandbox/XGMix

'Usage:
   global-ancestry-xgmix.R [-m <msp> -c <chromosome> -o <output>]
   global-ancestry-xgmix.R (-h | --help)

Options:
   -m input msp file
   -c Chromosome name
   -o output file name [default: xgmix-chrN-Q.csv]
 ]' -> doc


suppressPackageStartupMessages(library(docopt))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))

opts <- docopt(doc)

# functions ---------------------------------------------------------------


#' Read population code
#'
#' The first line is a comment line (msp_file), that specifies 
#' the order and encoding of populations,
#' eg: #Sub_population order/code: golden_retriever=0 
#' labrador_retriever=1 poodle poodle_small=2
#' @param msp_file path to <output_basename>.msp.tsv
#'
#' @return a tibble mapping population names to codes
#' @export
#'
#' @examples
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


#' MSP file to bed like format
#' 
#' Converts the msp to long format
#'
#' @param msp_file 
#'
#' @return
#' @export
#'
#' @examples
load_msp_as_bed <- function(msp_file) {
  
  msp <- 
    read_tsv(msp_file, skip = 1) %>% 
    rename(chm = `#chm`, n_snps = `n snps`)
  
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


# main call ---------------------------------------------------------------


msp_file <- opts$m
chrn <- opts$c
out_file <- opts$o

msp <- load_msp_as_bed(msp_file)
q_global <- 
  msp %>% 
  count(Individual, assignment) %>% 
  group_by(Individual) %>% 
  mutate(p = n / sum(n)) %>% 
  ungroup() %>% 
  mutate(
    chromosome = chrn
  )

# add 0% ancestry ---------------------------------------------------------
# the result in the computation above will not iclude cases in which
# and ancestry is 0%
# the next code will add those cases.

q_0s <- 
  crossing(
    tibble(Individual = unique(q_global$Individual)),
    tibble(assignment = unique(q_global$assignment))
  )


q_global <- 
  setdiff(
    q_0s, select(q_global, Individual, assignment)
  ) %>% 
  mutate(
    n = 0,
    p = 0
  ) %>% 
  bind_rows(q_global)

write_csv(q_global, out_file)


