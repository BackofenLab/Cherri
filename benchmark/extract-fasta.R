
library(tidyverse)
library(readxl)

# read data
dat <- 
  read_xlsx("benchmark-data.xlsx", sheet = "benchmark-data", na="NA") |> 
  drop_na(id)

# write fasta files
for (d in c("query","target")) {
    str_c(">",dat$id,".",str_sub(d,1,1),"\n",dat[[d]], collapse="\n") |> 
    cat(sep="\n", file = str_c(d,".fa"))
}

