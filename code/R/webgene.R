rm(list = ls())
library(mvtnorm)
library(tidyverse)
setwd("C:/Users/langzx/Desktop")
betas <- read.table("web_dat_HCT16.txt", sep = "\t", header = TRUE)
betas$ID_REF
betas$VALUE
write.csv(x = betas, file = "HCT16.csv", row.names = FALSE)

LOVO <- read.table("LoVo.txt", sep = "\t", header = TRUE)
MB231 <- read.table("MB231.txt", sep = "\t", header = TRUE)
DLD1 <- read.table("DLD1.txt", sep = "\t", header = TRUE)
MB453 <- read.table("MB453.txt", sep = "\t", header = TRUE)
Calu6 <- read.table("Calu6.txt", sep = "\t", header = TRUE)
write.csv(x = LOVO, file = "LOVO.csv", row.names = FALSE)
write.csv(x = MB231, file = "MB231.csv", row.names = FALSE)
write.csv(x = DLD1, file = "DLD1.csv", row.names = FALSE)
write.csv(x = MB453, file = "MB453.csv", row.names = FALSE)
write.csv(x = Calu6, file = "Calu6.csv", row.names = FALSE)
