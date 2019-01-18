#!/usr/bin/env Rscript

library("evobiR")
args = commandArgs(trailingOnly=TRUE)
file_list <- list.files(args[1])
og_list <- vector(mode="list", length=length(file_list))

i <- 0
for (file in file_list) {
	print("|||")
	i <- i + 1
	og <- strsplit(basename(file), "_")[[1]][2]
	og_num <- strsplit(basename(file), "_")[[1]][3]
	print(og)
	print(substr(og_num, 1, nchar(og_num)-6))
	og_list[[i]] <- CalcPopD(alignment=paste(args[1],file, sep="/"))
	print("===")
}