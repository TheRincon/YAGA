#!/usr/bin/env Rscript

library("evobir")

args = commandArgs(trailingOnly=TRUE)

file_list <- list.files(args[3])

og_list <- []

for (file in file_list) {
	og <- basename(file).
	og_list[] <- CalcPopD(alignment = system.file(file)
}