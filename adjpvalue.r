#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# setwd("/home/danielw1234/Desktop/Results/")
dat <- read.table(args[1], sep = "\t", header = TRUE)
dat <- dat[-1,]
j <- list()
h <- list()
display_line <- list()
fileConn<-file(args[2])
for (i in seq(1,nrow(dat))) {
  x <- matrix(c(dat[i,2], dat[i,3], dat[i,4], dat[i,5]), byrow = TRUE, 2, 2)
  display_line[i] <- paste(dat[i,2], dat[i,3], dat[i,4], dat[i,5])
  y <- fisher.test(x, alternative = "greater")
  h[i] <- y$p.value
}
for (p in seq(1, nrow(dat))) {
  b = dat[p,1]
  dat[p,1] = toString(b)
}
adjust <- p.adjust(h, method = "fdr", n = length(h))
file_text <- paste(dat$GOTERM, bon, display_line)
write(file_text, fileConn, append=TRUE)
close(fileConn)