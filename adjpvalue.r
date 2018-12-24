#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# setwd("/home/danielw1234/Desktop/Results/")
dat <- read.table(args[1], sep = "\t", header = TRUE)
dat <- dat[-1,]
j <- list()
h <- list()
display_line <- list()
# file.create(args[2])
# f <- file(args[2])
# open(f, "w+")
# f<- read.csv(args[2], sep="\t", header=TRUE)
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
GO_HEADER <- paste("GO_TERM	ADJPVALUE	ENRICHED_SET_AND_GO	ALL_GOTERMS_MINUS_ENRICHED	ENRICHED_NOT_IN_GO	NOT_IN_GO_NOT_IN_ENRICHED")
go_text <- paste(dat$GOTERM, adjust, display_line)
write(GO_HEADER, file="")
write(go_text, file="")