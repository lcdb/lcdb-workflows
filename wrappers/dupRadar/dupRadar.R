#!/usr/bin/env Rscript


library(dupRadar)
args   <- commandArgs(TRUE)
bam <- args[1]
gtf <- gsub("gtf=","",args[2])
stranded <- gsub("stranded=","",args[3])
paired   <- gsub("paired=","",args[4])
outfile   <- gsub("outfile=","",args[5])
threads  <- as.integer(gsub("threads=","",args[6]))

if(length(args) != 6) {
  stop (paste0("Usage: ./dupRadar.sh <file.bam> <genes.gtf> ",
               "<stranded=[no|yes|reverse]> paired=[yes|no] ",
               "outfile=./ threads=1"))
}

if(!file.exists(bam)) {
  stop(paste("File",bam,"does NOT exist"))
}

if(!file.exists(gtf)) {
  stop(paste("File",gtf,"does NOT exist"))
}


if(is.na(stranded) | !(grepl("no|yes|reverse",stranded))) {
  stop("Stranded has to be no|yes|reverse")
}

if(is.na(paired) | !(grepl("no|yes",paired))) {
  stop("Paired has to be no|yes")
}

if(is.na(threads)) {
  stop("Threads has to be an integer number")
}

stranded <- if(stranded == "no") 0 else if(stranded == "yes") 1 else 2


## calculate duplication rate matrix
dm <- analyzeDuprates(bam,
                      gtf,
                      stranded,
                      (paired == "yes"),
                      threads)

## duprate vs. expression smooth scatter
png(file=outfile, width=1000, height=1000)
duprateExpDensPlot(dm, main=basename(bam))
dev.off()
