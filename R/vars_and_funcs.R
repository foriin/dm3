library(GenomicRanges)
library(dplyr)

# Variables and data

load("genes.RData")

all.chroms <-  c("2L", "2LHet", "2R", "2RHet", "3L", "3LHet",
                 "3R", "3RHet", "4", "X", "XHet", "YHet")
euc.chroms <- c("2L", "2R", "3L", "3R", "X")
euc.gr <- GRanges(
  seqnames = Rle(c("2L", "2R", "3L", "3R", "X")),
  ranges = IRanges(
    start = c(1, 1600000, 1, 1, 1),
    end = c(22000000, 21146700, 22900000, 27900000, 22422800)
  )
)

# Functions

bedTools.shuffle.jac <- function(bed.1, bed.2, shuf = F, opt.string="-chrom"){

  bed.file.1 <- tempfile()
  bed.file.2 <- tempfile()

  shuf.1 <- tempfile()
  shuf.2 <- tempfile()

  jac <- tempfile()

  options(scipen = 99)

  write.table(bed.1, file = bed.file.1, quote = F, sep = "\t", col.names = F, row.names = F)
  write.table(bed.2, file = bed.file.2, quote = F, sep = "\t", col.names = F, row.names = F)
  if (shuf){
    command = paste("bedtools shuffle -i", bed.file.1,
                    "-g /home/artem/IMG/data/dmel/Genome/dm3.genome", opt.string, "|",
                    "sort -k1,1 -k2,2n - >", shuf.1, ";",
                    "bedtools shuffle -i", bed.file.2,
                    "-g /home/artem/IMG/data/dmel/Genome/dm3.genome", opt.string, "|",
                    "sort -k1,1 -k2,2n - >", shuf.2, ";",
                    "bedtools jaccard -a", shuf.1, "-b", shuf.2, ">", jac)
    # cat(command, "\n")
  } else {
    command = paste("bedtools jaccard -a", bed.file.1, "-b", bed.file.2, ">", jac)
  }
  try(system(command))

  res=read.table(jac, header = T)
  unlink(bed.file.1); unlink(bed.file.2); unlink(shuf.2); unlink(shuf.1); unlink(jac)
  return(res$jaccard)
}

df.from.GRanges <- function(gr){
  data.frame(
    chr = as.vector(seqnames(gr)),
    start = start(gr),
    end = end(gr)
  )
}

delete.het <- function(data, r6 = FALSE){
  if (!r6){
    euc.coords <- GRanges(
      seqnames = Rle(c("2L", "2R", "3L", "3R", "X")),
      ranges = IRanges(
        start = c(1, 1600000, 1, 1, 1),
        end = c(22000000, 21146700, 22900000, 27900000, 22422800)
      )
    )
  } else {
    euc.coords <- GRanges(
      seqnames = Rle(c("2L", "2R", "3L", "3R", "X")),
      ranges = IRanges(
        start = c(1, 5712495, 1, 4174279, 103614),
        end = c(22000000, 25259177, 22906900, 32074278, 23020964)
      )
    )
  }

  data.coords <- GRanges(
    seqnames = Rle(data$chr),
    ranges = IRanges(
      start = data$start,
      end = data$end
    )
  )

  data.x.euc <- subsetByOverlaps(data.coords, euc.coords)
  return(merge(data,
               data.frame(chr = seqnames(data.x.euc),
                          start = start(data.x.euc),
                          end = end(data.x.euc)),
               by = c("chr", "start", "end")) %>% arrange(chr, start))
}


delete.het.gr <- function(data, r6 = FALSE){
  if (!r6){
    euc.coords <- GRanges(
      seqnames = Rle(c("2L", "2R", "3L", "3R", "X")),
      ranges = IRanges(
        start = c(1, 1600000, 1, 1, 1),
        end = c(22000000, 21146700, 22900000, 27900000, 22422800)
      )
    )
  } else {
    euc.coords <- GRanges(
      seqnames = Rle(c("2L", "2R", "3L", "3R", "X")),
      ranges = IRanges(
        start = c(1, 5712495, 1, 4174279, 103614),
        end = c(22000000, 25259177, 22906900, 32074278, 23020964)
      )
    )
  }

  subsetByOverlaps(data, euc.coords)

}


