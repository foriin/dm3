library(GenomicRanges)
library(dplyr)

# Variables and data

# dm3.genes <- fread("~/IMG/data/dmel/gene_list/Drosophila_melanogaster.BDGP5.78.full.genes.gtf") %>%
#   select(c(1, 4, 5, 7, 9)) %>%
#   setNames(c("chr", "start", "end", "strand", "attr")) %>%
#   mutate(FBgn = sub('.*gene_id "(FBgn[0-9]+)";.*', '\\1', attr), gene_id = sub('.*gene_name "([^;]+)";.*', '\\1', attr),
#          TSS = ifelse(strand == "+", start, end)) %>%
#   select(-attr)
#
# dm3.genes <- dm3.genes %>% mutate(chr = sub("chr", "", chr))
# dm3.genes.gr <- makeGRangesFromDataFrame(dm3.genes, keep.extra.columns = T)
# dm3.tss.gr <- makeGRangesFromDataFrame(dm3.genes %>% select(-start, -end),
#                                        start.field = "tss",
#                                        end.field = "tss",
#                                        keep.extra.columns = T)
# save(dm3.genes, dm3.genes.gr, dm3.tss.gr, file = "genes.RData")


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

chr.lengths <- c(
  23011544,
  368872,
  21146708,
  3288761,
  24543557,
  2555491,
  27905053,
  2517507,
  1351857,
  22422827,
  204112,
  347038
)
names(chr.lengths) <- all.chroms

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

add.chr <- function(obj, rev = F){
  # Function to add or remove "chr" from chromosome names either in
  # data frames with designated "chr" column (chr being a name of a column)
  # or in GRanges objects
  if (typeof(obj) == "list"){
    if (rev) {
      obj$chr <-  gsub("chr", "", obj$chr)
      return(obj)
    } else{
      obj$chr <-  paste0("chr", obj$chr)
      return(obj)
    }
  } else if (typeof(obj) == "S4"){
    if (rev){
      seqlevels(obj) <-  gsub("chr", "", seqlevels(obj))
      return(obj)
    } else {
      seqlevels(obj) <- paste0("chr", seqlevels(obj))
      return(obj)
    }
  } else {
    stop("Wrong type of object. Must be either data frame or GRanges object")
  }
}


