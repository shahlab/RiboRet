#' Remove contaminating reads in a parallel fashion.
#'
#' This function removes contaminating reads, generally rRNA/tRNA reads, but
#' can be configured to remove reads mapping to wherever the user chooses.
#' It uses a user provided GFF to remove reads that map to user designated
#' regions.
#'
#' @param bamList A named list of bam files
#' @param gffLoc The path to a gff
#' @param contaminants A character vector containing the thing in gffCol you
#' want to remove, for example, rRNA. This is searched for in the gff "type"
#' column to make a gff that defines the contaminating reads. This is currently
#' not flexible because it requires a "type" column.
#' @return A list of bams depleted of contaminating reads
#' @export
#'
#'
remove_contam_reads <- function(bamList, gffLoc, contaminants) {
  # read in the gff
  gff <- rtracklayer::readGFFAsGRanges(gffLoc)

  # make a rtRNA gff
  rtrna.gff <- gff[gff$type %in% contaminants]

  # determine possible threads, maximally 4
  if (length(bamList) > 4) {
    threads <- 4
  } else {
    threads <- length(bamList)
  }

  # for each bam, remove the rtRNA reads
  parallel::mclapply(bamList, function(bam) {
    # define the set of reads to remove from the BAM
    reads.to.remove <-
      S4Vectors::queryHits(IRanges::findOverlaps(bam, rtrna.gff))

    #remove them
    return(bam[-reads.to.remove])
  }, mc.cores = threads)
}
