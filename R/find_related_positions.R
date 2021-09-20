#' Finds genomic positions that have many reads mapping to them.
#'
#' This function helps the user find positions that have many reads mapping to
#' them across all their samples under the logic that positions with many reads
#' mapping to them in at least one sample (though not necessarily others) could
#' be important.
#'
#' @param bamDf A data frame produced by the tally_reads function
#' @param maxDiff an integer, the maximum distance in number of bases to
#' consider something as "close". If a read occurs at position 1 and another at
#' position 4, these reads will be in the same "group" if maxDiff >= 3 but in
#' separate groups if maxDist <= 2
#' @param gffLoc the location of a gff
#' @returns a GRanges object that shows positions that have many reads at them
#' across all samples.
#' @export
#'
#'
find_related_positions <- function(bamDf, maxDiff, gffLoc){
  # read in the gff
  gff <- rtracklayer::readGFFAsGRanges(gffLoc)

  # find all the unique positions across the bams
  unique.pos.df <- bamDf |>
    dplyr::select(position) |>
    dplyr::arrange(position) |>
    dplyr::distinct()

  # a new col that will find the diff between a row and prev row
  downshifted.row <- c(NA, unique.pos.df$position[1:(nrow(unique.pos.df) - 1)])

  # find the differences
  bam.with.downshited.col <- unique.pos.df |>
    dplyr::mutate(pos2 = downshifted.row,
                  difference = position - pos2)

  # define some variables for the loop
  unified.group.list <- list()
  unified.group.list[[1]] <- bam.with.downshited.col$position[1]
  cc <- 1

  # for each row in the unified bam
  for(i in 2:nrow(bam.with.downshited.col)){
    # if the difference between a number and the next number is less than maxDiff
    if (bam.with.downshited.col$difference[i] <= maxDiff){
      # add that number to the current group
      unified.group.list[[cc]] <- c(unified.group.list[[cc]], bam.with.downshited.col$position[i])
    } else {
      # otherwise, add it to the next group
      cc <- cc + 1

      unified.group.list[[cc]] <- bam.with.downshited.col$position[i]
    }
  }

  # find the set of ranges from all the bams
  unified.list.range <- sapply(unified.group.list, function(x){c(min(x), max(x))})

  # make a GRanges object from that
  unified.grange <- GenomicRanges::GRanges(seqnames = GenomicAlignments::seqnames(gff)[1],
                                           ranges = IRanges::IRanges(start = unified.list.range[1,],
                                                                     end = unified.list.range[2,]),
                                           strand = "+")

  return(unified.grange)
}
