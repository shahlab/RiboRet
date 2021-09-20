#' A function that counts the reads at specific positions in each sample
#'
#' This function is intended to use the GRanges output from
#' find_related_positions to count the number of reads that occur in groups of
#' positions that looked interesting across the samples but can use any GRanges
#' object you're interested in, for example, a start codon only GFF. This
#' function defaults to 4 threads.
#'
#' @param bamList A list of bams produced by parallel_read_bams, or one cleaned
#' of contaminating reads by remove_contam_reads.
#' @param gRangesToLookFor a GRanges object. Reads corresponding to the ranges
#' in this object will be found from each of the bam files
#' @returns Returns a data frame that shows the number of reads from each file
#' mapping within the ranges specified in gRangesToLookFor
#' @export
#'
count_reads_at_groups <- function(bamList, gRangesToLookFor) {
  # for each bam
  parallel::mclapply(bamList, function(bam) {
    gRangesToLookFor[S4Vectors::subjectHits(IRanges::findOverlaps(bam, gRangesToLookFor))] |> # find the reads that fall in the specified range
      tibble::as_tibble() |> # convert to df
      dplyr::group_by(start, end) |> # for each start and end, i.e. each range
      dplyr::tally() |> # count the number of reads there
      dplyr::ungroup()
  }, mc.cores = 4) |>
    dplyr::bind_rows(.id = "sample") |> # bind to single df
    dplyr::mutate(range_size = end - start) # find the size of the range
}
