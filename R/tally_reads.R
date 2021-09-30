#' Count reads at positions given an offset
#'
#' This function counts the number of reads with 3' ends at a given position and
#' separates them by read lengths. This function defaults to using 4 threads. IT
#' CURRENTLY ONLY CONSIDERS POSITIVE STRAND READS.
#'
#' @param bamList A named list of bam files, presumably from the remove_contam_reads function
#' @param offset The offset from the 3' end, this is applied to all reads of all lengths
#' @returns A single data frame with counts per position per read length for all samples
#' @export
#'
tally_reads <- function(bamList, offset, ncores = 4){
  parallel::mclapply(bamList, function(bam){      # for each bam
    bam |>                                        #
      tibble::as_tibble() |>                      # convert to df
      dplyr::mutate(position = ifelse(strand == "+", end - offset, start + offset)) |>
      dplyr::group_by(strand, qwidth, position) |>            # for each end and read length
      dplyr::tally() |>                          # count reads that occur there
      dplyr::ungroup() |>                        #
      dplyr::mutate(rpm = n*1e6/sum(n))
  }, mc.cores = ncores) |>                        # number of cores
    dplyr::bind_rows(.id = "sample")              # combine to one df
}
