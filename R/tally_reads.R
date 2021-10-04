#' Count reads at positions given an offset
#'
#' This function counts the number of reads with 3' ends at a given position and
#' separates them by read lengths.
#'
#' @param bamList A named list of bam files, presumably from the remove_contam_reads function
#' @param offset A data frame of columns qwidth, offset, specifying the read length specific offsets
#' @returns A single data frame with counts per position per read length for all samples
#' @export
#'
tally_reads <- function(bamList, offset, ncores){
  parallel::mclapply(bamList, function(bam){        # for each bam
    bam |>                                          #
      tibble::as_tibble() |>                        # convert to df
      dplyr::left_join(., offset, by = "qwidth") |> # join the offset data frame
      dplyr::mutate(position = ifelse(strand == "+", end - offset, start + offset)) |>
      dplyr::group_by(strand, qwidth, position) |>  # for each end and read length
      dplyr::tally() |>                             # count reads that occur there
      dplyr::ungroup() |>                           #
      dplyr::mutate(rpm = n*1e6/sum(n))             # normalize reads
  }, mc.cores = ncores) |>                          # number of cores
    dplyr::bind_rows(.id = "sample")                # combine to one df
}
