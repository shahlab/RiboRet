#' Read in BAM files in a parallel manner.
#'
#' This function simply reads BAM files in using parallel's mclapply and
#' GenomicAlignments readGAlignments. It results in a named list of bam files. 
#' By default, this reads all BAMs it can find, if you want to use a subset of
#' BAMs from a folder, you make your own named list of BAM files.
#'
#' @param dirWithBams path to a directory with bam files in it, it will only.
#' attempt to read things that end in .bam.
#' @param ncores the number of cores to use, typically 1 per bam.
#' @return A list of bams read in using readGAlignments.
#' @export
#' 
#' 
parallel_read_bams <- function(dirWithBams, ncores) {
  # find the bams
  bam.locs <-
    dir(
      dirWithBams,
      pattern = ".bam",
      full.names = TRUE,
      recursive = TRUE
    )
  
  # name the bams with themselves, but only the important part
  bam.names <- sapply(bam.locs, function(bam) {
    stringr::str_split(bam, "/") |>
      unlist() |>
      tail(1) |>
      stringr::str_remove(".bam")
  })
  
  # set names
  names(bam.locs) <- bam.names
  
  # read each bam in
  return(
    parallel::mclapply(bam.locs, GenomicAlignments::readGAlignments, mc.cores = ncores)
  )
}