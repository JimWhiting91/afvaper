#' Transform genome windows/co-ordinates to df of chr, start, pos
#'
#' Takes as input genome windows/co-ordinates in the format of chrX:start-end and splits to 3 column matrix. The function is just a wrapper for tidyr::separate that assumes the chr:start-end format 
#'
#' @param genome_windows character string or character vector of genome windows/co-ordinates in the format of chr:start-end
#'
#' @import tidyr
#'
#' @return A data.frame of chr, start, end identifiers for genome windows
#' @export
window2pos_df <- function(genome_windows){
  tidyr::separate(data.frame(genome_windows),col = "genome_windows",into=c("chr","start","end"),sep = ":|-")
}
  