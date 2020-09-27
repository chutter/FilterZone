#' @title filterSummary
#'
#' @description Summarizes filtered datasets when provided with filter parameters
#'
#' @param alignment.data Alignment summary stats calculcated from summarizeAlignments
#'
#' @param alignment.folder The alignment folder from which the stats were calculated from in alignment.data
#'
#' @param dataset.name The name of your dataset, where all filtered datasets will be placed in this folder
#'
#' @param file.out if you wish to save to file, provide a file name for the summary
#'
#' @param length.filters Your selected length filters as a vector of values for the alignment length
#'
#' @param sample.filters Your selected sampling fraction filters as a vector of values between 0-1
#'
#' @param prop.pis.filters Your selected parsimony informatives sites filter as a vector of values between 0-1
#'
#' @param count.pis.filters Your selected parsimony informatives sites filter as a vector of base pair counts
#'
#' @return Returns a summary of your selected filters applied to the dataset
#'
#' @examples
#'
#' your.tree = ape::read.tree(file = "file-path-to-tree.tre")
#' astral.data = astralPlane(astral.tree = your.tree,
#'                           outgroups = c("species_one", "species_two"),
#'                           tip.length = 1)
#'
#'
#' @export

filterSummary = function(alignment.data = NULL,
                         alignment.folder = NULL,
                         dataset.name = NULL,
                         file.out = NULL,
                         length.filters = c(0),
                         sample.filters = c(0),
                         prop.pis.filters = c(0),
                         count.pis.filters = c(0)) {

  #Read in alignment data and set up
  if (length(alignment.data) == 1) {
    alignment.stats = data.table::fread(alignment.data, header = T)
  } else {
    alignment.stats = alignment.data
  }

  filter.summary = c()
  if (length(filter.length) != 1){
    filter.summary = rbind(filter.summary,
                           filterStats(data = alignment.stats,
                                       filter.name = "alignment_length",
                                       filter.values = filter.length,
                                       align.dataset = dataset.name) )
  }#end length

  #Sampling
  if (length(filter.sample) != 1){
    filter.summary = rbind(filter.summary,
                           filterStats(data = alignment.stats,
                                       filter.name = "proportion_samples",
                                       filter.values = filter.sample,
                                       align.dataset = dataset.name) )
  }#end sample

  #Count
  if (length(filter.count.pis) != 1){
    filter.summary = rbind(filter.summary,
                           filterStats(data = alignment.stats,
                                       filter.name = "count_pis",
                                       filter.values = filter.count.pis,
                                       align.dataset = dataset.name) )
  }#end sample

  #filter.prop.pis
  if (length(filter.prop.pis) != 1){
    filter.summary = rbind(filter.summary,
                           filterStats(data = alignment.stats,
                                       filter.name = "proportion_pis",
                                       filter.values = filter.prop.pis,
                                       align.dataset = dataset.name) )
  }#end sample

  if (is.null(file.out) != T) {
    #Saves the alignment stats
    write.csv(filter.summary, file = paste0(file.out, ".csv"), row.names = F)
  }#end if

  return(filter.summary)

}#end function

