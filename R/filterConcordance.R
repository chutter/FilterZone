#' @title filterConcordance
#'
#' @description Function for obtaining filtered concordance factor data
#'
#' @param input.dir summary data file from filterSummary
#'
#' @param clade.list a named list of clades of interest to test for concordance factors
#'
#' @param outgroups outgroups to root the tree
#'
#' @param all.data select T to obtain all the node data or F to obtain only nodes from clade list
#'
#' @return data.frame with concordance factor data for each filtered replicate
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


# Filter concordance function
filterConcordance = function(input.dir = NULL,
                             clade.list = NULL,
                             outgroups = NULL,
                             all.data = TRUE) {

  #Debug
  # input.dir = "concordance-factors"
  # clade.list = taxa.set
  # outgroups  = outgroup.taxa
  # all.data = TRUE

  cf.files = list.files(input.dir)
  cf.files = unique(gsub("\\..*", "", cf.files))

  all.data = c()
  for (x in 1:length(cf.files)){

    #Reads in the plane data
    plane.data = AstralPlane::createAstralPlaneCF(cf.file.name = paste0("concordance-factors/", cf.files[x]),
                                                  outgroups = outgroups,
                                                  tip.length = 1)

    #Merges the CF and node data
    merge.data = merge(plane.data@concordanceFactorData,
                       plane.data@nodeData,
                       by = "node")

    #Sets up new columns
    merge.data = cbind(monophyletic = NA, merge.data)
    merge.data = cbind(clade = NA, merge.data)
    merge.data = cbind(dataset = gsub(".*/", "", cf.files[x]), merge.data)

    for (y in 1:length(taxa.set)){

      #Overall target tree
      mrca.node = ape::getMRCA(plane.data@phylo, taxa.set[[y]])

      if (ape::is.monophyletic(plane.data@phylo,  taxa.set[[y]]) == T){
        mono.clade = TRUE
      } else{
        mono.clade = FALSE
      }#end else

      #What is the goal here?
      merge.data[merge.data$node %in% mrca.node,]$clade = names(taxa.set)[y]
      merge.data[merge.data$node %in% mrca.node,]$monophyletic = mono.clade

      #Obtain other close nodes
      temp.a = plane.data@edgeData
      node.data1 = temp.a[temp.a$node1 %in% mrca.node,]
      node.data2 = temp.a[temp.a$node2 %in% mrca.node,]
      node.data = rbind(node.data1, node.data2)
      more.nodes = unique(c(node.data$node1, node.data$node2))
      merge.data[merge.data$node %in% more.nodes,]$clade = names(taxa.set)[y]
    }#y loop

    all.data = rbind(all.data, merge.data)

  }#end x loop

  #Returns all the data
  all.data$node = as.numeric(all.data$node)
  return(all.data)

}#end function



