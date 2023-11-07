#### f14) DR() -----------------------------------------------------------------
# Calculates the DR based on Jetz et al., 2012 (?). This algorithm is able to calculate
# capable to calculate DR for sliced trees, distinguishing "void nodes" from a cutted tree.

# tree = a phylogenetic tree of the class "phylo" from ape package.

# Creating my DR function
#' Title
#'
#' @param tree
#'
#' @return
#' @export
#'
#' @examples

DR <- function(tree){

  # Capturing nodes information
  df <- as.data.frame(nodes_config(tree)$node_matrix)

  drs <- sapply(1:ncol(df), function(x){  #  x <- 8

    # Which nodes give origin to this specie
    nodes <- as.numeric(row.names(df[which(df[, x] == 1),]))

    # Which lines from the edge matrix it occupies?
    lines <- which(tree$edge[, 1] %in% nodes &             # Which nodes are in my spliting side from matrix
                     tree$edge[, 2] %in% c(x, nodes))

    # Save the branch lengths into a vector
    vec <- tree$edge.length[lines]

    # If there is an node of length 0, it is the species node conserved from the slice,
    # which could not be accounted for the DR calculation. Thus, must be removed.
    if(length(which(vec == 0)) > 0){
      vec <- vec[-c(which(vec == 0))]
    }

    # Calculate the DR
    dr_spp <- (sum(vec * sapply(length(vec):1, function(z){return(1/(2^(z-1)))})))^-1
    return(dr_spp)
  })

  # Saving it into a data.frame
  drs <- data.frame(Species = tree$tip.label, DR = drs)
  # Return the DRs of all tree species
  return(drs)
}
