#### f1) nodes_config() ----------------------------------------------------------
# Capture nodes and branches configurations.

# Arguments:
# tree = a phylogenetic tree of the phylo class;

## Capturing nodes and tree information
#' Title
#'
#' @param tree
#'
#' @return
#' @export
#'
#' @examples

nodes_config <- function(tree){

  # if the tree is ultrametric, make the evaluations
  if(ape::is.ultrametric(tree) == TRUE){
    ## Capturing nodes distances and configurations
    matrix_nodes <- ape::dist.nodes(tree)
    # Putting it into my tree
    tree$config <- data.frame(NodeBegin = tree$edge[, 1],
                              NodeEnd = tree$edge[, 2],
                              NodeLength = tree$edge.length,
                              YearBegin = matrix_nodes[tree$edge[, 2], tree$edge[1, 1]] - tree$edge.length,
                              YearEnd = matrix_nodes[tree$edge[, 2], tree$edge[1, 1]])

    ## Creating a node path matrix
    node_mat <- matrix(0, ncol = length(tree$tip.label), nrow = tree$Nnode)
    # Naming its columns as species and its rows as nodes
    row.names(node_mat) <- unique(tree$edge[,1])
    colnames(node_mat) <- tree$tip.label

    # Creating a presence-abscense nodes-matrix
    paths <- ape::nodepath(tree)
    for(i in 1:length(tree$tip.label)){
      node_mat[which(row.names(node_mat) %in% paths[[i]]), i] <- 1
    }
    # Adding it to the tree
    tree$node_matrix <- node_mat

    # Calculating the total tree depth
    tree$tree_depth <- matrix_nodes[tree$edge[1, 1], 1] # distance from root to the species 1 (Just for ultrametric trees)

    return(tree)
  } else {

    stop("The inputted phylogenetic tree is not ultrametric.")
  }
}
