#### f2) squeeze_tips() ----------------------------------------------------------
# Make a phylogenetic slice following a tipwards orientation (tips -> roots);
# How many years you want to collapse tipwardly?
# How many PD you want to collapse tipwardly?

# Arguments:
# tree = a phylogenetic tree of the phylo class;
# time = numeric temporal threshold at which the tips will be squeezed (pd or time threshold);
# criteria = temporal criteria to make the slices, it can be "my" (DEFAULT) or "pd";
# dropNodes = remove those void nodes without information (TRUE or FALSE).

#' Title
#'
#' @param tree
#' @param time
#' @param criteria
#' @param dropNodes
#'
#' @return
#' @export
#'
#' @examples

squeeze_tips <- function(tree, time, criteria = "my", dropNodes = FALSE){

  ## CUTTING THE PHYLOGENY TIPWARDLY (FROM TIPS -> ROOTS) ##

  # Getting the nodes configurations
  tree <- nodes_config(tree)

  # If the criteria for cutting is based on time:
  if(criteria == "my"){

    if(time < max(tree$config$YearEnd)){
      time <- max(tree$config$YearEnd) - time # capture the tree length

      ## Cutting the phylogeny TIPWARDLY (TIPS -> ROOT):
      # Which nodes ends after j but are SMALLER or bigger than my thresholds?
      bigger <- which((tree$config$YearEnd[which(tree$config$YearEnd > time)] - time) >= tree$edge.length[which(tree$config$YearEnd > time)])
      smaller <- which((tree$config$YearEnd[which(tree$config$YearEnd > time)] - time) < tree$edge.length[which(tree$config$YearEnd > time)])

      # For those smaller, turn their edge.length 0 (ZERO)
      tree$edge.length[which(tree$config$YearEnd > time)][bigger] <- 0

      # Those values bigger than my threshold when subtracting j,
      # remove its remaining difference to reach the threshold
      tree$edge.length[which(tree$config$YearEnd > time)][smaller] <-
        tree$edge.length[which(tree$config$YearEnd > time)][smaller] - (tree$config$YearEnd[which(tree$config$YearEnd > time)][smaller] - time)
    } else {
      stop("The threshold inputted is bigger than the available by the phylogeny")
    }
  }

  # If the criteria for cutting is based on PD:
  if(criteria == "pd"){

    if(time < sum(tree$edge.length)){
      ## How many PD we have at each nodes spliting?
      # Separating only the nodes inital information
      nodes <- tree$config[!duplicated(tree$config$NodeBegin), c(1, 2, 4)]

      ## Making a matrix with these informations
      df <- data.frame(time = sort(unique(nodes$YearBegin)), # Time which init the node
                       nBranch = 2:(length(unique(nodes$YearBegin)) + 1), # Number of branches
                       timeLength = c(sort(unique(nodes$YearBegin)), max(tree$config$YearEnd))[-1] - # time length of the node
                         sort(unique(nodes$YearBegin)))

      # Calculating the cumulative TIME and PD per node split
      df$cumulativeTIME <- cumsum(df$timeLength)
      df$cumulativePD <- cumsum(df$nBranch * df$timeLength)

      # Separating the node information
      node_info <- df[which(df$cumulativePD >= time)[1],]
      # How much of PD i will need to remove for each node
      rmv <- (node_info[, 5] - time)/node_info[, 2]
      # How much time i will need to remove? (based on the PD)
      time <- node_info[, 4] - rmv

      ## Cutting the phylogeny TIPWARDLY (TIPS -> ROOT):
      # Which nodes ends after j but are SMALLER or bigger than my thresholds?
      bigger <- which((tree$config$YearEnd[which(tree$config$YearEnd > time)] - time) >= tree$edge.length[which(tree$config$YearEnd > time)])
      smaller <- which((tree$config$YearEnd[which(tree$config$YearEnd > time)] - time) < tree$edge.length[which(tree$config$YearEnd > time)])

      # For those smaller, turn their edge.length 0 (ZERO)
      tree$edge.length[which(tree$config$YearEnd > time)][bigger] <- 0

      # Those values bigger than my threshold when subtracting j,
      # remove its remaining difference to reach the threshold
      tree$edge.length[which(tree$config$YearEnd > time)][smaller] <-
        tree$edge.length[which(tree$config$YearEnd > time)][smaller]-(tree$config$YearEnd[which(tree$config$YearEnd > time)][smaller] - time)
    } else {
      stop("The threshold inputted is bigger than the available by the phylogeny")
    }
  }

  ## If there are some void nodes/edges inside our algorithm, remove them
  if(dropNodes == TRUE){

    # Which are our void nodes with 0 length?
    rm_nd <- which(!(tree$edge[,2] %in% c(1:length(tree$tip.label))) & tree$edge.length == 0)

    if(length(rm_nd) > 0){
      # Correcting their values
      for (i in 1:length(rm_nd)) { # i <- 1
        # which node comes after our node
        val <- tree$edge[rm_nd[i], 2]
        # all places with this value, need to be turned into its previous node
        tree$edge[which(tree$edge[,1] == val), 1] <- tree$edge[rm_nd[i], 1]
      }

      # Them, removing those edges and edgelenghts that are 0 and are node-to-node
      tree$edge <- tree$edge[-c(rm_nd), ]               # Removing from the edge matrix
      tree$edge.length <- tree$edge.length[-c(rm_nd)]   # Removing from the edge length vectors
      tree$Nnode <- length(unique(tree$edge[,1]))       # Adding the new number of nodes

      # Unique tree nodes
      oldnodes <- sort(unique(tree$edge[,1]))

      # Renaming the tree nodes
      for (i in 1:length(oldnodes)) { # i <- 1
        tree$edge[which(tree$edge[,1] == oldnodes[i]), 1] <- length(tree$tip.label) + i
        tree$edge[which(tree$edge[,2] == oldnodes[i]), 2] <- length(tree$tip.label) + i
      }
    }
  }

  return(tree)
}
