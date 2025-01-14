#### f3) squeeze_root() ----------------------------------------------------------
# Make a phylogenetic slice following a rootwards orientation (roots -> tips);

# Arguments:
# tree = a phylogenetic tree of the phylo class;
# time = numeric temporal threshold at which the root will be squeezed (pd or time threshold);
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

squeeze_root <- function(tree, time, criteria = "my", dropNodes = FALSE){

  ## CUTTING THE PHYLOGENY ROOTWARDLY (FROM ROOTS -> TIPS) ##

  # Getting the nodes configurations
  tree <- nodes_config(tree)

  # If the criteria for cutting is based on time:
  if(criteria == "my"){

    if(time < max(tree$config$YearEnd)){
      # Capturing the tree length (cutting end)
      j <- max(tree$config$YearEnd)

      if(time > 0){

        ### Cutting the phylogeny ROOTWARDLY (ROOT -> TIPS):
        # The threshold occupy different nodes with different requirements for node slicing?
        # (ex: several nodes starting before and after a given threshold)
        if(length(unique(c(which(sort(unique(tree$config$YearEnd)) >= j)[1], which(sort(unique(tree$config$YearEnd)) >= (j-time))[1]))) > 1){ # + de 1 TRUE

          ## Correcting the length of those nodes that ends inside the given interval,
          # but start before the threshold established in (j - n):
          tree$edge.length[tree$config$YearEnd >= (j-time) & tree$config$YearEnd <= (j) & tree$config$YearBegin < (j-time)] <-
            (tree$config$YearEnd[tree$config$YearEnd >= (j-time) & tree$config$YearEnd <= (j) & tree$config$YearBegin < (j-time)] - (j-time))

          ## All nodes ending before my threshold (j-n), turn into zero their lengths
          tree$edge.length[tree$config$YearEnd < (j-time)] <- 0

          ## All nodes originating befores my threshold (j - n) and ending after my threshold (j)
          # (an entire branch length within the interval), assign (n) to its edge.length
          tree$edge.length[tree$config$YearBegin < (j-time) & tree$config$YearEnd >= j] <- time
        }

        ## If all remaining nodes had its edges bigger than the interval:
        if(length(unique(c(which(sort(unique(tree$config$YearEnd)) >= j)[1], which(sort(unique(tree$config$YearEnd)) >= (j-time))[1]))) == 1){ # ! de 0, porque ele n?o foi cortado; o bra?o ?ntegro (sem recorte) que o valor de recorte

          ## Turn into 0 those edges ending before my threshold
          tree$edge.length[which(tree$config$YearEnd < (j))] <- 0

          ## Those remaining nodes, turn into n their edge.lengths
          tree$edge.length[tree$edge.length > 0] <- time
        }
      }
    } else {
      stop("The threshold inputted is bigger than the available by the phylogeny")
    }
  }


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
      rmv <- (node_info[, 5] - time)/node_info[, 2]             # time <- 5
      # How much time i will need to remove? (based on the PD)
      time <- node_info[, 4] - rmv # it is where the node begin

      # cutting ending
      j <- max(tree$config$YearEnd)

      # Cutting
      if(time > 0){

        ### Cutting the phylogeny ROOTWARDLY (ROOT -> TIPS):
        # The threshold occupy different nodes with different requirements for node slicing?
        # (ex: several nodes starting before and after a given threshold)
        if(length(unique(c(which(sort(unique(tree$config$YearEnd)) >= j)[1], which(sort(unique(tree$config$YearEnd)) >= time)[1]))) > 1){ # + de 1 TRUE

          ## Correcting the length of those nodes that ends inside the given interval,
          # but start before the threshold established in (j - n):
          tree$edge.length[tree$config$YearEnd >= time & tree$config$YearEnd <= (j) & tree$config$YearBegin < time] <-
            (tree$config$YearEnd[tree$config$YearEnd >= time & tree$config$YearEnd <= (j) & tree$config$YearBegin < time] - time)

          ## All nodes ending before my threshold (j-n), turn into zero their lengths
          tree$edge.length[tree$config$YearEnd < time] <- 0

          ## All nodes originating befores my threshold (j - n) and ending after my threshold (j)
          # (an entire branch length within the interval), assign (n) to its edge.length
          tree$edge.length[tree$config$YearBegin < time & tree$config$YearEnd >= j] <- j-time
        }

        ## If all remaining nodes had its edges bigger than the interval:
        if(length(unique(c(which(sort(unique(tree$config$YearEnd)) >= j)[1], which(sort(unique(tree$config$YearEnd)) >= time)[1]))) == 1){ # ! de 0, porque ele n?o foi cortado; o bra?o ?ntegro (sem recorte) que o valor de recorte

          ## Turn into 0 those edges ending before my threshold
          tree$edge.length[which(tree$config$YearEnd < (j))] <- 0

          ## Those remaining nodes, turn into n their edge.lengths
          tree$edge.length[tree$edge.length > 0] <- j-time
        }
      }
    } else {
      stop("The threshold inputted is bigger than the available by the phylogeny")
    }
  }

  ## If there are some void nodes/edges inside our algorithm, remove them
  if(dropNodes == TRUE){

    # Which are our void nodes with 0 length?
    rm_nd <- which(!(tree$edge[,2] %in% c(1:length(tree$tip.label))) & tree$edge.length == 0) # tree <- teste

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

  return(tree) # Subtrai o comprimento observado de cada n? pela diferen?a
}
