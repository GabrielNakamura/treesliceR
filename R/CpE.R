#### f7) CpE() -----------------------------------------------------------------

# Calculates the rate of Cumulation of phylogenetic Endemism.

# tree = a phylogenetic tree of the phylo class;
# n = number of temporal slices (method = 1) or time interval;
# mat = a complete presence absence matrix of all studied species and sites;
# criteria = temporal criteria for slices ("pd" or "my");
# pEO = numeric proportion to define the temporal origin that the phylogenetic diversity
# starts to accumulate on a given site. The default is 5%;
# ncor = number of cores the user wants to parallelize.


#' Title
#'
#' @param tree
#' @param n
#' @param mat
#' @param criteria
#' @param pEO
#' @param ncor
#'
#' @return
#' @export
#'
#' @examples

CpE <- function(tree, n, mat, criteria = "my", pEO = 5, ncor = 0){

  ## Cleaning the phylogeny (if necessary) and cutting it into pieces ----------

  # If a doesnt have a matrix, but a vector of species,
  # transform it into a species matrix
  if((is.matrix(mat) | is.data.frame(mat)) == FALSE){
    mat <- t(as.matrix(mat))
  }

  # Checking if there is species in the matrix without presence and removing them
  spps_pa <- colSums(mat)
  if(sum(spps_pa == 0) > 0){

    if(nrow(mat) == 1){
      # Removing those species
      mat <- t(as.matrix(mat[, -(which(spps_pa == 0))]))
      warning("Removing the species in presence abscence matrix without any occurrence")
    } else {
      # Removing those species
      mat <- mat[, -(which(spps_pa == 0))]
      warning("Removing the species in presence abscence matrix without any occurrence")
    }
  }

  # Dropping the lineage tips that arent in my spp matrix
  if(all(tree$tip.label %in% colnames(mat)) == FALSE){
    tree <- ape::keep.tip(tree, intersect(tree$tip.label, colnames(mat)))
    warning("Removing tips from phylogeny that are absent on species matrix")
  }

  # Cutting the phylogenetic tree into equal width slices
  branch_pieces <- phylo_pieces(tree, n, criteria = criteria,
                                timeSteps = TRUE, returnTree = TRUE)

  # Separating the time steps from the phylogenetic pieces
  age <- branch_pieces[[2]][length(branch_pieces[[2]]):1]
  tree <- branch_pieces[[3]]
  branch_pieces <- branch_pieces[[1]]
  i <- NULL

  ## Calculating the branch range size within each node branch -----------------

  # If there is only one site, the range of all nodes equals 1
  if(nrow(mat) == 1){
    r_sizes <- rep(1, length(tree$edge[,2]))
  } else {
    r_sizes <- sapply(tree$edge[,2], function(x){
      # If its a tip branch
      if(x < tree$edge[1, 1]){
        # There is only one species within the node, calculate and save its range
        deno <- sum(mat[, tree$tip.label[x]])
      } else {
        # If its a node branch
        # Which species share those nodes
        spps_node <- names(which(tree$node_matrix[as.character(x), ] > 0))  # node 1    x <- 689

        # Capturing the range size denominator from the node
        if(length(spps_node) == 1){
          # If there is only one species within the node
          deno <- sum(mat[, spps_node])
        } else { # If there is more species
          # Which is the shared range size of the species sharing this node
          deno <- sum(rowSums(mat[, spps_node]) > 0)
        }
      }
      return(deno)
    })
  }


  ## Calculating assemblages CpE, PE and pEO -----------------------------------
  # The user wants to use more CPU cores?
  if(ncor > 0){
    # Register the number of desired clusters
    doParallel::registerDoParallel(ncor)
    # Loop and capture the values for each assemblage
    CPE <- foreach::foreach(i = 1:nrow(mat), .combine = rbind) %dopar% {

      # Which species are within this assemblage
      tips <- which(mat[i,] > 0) # (spp numbers follows same position on node matrix)  i <- 1000

      if(length(tips) == 0) {
        CPErate <- c(NA, NA, NA)
        return(CPErate)

      } else {
        # Obtaining the node matrix for those species
        if(length(tips) == 1){
          #  if there is a single spp, which nodes give origin to my species
          nodes <- as.numeric(names(which(rowSums(as.data.frame(tree$node_matrix[, tips])) > 0)))
        }

        if(length(tips) > 1){
          # Which nodes give origin to my species
          nodes <- as.numeric(names(which(rowSums(tree$node_matrix[, tips]) > 0)))
        }

        # Obtaining the tips and node positions
        positions <- which(tree$edge[,2] %in% c(tips, nodes))

        # Calculating the relative PE on each phylo slice
        CPE <- sapply(branch_pieces, function(x){     # x <- branch_pieces[[1800]]
          # Calculating the PE stored on each tree slice
          return(sum(x$edge.length[positions]/r_sizes[positions]))
        })

        # Saving the CPE
        CPE <- cumsum(CPE)/sum(CPE)  #  CPE[1809]

        # Fitting the non-linear model to obtain the CpB-rate:
        if(is.na(CPE[1]) == FALSE){ # Checking if there is any NA emerging due 0/0 fraction
          nonlinear_cum <- stats::nls(CPE ~ exp(-r*age),
                                      start = list(r = 0.2),
                                      control = stats::nls.control(maxiter = 1000))   # plot(CPE ~ age)
          # Saving the exponentital coefficient
          CPErate <- stats::coef(nonlinear_cum)
        } else {
          # If it cant be calculated, return NA
          CPErate <- NA
        }

        # Obtaining the PE
        PE <- sum(tree$edge.length[positions]/r_sizes[positions])

        # Creating a vec output containing the CpE-rate, pEO
        CPErate <- c(CPErate, PE, ((-1)*(log(pEO/100)/CPErate)))
        return(CPErate)
      }
    }
    # stop the clusters after running the algorithm
    doParallel::stopImplicitCluster()
  }

  # The user do not set clusters for being used
  if(ncor == 0){
    # Loop and capture the values for each assemblage
    CPE <- foreach::foreach(i = 1:nrow(mat), .combine = rbind) %do% {

      # Which species are within this assemblage
      tips <- which(mat[i,] > 0) # (spp numbers follows same position on node matrix)  i <- 1

      if(length(tips) == 0) {
        CPErate <- c(NA, NA, NA)
        return(CPErate)

      } else {
        # Obtaining the node matrix for those species
        if(length(tips) == 1){
          #  if there is a single spp, which nodes give origin to my species
          nodes <- as.numeric(names(which(rowSums(as.data.frame(tree$node_matrix[, tips])) > 0)))
        }

        if(length(tips) > 1){
          # Which nodes give origin to my species
          nodes <- as.numeric(names(which(rowSums(tree$node_matrix[, tips]) > 0)))
        }

        # Obtaining the tips and node positions
        positions <- which(tree$edge[,2] %in% c(tips, nodes))

        # Calculating the relative PE on each phylo slice
        CPE <- sapply(branch_pieces, function(x){     # x <- branch_pieces[[1800]]
          # Calculating the PE stored on each tree slice
          return(sum(x$edge.length[positions]/r_sizes[positions]))
        })

        # Saving the CPE
        CPE <- cumsum(CPE)/sum(CPE)  #  CPE[1809]

        # Fitting the non-linear model to obtain the CpB-rate:
        if(is.na(CPE[1]) == FALSE){ # Checking if there is any NA emerging due 0/0 fraction
          nonlinear_cum <- stats::nls(CPE ~ exp(-r*age),
                                      start = list(r = 0.2),
                                      control = stats::nls.control(maxiter = 1000))   # plot(CPE ~ age)
          # Saving the exponentital coefficient
          CPErate <- stats::coef(nonlinear_cum)
        } else {
          # If it cant be calculated, return NA
          CPErate <- NA
        }

        # Obtaining the PE
        PE <- sum(tree$edge.length[positions]/r_sizes[positions])

        # Creating a vec output containing the CpE-rate, pEO
        CPErate <- c(CPErate, PE, ((-1)*(log(pEO/100)/CPErate)))
        return(CPErate)
      }
    }
  }

  # Renaming the columns and rows, and returning them as output
  if(nrow(mat) > 1){
    colnames(CPE) <- c("CpE", "PE", "pEO")
    rownames(CPE) <- 1:nrow(mat)
  } else {
    names(CPE) <- c("CpE", "PE", "pEO")
  }

  return(CPE)
}
