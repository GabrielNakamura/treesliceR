#### f6) CpD() -----------------------------------------------------------------

# Calculates the rate of Cumulation of phylogenetic B-diversity.

# tree = a phylogenetic tree of the phylo class;
# n = number of temporal slices (method = 1) or time interval;
# mat = a complete presence absence matrix of all studied species and sites;
# criteria = temporal criteria for slices ("pd" or "my");
# pDO = numeric proportion to define the temporal origin that the phylogenetic
# diversity starts to accumulate on a given site. The default is 5%;
# ncor = number of cores the user wants to parallelize.


CpD <- function(tree, n, mat, criteria = "my", pDO = 5, ncor = 0){

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

  ## Calculating assemblages CpD, PD and pDO -----------------------------------
  # The user wants to use more CPU cores?
  if(ncor > 0){
    # Register the number of desired clusters
    doParallel::registerDoParallel(ncor)
    # Loop and capture the values for each assemblage
    CPD <- foreach::foreach(i = 1:nrow(mat), .combine = rbind) %dopar% {

      # Which species are within this assemblage
      tips <- which(mat[i,] > 0) # (spp numbers follows same position on node matrix)  i <- 1000

      if(length(tips) == 0) {
        CPDrate <- c(NA, NA, NA)
        return(CPDrate)

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
        CPD <- sapply(branch_pieces, function(x){     # x <- branch_pieces[[1800]]
          # Calculating the PE stored on each tree slice
          return(sum(x$edge.length[positions]))
        })

        # Saving the CPD
        CPD <- cumsum(CPD)/sum(CPD)  #  CPD[1809]

        # Fitting the non-linear model to obtain the CpB-rate:
        if(is.na(CPD[1]) == FALSE){ # Checking if there is any NA emerging due 0/0 fraction
          nonlinear_cum <- stats::nls(CPD ~ exp(-r*age),
                                      start = list(r = 0.2),
                                      control = stats::nls.control(maxiter = 1000))   # plot(CPE ~ age)
          # Saving the exponentital coefficient
          CPDrate <- stats::coef(nonlinear_cum)
        } else {
          # If it cant be calculated, return NA
          CPDrate <- NA
        }

        # Obtaining the PE
        PD <- sum(tree$edge.length[positions])

        # Creating a vec output containing the CpE-rate, pEO
        CPDrate <- c(CPDrate, PD, ((-1)*(log(pDO/100)/CPDrate)))
        return(CPDrate)
      }
    }
    # stop the clusters after running the algorithm
    doParallel::stopImplicitCluster()
  }

  # The user do not set clusters for being used
  if(ncor == 0){
    # Loop and capture the values for each assemblage
    CPD <- foreach::foreach(i = 1:nrow(mat), .combine = rbind) %do% {

      # Which species are within this assemblage
      tips <- which(mat[i,] > 0) # (spp numbers follows same position on node matrix)  i <- 1000

      if(length(tips) == 0) {
        CPDrate <- c(NA, NA, NA)
        return(CPDrate)

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
        CPD <- sapply(branch_pieces, function(x){     # x <- branch_pieces[[1800]]
          # Calculating the PE stored on each tree slice
          return(sum(x$edge.length[positions]))
        })

        # Saving the CPD
        CPD <- cumsum(CPD)/sum(CPD)  #  CPD[1809]

        # Fitting the non-linear model to obtain the CpB-rate:
        if(is.na(CPD[1]) == FALSE){ # Checking if there is any NA emerging due 0/0 fraction
          nonlinear_cum <- stats::nls(CPD ~ exp(-r*age),
                                      start = list(r = 0.2),
                                      control = stats::nls.control(maxiter = 1000))   # plot(CPE ~ age)
          # Saving the exponentital coefficient
          CPDrate <- stats::coef(nonlinear_cum)
        } else {
          # If it cant be calculated, return NA
          CPDrate <- NA
        }

        # Obtaining the PE
        PD <- sum(tree$edge.length[positions])

        # Creating a vec output containing the CpE-rate, pEO
        CPDrate <- c(CPDrate, PD, ((-1)*(log(pDO/100)/CPDrate)))
        return(CPDrate)
      }
    }
  }

  # Renaming the columns and rows, and returning them as output
  if(nrow(mat) > 1){
    colnames(CPD) <- c("CpD", "PD", "pDO")
    rownames(CPD) <- 1:nrow(mat)
  } else {
    names(CPD) <- c("CpD", "PD", "pDO")
  }

  return(CPD)
}
