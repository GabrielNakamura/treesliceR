#### f11) CpR_sensitivity() -----------------------------------------------------

### Calculates the sensitivity of the CpB-rates based on this index calculated
# through different time-slices.

# tree = a phylogenetic tree of the phylo class;
# vec = a numeric vector containing a serie of number of slices;
# mat = a complete presence absence matrix of all studied species and sites;
# asb = list with assemblages and its adjacent cells;
# rate = need to provide the CpRate desired (i.e., "CpD", "CpE", "CpB_RW", "CpB");
# samp = number of sites to be sampled to make the sensitivity analysis;
# comp = component of the beta-diversity if the user wants the CpB (Deafault is
# "sorensen", but "turnover" and "nestedness" can be set in);
# method = method for calculation of the CpB or CpB_RW, which can be "pairwise"
# or "multisite" (which is the default);
# criteria = temporal criteria for making phylogenetic slices ("pd" or "my");
# ncor = number of cores the user wants to parallelize.


CpR_sensitivity <- function(tree, vec, mat = NULL, asb = NULL, rate = NULL, samp = 0,
                            comp = "sorensen", method = "multisite",
                            criteria = "my", ncor = 0){

  # The user provided the rate to be tested?
  if(is.null(rate) == TRUE){
    stop("The desired rate must be provided on rate argument")
  } else {
    # The user provided the number of sites to be sampled?
    if(samp == 0){
      stop("The number of sites for running the sensitivity analysis was not provided")
    } else {

      # Creating an empty "n"
      n <- NULL

      ## Cumulative Phylogenetic Diversity (CpD) -------------------------------

      if(rate == "CpD"){
        # If the samp inputed is bigger than the number of sites, return a error
        if(samp > nrow(mat)){
          stop("The number of site samples inputted in bigger than the sites available")
        } else {
          # Sampling random sites to run the algorithm
          sites <- sample(1:nrow(mat), samp)

          if(length(sites) == 1){
            # Running the sensitivity algorithm
            sensitivity <- foreach::foreach(n = vec, .combine = cbind) %do% {
              output <- CpD(tree, n, mat[sites,], criteria = criteria, ncor = ncor)[1]
              return(output)
            }
            # Renaming the sensitivity data.frame object
            colnames(sensitivity) <- as.character(vec)
          } else {
            # Running the sensitivity algorithm
            sensitivity <- foreach::foreach(n = vec, .combine = cbind) %do% {
              output <- CpD(tree, n, mat[sites,], criteria = criteria, ncor = ncor)[,1]
              return(output)
            }
            # Renaming the sensitivity obj
            colnames(sensitivity) <- as.character(vec)
            rownames(sensitivity) <- 1:nrow(sensitivity)
          }
        }
      }

      ## Cumulative Phylogenetic Endemism (CpE) --------------------------------

      if(rate == "CpE"){
        # If the samp inputed is bigger than the number of sites, return a error
        if(samp > nrow(mat)){
          stop("The number of site samples inputted in bigger than the sites available")
        } else {
          # Sampling random sites to run the algorithm
          sites <- sample(1:nrow(mat), samp)

          if(length(sites) == 1){
            # Running the sensitivity algorithm
            sensitivity <- foreach::foreach(n = vec, .combine = cbind) %do% {
              output <- CpE(tree, n, mat[sites,], criteria = criteria, ncor = ncor)[1]
              return(output)
            }
            # Renaming the sensitivity data.frame object
            colnames(sensitivity) <- as.character(vec)

          } else {
            # Running the sensitivity algorithm
            sensitivity <- foreach::foreach(n = vec, .combine = cbind) %do% {
              output <- CpE(tree, n, mat[sites,], criteria = criteria, ncor = ncor)[,1]
              return(output)
            }
            # Renaming the sensitivity obj
            colnames(sensitivity) <- as.character(vec)
            rownames(sensitivity) <- 1:nrow(sensitivity)
          }
        }
      }

      ## Cumulative Phylogenetic B-Diversity (CpB) -----------------------------

      if(rate == "CpB"){

        # Checking if a single matrix or a list of matrixes was provided
        if(sum(class(asb) != "list") >= 1){
          asb <- list(asb)  # Transforming the list into a single matrix
        }

        # If the samp inputed is bigger than the number of sites, return a error
        if(samp > length(asb)){
          stop("The number of site samples inputted in bigger than the sites available")
        } else {
          # Sampling random sites to run the algorithm
          sites <- sample(1:length(asb), samp)

          if(length(sites) == 1){
            # Running the sensitivity algorithm
            sensitivity <- foreach::foreach(n = vec, .combine = cbind) %do% {
              output <- CpB(tree, n, asb[sites], method = method, comp = comp,
                            criteria = criteria, ncor = ncor)[1]
              return(output)
            }
            # Renaming the sensitivity data.frame object
            colnames(sensitivity) <- as.character(vec)

          } else {
            # Running the sensitivity algorithm
            sensitivity <- foreach::foreach(n = vec, .combine = cbind) %do% {
              output <- CpB(tree, n, asb[sites], method = method, comp = comp,
                            criteria = criteria, ncor = ncor)[,1]
              return(output)
            }
            # Renaming the sensitivity obj
            colnames(sensitivity) <- as.character(vec)
            rownames(sensitivity) <- 1:nrow(sensitivity)
          }
        }
      }

      ## Cumulative Phylogenetic B-Diversity Range-Weighted (CpB RW) -----------

      if(rate == "CpB_RW"){

        # Checking if a single matrix or a list of matrixes was provided
        if(sum(class(asb) != "list") >= 1){
          asb <- list(asb)  # Transforming the list into a single matrix
        }

        # If the samp inputed is bigger than the number of sites, return a error
        if(samp > length(asb)){
          stop("The number of site samples inputted in bigger than the sites available")
        } else {
          # Sampling random sites to run the algorithm
          sites <- sample(1:length(asb), samp)

          if(length(sites) == 1){
            # Running the sensitivity algorithm
            sensitivity <- foreach::foreach(n = vec, .combine = cbind) %do% {
              output <- CpB_RW(tree, n, mat, asb[sites], method = method,
                               criteria = criteria, ncor = ncor)[1]
              return(output)
            }
            # Renaming the sensitivity data.frame object
            colnames(sensitivity) <- as.character(vec)

            } else {
              # Running the sensitivity algorithm
              sensitivity <- foreach::foreach(n = vec, .combine = cbind) %do% {
                output <- CpB_RW(tree, n, mat, asb[sites], method = method,
                                 criteria = criteria, ncor = ncor)[,1]    ## PRECISO TER UMA CONDIÇÃO POR AQUI PARA IDENTIFICAR SE EU TENHO SÓ UM SITE
                return(output)
              }
              # Renaming the sensitivity obj
              colnames(sensitivity) <- as.character(vec)
              rownames(sensitivity) <- 1:nrow(sensitivity)
          }
        }
      }

      # Returning it
      return(sensitivity)
    }
  }
}
