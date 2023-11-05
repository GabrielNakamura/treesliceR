#### f4) squeeze_int() ------------------------------------------------------
# Retrieve a temporal phylogenetic interval from an inputed phylogenetic tree

# Arguments:
# tree = a phylogenetic tree of the phylo class;
# from = numeric temporal threshold at which the root will be squeezed;
# to = numeric temporal threshold at which the tips will be squeezed;
# invert = invert the orientation of the slice. Instead of capturing the interval
# it removes the defined interval and return a root and a tip slice;
# criteria = temporal criteria to make the slices, it can be "my" (DEFAULT) or "pd";
# dropNodes = remove those void nodes without information (TRUE or FALSE).

# Make a phylogenetic slice in both orientation to return a single phylo piece;
#' Title
#'
#' @param tree
#' @param from
#' @param to
#' @param invert
#' @param criteria
#' @param dropNodes
#'
#' @return
#' @export
#'
#' @examples

squeeze_int <- function(tree, from, to, invert = FALSE, criteria = "my", dropNodes = FALSE){

  # The used want a phylogenetic interval?
  if(invert == FALSE){
    # The order for making the slices are important depending on the criteria used
    if(criteria == "my"){
      if(from > to){
        tree <- squeeze_root(tree, from, criteria = criteria, dropNodes = dropNodes)
        tree <- squeeze_tips(tree, to, criteria = criteria, dropNodes = dropNodes)
      } else {
        stop("The thresholds set in arguments [from] and [to] are incompatible")
      }

    }
    if(criteria == "pd"){
      if(from < to){
        tree <- squeeze_tips(tree, to, criteria = criteria, dropNodes = dropNodes)
        tree <- squeeze_root(tree, from, criteria = criteria, dropNodes = dropNodes)
      } else {
        stop("The thresholds set in arguments [from] and [to] are incompatible")
      }

    }
    return(tree)
  }

  # or he wants to remove a phylogenetic interval?
  if(invert == TRUE){

    if(criteria == "my"){
      if(from < to){
        stop("The thresholds set in arguments [from] and [to] are incompatible")
      }
    } else {
      tree1 <- squeeze_root(tree, to, criteria = criteria, dropNodes = dropNodes)
      tree2 <- squeeze_tips(tree, from, criteria = criteria, dropNodes = dropNodes)
    }

    if(criteria == "pd"){
      if(from > to){
        stop("The thresholds set in arguments [from] and [to] are incompatible")
      }
    } else {
      tree1 <- squeeze_root(tree, to, criteria = criteria, dropNodes = dropNodes)
      tree2 <- squeeze_tips(tree, from, criteria = criteria, dropNodes = dropNodes)
    }

    return(list(tree1, tree2))
  }
}
