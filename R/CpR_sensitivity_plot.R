#### f12) CpR_sensitivity_plot() -----------------------------------------------------
##  Make a plot of the sensitivity analysis output obtained from the function CpB_sensitivity().

# sst_output = A dataframe outputed from the CpB_sensitivity();
# stc = a character vector defining a descriptive statistic to make the plot (default is "mean").
# rate = need to provide the CpRate desired (i.e., "CpD", "CpE", "CpB_RW", "CpB");

#' Title
#'
#' @param sst_output
#' @param rate
#' @param stc
#'
#' @return
#' @export
#'
#' @examples

CpR_sensitivity_plot <- function(sst_output, rate = NULL, stc = "mean"){

  if(is.null(rate) == TRUE){
    stop("The user must inform the rate inputted to display it graphically")
  } else {
    if(rate == "CpD"){
      lab_R <- bquote(.(stc) ~  CpD["rate"])
    }
    if(rate == "CpE"){
      lab_R <- bquote(.(stc) ~  CpE["rate"])
    }
    if(rate == "CpB"){
      lab_R <- bquote(.(stc) ~  CpB["rate"])
    }
    if(rate == "CpB_RW"){
      lab_R <- bquote(.(stc) ~  CpB["rate"])
    }
  }

  # Removing NA's from the outputs
  if(sum(is.na(sst_output)) > 0){
    # Finding the rows contaning NA's and removing them
    sst_output <- sst_output[-(which(rowSums(is.na(sst_output)) > 0)),]
  }

  # Creating a df with summarizing statistic informations
  val <- apply(sst_output, 2, stc)
  vec <- as.numeric(colnames(sst_output))
  df <- data.frame(val, vec)

  # Plotting those informations
  return(ggplot2::ggplot(df, ggplot2::aes(x = vec)) +
           ggplot2::geom_line(ggplot2::aes(y = val), colour = "grey50", size = 1.2) +
           ggplot2::geom_point(ggplot2::aes(y = val), col = "black", fill = "#000000", size = 2, shape = 21) +
           ggplot2::labs(x = "Number of phylogenetic slices", y = lab_R) +
           ggplot2::scale_x_continuous(breaks = vec,
                                       labels = vec) +
           ggplot2::theme_classic() +
           ggplot2::theme(axis.text.y = ggplot2::element_text(colour = "black", size=9),
                          axis.text.x = ggplot2::element_text(colour = "black", size=9),
                          axis.title = ggplot2::element_text(size=14),
                          axis.ticks = ggplot2::element_line(colour = "black")))
}
