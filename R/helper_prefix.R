#' Helper function
#' @description
#' This function creates the plink out prefix for clumping, scoring and EBayes
#' scoring. Not exported
#' @param out_prefix  description
#' @param ebayes Default FALSE
#' @keywords plink prefix

helper_prefix <- function(out_prefix,
                          ebayes=FALSE){
  if (is.null(out_prefix)) {
    prefix <- paste0("")
    eb_prefix <- paste0("eb_")

    if (ebayes) {
      return(eb_prefix)
    }else{
      return(prefix)
    }

  } else {
    prefix <- paste0(out_prefix, "_")
    eb_prefix <- paste0(out_prefix, "_eb_")

    if (ebayes) {
      return(eb_prefix)
    }else{
      return(prefix)
    }
  }
}
