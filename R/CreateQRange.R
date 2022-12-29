#' Create q_range object for PreparePlinkFile()
#' @description
#' Calls CreateQRange() with either the default pthres or user input values
#' Not exported
#' @param x pthres
#' @return CreateQRange() object
#'
CreateQRange <- function(x)
{
  print("executing CreateQRange()... ")
  if (is.null(x)) {
    pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,
                5E-02,5E-01,1.0)
    print(paste0("pthres is NULL...using default values "))
  } else {
    pthres <- x
  }

  n_pthres <- length(pthres)
  q_range <- data.frame(filename = rep("p_tar",n_pthres),
                        small_P = rep(0,n_pthres),
                        max_P = rep(1,n_pthres))

  temp <- 1

  for(k2 in 1:length(pthres)){

    q_range[temp,1] = paste0("p_tar_",k2)
    q_range[temp,3] = pthres[k2]
    temp = temp+1
  }
  q_range = q_range[1:(temp-1),]
  print("CreateQRange() complete... ")
  return(q_range)
}
