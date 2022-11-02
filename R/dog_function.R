#' A Dog Function
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords dogs
#' @export
#' @examples
#' dog_function()

dog_function <- function(love=TRUE){
  if(love==TRUE){
    print("Liliuokalani rules!")
  }
  else {
    print("I am not a cool person.")
  }
}
