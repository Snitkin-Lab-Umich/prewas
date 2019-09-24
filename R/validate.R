#' Check if the input to the function is a number.
#'
#' @param num Number.
#'
#' @export
check_is_number <- function(num){
  if (class(num) %in% c("numeric")) {
    stop("Input must be a numeric")
  }
}


check_inputs <- function(inputs){


}

read_gff <- function(gff_path){
  gff <- read.table(gff_path,
                    sep = '\t',
                    header = FALSE,
                    stringsAsFactors = FALSE,
                    quote = "",
                    comment.char = '#'
  )
  return(gff)
}
