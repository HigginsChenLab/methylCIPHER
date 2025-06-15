#' trafo
#'
#' @param x A vector of sample ages
#' @param adult.age Age boundary of adulthood--set to 20
#'
#' @return a vector of transformed ages
#' @export
trafo= function(x,adult.age= 20) {
  x=(x+1)/(1+adult.age); y=ifelse(x<=1, log( x),x-1);y
}

#' Developmental Age Transformation
#'
#' @param x A vector of sample ages
#' @param adult.age The age considered to be the cutoff for adulthood
#'
#' @return transformed age prediction
#' @export
anti.trafo= function(x,adult.age=20) {
  ifelse(x<0, (1+adult.age)*exp(x)-1, (1+adult.age)*x+adult.age)
}
