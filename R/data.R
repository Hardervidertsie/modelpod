

#' Typical dose response GFP reporter data
#'
#' Data from a recent publication
#'
#'
#' @format A data frame with 911 rows and 4 variables:
#' \describe{
#'   \item{plateID}{plate id functions as replicate id in this case, 3 unique ids}
#'   \item{treatment}{set of treatments, with DMSO as negative control}
#'   \item{dose_uM}{concentrations}
#'   \item{GFP_int}{GFP intensity levels, fraction of GFP positive cells}
#' }
#' @source \url{DDS LACDR LU}
"test_data"


#' Typical time lapse dose response GFP reporter data
#'
#' Data from a recent publication
#'
#'
#' @format A data frame with 29852 rows and 9 variables:
#' \describe{
#'   \item{plateID}{plate id functions as replicate id in this case, 3 unique ids}
#'   \item{treatment}{set of treatments, with DMSO as negative control}
#'   \item{dose_uM}{concentrations}
#'   \item{GFP_int}{response variable, GFP intensity levels}
#'   \item{cell_line}{srxn1 cell line}
#'   \item{replID}{replicate id}
#'   \item{timeID}{time id}
#'   \item{cell_count}{response variable, number of cells}
#'   \item{GFP_2m}{response variable, fraction of GFP positive cells}
#' }
#' @source \url{DDS LACDR LU}
"summary_data"