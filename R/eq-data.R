#' Earthquake data
#'
#' This data set contains simulated earthquakes, providing date and time,
#' location, and magnitude of each event.
#'
#' @docType data
#'
#' @usage data(eq)
#'
#' @format An object of class \code{"cross"}; see \code{\link[qtl]{read.cross}}.
#'
#' @keywords datasets
#'
# #' @references Moore et al. (2013) Genetics 195:1077-1086
# #' (\href{https://www.ncbi.nlm.nih.gov/pubmed/23979570}{PubMed})
# #'
# #' @source \href{https://phenome.jax.org/projects/Moore1b}{QTL Archive}
# #'
# #' @examples
# #' data(grav)
# #' times <- attr(grav, "time")
# #' phe <- grav$pheno
# #' \donttest{iplotCurves(phe, times)}
"eq"
