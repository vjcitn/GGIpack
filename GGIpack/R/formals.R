# any ggiResource should explicitly have seqnames, score,
# start, end, space (Genome tag, like hg19)

#' manage a GGI resource which must have start, end, seqnames, space, score
#' it also can have a tbl element
#' @export
setClass("ggiResource", slots=c(start="numeric",
 end="numeric", seqnames="ANY", space="ANY",
 score="numeric", tbl="ANY"))

setMethod("show", "ggiResource", function(object) {
 callNextMethod()
})

#' ABRIGresource is tailored to ABRIG data resources
#' @export
setClass("ABRIGresource", contains="ggiResource")

setMethod("show", "ABRIGresource", function(object) {
 cat("CDNM ggiResource\n")
 print(slot(object, "tbl"))
})


#' GTEx resource is tallored to the  wholeblpl05.parquet and lungpl05.parquet  resources.
#' @export
setClass("GTExresource", contains="ggiResource")

#' GTEx resource is tallored to the  wholeblpl05.parquet and lungpl05.parquet  resources.
#' @export
setClass("GTExresource", contains="ggiResource")

#' constructor for GTEx examples
#' @param con DBI connection
#' @param space character(1) must indicate build
#' @param pfile character(1) path to parquet file
#' @export
GTExresource = function (con, space = "hg19", pfile) 
{
    pp = sprintf("read_parquet(%s)", sQuote(pfile))
    tb = tbl(con, pp)
    ans = dplyr::mutate(tb, score = pvalue, seqnames = chromosome)
    new("GTExresource", space = space, tbl = ans)
}

