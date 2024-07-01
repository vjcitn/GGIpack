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
