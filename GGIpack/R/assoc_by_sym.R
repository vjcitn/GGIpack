#' produce data.frame for GTExresource content filtered for a specific gene
#' @param gtexres GTExresource instance
#' @param sym character(1) gene symbol
#' @examples
#' lu = ggi_gtex_cache("lungpl05.parquet")
#' con = DBI::dbConnect(duckdb::duckdb())
#' lures = GTExresource(con, tisstag="lung", pfile=lu)
#' aa = assoc_by_sym( lures )
#' DBI::dbDisconnect(con)
#' @export
assoc_by_sym = function( gtexres, sym = "ORMDL3" ) {
 # assumes number of associations for a gene is manageable as data.frame
 stopifnot(is(gtexres, "GTExresource"))
 data("gloc_hg19", package="GGIpack")
 mid = which(gloc_hg19$symbol == sym)
 stopifnot(length(mid)>=1)
 mid = names(gloc_hg19[mid[1]]) # should message if length(mid)>1
 tmp = gtexres@tbl |> dplyr::filter(molecular_trait_id == mid) |> as.data.frame()
 dplyr::mutate(tmp, genesym=sym)
}
