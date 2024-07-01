#' define constructor for ABRIGresource 
#' @import dplyr
#' @param con is a DBI connection (typically duckdb)
#' @param tissue character(1)
#' @param space character(1) e.g., "hg19"
#' @examples
#' con = DBI::dbConnect(duckdb())
#' ll = ABRIGresource( con, "BAL" )
#' @return ABRIGresource instance for the identified tissue
#' @export
ABRIGresource = function(con, tissue, space="hg19") {
   ttypes = c("BAL", "BroncEpiBrush", "CD4Stim", "CD4Unstim",
                  "AlvMacphage", "PaxRNA")
   stopifnot(tissue %in% ttypes)
   ans = dplyr::tbl(con, pfiles[tissue]) |> mutate(score=FDR, seqnames=CHR)
   new("ABRIGresource", space=space, tbl=ans)
}

#' filter a GGI (ABRIGresource) instance by range
#' @param res ABRIGresource instance
#' @param ggr GenomicRanges instance
#' @param tag character(1) value in `ggr_field` used for filtering, e.g., a gene
#' @param radius numeric(1) flanking region size
#' @param ggr_field character(1) metadata element in mcols(ggr) for filtering
#' @export
filterByRange = function(res, ggr, tag, radius=1e5, ggr_field="gene_name") {
  stopifnot(inherits(res, "ABRIGresource"))
  ok = ggr[which(mcols(ggr)[[ggr_field]] == tag)]
  stopifnot(width(ok)>0)
  anac = function(x) as.numeric(as.character(x)) # for Rle
  ans = slot(res, "tbl") |> dplyr::filter(CHR == local(anac(seqnames(ok)[1])), BP >= local(IRanges::start(ok)-radius), BP<= local(IRanges::end(ok)+radius))
  new("ABRIGresource", space = slot(res, "space"), tbl=ans)
}
