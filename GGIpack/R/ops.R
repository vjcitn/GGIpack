#' define constructor for ABRIGresource 
#' @import dplyr
#' @param con is a DBI connection (typically duckdb)
#' @param tissue character(1)
#' @param space character(1) e.g., "hg19"
#' @param pfiles list of absoulte paths to the data for each tissue. 
#' @examples
#' con = DBI::dbConnect(duckdb::duckdb())
#' ll = ABRIGresource( con, "BAL" , pfiles= ABRIGparquet_paths())
#' print(ll)
#' DBI::dbDisconnect(con)
#' @return ABRIGresource instance for the identified tissue
#' @export
ABRIGresource = function(con, tissue, space="hg19", pfiles) {
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
#' @examples
#' con = DBI::dbConnect(duckdb::duckdb())
#' ll = ABRIGresource( con, "BAL" , pfiles= ABRIGparquet_paths())
#' utils::data("gloc_hg19", package = "GGIpack")
#' kk <- filterByRange(ll, gloc_hg19, "DSP", ggr_field="gene_name")
#' print(kk)
#' DBI::dbDisconnect(con)
#' @return filterByRange instance for the identifies gene_name or gene_id 
#' @export
filterByRange = function(res, ggr, tag, radius=1e5, ggr_field="gene_name") {
  stopifnot(inherits(res, "ABRIGresource"))
  ok = ggr[which(GenomicRanges::mcols(ggr)[[ggr_field]] == tag)]
  stopifnot(width(ok)>0)
  anac = function(x) as.numeric(as.character(x)) # for Rle
  ans = methods::slot(res, "tbl") |> dplyr::filter(CHR == local(anac(GenomicRanges::seqnames(ok)[1])), BP >= local(IRanges::start(ok)-radius), BP<= local(IRanges::end(ok)+radius))
  new("ABRIGresource", space = methods::slot(res, "space"), tbl=ans)
}

#' Takes in a path to a data table and coherses the data into the proper format 
#' @param path The complete path to the DN8 file.
#' @examples
#' path = system.file("extdata/Alveolar_Macrophages_IS.MICA:ILMN_3241692.CAU.meta", package="GGIpack")
#' head(checkData(path = path))
#' @return a data set that can be used to graph in JBrowseR.
#' @export
checkData = function(path){
    df = utils::read.table(path, header = TRUE)
    needCols = c("BP", "P", "CHR")
    stopifnot(needCols %in% colnames(df))
    data =dplyr::mutate(df, start = BP, end =BP)
    return(data)
  }
}

