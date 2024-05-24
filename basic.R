library(methods)
library(duckdb)
library(dplyr)

setClass("ggiResource", slots=c(start="numeric",
 end="numeric", seqnames="ANY", space="ANY",
 score="numeric", tbl="ANY"))

setMethod("show", "ggiResource", function(object) {
 callNextMethod()
})

setClass("ABRIGresource", contains="ggiResource")

setMethod("show", "ABRIGresource", function(object) {
 cat("CDNM ggiResource\n")
 print(slot(object, "tbl"))
})

library(EnsDb.Hsapiens.v75)

pfiles =
c("/udd/remcr/abrig/ca_ba.parquet", "/udd/remcr/abrig/ca_brep.parquet",
"/udd/remcr/abrig/cc_s.parquet", "/udd/remcr/abrig/cc_unstim.parquet",
"/udd/remcr/abrig/ch_am.parquet", "/udd/remcr/abrig/comb.parquet")

names(pfiles) = c("BAL", "BroncEpiBrush", "CD4Stim", "CD4Unstim",
                  "AlvMacphage", "PaxRNA")


if (!exists("gloc_hg19")) gloc_hg19 = genes(EnsDb.Hsapiens.v75)

ABRIGresource = function(con, tissue) {
   ttypes = c("BAL", "BroncEpiBrush", "CD4Stim", "CD4Unstim",
                  "AlvMacphage", "PaxRNA")
   stopifnot(tissue %in% ttypes)
   ans = dplyr::tbl(con, pfiles[tissue])
   new("ABRIGresource", tbl=ans)
}

filterByHGNC = function(res, ggr, sym, radius=1e5, field="gene_name") {
  stopifnot(inherits(res, "ABRIGresource"))
  ok = ggr[which(mcols(ggr)[[field]] == sym)]
  stopifnot(width(ok)>0)
  anac = function(x) as.numeric(as.character(x)) # for Rle
  ans = slot(res, "tbl") |> dplyr::filter(CHR == local(anac(seqnames(ok)[1])), BP >= IRanges::start(ok)-radius, BP<= IRanges::end(ok)+radius)
  new("ABRIGresource", tbl=ans)
}

kk <- filterByHGNC(ll, gloc_hg19, "DSP")
slot(kk, "tbl") |> group_by(gene,CHR,BP,SNP) |> summarise(n=n())

