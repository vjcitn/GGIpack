library(methods)
library(duckdb)
library(dplyr)

setClass("ggiResource", slots=c(start="numeric",
 end="numeric", seqnames="ANY", space="ANY",
 score="numeric"))

setClass("ABRIGresource", contains="ggiResource")

library(EnsDb.Hsapiens.v75)

pfiles =
c("/udd/remcr/abrig/ca_ba.parquet", "/udd/remcr/abrig/ca_brep.parquet",
"/udd/remcr/abrig/cc_s.parquet", "/udd/remcr/abrig/cc_unstim.parquet",
"/udd/remcr/abrig/ch_am.parquet", "/udd/remcr/abrig/comb.parquet")

names(pfiles) = c("BAL", "BroncEpiBrush", "CD4Stim", "CD4Unstim",
                  "AlvMacphage", "PaxRNA")


gloc_hg19 = genes(EnsDb.Hsapiens.v75)

ABRIGresource = function(con, tissue) {
   ttypes = c("BAL", "BroncEpiBrush", "CD4Stim", "CD4Unstim",
                  "AlvMacphage", "PaxRNA")
   stopifnot(tissue %in% ttypes)
   dplyr::tbl(con, pfiles[tissue])
}
