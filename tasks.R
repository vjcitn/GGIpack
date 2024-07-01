library(GGIpack)

library(EnsDb.Hsapiens.v75)

pfiles = ABRIGparquet_paths()

if (!exists("gloc_hg19")) gloc_hg19 = genes(EnsDb.Hsapiens.v75)

con = DBI::dbConnect(duckdb::duckdb())
ll = ABRIGresource( con, "BAL" )
print(ll)

kk <- filterByRange(ll, gloc_hg19, "DSP", ggr_field="gene_name")
kk2 <- filterByRange(ll, gloc_hg19, "ENSG00000096696", ggr_field="gene_id")
kk@tbl
slot(kk, "tbl") |> dplyr::group_by(gene,CHR,BP,SNP) |> dplyr::summarise(n=n())
slot(kk, "space")
kk2

DBI::dbDisconnect(con)
