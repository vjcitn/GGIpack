.libPaths("/root/GGIpack/R/")
utils::data("gloc_hg19", package = "GGIpack")
con = DBI::dbConnect(duckdb::duckdb())
tinyapp2(con, gloc_hg19)
