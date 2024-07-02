#' demo app
#' @import shiny
#' @import dplyr
#' @import duckdb
#' @import EnsDb.Hsapiens.v75
#' @param con a DBI connection
#' @param genelocs a GRanges instance with gene addresses
#' @export
tinyapp = function(con, genelocs = genes(EnsDb.Hsapiens.v75) ) {
 pfiles <<- ABRIGparquet_paths()
 ui = fluidPage(
  sidebarLayout(
   sidebarPanel(
    helpText("GGI demo"),
    radioButtons("tiss", "tissue", names(pfiles)),
    selectInput("gene", "gene", c("ORMDL3", "DSP", "AGL", "DBT", "SASS6"))
    ), 
   mainPanel(
    verbatimTextOutput("stuff")
    )
  )
 )   # also need tabs, about etc.
 
 server = function(input, output) {
  output$stuff = renderPrint({
   mygene = input$gene
   mytiss = input$tiss
   newres = ABRIGresource(con, input$tiss, pfiles=pfiles)
   kk <- filterByRange(newres, genelocs, mygene, ggr_field="gene_name")
   kk@tbl
  })
 }
 
 runApp(list(ui=ui, server=server))
}
