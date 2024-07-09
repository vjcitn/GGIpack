#' demo app
#' @import shiny
#' @import dplyr
#' @import duckdb
#' @import EnsDb.Hsapiens.v75
#' @import ensembldb
#' @param con a DBI connection
#' @param genelocs a GRanges instance with gene addresses
#' @note Very specialized, just has a few genes, uses specific
#' field from genelocs argument.
#' @examples
#' gloc_hg19 = ensembldb::genes(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75)
#' con = DBI::dbConnect(duckdb::duckdb())
#' tinyapp(con, gloc_hg19)
#' @export
tinyapp = function(con, genelocs = ensembldb::genes(EnsDb.Hsapiens.v75) ) {
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

#' demo app 2
#' @import shiny
#' @param con a DBI connection
#' @param genelocs a GRanges instance with gene addresses
#' @note Very specialized, just has a few genes, uses specific
#' field from genelocs argument.
#' @examples
#' gloc_hg19 = ensembldb::genes(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75)
#' con = DBI::dbConnect(duckdb::duckdb())
#' tinyapp2(con, gloc_hg19)
#' @export
tinyapp2 = function(con, genelocs) {
 pfiles <<- ABRIGparquet_paths()
 ui = fluidPage(
  sidebarLayout(
   sidebarPanel(
    helpText("GGI demo"),
    selectInput("gene", "gene", c("ORMDL3", "DSP", "AGL", "DBT", "SASS6"))
    ), 
   mainPanel(
    #verbatimTextOutput("stuff")
    uiOutput("alltabs")
    )
  )
 )   # also need tabs, about etc.
 
 server = function(input, output) {
  output$stuff = renderPrint({
   mygene = input$gene
   mytiss = input$tiss
   newres = ABRIGresource(con, input$tiss, pfiles = pfiles)
   kk <- filterByRange(newres, genelocs, mygene, ggr_field="gene_name")
   kk@tbl
  })
#
# highly repetitious, use reactive better or build a list, possibly
# parallelized
#
  output$BALstuff = renderPrint({
   mygene = input$gene 
   BALres = ABRIGresource(con, "BAL", pfiles = pfiles)
   kk <- filterByRange(BALres, genelocs, mygene, ggr_field="gene_name")
   kk@tbl
   })
  output$BEBstuff = renderPrint({
   mygene = input$gene 
   BEBres = ABRIGresource(con, "BroncEpiBrush", pfiles = pfiles)
   kk <- filterByRange(BEBres, genelocs, mygene, ggr_field="gene_name")
   kk@tbl
   })
  output$CD4stim = renderPrint({
   mygene = input$gene 
   BEBres = ABRIGresource(con, "CD4Stim", pfiles = pfiles)
   kk <- filterByRange(BEBres, genelocs, mygene, ggr_field="gene_name")
   kk@tbl
   })
   output$CD4Unstim = renderPrint({
    mygene = input$gene 
    BEBres = ABRIGresource(con, "CD4Unstim", pfiles = pfiles)
    kk <- filterByRange(BEBres, genelocs, mygene, ggr_field="gene_name")
    kk@tbl
  })
   output$AlvMacphage = renderPrint({
     mygene = input$gene 
     BEBres = ABRIGresource(con, "AlvMacphage", pfiles = pfiles)
     kk <- filterByRange(BEBres, genelocs, mygene, ggr_field="gene_name")
     kk@tbl
   })
   output$PaxRNA = renderPrint({
     mygene = input$gene 
     BEBres = ABRIGresource(con, "PaxRNA", pfiles = pfiles)
     kk <- filterByRange(BEBres, genelocs, mygene, ggr_field="gene_name")
     kk@tbl
   })
  output$alltabs = renderUI({
   tabsetPanel(
    tabPanel("BAL", helpText("A"), verbatimTextOutput("BALstuff")),
    tabPanel("BronchEpiBrush", verbatimTextOutput("BEBstuff")),
    tabPanel("CD4stim", verbatimTextOutput("CD4stim")),
    tabPanel("CD4Unstim", verbatimTextOutput("CD4Unstim")),
    tabPanel("AlvMacphage", verbatimTextOutput("AlvMacphage")),
    tabPanel("PaxRNA", verbatimTextOutput("PaxRNA"))
    )
   })
 }
 
 runApp(list(ui=ui, server=server))
}
