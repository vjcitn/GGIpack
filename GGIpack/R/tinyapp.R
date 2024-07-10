
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
  
  allrefs = reactive({
   ttypes = c("BAL", "BroncEpiBrush", "CD4Stim", "CD4Unstim",
                  "AlvMacphage", "PaxRNA")
   mygene = input$gene
   allres = lapply(ttypes, function(x) ABRIGresource(con, x, pfiles=pfiles))
   allfilt = lapply(allres, function(x) filterByRange(x,
         genelocs, mygene, ggr_field="gene_name"))
   names(allfilt) = ttypes
   allfilt
   })
  dorounds = function(mydf) {
   mydf$P = round(mydf$P, 3)
   mydf$SE = round(mydf$SE, 3)
   mydf$MAF = round(mydf$MAF, 3)
   mydf$BETA = round(mydf$BETA, 3)
   mydf |> dplyr::select(-score, -seqnames)
   }
  output$BALstuff = DT::renderDataTable({
   refs = allrefs()
   refs[["BAL"]]@tbl |> dplyr::arrange(FDR) |> as.data.frame() |> dorounds()
   })
  output$BEBstuff = DT::renderDataTable({
   refs = allrefs()
   refs[["BroncEpiBrush"]]@tbl |> dplyr::arrange(FDR) |> as.data.frame() |> dorounds()
   })
  output$CD4stim = DT::renderDataTable({
   refs = allrefs()
   refs[["CD4Stim"]]@tbl |> dplyr::arrange(FDR) |> as.data.frame() |> dorounds()
   })
   output$CD4Unstim = DT::renderDataTable({
   refs = allrefs()
   refs[["CD4Unstim"]]@tbl |> dplyr::arrange(FDR) |> as.data.frame() |> dorounds()
   })
   output$AlvMacphage = DT::renderDataTable({
   refs = allrefs()
   refs[["AlvMacphage"]]@tbl |> dplyr::arrange(FDR) |> as.data.frame() |> dorounds()
   })
   output$PaxRNA = DT::renderDataTable({
   refs = allrefs()
   refs[["PaxRNA"]]@tbl |> dplyr::arrange(FDR) |> as.data.frame() |> dorounds()
   })
  output$alltabs = renderUI({
   tabsetPanel(
    tabPanel("BAL", helpText("A"), DT::dataTableOutput("BALstuff")),
    tabPanel("BronchEpiBrush", DT::dataTableOutput("BEBstuff")),
    tabPanel("CD4stim", DT::dataTableOutput("CD4stim")),
    tabPanel("CD4Unstim", DT::dataTableOutput("CD4Unstim")),
    tabPanel("AlvMacphage", DT::dataTableOutput("AlvMacphage")),
    tabPanel("PaxRNA", DT::dataTableOutput("PaxRNA"))
    )
   })
 }
 
 runApp(list(ui=ui, server=server))
}
