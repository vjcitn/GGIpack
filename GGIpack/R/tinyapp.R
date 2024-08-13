#' demo app 2
#' @rawNamespace import(shiny, except=c(dataTableOutput, renderDataTable))
#' @rawNamespace import(GenomicRanges, except=c(intersect, union, setdiff))
#' @import DT
#' @param con a DBI connection
#' @param genelocs a GRanges instance with gene addresses
#' @note Very specialized, just has a few genes, uses specific
#' field from genelocs argument.
#' @examples
#' if (interactive()){
#' utils::data("gloc_hg19", package = "GGIpack")
#' con = DBI::dbConnect(duckdb::duckdb())
#' tinyapp2(con, gloc_hg19)
#' }
#' @export
tinyapp2 = function(con, genelocs) {
 pfiles <<- ABRIGparquet_paths()
 utils::data("geneNames", package = "GGIpack")
 ui = fluidPage(
  sidebarLayout(
   sidebarPanel(
    helpText("GGI demo"),
    selectizeInput("gene", "gene", geneNames)
    ), 
   mainPanel(
    uiOutput("alltabs")
    ) #mainPanel
  ) #sidebarLayout
 )   # also need tabs, about etc.
 
 server = function(input, output, session) {
  updateSelectizeInput(session, "gene", choices =sort(geneNames),server = TRUE)
   
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
    req(input$gene)
   ttypes = c("BAL", "BroncEpiBrush", "CD4Stim", "CD4Unstim",
                  "AlvMacphage", "PaxRNA")
   mygene = input$gene
   allres = lapply(ttypes, function(x) ABRIGresource(con, x, pfiles=pfiles))
   allfilt = lapply(allres, function(x) filterByRange(x,
         genelocs, mygene, ggr_field="gene_name"))
   names(allfilt) = ttypes
   allfilt
   })#allrefs
  

  output$aboutPage= renderPrint({
    br(),
          p(#sprintf(
            "GGIpack tinyapp2 version %s.  This app uses parquet files made 
                 from the abrig  data release on 05/15/2023. This new data release merges the 
                 population data for each of the cell types together thus there is no more choice 
                 for population since there is no way to separate the data anymore. The original
                 data can be found in the following  path on the Nantucket server."#,
            #packageVersion("GGIpack") ) 
          ),
          br(),
          p("/proj/regeps/regep00/studies/ABRIG/analyses/reahs/cis_eqtl_matrixEqtl.Release.15.05.23/"),
          br(),
          p("The parquet files were made using the package duckdb . In short these files
                   were made by merging all the 23 chromosome files for each of the cell types 
                   into one file. The table below shows the name of the cell type, the name 
                   that it is called in the app, and the name of the actual file."),
          br(),
          #dataTableOutput("parquet"),
          p("The files can be found in the Nantucket server under."),
          br(),
          p("/udd/remcr/abrig/"),
          br(),
          p("The code that was used to create these files can be found in the following 
                   Changit Repository."),
          br(),
          p("https://changit.bwh.harvard.edu/remcr/abrigResources")

})#stuff
   
           
  dorounds = function(mydf) {
   mydf$P = formatC(mydf$P, format = "e", digits= 3)
   mydf$SE = round(mydf$SE, 3)
   mydf$MAF = round(mydf$MAF, 3)
   mydf$BETA = round(mydf$BETA, 3)
   mydf$FDR= round(mydf$FDR, 3)
   mydf$statistic= round(mydf$statistic, 3)
   mydf |> dplyr::select(-score, -seqnames, -SNP)
   }#dorounds
  output$BALstuff = DT::renderDataTable({
   refs = allrefs()
   refs[["BAL"]]@tbl |> dplyr::arrange(FDR) |> as.data.frame() |> dorounds()
   })#BALstuff
  output$BEBstuff = DT::renderDataTable({
   refs = allrefs()
   refs[["BroncEpiBrush"]]@tbl |> dplyr::arrange(FDR) |> as.data.frame() |> dorounds()
   })#BEBstuff
  output$CD4stim = DT::renderDataTable({
   refs = allrefs()
   refs[["CD4Stim"]]@tbl |> dplyr::arrange(FDR) |> as.data.frame() |> dorounds()
   }) #CD4Stim
   output$CD4Unstim = DT::renderDataTable({
   refs = allrefs()
   refs[["CD4Unstim"]]@tbl |> dplyr::arrange(FDR) |> as.data.frame() |> dorounds()
   })#CD4Unstim
   output$AlvMacphage = DT::renderDataTable({
   refs = allrefs()
   refs[["AlvMacphage"]]@tbl |> dplyr::arrange(FDR) |> as.data.frame() |> dorounds()
   })#AlvMacphage
   output$PaxRNA = DT::renderDataTable({
   refs = allrefs()
   refs[["PaxRNA"]]@tbl |> dplyr::arrange(FDR) |> as.data.frame() |> dorounds()
   }) #paxRNA
  output$alltabs = renderUI({
   tabsetPanel(
    tabPanel("BAL",  DT::dataTableOutput("BALstuff")),
    tabPanel("BronchEpiBrush", DT::dataTableOutput("BEBstuff")),
    tabPanel("CD4stim", DT::dataTableOutput("CD4stim")),
    tabPanel("CD4Unstim", DT::dataTableOutput("CD4Unstim")),
    tabPanel("AlvMacphage", DT::dataTableOutput("AlvMacphage")),
    tabPanel("PaxRNA", DT::dataTableOutput("PaxRNA")),
    tabPanel("about", "About", verbatimTextOutput("aboutPage"))
    ) #tabsetPanel
   }) #renderUI
  
 }#end of server
 runApp(list(ui=ui, server=server))
 DBI::dbDisconnect(con)
} #end of tinyapp2
