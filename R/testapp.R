#' testapp app 
#' @rawNamespace import(shiny, except=c(dataTableOutput, renderDataTable))
#' @rawNamespace import(GenomicRanges, except=c(intersect, union, setdiff))
#' @import igvShiny
#' @import DT
#' @examples
#' if (interactive()){
#' testapp()
#' }
#' @export
testapp = function() {
  data("gloc_hg19", package="GGIpack")
  utils::data("ensg", package = "GGIpack")
  ensg = ensg[order(names(ensg))]
  
  Sys.setenv("GGI_PARQUET_FOLDER"= tempdir())
  
  genomeVersion = "hg19"
  
  # set up data resources
  con = DBI::dbConnect(duckdb::duckdb())
  lungpa = try(ggi_gtex_cache("lungpl05.parquet"))
  if (inherits(lungpa, "try-error") | nchar(lungpa)==0)
    lungpa = file.path(Sys.getenv("GGI_PARQUET_FOLDER"), "lungpl05.parquet")
  lungres = GTExresource(con, tisstag="lung", pfile=lungpa)
  wbpa = try(ggi_gtex_cache("wholeblpl05.parquet"))
  if (inherits(wbpa, "try-error") | nchar(wbpa)==0)
    wbpa = file.path(Sys.getenv("GGI_PARQUET_FOLDER"), "wholeblpl05.parquet")
  whblres = GTExresource(con, tisstag="wholebl", pfile=wbpa)
  resl = list(lung=lungres, wholebl=whblres)
  fvec = names(resl[[1]]@tbl |> head(2) |> as.data.frame())
  
  #setup needed for the igvShiny table.
  #-----------------------------------------------------------------------------------
  if(!dir.exists("tracks"))
    dir.create("tracks")
  shiny::addResourcePath("tracks", "tracks")
  #-------------------------------------------------------------------------------------
  
  # simple UI based on selections
  
  ui = fluidPage(
    sidebarLayout(
      sidebarPanel(
        helpText("using gtex eqtl data"),
 actionButton("addGwasTrackButton", "Add GWAS Track"),
        checkboxGroupInput("respicks", "resources",
                           choices=names(resl), selected=names(resl)[1]),
        numericInput("nrecs", "nrecs", min=5, max=100000, value=10000), 
        radioButtons("focus", "focus", choices=c("chr", "gene", "rsid")),
        conditionalPanel(
          condition = "input.focus == 'chr'",
          radioButtons("chr", "chr", choices=1:22, selected=1, inline=TRUE)
        ),
        conditionalPanel(
          condition = "input.focus == 'gene'",
          selectInput("gene", "gene", choices=NULL)
        ),
        conditionalPanel(
          condition = "input.focus == 'rsid'",
          textInput("snp", "snp")
        ),
        actionButton("graphGWAS", "Graph GWAS"), 
        actionButton("stop", "stop app"),
        width=2
      ),
      mainPanel(
        tabPanel("igvShiny",igvShiny::igvShinyOutput("igvShiny_0")),
        uiOutput("all")
      )
    )
  )
  
  server = function(input, output, session) {

f <- system.file(package="igvShiny", "extdata", "gwas.RData")
stopifnot(file.exists(f))
tbl.gwas <- get(load(f))

    updateSelectizeInput(session, 'gene', choices = ensg, server = TRUE)
    observeEvent(input$stop, {
      DBI::dbDisconnect(con)
      stopApp()
    })#updateSelectizeInput
    

    
server2 = function(input, output, session) {
           

  observeEvent(input$addGwasTrackButton, {
    sprintf("---- addGWASTrack")
    sprintf("current working directory: %s", getwd())
    showGenomicRegion(session, id="igvShiny_0", "chr19:45,248,108-45,564,645")
    loadGwasTrack(session, id="igvShiny_0", trackName="gwas", tbl=tbl.gwas, deleteTracksOfSameName=FALSE)
  }) 
 output$igvShiny_0 <- renderIgvShiny({
    cat("--- starting renderIgvShiny\n");
    genomeOptions <- parseAndValidateGenomeSpec(genomeName="hg38",  initialLocus=loci[7])
    x <- igvShiny(genomeOptions,
                  displayMode="SQUISHED",
                  tracks=list()
    )
    cat("--- ending renderIgvShiny\n");
    return(x)
  })
    
}  
    
  observeEvent(input$addGwasTrackButton, {
    sprintf("---- addGWASTrack")
    sprintf("current working directory: %s", getwd())
    showGenomicRegion(session, id="igvShiny_0", "chr19:45,248,108-45,564,645")
    loadGwasTrack(session, id="igvShiny_0", trackName="gwas", tbl=tbl.gwas, deleteTracksOfSameName=FALSE)
  }) 
    
    ##igvshiny setup the baseline
    output$igvShiny_0 = igvShiny::renderIgvShiny({
      genomeOptions <- igvShiny::parseAndValidateGenomeSpec(genomeName = genomeVersion, initialLocus = "all")
      igvShiny(genomeOptions)
    })#igvShiny_0
    
    
    observeEvent(input$chr ,{
      print("I made it to observeEvent")
      print("input chr")
      print(input$chr )
      
      
      
      dat = resl[["lung"]]@tbl |> dplyr::filter(seqnames == as.character(local(6))) |>
        head(10000) |> as.data.frame() 
      gwasTrack = makeGWASTrack( dat = dat)
      print("made it to input$focus == chr")
      display(gwasTrack, session, id = "igvShiny_0")
      
      
    })#observeEvent
    
    
    
    z = lapply(names(resl), function(x) {
      output[[x]] = 
        DT::renderDataTable({
          if (input$focus == "chr")
            resl[[x]]@tbl |> dplyr::filter(seqnames == as.character(local(input$chr))) |>
            dplyr::arrange(score) |>
            head(input$nrecs) |> as.data.frame() |> dorounds()
          else if (input$focus == "gene") 
            resl[[x]]@tbl |> dplyr::filter(molecular_trait_id == as.character(local(input$gene))) |>
            as.data.frame() |> dorounds()
          else if (input$focus == "rsid") 
            dat = resl[[x]]@tbl |> dplyr::filter(rsid == as.character(local(input$snp))) |>
              as.data.frame() |> dorounds()
        })#renderDataTable
    }#function
    )#lapply
    
    
    
    # communicate selected components to UI
    output$all = renderUI({
      o = lapply(c(input$respicks, "viz"), function(x) {
        if (x != "viz") {
          tabPanel(x, DT::dataTableOutput(x))
          #tabPanel("igvShiny",igvShiny::igvShinyOutput("igvShiny_0"))
        }
        
        
      })
      do.call(tabsetPanel, o)
    }) #renderUI
  }
  
  runApp(list(ui=ui, server=server))
}

