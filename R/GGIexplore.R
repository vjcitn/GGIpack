#' GGIexplore app 
#' @rawNamespace import(shiny, except=c(dataTableOutput, renderDataTable))
#' @rawNamespace import(GenomicRanges, except=c(intersect, union, setdiff))
#' @import igvShiny
#' @import DT
#' @examples
#' if (interactive()){
#' GGIexplore()
#' }
#' @export
GGIexplore = function() {
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
      actionButton("stop", "stop app"),
      width=2
    ),
    mainPanel(
      uiOutput("all")
    )
  )
)

server = function(input, output, session) {
  updateSelectizeInput(session, 'gene', choices = ensg, server = TRUE)
  observeEvent(input$stop, {
    DBI::dbDisconnect(con)
    stopApp()
  })#updateSelectizeInput
  
  
  
  
  ##igvshiny setup the baseline
  output$igvShiny_0 = igvShiny::renderIgvShiny({
    genomeOptions <- igvShiny::parseAndValidateGenomeSpec(genomeName = genomeVersion, initialLocus = "all")
    igvShiny(genomeOptions)
  })#igvShiny_0
  
  
  observeEvent(list(input$respicks, input$focus, input$nrecs) ,{
    print("I made it to observeEvent")
    print("input respick")
    print(input$respicks )
    print("input focus")
    print(input$focus )
    print("input nrecs")
    print(input$nrecs )
    par(mfrow=c(length(input$respicks), 1))
    for (x in input$respicks) {
      if (input$focus == "chr"){
        dat = resl[[x]]@tbl |> dplyr::filter(seqnames == as.character(local(input$chr))) |>
          head(input$nrecs) |> as.data.frame() 
        gwasTrack = makeGWASTrack( dat = dat)
        #print(head(GWASTrack))
        display(gwasTrack, session, id = "igvShiny_0")
        
      }
      else if (input$focus == "gene") {
        dat = resl[[x]]@tbl |> dplyr::filter(molecular_trait_id == as.character(local(input$gene))) |>
          as.data.frame() 
        gwasTrack = makeGWASTrack( dat = dat)
        #print(head(GWASTrack))
        display(gwasTrack, session, id = "igvShiny_0")
        
      }
      else if (input$focus == "rsid") {
        dat = resl[[x]]@tbl |> dplyr::filter(rsid == as.character(local(input$snp))) |>
          as.data.frame() 
        gwasTrack = makeGWASTrack( dat = dat)
        #print(head(GWASTrack))
        display(gwasTrack, session, id = "igvShiny_0")
        
      }
      #gwasTrack = makeGWASTrack( dat = dat)
      #print(head(GWASTrack))
      #igvShiny::display(gwasTrack, session, id = "igvShiny_0")
    }#for loop
    print("made it to the end")
    #I need to make a name for the track 
    gwasTrack = makeGWASTrack( dat = dat)
    #print(head(GWASTrack))
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
  
  
  output$theplot = renderPlot({ 
    par(mfrow=c(length(input$respicks), 1))
    for (x in input$respicks) {
      if (input$focus == "chr")
        dat = resl[[x]]@tbl |> dplyr::filter(seqnames == as.character(local(input$chr))) |>
          head(input$nrecs) |> as.data.frame() 
      else if (input$focus == "gene") 
        dat = resl[[x]]@tbl |> dplyr::filter(molecular_trait_id == as.character(local(input$gene))) |>
          as.data.frame() 
      else if (input$focus == "rsid") 
        dat = resl[[x]]@tbl |> dplyr::filter(rsid == as.character(local(input$snp))) |>
          as.data.frame() 
      plot(dat$start, -log10(dat$score), main=x)
    }
  }, height=800L)
  
  
  # communicate selected components to UI
  output$all = renderUI({
    o = lapply(c(input$respicks, "viz"), function(x) {
      if (x != "viz") {
        tabPanel(x, DT::dataTableOutput(x))
        #tabPanel("igvShiny",igvShiny::igvShinyOutput("igvShiny_0"))
        }
      #else tabPanel("viz", shiny::plotOutput("theplot"))
      else {
        tabPanel("igvShiny",igvShiny::igvShinyOutput("igvShiny_0"))
        #igvShiny::display(input$gwasgraph, session, id = "igvShiny_0")
      }
    })
    do.call(tabsetPanel, o)
  }) #renderUI
}

runApp(list(ui=ui, server=server))
}

