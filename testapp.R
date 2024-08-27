

testapp = function() {
  patquetTableLoc <-system.file("extdata","parquetDataTable.csv", package = "GGIpack" )
  patquetTable <- utils::read.csv(patquetTableLoc) 
  genomeVersion <- "hg19"
  utils::data("geneNames", package = "GGIpack")
  #-----------------------------------------------------------------------------------
  if(!dir.exists("tracks"))
    dir.create("tracks")
  shiny::addResourcePath("tracks", "tracks")
  #-------------------------------------------------------------------------------------
  
  
  ui = fluidPage(
    sidebarLayout(
      sidebarPanel(
        helpText("GGI demo"),
        selectizeInput("gene", "gene", geneNames) #,
        #actionButton("zoomButton", label = "Zoom")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel(
      DT::DTOutput("parquet")
        ) ,# parquet tabpanel
      tabPanel(igvShiny::igvShinyOutput("igvShiny_0")
               #shinyFeedback::useShinyFeedback()
      )#igvshiny tabset
      ) #tabsetpanel
    ) #mainpanel
    ) #sidebarlayout 
  )#fluidpage
    
    server = function(input, output, session) {  
      
  updateSelectizeInput(session, "gene", choices =sort(geneNames),server = TRUE)
  output$parquet<- DT::renderDT({
    data.frame(patquetTable)
  })#parquet
  
  output$igvShiny_0 = igvShiny::renderIgvShiny({
    genomeOptions <- igvShiny::parseAndValidateGenomeSpec(genomeName = genomeVersion, initialLocus = "all")
    igvShiny(genomeOptions)
  })#igvShiny_0
  
    }# server
    
    runApp(list(ui=ui, server=server))
}









