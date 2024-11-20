# derived from igvShiny demo file by Vince 20 Nov 2024

#' app with a button to add specific GWAS track
#' @examples
#' if (interactive()) {
#'   simpleIGV()
#' }
#' @export
simpleIGV = function() {
# constants
  
  genomes <- c("hg38", "hg19", "mm10", "tair10", "rhos")
  loci <- c("chr5:88,466,402-89,135,305",  "chr1:7,426,231-7,453,241", "MEF2C", "Mef2c",
            "1:7,432,931-7,440,395", "NC_007494.2:370,757-378,078",
            "chr1:6,575,383-8,304,088")
                       
# an example data.frame with specific columns and data types

f <- system.file(package="igvShiny", "extdata", "gwas.RData")
stopifnot(file.exists(f))
tbl.gwas <- get(load(f))

# a stripped-down UI

ui = shinyUI(fluidPage(

  sidebarLayout(
    sidebarPanel(
      helpText("testing"),
      
      br(),
      actionButton("addGwasTrackButton", "Add GWAS Track"),
      
#      div(style="background-color: white; width: 200px; height:30px; padding-left: 5px;
#                   margin-top: 10px; border: 1px solid blue;",
#          htmlOutput("chromLocDisplay")),
#      hr(),
      width=2
    ),
    mainPanel(
      igvShinyOutput('igvShiny_0'),
      width=10
    )
  ) # sidebarLayout
))

# a stripped down server

server = function(input, output, session) {
                       

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
runApp(list(ui=ui, server=server))

}
