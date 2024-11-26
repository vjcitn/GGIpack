# derived from igvShiny demo file by Vince 20 Nov 2024
#' app with a button to add specific GWAS track'
#' @rawNamespace import(shiny, except=c(dataTableOutput, renderDataTable))
#' @rawNamespace import(GenomicRanges, except=c(intersect, union, setdiff))
#' @import igvShiny
#' @import DT
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
                       
# add some code that can get and print a prespecified subtable of
# lung gtex

fo = ggi_gtex_cache("lungpl05.parquet")
lungpa = fo
con = DBI::dbConnect(duckdb::duckdb())
lu = GTExresource(con, tisstag="lung", pfile=lungpa)
mydat = lu@tbl |> head(50) |> as.data.frame()

#
# rename the columns of mydat appropriately so that they agree with tbl.gwas
#

# then assign this new table to "tbl.gwas"


  observeEvent(input$addGwasTrackButton, {
    #this sprintf function do not work but the code does make it here. 
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
