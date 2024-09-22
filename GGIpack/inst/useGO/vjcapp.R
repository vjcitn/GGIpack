library(GGIpack)
library(shiny)
data("gloc_hg19", package="GGIpack")

# set up data resources
con = DBI::dbConnect(duckdb::duckdb())
lungres = GTExresource(con, tisstag="lung", pfile=file.path(Sys.getenv("GGI_PARQUET_FOLDER"),
    "lungpl05.parquet"))
whblres = GTExresource(con, tisstag="wholebl", pfile=file.path(Sys.getenv("GGI_PARQUET_FOLDER"),
    "wholeblpl05.parquet"))
resl = list(lung=lungres, wholebl=whblres)
fvec = names(resl[[1]]@tbl |> head(2) |> as.data.frame())

# simple UI based on selections

ui = fluidPage(
 sidebarLayout(
  sidebarPanel(
   helpText("using gtex eqtl data"),
   checkboxGroupInput("respicks", "resources",
        choices=names(resl), selected=names(resl)[1]),
   numericInput("nrecs", "nrecs", min=5, max=1000, value=10), 
   radioButtons("chr", "chr", choices=1:22, selected=1, inline=TRUE),
   actionButton("stop", "stop app"),
   width=2
   ),
  mainPanel(
   uiOutput("all")
    )
   )
  )

server = function(input, output) {
  observeEvent(input$stop, {
    DBI::dbDisconnect(con)
    stopApp()
    })
# prepare output components
# for now very simple processing of tables, later, perform filtering
# based on gene selection, must be reactive (would not use 'head()' but 
# a filter based on locations
  z = lapply(names(resl), function(x) {
  print(x)
  output[[x]] = 
     DT::renderDataTable(resl[[x]]@tbl |> dplyr::filter(seqnames == as.character(local(input$chr))) |>
                head(input$nrecs) |> as.data.frame() ) # |> DT::datatable())
  })
  output$theplot = renderPlot({ plot(1,1) })
# communicate selected components to UI
  output$all = renderUI({
   o = lapply(c(input$respicks, "viz"), function(x) {
            if (x != "viz") tabPanel(x, DT::dataTableOutput(x))
            else tabPanel("viz", shiny::plotOutput("theplot"))
            })
   do.call(tabsetPanel, o)
  })
 }

runApp(list(ui=ui, server=server))
