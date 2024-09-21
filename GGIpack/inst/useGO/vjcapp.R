library(GGIpack)
library(shiny)
data("gloc_hg19", package="GGIpack")

# set up data resources
con = DBI::dbConnect(duckdb::duckdb())
lungres = GTExresource(con, pfile=file.path(Sys.getenv("GGI_PARQUET_FOLDER"),
    "lungpl05.parquet"))
whblres = GTExresource(con, pfile=file.path(Sys.getenv("GGI_PARQUET_FOLDER"),
    "wholeblpl05.parquet"))
resl = list(lung=lungres, wholebl=whblres)

# simple UI based on selections

ui = fluidPage(
 sidebarLayout(
  sidebarPanel(
   helpText("using gtex eqtl data"),
   checkboxGroupInput("respicks", "resources",
        choices=names(resl), selected=names(resl)[1])
   ),
  mainPanel(
   uiOutput("all")
    )
   )
  )

server = function(input, output) {
# prepare output components
# for now very simple processing of tables, later, perform filtering
# based on gene selection, must be reactive (would not use 'head()' but 
# a filter based on locations
  z = lapply(names(resl), function(x) {
  print(x)
  output[[x]] = 
     DT::renderDataTable(resl[[x]]@tbl |> head() |> as.data.frame() ) # |> DT::datatable())
  })
# communicate selected components to UI
  output$all = renderUI({
   o = lapply(input$respicks, function(x) tabPanel(x, DT::dataTableOutput(x))) #output[[x]])
   do.call(tabsetPanel, o)
  })
 }

runApp(list(ui=ui, server=server))
