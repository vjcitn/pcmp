  library(SingleCellExperiment)
  library(d3scatter)
  library(crosstalk)
  sce = pcmp::sce300xx
  rd = reducedDims(sce)
  nrd = names(rd)
  ncomps = vapply(rd, ncol, numeric(1))
  stopifnot(all(ncomps == ncomps[1]))
  ncomp <- ncomps[1]

ui <- fluidPage(
 sidebarPanel(width=2,
   helpText("abc"),
   selectInput("meth1", "method1", c("PC", "UM", "TS"), "PC"),
   selectInput("meth2", "method1", c("PC", "UM", "TS"), "TS"),
   numericInput("topx", "top x", 1, min=1, max=ncomp-1, step=1),
   numericInput("topy", "top y", 2, min=2, max=ncomp, step=1),
   numericInput("botx", "bot x", 1, min=1, max=ncomp-1, step=1),
   numericInput("boty", "bot y", 2, min=2, max=ncomp, step=1),
   actionButton("btnSend", "Stop app")
   ),
 mainPanel(
  tabsetPanel(
  tabPanel("by donor",
   fluidRow(
     column(6, d3scatterOutput("scatter1")),
     column(6, d3scatterOutput("scatter2"))
     ),
    fluidRow(
     column(6, d3scatterOutput("scatter3")),
     column(6, d3scatterOutput("scatter4"))
     )
   ), # end panel
  tabPanel("by region",
   fluidRow(
     column(6, d3scatterOutput("scatter5")),
     column(6, d3scatterOutput("scatter6"))
     ),
    fluidRow(
     column(6, d3scatterOutput("scatter7")),
     column(6, d3scatterOutput("scatter8"))
     )
   ), # end panel
  tabPanel("by subregion",
   fluidRow(
     column(6, d3scatterOutput("scatter9")),
     column(6, d3scatterOutput("scatter10"))
     ),
    fluidRow(
     column(6, d3scatterOutput("scatter11")),
     column(6, d3scatterOutput("scatter12"))
     )
   ) # end panel
  )
 )
 )
 
