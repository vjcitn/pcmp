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
   helpText("pcmp is a prototype of using crosstalk to survey different dimension reduction procedures for SingleCellExperiment data.  Different methods are used in the left and right columns, and different projection components can are used in the top and bottom rows, as selected using the method/top/bot controls below.  See the 'about' tab for more information."),
   selectInput("pickedStrat", "stratby", discreteColdVars(sce), "donor_id"),
   selectInput("meth1", "method 1", c("PC", "UM", "TS"), "PC"),
   selectInput("meth2", "method 2", c("PC", "UM", "TS"), "TS"),
   numericInput("topx", "top x", 1, min=1, max=ncomp-1, step=1),
   numericInput("topy", "top y", 2, min=2, max=ncomp, step=1),
   numericInput("botx", "bot x", 2, min=1, max=ncomp-1, step=1),
   numericInput("boty", "bot y", 3, min=2, max=ncomp, step=1),
   actionButton("btnSend", "Stop app")
   ),
 mainPanel(
  tabsetPanel(
  tabPanel("scatter",
   fluidRow(
     column(6, d3scatterOutput("scatter1")),
     column(6, d3scatterOutput("scatter2"))
     ),
    fluidRow(
     column(6, d3scatterOutput("scatter3")),
     column(6, d3scatterOutput("scatter4"))
     )
   ), # end panel
  tabPanel("about",
   helpText("pcmp is a prototype of using crosstalk to survey different dimension reduction procedures for SingleCellExperiment data.  Different methods are used in the left and right columns, and different projection components can are used in the top and bottom rows, as selected using the method/top/bot controls below.  The ColorBy button will recolor points according to discrete covariates in the colData of the input object."),
   helpText("The initial example uses the Allen Brain Atlas RNA-seq data on anterior cingulate cortex (ACC) and primary visual cortex (VIS) brain regions.  Strata were formed using donor (3 levels) and region (2 levels) and 300 cells were sampled at random in each stratum.")
 )
)
)
) 
