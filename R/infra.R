pkgVersion = function() as.character(read.dcf(system.file("DESCRIPTION", package="pcmp"))[,"Version"])

#' find discrete variables in colData
#' @import SingleCellExperiment crosstalk d3scatter dplyr shiny
#' @import Rtsne irlba umap methods
#' @importFrom S4Vectors SimpleList
#' @param se instance of SummarizedExperiment; SingleCellExperiment also works
#' @param maxnuv numeric(1) largest number of unique values that may be present in a variable that will be regarded as 'discrete'
#' @return character vector of variable names satisfying discreteness condition
#' @examples
#' discreteColdVars  # need example se
#' @export
discreteColdVars = function(se, maxnuv=25) {
 nuv = sapply(colData(se), function(x) length(unique(x)))
 cdn = names(nuv)
 cdn[ which(nuv>1 & nuv<=maxnuv) ]
}
 

#' create vector of indices constituting a stratified sample using a discrete variable to partition strata
#' @param strata atomic vector
#' @param nperstrat size of sample to be retrieved
#' @note fails if any stratum has size less than 'nperstrat'
#' @return numeric vector of indices
#' @examples
#' mytags = c(1,1,1,2,2,3,1,1,1,2,2,3,3,3,3,2,2,2)
#' mydat = 1:length(mytags)
#' names(mydat) = mytags
#' inds = stratsampInds(mytags, nperstrat=2)
#' mydat[inds]
#' @export
stratsampInds = function(strata, nperstrat=300) {
  levs = unique(strata)
  inds = lapply(levs, function(x) which(strata == x))
  ns = vapply(inds, length, numeric(1))
  stopifnot(all(ns>nperstrat))
  samp = lapply(inds, function(x) sample(x, size=nperstrat))
  sort(unlist(samp))
}

#' schematize a sequence of projection methods
#' @param warmdim numeric(1) defaults to 15 -- defines the dimension
#' reduction using irlba()$u for 'warming up' prior to UMAP or t-SNE runs
#' @return a list with two components, 'projectors'
#' and 'retrievers'.  'projectors' is a named list of
#'  functions with parameters 'x' and 'ncomp',
#' where the intended usage has 'x' inheriting from SingleCellExperiment, 
#' and 'ncomp' numeric defining the number of projection components
#' to be retained.  'retrievers' is a named list of functions
#' with parameter 'x', which extracts the ncomp x nsamp
#' projection for each method.
#' @note The default sequence is irlba::irlba, umap::umap,
#' Rtsne::Rtsne; see the note on `warmdim` parameter.
#' Prior to pcmp version 0.3.0, TSNE and UMAP projections were
#' computed 'in one shot' with high target numbers of dimensions.
#' Beginning in 0.3.0, TSNE and UMAP projections are computed
#' targeting 2 and 3 dimensions separately for TSNE, and 2-4
#' dimensions for UMAP.
#' @examples
#' defaultProjectors
#' @export
defaultProjectors = function(warmdim=15) list(
  projectors=list(
   PPCA=function(x) {
      if (isTRUE(options()$verbose))
         message("starting irlba at ", date())
      irlba::irlba(x, nv=4)
      },
   UMAP2=function(x) {
      if (isTRUE(options()$verbose))
         message("starting umap at ", date())
      basic = irlba::irlba(x, nv=warmdim)$u
      umap::umap(basic, n_components=2)
      },
   UMAP3=function(x) {
      if (isTRUE(options()$verbose))
         message("starting umap at ", date())
      basic = irlba::irlba(x, nv=warmdim)$u
      umap::umap(basic, n_components=3)
      },
   UMAP4=function(x) {
      if (isTRUE(options()$verbose))
         message("starting umap at ", date())
      basic = irlba::irlba(x, nv=warmdim)$u
      umap::umap(basic, n_components=4)
      },
   TSNE2=function(x) {
      if (isTRUE(options()$verbose))
         message("starting Rtsne at ", date())
      basic = irlba::irlba(x, nv=warmdim)$u
      Rtsne::Rtsne(basic, dims=2)
     },
   TSNE3=function(x) {
      if (isTRUE(options()$verbose))
         message("starting Rtsne at ", date())
      basic = irlba::irlba(x, nv=warmdim)$u
      Rtsne::Rtsne(basic, dims=3)
     }
    ),
  retrievers = list(
   PPCA = function(x) x$u,
   TSNE2 = function(x) x$Y,
   TSNE3 = function(x) x$Y,
   UMAP2 = function(x) x$layout,
   UMAP3 = function(x) x$layout,
   UMAP4 = function(x) x$layout
    )
   )

.checkDefaultProjectors = function(dp) {
  stopifnot(all(names(dp) == c("projectors", "retrievers")))
  prn = names(dp$projectors)
  ren = names(dp$retrievers)
  stopifnot(all.equal(prn,ren))
}

#' add a series of projectors to transformed assay data
#' in a SingleCellExperiment
#' @param sce SingleCellExperiment instance
#' @param tx function that will transform assay quantifications,
#' defaults to `function(x) log(x+1)`
#' @param projectors a list of list of functions like that created
#' @param assayind numeric(1) second argument ('i') to assay(), defaults to 1
#' @param \dots not used
#' by `defaultProjectors()`
#' @examples
#' if (interactive()) {
#' sce = sce300xx
#' reducedDims(sce) = NULL
#' ov = options()$verbose
#' options(verbose=TRUE)
#' sce = addProjections(sce, proj=defaultProjectors())
#' sce
#' options(verbose=ov) # restore
#' }
#' @export
addProjections = function(sce, tx=function(x) log(x+1),
   projectors, assayind=1, ...) {
  stopifnot(is(sce, "SingleCellExperiment"))
  matToUse = t(tx(assay(sce, i=assayind)))
  prs = lapply(projectors$projectors, function(p) 
      p(x=matToUse))
  #ans = lapply(seq_len(length(prs)), function(i) 
  #  projectors$retrievers[[i]](prs[[i]])
  #  )
  pn = names(prs)
  ans = lapply(pn, function(x) projectors$retrievers[[x]](prs[[x]]))
  names(ans) = names(projectors$retrievers)
  ini = reducedDims(sce)
  if (length(ini) == 0) reducedDims(sce) = SimpleList(ans)
  else reducedDims(sce) = c(ini, SimpleList(ans))
  sce
}

#' Define a shiny app for comparing projections
#' @importFrom limma lmFit eBayes topTable
#' @importFrom stats model.matrix
#' @param sce SingleCellExperiment instance with reducedDims populated
#' with various candidates for projection from assay data --
#' @param assayind numeric(1) 'i' argument to assay(sce, i), defaults to 1, for extracting expression values
#' @note In sce300xx, names(reducedDims) are 'PC', 'UM', 'TS'
#' for projections generated by
#' irlba::irlba, umap::umap, and Rtsne::Rtsne respectively.  Function
#' will delete .pcmpTab, .pcmpSelNum, .pcmpSelCells from .GlobalEnv
#' if found, so that it can assemble information on selections.
#' Note that tabulation of DE genes can be very slow if
#' the input SingleCellExperiment has more than 10k rows.
#' @return will return an instance of PcmpSels
#' @examples
#' if (interactive()) {
#'  message("be sure to use example(pcmpApp, ask=FALSE)")
#'  lit = pcmp::sce300xx
#'  sds = rowSds(assay(lit))
#'  top1k = order(sds,decreasing=TRUE)[1:1000]
#'  lit = lit[top1k,]
#'  chk = pcmpApp(lit) # fires up browser, if ask != FALSE you must
#'  # assent to each update 
#'  try(chk)  # may error if no selections
#'  try(head(geneTable(chk),3))  # likewise
#' }
#' @export
pcmpApp = function(sce) {
 if (nrow(sce)>5000) message("note that performance is greatly enhanced by filtering the feature set down to 5k or so.")
 stores = c(".pcmpTab", ".pcmpSelNum", ".pcmpSelCells")
 sapply(stores, function(x) if(x %in% ls(.GlobalEnv, all.names=TRUE))
    rm(list=x, envir=.GlobalEnv))

  sce$.cellid = colnames(sce)
  rd = reducedDims(sce)
  nrd = names(rd)
  ncomps = vapply(rd, ncol, numeric(1))
  names(ncomps) = nrd # projection-specific numbers of components
  ncomp = 3 # MUST FIX!
#  stopifnot(all(ncomps == ncomps[1]))
#  ncomp <- ncomps[1]
  discv = discreteColdVars(sce)

# uiMaker = function(sce) {
#
#  rd = reducedDims(sce)
#  nrd = names(rd)
#  ncomps = vapply(rd, ncol, numeric(1))
#  stopifnot(all(ncomps == ncomps[1]))
#  ncomp <- ncomps[1]
#  discv = discreteColdVars(sce)
#
pkgVersion = function() as.character(read.dcf(system.file("DESCRIPTION", package="pcmp"))[,"Version"])
UI =   fluidPage(
   sidebarPanel(width=3,
    fluidRow(
     column(12,
       helpText(h3(sprintf("pcmp %s: crosstalk-based interactive graphics \
for dimension reduction in single-cell transcriptomics. \ 
See the 'about' tab for more information.", pkgVersion())))
       )
      ),
    fluidRow( 
     column(12,
       selectInput("pickedStrat", "stratby", discv, discv[1])
      )
     ),
    fluidRow( 
     column(12,
       helpText(h4("projection methods:"))
      )
     ),
    fluidRow(
      column(6,
       selectInput("meth1", "left", nrd, nrd[1])),
      column(6,
       selectInput("meth2", "right", nrd, nrd[2]))
     ),
    fluidRow( 
     column(12,
       helpText(h4("dimensions to use:"))
      )
     ),
#
# the following needs to be conditioned on selected projections
# but for now we are assuming we need no more than 3
#
    fluidRow(
      column(6,
       numericInput("topx", "top x", 1, min=1, max=ncomp-1, step=1)),
      column(6,
       numericInput("topy", "top y", 2, min=2, max=ncomp, step=1))
     ),
    fluidRow(
      column(6,
       numericInput("botx", "bot x", 2, min=1, max=ncomp-1, step=1)),
      column(6,
       numericInput("boty", "bot y", 3, min=2, max=ncomp, step=1))
     ),
    fluidRow(
      column(12,
       helpText(h4("downloads:")))
     ),
    fluidRow(
      column(6,
       downloadButton("downloadData", "DE genes")),
      column(6,
       downloadButton("downloadData2", "cellSets"))
     ),
    fluidRow(
      column(12,
       helpText(h4("to conclude:")))
     ),
    fluidRow(
      column(12,
       actionButton("btnSend", "Stop app"))
     )
     ),
   mainPanel(
    tabsetPanel(
    tabPanel("scatter",
     helpText("Drag over points to make a selection, then use the selTable tab to generate DE signature; selections accumulate as this process is repeated."),
     fluidRow(
       column(6, d3scatterOutput("scatter1")),
       column(6, d3scatterOutput("scatter2"))
       ),
      fluidRow(
       column(6, d3scatterOutput("scatter3")),
       column(6, d3scatterOutput("scatter4"))
       )
     ), # end panel
    tabPanel("selTable",
     DT::dataTableOutput("summary")
     ),
    tabPanel("accum",
       helpText("Method and dimensions taken from bottom left panel"),
        fluidRow(column(12,
          plotOutput("accum")
          )),
        fluidRow(column(12,
          helpText("up to four biplots displayed here; use biplotSel on PcmpSel objects to work with additional selections")
          )),
        fluidRow(column(12,
          plotOutput("accum2", height="800px")
          ))
        ),
    tabPanel("about",
     helpText(h3("pcmp demonstrates crosstalk-based interactive graphics for surveying different dimension reduction procedures for data in SingleCellExperiment containers.  The reducedDims component must be populated with several reductions, each including at least 4 dimensions.   Different methods are used in the left and right columns, and different projection components can are used in the top and bottom rows, as selected using the method/top/bot controls.  The 'stratby' button will recolor points according to discrete covariates in the colData of the input object.")),
     helpText(h3("current input data structure:")),
     verbatimTextOutput("scedump"),
     helpText(h3("metadata strings:")),
     verbatimTextOutput("scedump2"),
     helpText(h3("pcmpApp can be demonstrated with the object pcmp::sce300xx, an extract from the Allen Brain Atlas RNA-seq data on anterior cingulate cortex (ACC) and primary visual cortex (VIS) brain regions.  Strata were formed using donor (3 levels) and region (2 levels) and 300 cells were sampled at random in each stratum.  The murLung3k app at vjcitn.shinyapps.io uses an extract from the Tabula Muris project data focused on a collection of cells from the mouse lung; the uncorrected data are available in the github repo vjcitn/pcmpshin in data/mouse3kf.rda."))
   )
  )
  )
  ) 
# end UI


# server for pcmp app
# defines data flow for a pair of projection types, each shown
# in two view based on different choices of dimensions

#basicServer <- function(sce) function(input, output, session) {
SERVER <- function(input, output, session) {
  requireNamespace("limma")
#
# add colnames as a column in colData -- uses .cellid field silently
#
#  sce$.cellid = colnames(sce)
#
# create a named data.frame combining the reducedDims with the colData
#
  rd = reducedDims(sce)
  nmeth = length(rd) # list of matrices of projected data
  methnames = names(rd)
  nrd = names(rd)
  ncomps = vapply(rd, ncol, numeric(1))
#  stopifnot(all(ncomps == ncomps[1]))  # requires balanced representation 
    # of all projections
  ncomp <- ncomps[1]
  ncomp = 3 # MUST FIX !
  indf = data.frame(do.call(cbind, as.list(rd)))
  cn = paste0(nrd[1], 1:ncomp)
  if (length(nrd) > 1) {
      for (j in 2:length(nrd)) cn = c(cn, paste0(nrd[j], 1:ncomp))
      }
  colnames(indf) = cn
  indf <- as.data.frame(cbind(indf, colData(sce)))
  
#
# build the formulas needed for d3scatter
#
  fmlist = lapply(methnames, function(x) list())
  names(fmlist) = methnames
  for (i in 1:nmeth) {
   curtags = paste0(methnames[i], 1:ncomp)
   fmlist[[i]] = lapply(curtags, function(x) as.formula(c("~", x)))
   names(fmlist[[i]]) = curtags
  }
  
#
# build the shared data
#
  enhDf = reactive({
   indf$strat = colData(sce)[[input$pickedStrat]]
   indf$key = 1:nrow(indf)
   indf
   })  

  shared_dat <- SharedData$new(enhDf) #enhDf, key=~key)

#
# set up reactive download entities: table of limma results, table of selected cells with selection sequence number
#
    output$downloadData <- downloadHandler(
       filename = function() {
         paste('data-', Sys.Date(), '.csv', sep='')
       },
       content = function(con) {
         write.csv(.GlobalEnv$.pcmpTab, con)
       }
     )
    output$downloadData2 <- downloadHandler(
       filename = function() {
         paste('data-', Sys.Date(), '.csv', sep='')
       },
       content = function(con) {
         dat = .GlobalEnv$.pcmpSelCells
         nsel = length(dat)
         selind = rep(1:nsel,sapply(dat,length))
         ans = data.frame(group=selind, cellid=unlist(dat))
         write.csv(ans, con)
       }
     )

#
# for the 'about' tab, show the SCE in use and some metadata
#

    output$scedump = renderPrint({
        print(sce)
    })
    output$scedump2 = renderPrint({
        print(metadata(sce)[c("note", "origin")])
    })


#
# produce the panels
#
  output$scatter1 <- renderD3scatter({
    methx = paste0(input$meth1, input$topx)
    methy = paste0(input$meth1, input$topy)
    d3scatter(shared_dat, fmlist[[input$meth1]][[methx]], 
            fmlist[[input$meth1]][[methy]], ~strat, width = "100%")
  })
  output$scatter2 <- renderD3scatter({
    methx = paste0(input$meth2, input$topx)
    methy = paste0(input$meth2, input$topy)
    d3scatter(shared_dat, fmlist[[input$meth2]][[methx]], 
            fmlist[[input$meth2]][[methy]], ~strat, width = "100%")
  })
  output$scatter3 <- renderD3scatter({
    methx = paste0(input$meth1, input$botx)
    methy = paste0(input$meth1, input$boty)
    d3scatter(shared_dat, fmlist[[input$meth1]][[methx]], 
            fmlist[[input$meth1]][[methy]], ~strat, width = "100%")
  })
  output$scatter4 <- renderD3scatter({
    methx = paste0(input$meth2, input$botx)
    methy = paste0(input$meth2, input$boty)
    d3scatter(shared_dat, fmlist[[input$meth2]][[methx]], 
            fmlist[[input$meth2]][[methy]], ~strat, width = "100%")
  })
#  output$try3d <- renderScatterplotThree({
#    colors = palette(rainbow(30))[ as.numeric(
#                      factor(colData(sce)[[input$pickedStrat]])) ]
#    scatterplot3js(PC1, PC2, PC3, crosstalk=shared_dat, brush=TRUE,
#       color=colors)
#    })

#
# collect the information on selections so far
#
 output$accum = renderPlot({
 invalidateLater(2500)
 ans = list(cells = .GlobalEnv$.pcmpSelCells, limmaTab=.GlobalEnv$.pcmpTab)
 tmp = new("PcmpSels", cellSets=ans$cells, geneTable=ans$limmaTab)
 replay(sce, tmp, input$meth1, input$botx, input$boty) 
 })

 output$accum2 = renderPlot({
 invalidateLater(2500)
 ans = list(cells = .GlobalEnv$.pcmpSelCells, limmaTab=.GlobalEnv$.pcmpTab)
 npl = min(c(length(ans[["cells"]]), 4))
 tmp = new("PcmpSels", cellSets=ans$cells, geneTable=ans$limmaTab)
 #replay(sce, tmp, input$meth1, input$botx, input$boty) 
 #opar = par(no.readonly=TRUE)
 #par(mfrow=c(2,2), mar=c(3,3,1,1))
 mym = matrix(1:4,byrow=T,nc=2)
 layout(mym, widths=c(2,2))
 for (i in 1:npl) biplotSel(sce, tmp, which=i, main=paste("selection", i))
 #par(opar)
 })
    
#
# very rudimentary approach to acquiring a signature of a selected group of cells
# use limma on log-transformed counts comparing selected to non-selected
# could do something to balance sample sizes ...
#

output$summary <- DT::renderDataTable({
    df <- shared_dat$data(withSelection = TRUE) %>%
      filter(selected_ | is.na(selected_)) %>%
      mutate(selected_ = NULL)
    sel=rep(0, ncol(sce))
    names(sel) = colnames(sce)
    sel[df$.cellid] = 1
    mm = stats::model.matrix(~sel, data=data.frame(sel=sel))
   showNotification(paste("starting table processing", date()), id="limnote")
    X = log(assay(sce)+1)
    f1 = lmFit(X, mm)
    ef1 = eBayes(f1)
    options(digits=3)

    tt = topTable(ef1, 2, n=20)
    tt$featid = rownames(tt)
    if (!(".pcmpSelNum" %in% ls(.GlobalEnv, all.names=TRUE))) assign(".pcmpSelNum", 1, .GlobalEnv)
      else assign(".pcmpSelNum", .GlobalEnv$.pcmpSelNum + 1, .GlobalEnv)
    if (!(".pcmpSelCells" %in% ls(.GlobalEnv, all.names=TRUE))) assign(".pcmpSelCells", list(df$.cellid), .GlobalEnv)
      else assign(".pcmpSelCells", c(.GlobalEnv$.pcmpSelCells, list(df$.cellid)), .GlobalEnv)
    tt = cbind(tt, selnum=.GlobalEnv$.pcmpSelNum[1])
    if (!(".pcmpTab" %in% ls(.GlobalEnv, all.names=TRUE))) assign(".pcmpTab", tt, .GlobalEnv)
      else assign(".pcmpTab", rbind(.GlobalEnv$.pcmpTab, tt, make.row.names=FALSE), .GlobalEnv)
    ans = DT::formatRound(DT::datatable(tt), 1:7, digits=3)
   removeNotification(id="limnote")
    ans
  })

#
# prepare stop button
#

   observe({
                    if(input$btnSend > 0)
                        isolate({
                           stopApp(returnValue=0)
                        })  
           })  
  } # end server
 tmp = runApp(list(ui=UI, server=SERVER))
#
 ans = list(cells = .GlobalEnv$.pcmpSelCells, limmaTab=.GlobalEnv$.pcmpTab)
 new("PcmpSels", cellSets=ans$cells, geneTable=ans$limmaTab)
}

rdprops = function(sce) {
 rd = reducedDims(sce)
 rdnames = names(rd)
 rdncols = sapply(rd, ncol)
 list(names=rdnames, ncols=rdncols)
}

newApp = function(sce, cellIdTag = ".cellid", selectionPrefix="sel_",
   myenv = new.env()) {

##
## general annotation tasks
##
discv = discreteColdVars(sce) # find discrete vbls suitable for coloring
colData(sce)[[cellIdTag]] = colnames(sce) # add cell identifier in colData

newUI = fluidPage(
 sidebarPanel(width=3,
   fluidRow(
     column(12,
       helpText(h3(sprintf("pcmp %s: crosstalk-based interactive graphics \
for dimension reduction in single-cell transcriptomics. \ 
See the 'about' tab for more information.", pkgVersion())))
       )
      ),
    fluidRow( 
     column(12,
       selectInput("pickedStrat", "stratby", discv, discv[1])
      )
     ),
  helpText("demo"),
  uiOutput("rdpicker"),
  uiOutput("dimpickerTop"),
  uiOutput("dimpickerBot"),
    fluidRow(
      column(12,
       helpText(h4("downloads:")))
     ),
    fluidRow(
      column(6,
       downloadButton("downloadData", "DE genes")),
      column(6,
       downloadButton("downloadData2", "cellSets"))
     ),
    fluidRow(
      column(12,
       helpText(h4("to conclude:")))
     ),
    fluidRow(
      column(12,
       actionButton("btnSend", "Stop app"))
    )
  ), # end sidebarPanel
   mainPanel(
    tabsetPanel(
    tabPanel("scatter",
     helpText("Drag over points to make a selection, then use the selTable tab to generate DE signature; selections accumulate as this process is repeated."),
     fluidRow(
       column(6, d3scatterOutput("scatter1")),
       column(6, d3scatterOutput("scatter2"))
       ),
      fluidRow(
       column(6, d3scatterOutput("scatter3")),
       column(6, d3scatterOutput("scatter4"))
       )
     ), # end panel
    tabPanel("selTable",
     DT::dataTableOutput("summary")
     ),
    tabPanel("accum",
       helpText("Method and dimensions taken from bottom left panel"),
        fluidRow(column(12,
          plotOutput("accum")
          )),
        fluidRow(column(12,
          helpText("up to four biplots displayed here; use biplotSel on PcmpSel objects to work with additional selections")
          )),
        fluidRow(column(12,
          plotOutput("accum2", height="800px")
          ))
        ),
    tabPanel("about",
     helpText(h3("pcmp demonstrates crosstalk-based interactive graphics for surveying different dimension reduction procedures for data in SingleCellExperiment containers.  The reducedDims component must be populated with several reductions, each including at least 4 dimensions.   Different methods are used in the left and right columns, and different projection components can are used in the top and bottom rows, as selected using the method/top/bot controls.  The 'stratby' button will recolor points according to discrete covariates in the colData of the input object.")),
     helpText(h3("current input data structure:")),
     verbatimTextOutput("scedump"),
     helpText(h3("metadata strings:")),
     verbatimTextOutput("scedump2"),
     helpText(h3("pcmpApp can be demonstrated with the object pcmp::sce300xx, an extract from the Allen Brain Atlas RNA-seq data on anterior cingulate cortex (ACC) and primary visual cortex (VIS) brain regions.  Strata were formed using donor (3 levels) and region (2 levels) and 300 cells were sampled at random in each stratum.  The murLung3k app at vjcitn.shinyapps.io uses an extract from the Tabula Muris project data focused on a collection of cells from the mouse lung; the uncorrected data are available in the github repo vjcitn/pcmpshin in data/mouse3kf.rda."))
   )
  )
  )
#  ) 
# mainPanel(
#  helpText("demo2"),
#  textOutput("picks"),
#  plotOutput("demo3")
#  )  # end mainPanel
 ) # end fluidPage

newserver = function(input, output, session) {
#
# rdprops acquires names and dimensions of reducedDims component
#
 props = rdprops(sce)
 nc = props$ncols
 methnames = props$names
#
# indf will collect all reduced dimensions and cbind with colData
#
 indf = data.frame(do.call(cbind, as.list(reducedDims(sce))))
 print(nc)
 print(methnames)
 alln = unlist(lapply(1:length(nc), function(i) {
   paste0(methnames[i], "_", 1:nc[i])
   }))
 colnames(indf) = alln
 indf <- as.data.frame(cbind(indf, colData(sce)))
##
## data preparations for crosstalk
##
#
# build the shared data
#
  enhDf = reactive({
   indf$strat = colData(sce)[[input$pickedStrat]]
   indf$key = 1:nrow(indf)
   indf
   })
  shared_dat <- SharedData$new(enhDf) #enhDf, key=~key)
#
# given methods leftRD, rightRD, and topx/y, botx/y create
# the formulas needed by d3scatter
#
  reactiveFmlas = reactive({
    validate(
     need(input$rightRD, "select top right x coord")
     )
    validate(
     need(input$TopX, "select top left x coord")
     )
   list(
    upleftX = as.formula(c("~", paste0(input$leftRD, "_", input$TopX))),
    upleftY = as.formula(c("~", paste0(input$leftRD, "_", input$TopY))),
    uprightX = as.formula(c("~", paste0(input$rightRD, "_", input$TopX))),
    uprightY = as.formula(c("~", paste0(input$rightRD, "_", input$TopY))),
    lowleftX = as.formula(c("~", paste0(input$leftRD, "_", input$BotX))),
    lowleftY = as.formula(c("~", paste0(input$leftRD, "_", input$BotY))),
    lowrightX = as.formula(c("~", paste0(input$rightRD, "_", input$BotX))),
    lowrightY = as.formula(c("~", paste0(input$rightRD, "_", input$BotY)))
     ) })

##
## basic UI control rendering tasks
##
 output$rdpicker = renderUI( {
  fluidRow(
   column(6,
    selectInput("leftRD", "left", props$names, props$names[1])), # leftRD: left reduced dim
   column(6,
    selectInput("rightRD", "right", props$names, props$names[2]))
   )
   } )
  output$dimpickerTop = renderUI( {
    validate(
     need(input$rightRD, "select top right x coord")
     )
    allNcols = props$ncols
    activeNcols = allNcols[c(input$rightRD, input$leftRD)]
    minx = min(activeNcols)
   fluidRow(
    column(6,
     numericInput("TopX", "topX", value=1, min=1, max=minx-1)),
    column(6,
     numericInput("TopY", "topY", value=2, min=2, max=minx))
    ) # end fluidRow
   } )
  output$dimpickerBot = renderUI( {
    validate(
     need(input$rightRD, "select top right x coord")
     )
    allNcols = props$ncols
    activeNcols = allNcols[c(input$rightRD, input$leftRD)]
    minx = min(activeNcols)
   fluidRow(
    column(6,
     numericInput("BotX", "botX", value=1, min=1, max=minx-1)),
    column(6,
     numericInput("BotY", "botY", value=2, min=2, max=minx))
    ) # end fluidRow
   } )
  output$picks = renderText({
    validate(
     need(input$rightRD, "select top right x coord")
     )
   paste(input$leftRD, input$rightRD)
   })
##
## end basic UI control rendering
##
  output$demo3 = renderPlot({
    validate(
     need(input$rightRD, "select top right x coord")
     )
    validate(
     need(input$leftRD, "select top left x coord")
     )
    validate(
     need(input$TopX, "select top left x coord")
     )
   par(mfrow=c(2,2))
   rd = reducedDims(sce)
   plot(rd[[input$leftRD]][ ,input$TopX ],
        rd[[input$leftRD]][ ,input$TopY ], col=sce$PIMMquart)
   for (i in 2:4) plot(1,1)
   })
  output$scatter1 <- renderD3scatter({
    validate(
     need(input$leftRD, "select top left x coord")
     )
    validate(
     need(input$TopX, "select top left x coord")
     )
    d3scatter(shared_dat, reactiveFmlas()$upleftX, reactiveFmlas()$upleftY,
            ~strat, width = "100%")
  })
  output$scatter2 <- renderD3scatter({
    validate(
     need(input$leftRD, "select top left x coord")
     )
    validate(
     need(input$TopX, "select top left x coord")
     )
    d3scatter(shared_dat, reactiveFmlas()$uprightX, reactiveFmlas()$uprightY,
            ~strat, width = "100%")
  })
  output$scatter3 <- renderD3scatter({
    validate(
     need(input$leftRD, "select top left x coord")
     )
    validate(
     need(input$TopX, "select top left x coord")
     )
    validate(
     need(input$BotX, "select top left x coord")
     )
    d3scatter(shared_dat, reactiveFmlas()$lowleftX, reactiveFmlas()$lowleftY,
            ~strat, width = "100%")
  })
  output$scatter4 <- renderD3scatter({
    validate(
     need(input$leftRD, "select top left x coord")
     )
    validate(
     need(input$TopX, "select top left x coord")
     )
    validate(
     need(input$BotX, "select bot left x coord")
     )
    d3scatter(shared_dat, reactiveFmlas()$lowrightX, reactiveFmlas()$lowrightY,
            ~strat, width = "100%")
  })
#
# operations for the 'about' tab
#
    output$scedump = renderPrint({
        print(sce)
    })  
    output$scedump2 = renderPrint({
        print(metadata(sce))
    })  
#
# collect selection reactively
#
    getSels = reactive({
     df <- shared_dat$data(withSelection = TRUE) %>%
      filter(selected_ | is.na(selected_)) %>%
      mutate(selected_ = NULL)
     sel=rep(0, ncol(sce))
     names(sel) = colnames(sce)
     sel[df$.cellid] = 1
  # use the environment to retain selection number
     if (is.null(myenv$selnum)) myenv$selnum = 1
      else myenv$selnum = myenv$selnum + 1
     colData(sce) <<- cbind(colData(sce), sel)
     names(colData(sce))[ncol(colData(sce))] <<- paste0(
              selectionPrefix, myenv$selnum)
     sce
     })
#
# collect selection and run limma
#
output$summary <- DT::renderDataTable({
   sce <<- getSels()
   print(names(colData(sce)))
   print(table(sce$sel_1))
#    df <- shared_dat$data(withSelection = TRUE) %>%
#      filter(selected_ | is.na(selected_)) %>%
#      mutate(selected_ = NULL)
#    sel=rep(0, ncol(sce))
#    names(sel) = colnames(sce)
#    sel[df$.cellid] = 1
    sel = colData(sce)[,ncol(colData(sce))] # use rightmost col for current selection
    mm = stats::model.matrix(~sel, data=data.frame(sel=sel))
   showNotification(paste("starting table processing", date()), id="limnote")
    X = log(assay(sce)+1)
    f1 = lmFit(X, mm)
    ef1 = eBayes(f1)
    options(digits=3)

    tt = topTable(ef1, 2, n=20)
    tt$featid = rownames(tt)
#    if (!(".pcmpSelNum" %in% ls(myenv, all.names=TRUE))) assign(".pcmpSelNum", 1, myenv)
#      else assign(".pcmpSelNum", myenv$.pcmpSelNum + 1, envir=myenv)
#    if (!(".pcmpSelCells" %in% ls(myenv, all.names=TRUE))) assign(".pcmpSelCells", list(df$.cellid), envir=myenv)
#      else assign(".pcmpSelCells", c(myenv$.pcmpSelCells, list(df$.cellid)), envir=myenv)
#    tt = cbind(tt, selnum=myenv$.pcmpSelNum[1])
#    if (!(".pcmpTab" %in% ls(myenv, all.names=TRUE))) assign(".pcmpTab", tt, envir=myenv)
#      else assign(".pcmpTab", rbind(myenv$.pcmpTab, tt, make.row.names=FALSE), envir=myenv)
    ans = DT::formatRound(DT::datatable(tt), 1:7, digits=3)
   removeNotification(id="limnote")
    ans
  })
#
# prepare stop button
#

   observe({
            if(input$btnSend > 0)
               isolate({
                 stopApp(returnValue=0)
                        })  
           })  


 } #end server

runApp(list(ui=newUI, server=newserver))
list(myenv, sce)
}
