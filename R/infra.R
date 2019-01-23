pkgVersion = function() as.character(read.dcf(system.file("DESCRIPTION", package="pcmp"))[,"Version"])

#' find discrete variables in colData
#' @import SingleCellExperiment crosstalk d3scatter dplyr shiny ggplot2
#' @import Rtsne irlba umap methods
#' @importFrom S4Vectors SimpleList
#' @importFrom SummarizedExperiment colData<-
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
  names(ans) = pn # names(projectors$retrievers)
  ini = reducedDims(sce)
  if (length(ini) == 0) reducedDims(sce) = SimpleList(ans)
  else reducedDims(sce) = c(ini, SimpleList(ans))
  sce
}

rdprops = function(sce) {
 rd = reducedDims(sce)
 rdnames = names(rd)
 rdncols = sapply(rd, ncol)
 list(names=rdnames, ncols=rdncols)
}

#' Define a shiny app for comparing projections
#' @importFrom limma lmFit eBayes topTable
#' @importFrom stats model.matrix
#' @importFrom grDevices dev.off
#' @importFrom graphics layout legend par
#' @importFrom utils tail write.csv
#' @param sce SingleCellExperiment instance with reducedDims populated
#' with various candidates for projection from assay data --
#' @param assayind numeric(1) 'i' argument to assay(sce, i), defaults to 1, for extracting expression values
#' @param cellIdTag character(1) defaulting to ".cellid", added to colData(sce) for tracking selections
#' @param selectionPrefix character(1) defaulting to "sel_", prefix for the name of a selection indicator variable added to colData(sce)
#' @param numDiscGenes numeric(1) used to limit number of differentially expressed genes tabulated per selection
#' @note In sce300xx, names(reducedDims) are 'PPCA', 'UMAP*', 'TSNE*'
#' for projections generated by
#' irlba::irlba, umap::umap, and Rtsne::Rtsne respectively.
#' Prior to pcmp version 0.4.0, the projections were computed and prepared
#' differently.  As Rtsne no longer allows dims > 3, we have separately
#' computed and stored t-SNE projections into 2- and 3- dimensional
#' spaces, labeled TSNE2 and TSNE3.  (In previous versions of pcmp, the 4-d t-SNE was computed
#' and marginal planar projections were provided.)
#' Note that tabulation of DE genes can be very slow if
#' the input SingleCellExperiment has more than 10k rows.
#' @return will return an updated SingleCellExperiment
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
#'  try(head(metadata(chk)$limmaTabs,3))  # likewise
#'  try(table(metadata(chk)$limmaTabs$featid))  # likewise
#' }
#' @export
pcmpApp = function(sce, assayind=1, cellIdTag = ".cellid", 
     selectionPrefix="sel_", numDiscGenes=20) {
##
## general annotation tasks
##
discv = discreteColdVars(sce) # find discrete vbls suitable for coloring
colData(sce)[[cellIdTag]] = colnames(sce) # add cell identifier in colData
myenv = new.env()  # can this be dropped?
cl = match.call()
metadata(sce)$params = allParamBindings(cl)
##
## define app UI
##
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
#    fluidRow(
#      column(6,
#       downloadButton("downloadData", "DE genes")),
#      column(6,
#       downloadButton("downloadData2", "cellSets"))
#     ),
    fluidRow(
      column(12,
       helpText(h4("return updated SCE:")))
     ),
    fluidRow(
      column(12,
       actionButton("btnSend", "Stop app"))
    )
  ), # end sidebarPanel
   mainPanel(
    tabsetPanel(
    tabPanel("scatter", id="tabs",
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
     DT::dataTableOutput("summary"), id="seller"
     ),
    tabPanel("accum", id="accer",
       helpText("Method and dimensions taken from bottom left panel;"),
       helpText("Toggle the projection method to update the plot"),
        fluidRow(column(12,
          plotOutput("accum")
          )),
        fluidRow(column(12,
          helpText("up to four biplots displayed here; use biplotSel on PcmpSel objects to work with additional selections")
          )),
        fluidRow(column(12,
#          plotOutput("accum2", height="800px")
          plotOutput("accum2", height="800px")
          ))
        ),
    tabPanel("about", id="about",
     helpText(h3("pcmp demonstrates crosstalk-based interactive graphics for surveying different dimension reduction procedures for data in SingleCellExperiment containers.  The reducedDims component must be populated with several reductions.   Different methods are used in the left and right columns, and different projection components can be displayed in the top and bottom rows, as selected using the method/top/bot controls.  The 'stratby' button will recolor points according to discrete covariates in the colData of the input object.")),
     helpText(h3("current input data structure:")),
     verbatimTextOutput("scedump"),
     helpText(h3("metadata strings:")),
     verbatimTextOutput("scedump2"),
     helpText(h3("pcmpApp can be demonstrated with the object pcmp::sce300xx, an extract from the Allen Brain Atlas RNA-seq data on anterior cingulate cortex (ACC) and primary visual cortex (VIS) brain regions.  Strata were formed using donor (3 levels) and region (2 levels) and 300 cells were sampled at random in each stratum.  The murLung3k app at vjcitn.shinyapps.io uses an extract from the Tabula Muris project data focused on a collection of cells from the mouse lung; the uncorrected data are available in the github repo vjcitn/pcmpshin in data/mouse3kf.rda."))
   )
  )
  )
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
    print(activeNcols)
    minx = min(activeNcols)
    print(minx)
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
    print(minx)
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
   sce <<- getSels()  # must update local SCE
   sel = colData(sce)[,ncol(colData(sce))] # use rightmost col for current selection
   mm = stats::model.matrix(~sel, data=data.frame(sel=sel))
   showNotification(paste("starting table processing", date()), id="limnote")
   X = log(assay(sce, assayind)+1)
   f1 = lmFit(X, mm)
   ef1 = eBayes(f1)
   options(digits=3)
   tt = topTable(ef1, 2, n=numDiscGenes)
   tt$featid = rownames(tt)
   ans = DT::formatRound(DT::datatable(tt), 1:7, digits=3)
   ntabs = length(metadata(sce)$limmaTabs) 
   if (ntabs == 0) {
      tt$selnum = 1
      metadata(sce)$limmaTabs <<- tt
      }
   else {
      lt = metadata(sce)$limmaTabs
      last = tail(lt$selnum,1)
      tt$selnum = last+1
      metadata(sce)$limmaTabs <<- rbind(metadata(sce)$limmaTabs, tt)
      }
   removeNotification(id="limnote")
   ans
  })

   accumDF = reactive({
     validate(need(input$leftRD, "select top left x"))
     rd = reducedDims(sce)
     ndf = data.frame(x = rd[[input$leftRD]][, input$TopX], 
              y=rd[[input$leftRD]][, input$TopY])
     ndf$grp = "unsel"
     cd = colData(sce)
     reneed = paste0("^", selectionPrefix)
     selinds = grep(reneed, names(cd))
     if (length(selinds)>0) {
         cur = 0
         for (i in selinds) {
           cur = cur+1
           ndf$grp[which(cd[,i]==1)] = paste0(selectionPrefix, cur)
           }
         }
     ndf
   })
   output$accum = renderPlot( {
     ndf = accumDF()
     ggplot(ndf, aes(x=x, y=y, colour=grp)) + geom_point() +
          guides(colour = guide_legend(override.aes = list(size=10),
             label.theme = element_text( size = 15),
             title.theme = element_text( size = 15)))
     })
   output$accum2 = renderPlot( {
       lt = metadata(sce)$limmaTabs
       validate(need(lt, "make a selection"))
       cd = colData(sce)
       cdselinds = grep(selectionPrefix, names(cd))
       selinds = cd[,cdselinds]
       nsels = length(unique(lt$selnum))
       mym = matrix(1:4,byrow=T,nc=2)
       layout(mym, widths=c(2,2))
       for (i in 1:min(c(nsels,4))) {
        curfeat = lt$featid[which(lt$selnum==i)]
        curass = t(log(assay(sce[curfeat, which(selinds[,i]==1)])+1))
        curpc = prcomp(curass)
        biplot(curpc, xlabs=rep(".", sum(selinds[,i]==1)), expand=.8, cex=c(2,1))
        }
      } )
#
# prepare stop button, which will recompute the selection map
# and create a rasterized image of the map, bound in the metadata slot
#
   observe({
            if(input$btnSend > 0)
               isolate({
                 stopApp(returnValue=0)
                        })  
           })  


 } #end server

runApp(list(ui=newUI, server=newserver))
sce
}


allParamBindings = function(cl) {
 # cl is a call
 stopifnot(inherits(cl, "call"))
 fname = as.character(as.list(cl)[[1]])
 defaults = as.list(get(fname))
 defaults = defaults[-length(defaults)]
 bound = as.list(cl)[-1]
 tmp = list(defaults=defaults, bound=bound)
 usedef = setdiff(names(tmp$defaults), names(tmp$bound))
 kp = names(tmp$bound)
 c(tmp$bound[kp], tmp$defaults[usedef])
}

#' plot selection map from output of pcmpApp
#' @param esce an extended SingleCellExperiment, with metadata extended using pcmpApp
#' @param proj character(1) name of a reducedDims component of esce
#' @param dim1 numeric(1) dimension to use as x axis
#' @param dim2 numeric(1) dimension to use as y axis
#' @param alpha numeric(1) defaults to .5 for transparency used in geom_point
#' @export
plotSelMap = function(esce, proj, dim1=1, dim2=2, alpha=.5) {
 stopifnot("params" %in% names(metadata(esce)))
 stopifnot(proj %in% names(reducedDims(esce)))
 selPref = metadata(esce)$params$selectionPrefix
 selinds = grep(selPref, names(colData(esce)),value=TRUE)
 seldat = colData(esce)[,selinds]
 grp = data.matrix(seldat) %*% (1:ncol(seldat))
 if (any(grp==0)) grp[which(grp==0)] = ncol(seldat)+1
 rd = reducedDims(esce)[[proj]]
 pldf = data.frame(x=rd[,dim1], y=rd[,dim2], grp=c(names(seldat), "unsel")[grp])
 names(pldf)[1:2] = paste0(proj, "_", c(dim1, dim2))
 ggplot(pldf, aes_string(x=names(pldf)[1], y=names(pldf[2]),
     colour="grp")) + geom_point(alpha=alpha)
}

getGroup = function (esce, gvname= "group" ) {
    stopifnot("params" %in% names(metadata(esce)))
    selPref = metadata(esce)$params$selectionPrefix
    selinds = grep(selPref, names(colData(esce)), value = TRUE)
    seldat = colData(esce)[, selinds]
    grp = data.matrix(seldat) %*% (1:ncol(seldat))
    if (any(grp == 0)) 
        grp[which(grp == 0)] = ncol(seldat) + 1
    grp
}

#SharedData$new(d, ~key) %>%
#  plot_ly(x = ~x, y = ~y) %>%
#  highlight("plotly_selected") %>%
#  layout(dragmode = "lasso")

#' simple biplot display of transformed expression data for a subset of cluster discriminating genes
#' @param sce a SingleCellExperiment as produced by pcmpApp
#' @param selnum numeric(1) the cluster to display
#' @param assaytx a function to be applied to the assay data
#' @param ntop numeric(1) number of genes to use in PCA
#' @param \dots additional parameters to biplot.prcomp, which should
#' exclude xlabs, expand, cex, and col which are hard coded.
#' @examples
#' biplotSel(pcmp::vigAccum, selnum=4, ntop=10)
#' @export
biplotSel = function(sce, selnum=1, 
 assaytx = function(x)log(x+1),  ntop=6, ...) {
 stopifnot("limmaTabs" %in% names(metadata(sce)))
 lt = metadata(sce)$limmaTabs
 stopifnot(selnum %in% unique(lt$selnum))
 genes = lt$featid[which(lt$selnum==selnum)]
 dat = assay(sce)[genes[1:ntop], which(getGroup(sce)==selnum)]
 biplot(prcomp(t(assaytx(dat))), 
    xlabs=rep(".", ncol(dat)), expand=.8, cex=c(2,1), col=c(1,2), ...)
}

