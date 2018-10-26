
server <- function(input, output, session) {
  library(pcmp)
  library(d3scatter)
  library(crosstalk)
# hardwired for shinyapps.io -- must keep code
# in sync with more flexible pcmpApp
  sce = pcmp::sce300xx 
  rd = reducedDims(sce)
  nrd = names(rd)
  ncomps = vapply(rd, ncol, numeric(1))
  stopifnot(all(ncomps == ncomps[1]))
  ncomp <- ncomps[1]
  indf = data.frame(do.call(cbind, as.list(rd)))
  cn = paste0(nrd[1], 1:ncomp)
  if (length(nrd) > 1) {
      for (j in 2:length(nrd)) cn = c(cn, paste0(nrd[j], 1:ncomp))
  }
  colnames(indf) = cn
  indf <- as.data.frame(cbind(indf, colData(sce)))
  
  PCtags = paste0("PC", 1:ncomp)
  UMtags = paste0("UM", 1:ncomp)
  TStags = paste0("TS", 1:ncomp)
  PCfmlas = lapply(PCtags, function(x) as.formula(c("~", x)))
  UMfmlas = lapply(UMtags, function(x) as.formula(c("~", x)))
  TSfmlas = lapply(TStags, function(x) as.formula(c("~", x)))
  names(PCfmlas) = PCtags
  names(UMfmlas) = UMtags
  names(TSfmlas) = TStags
  fmlist = list(PC=PCfmlas, UM=UMfmlas, TS=TSfmlas)

#  enhDf = reactive({
#   indf$strat = colData(se)[[input$str2use]]
#   indf$key = 1:nrow(indf)
#   print(head(indf))
#   indf
#   })
#

  enhDf = reactive({
   indf$strat = colData(sce)[[input$pickedStrat]]
   indf$key = 1:nrow(indf)
   indf
   })  


  shared_dat <- SharedData$new(enhDf) #enhDf, key=~key)

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

output$summary <- DT::renderDataTable({
    df <- shared_dat$data(withSelection = TRUE) %>%
      filter(selected_ | is.na(selected_)) %>%
      mutate(selected_ = NULL)
    sel=rep(0, ncol(sce))
    names(sel) = colnames(sce)
    sel[df$seq_name] = 1
    mm = model.matrix(~sel, data=data.frame(sel=sel))
    library(limma)
   showNotification(paste("starting table processing", date()), id="limnote")
    X = log(assay(sce)+1)
    f1 = lmFit(X, mm)
    ef1 = eBayes(f1)
  print(paste0("finish lmFit", date()))
    options(digits=3)

    tt = topTable(ef1, 2, n=20)
    if (!(".pcmpSelNum" %in% ls(.GlobalEnv, all=TRUE))) assign(".pcmpSelNum", 1, .GlobalEnv)
      else assign(".pcmpSelNum", .GlobalEnv$.pcmpSelNum + 1, .GlobalEnv)
    if (!(".pcmpSelCells" %in% ls(.GlobalEnv, all=TRUE))) assign(".pcmpSelCells", list(df$seq_name), .GlobalEnv)
      else assign(".pcmpSelCells", c(.GlobalEnv$.pcmpSelCells, list(df$seq_name)), .GlobalEnv)
    tt = cbind(tt, selnum=.GlobalEnv$.pcmpSelNum[1])
    if (!(".pcmpTab" %in% ls(.GlobalEnv, all=TRUE))) assign(".pcmpTab", tt, .GlobalEnv)
      else assign(".pcmpTab", rbind(.GlobalEnv$.pcmpTab, tt), .GlobalEnv)
    ans = DT::formatRound(DT::datatable(tt), 2:7, digits=3)
   removeNotification(id="limnote")
    ans
  })

   observe({
                    if(input$btnSend > 0)
                        isolate({
                           stopApp(returnValue=0)
                        })  
           })  


}

