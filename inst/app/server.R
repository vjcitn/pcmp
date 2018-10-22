
server <- function(input, output, session) {
  library(pcmp)
  library(d3scatter)
  library(crosstalk)

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

  shared_dat <- SharedData$new(indf) #enhDf, key=~key)

  output$scatter1 <- renderD3scatter({
    methx = paste0(input$meth1, input$topx)
    methy = paste0(input$meth1, input$topy)
    d3scatter(shared_dat, fmlist[[input$meth1]][[methx]], 
            fmlist[[input$meth1]][[methy]], ~donor_id, width = "100%")
  })
  output$scatter2 <- renderD3scatter({
    methx = paste0(input$meth2, input$topx)
    methy = paste0(input$meth2, input$topy)
    d3scatter(shared_dat, fmlist[[input$meth2]][[methx]], 
            fmlist[[input$meth2]][[methy]], ~donor_id, width = "100%")
  })
  output$scatter3 <- renderD3scatter({
    methx = paste0(input$meth1, input$botx)
    methy = paste0(input$meth1, input$boty)
    d3scatter(shared_dat, fmlist[[input$meth1]][[methx]], 
            fmlist[[input$meth1]][[methy]], ~donor_id, width = "100%")
  })
  output$scatter4 <- renderD3scatter({
    methx = paste0(input$meth2, input$botx)
    methy = paste0(input$meth2, input$boty)
    d3scatter(shared_dat, fmlist[[input$meth2]][[methx]], 
            fmlist[[input$meth2]][[methy]], ~donor_id, width = "100%")
  })
  output$scatter5 <- renderD3scatter({
    methx = paste0(input$meth1, input$topx)
    methy = paste0(input$meth1, input$topy)
    d3scatter(shared_dat, fmlist[[input$meth1]][[methx]], 
            fmlist[[input$meth1]][[methy]], ~brain_region, width = "100%")
  })
  output$scatter6 <- renderD3scatter({
    methx = paste0(input$meth2, input$topx)
    methy = paste0(input$meth2, input$topy)
    d3scatter(shared_dat, fmlist[[input$meth2]][[methx]], 
            fmlist[[input$meth2]][[methy]], ~brain_region, width = "100%")
  })
  output$scatter7 <- renderD3scatter({
    methx = paste0(input$meth1, input$botx)
    methy = paste0(input$meth1, input$boty)
    d3scatter(shared_dat, fmlist[[input$meth1]][[methx]], 
            fmlist[[input$meth1]][[methy]], ~brain_region, width = "100%")
  })
  output$scatter8 <- renderD3scatter({
    methx = paste0(input$meth2, input$botx)
    methy = paste0(input$meth2, input$boty)
    d3scatter(shared_dat, fmlist[[input$meth2]][[methx]], 
            fmlist[[input$meth2]][[methy]], ~brain_region, width = "100%")
  })
  output$scatter9 <- renderD3scatter({
    methx = paste0(input$meth1, input$topx)
    methy = paste0(input$meth1, input$topy)
    d3scatter(shared_dat, fmlist[[input$meth1]][[methx]], 
            fmlist[[input$meth1]][[methy]], ~brain_subregion, width = "100%")
  })
  output$scatter10 <- renderD3scatter({
    methx = paste0(input$meth2, input$topx)
    methy = paste0(input$meth2, input$topy)
    d3scatter(shared_dat, fmlist[[input$meth2]][[methx]], 
            fmlist[[input$meth2]][[methy]], ~brain_subregion, width = "100%")
  })
  output$scatter11 <- renderD3scatter({
    methx = paste0(input$meth1, input$botx)
    methy = paste0(input$meth1, input$boty)
    d3scatter(shared_dat, fmlist[[input$meth1]][[methx]], 
            fmlist[[input$meth1]][[methy]], ~brain_subregion, width = "100%")
  })
  output$scatter12 <- renderD3scatter({
    methx = paste0(input$meth2, input$botx)
    methy = paste0(input$meth2, input$boty)
    d3scatter(shared_dat, fmlist[[input$meth2]][[methx]], 
            fmlist[[input$meth2]][[methy]], ~brain_subregion, width = "100%")
  })

   observe({
                    if(input$btnSend > 0)
                        isolate({
                           metadata(se)$"abc" = "def"
                           stopApp(returnValue=se)
                        })  
           })  


}

