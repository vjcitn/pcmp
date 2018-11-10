library(testthat)
library(pcmp)

# [1] "acc4vis2"          "addPIMM"           "addProjections"   
# [4] "biplotSel"         "cellSets"          "defaultProjectors"
# [7] "discreteColdVars"  "geneTable"         "pcmpApp"          
#[10] "PcmpSels"          "replay"            "sce300xx"         
#[13] "stratsampInds"     "vigAccum"         

context("verify nonregression of serialized data, and test geneTable and cellSets functions")

test_that("acc4vis2 has expected content", {
  expect_true(is(acc4vis2, "PcmpSels"))
  expect_true(nrow(geneTable(acc4vis2))==120)
  expect_true(ncol(geneTable(acc4vis2))==8)
  expect_true(all.equal(head(cellSets(acc4vis2)[[1]]),
    c("SM-DAIFX_S79_E1-50", "SM-DAIFZ_S85_E1-50", "SM-DAIG2_S67_E1-50", 
      "SM-DAIG2_S93_E1-50", "SM-DAIG3_S78_E1-50", "SM-DAIGB_S93_E1-50")))
})

test_that("vigAccum has expected content", {
  expect_true(is(vigAccum, "PcmpSels"))
  expect_true(length(cellSets(vigAccum))==6)
})

test_that("state of defaultProjectors is known", {
  expect_true(all.equal(names(defaultProjectors()), c("projectors", 
     "retrievers")))
  expect_true(all.equal(names(defaultProjectors()[[1]]), c("PPCA", "UMAP", "TSNE")))
})

someRowSums = c(`3.8-1.2` = 1, `3.8-1.3` = 6, `3.8-1.4` = 0, `3.8-1.5` = 3, 
`5-HT3C2` = 449, A1BG = 6793, `A1BG-AS1` = 720, A1CF = 900, A2M = 19612, 
`A2M-AS1` = 3046, A2ML1 = 1496, A2MP1 = 465, A3GALT2 = 0, A4GALT = 3, 
A4GNT = 9, AA06 = 80, AAAS = 52244, AACS = 68461, AACSP1 = 2614, 
AADAC = 17)


test_that("sce300xx has expected content", {
  expect_true(is(sce300xx, "SingleCellExperiment"))
  expect_true(nrow(sce300xx)==50281)
  expect_true(ncol(sce300xx)==1800)
  expect_true(ncol(colData(sce300xx))==34)
  expect_true(all.equal(rowSums(assay(sce300xx[1:20,])), someRowSums))
})

# addPIMM requires two sces; fake it

test_that("addPIMM gets arithmetic right", {
  lit = sce300xx[1:20,1:30]
  expect_true(all(na.omit(addPIMM(lit,lit)$PIMM)==0.5))
})

context("check discreteColdVars")

coldv = c("donor_id", "sex", "age_days", "brain_region", "brain_subregion", 
"facs_date", "facs_sort_criteria", "strat", "PIMMquart")


test_that("finds right discrete vars", {
  expect_true(all.equal(discreteColdVars(sce300xx), coldv))
})
