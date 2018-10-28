library(data.table)
visx = fread("human_VISp_2018-10-04_exon-matrix.csv")
visex = as.matrix(visx[,-1])
gn = visx[,1,drop=TRUE][[1]]
rownames(visex) = gn
cd = read.csv("human_VISp_2018-10-04_samples-columns.csv",
  stringsAsFactors=FALSE)
rownames(cd) = cd[,1]
library(SummarizedExperiment)
VISexons = SummarizedExperiment(visex)
gfn = "human_VISp_2018-10-04_genes-rows.csv"
gr = read.csv(gfn, stringsAsFactors=FALSE)
rownames(gr) = gr[,1]
colData(VISexons) = DataFrame(cd)
rowData(VISexons) = DataFrame(gr)
assays(VISexons) = SimpleList(counts=visex)
metadata(VISexons)$date=date()
metadata(VISexons)$url = "http://celltypes.brain-map.org/api/v2/well_known_file_download/738610707"
