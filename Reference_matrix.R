library(dplyr)
library(tidyr)
library(scater)
library(scuttle)
library(ggplot2)
library(ggridges)
library(patchwork)
library(DropletUtils)
library(BiocParallel)
library(SingleCellExperiment)
bp <- SnowParam(10)
#source("utils.R")
set.seed(250112)

lys <- readRDS("C:/CosMx/Synovial data/sce-lv2.rds")
new_h5_path <- "C:/CosMx/Synovial data/sce-raw/assays.h5"  

lab <- list(
  imm=c(a="T/NK cell", b="B cell"),
  mye=c(a="mac1", b="mac2", c="mac3", e="mac4", 
        d="nphil", f="nphil"),
  str=c(
    a="fib1", b="fib2", d="fib3", 
    c="endo1", g="endo2", h="endo3", 
    e="blood", f="mural/peri"))

for (sub in names(lab)) {
  idx <- match(lys[[sub]]$kid_lv2, names(lab[[sub]]))
  lys[[sub]]$lab <- lab[[sub]][idx]
}
lys <- lapply(lys, \(.) { rowData(.)$sel <- NULL; .})

# Update the HDF5-backed assays
for (sub in names(lys)) {
  for (assay_name in assayNames(lys[[sub]])) {
    lys[[sub]]@assays@data@listData[["counts"]]@seed@seed@seed@filepath = new_h5_path
    lys[[sub]]@assays@data@listData[["logcounts"]]@seed@seed@seed@seed@seed@seed@seed@filepath = new_h5_path
  }
}

sce <- do.call(cbind, lys); sce <- sce[, !is.na(sce$lab)]
sce <- sce[, !is.na(sce$lab)]

ref_matrix <- aggregateAcrossCells(
  sce, ids = sce$lab, 
  statistics = "mean", use.assay.type = "logcounts"
)
