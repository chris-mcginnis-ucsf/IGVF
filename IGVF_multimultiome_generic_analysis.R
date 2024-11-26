######################################################
## Generic IGVF MULTI-multiome preprocessing script ##
## Chris McGinnis, PhD -- Satpathy Lab, Stanford #####
## November 26th 2024 ################################
######################################################

## initialize workspace
library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(ggplot2)
library(rhdf5)
library(deMULTIplex2)

## download from github: https://github.com/chris-mcginnis-ucsf/IGVF
load('~/bar.ref.Robj')
load('~/rna_whitelist.Robj')
load('~/atac_whitelist.Robj')

## set file paths -- example for experiment with 2 10x lanes
# from cellranger-atac count output
bed_path <- c('lane1_peaks.bed',
              'lane2_peaks.bed')
singlecell_path <- c('lane1_singlecell.csv',
                     'lane2_singlecell.csv')
frag_path <- c('lane1_fragments.tsv.gz',
               'lane2_fragments.tsv.gz')
# from cellranger count output
rna_path <- c('lane1_raw_feature_bc_matrix.h5',
              'lane2_raw_feature_bc_matrix.h5')
# multi-seq fastqs
multi_r1_path <- c('lane1_R1.fastq.gz',
                   'lane2_R1.fastq.gz')
multi_r2_path <- c('lane1_R2.fastq.gz',
                   'lane2_R2.fastq.gz')

## read in peak sets, convert to granges, unify, parse
plan("multisession", workers = 30)
options(future.globals.maxSize = 60000 * 424^2) 

peaks_lane1 <- read.table(file = bed_path[1], col.names = c("chr", "start", "end"))
gr_lane1 <- makeGRangesFromDataFrame(peaks_lane1)

peaks_lane2 <- read.table(file = bed_path[2], col.names = c("chr", "start", "end"))
gr_lane2 <- makeGRangesFromDataFrame(peaks_lane2)

peaks <- reduce(x = c(gr_lane1, gr_lane2))
peakwidths <- width(peaks)
peaks <- peaks[peakwidths < 4000 & peakwidths > 20]

## load metadata, identify nFrag threshold
meta_lane1 <- read.table(file = singlecell_path[1], stringsAsFactors = FALSE, sep = ",", header = TRUE, row.names = 1)[-1, ]
meta_lane2 <- read.table(file = singlecell_path[2], stringsAsFactors = FALSE, sep = ",", header = TRUE, row.names = 1)[-1, ]

meta_lane1 <- meta_lane1[which(meta_lane1$passed_filters >= 100), ] # set very low threshold at first
meta_lane2 <- meta_lane2[which(meta_lane2$passed_filters >= 100), ]
meta_lane1[,'lane'] <- 'lane1'
meta_lane2[,'lane'] <- 'lane2'

temp <- rbind(meta_lane1, meta_lane2)
ggplot(temp, aes(x = lane, y = log10(passed_filters), fill = lane)) + geom_violin() + theme_classic() + geom_hline(yintercept = 3) + ylim(c(2,6))

meta_lane1 <- meta_lane1[meta_lane1$passed_filters > 10^3, ] # visually identify threshold from violin plot
meta_lane2 <- meta_lane2[meta_lane2$passed_filters > 10^3, ]

## create fragment object and peak martrix
frag_lane1 <- CreateFragmentObject(path = frag_path[1], cells = rownames(meta_lane1))
peakmat_lane1 <- FeatureMatrix(fragments = frag_lane1, features = peaks, cells = rownames(meta_lane1), process_n = 5e3)

frag_lane2 <- CreateFragmentObject(path = frag_path[2], cells = rownames(meta_lane2))
peakmat_lane2 <- FeatureMatrix(fragments = frag_lane2, features = peaks, cells = rownames(meta_lane2), process_n = 5e3)

## read in rna data, tranlsate to ATAC cellIDs, subset on cells in peakmat, remove cells with < 100 UMI, remove genes w/ < 5 UMI, unite peakmat/rna cell IDs
rna_lane1 <- Read10X_h5(file = rna_path[1], name = 'group')
colnames(rna_lane1) <- paste0(atac_wlist, '-1')
rna_lane1 <- rna_lane1[ , colnames(peakmat_lane1)] 
rna_lane1 <- rna_lane1[ , which(colSums(rna_lane1) >= 100)]

rna_lane2 <- h5read(file = rna_path[2], name = 'group')
colnames(rna_lane2) <- paste0(atac_wlist, '-1')
rna_lane2 <- rna_lane2[ , colnames(peakmat_lane2)] 
rna_lane2 <- rna_lane2[ , which(colSums(rna_lane2) >= 100)]

temp <- cbind(rna_lane1, rna_lane2)
ind <- which(rowSums(temp) >= 5)
rna_lane1 <- rna_lane1[ind, ]
rna_lane2 <- rna_lane2[ind, ]
peakmat_lane1 <- peakmat_lane1[ , colnames(rna_lane1)]
peakmat_lane2 <- peakmat_lane2[ , colnames(rna_lane2)]
frag_lane1 <- subset(frag_lane1, cells = colnames(peakmat_lane1))
frag_lane2 <- subset(frag_lane2, cells = colnames(peakmat_lane2))

## create and merge seurat objects, preprocess
seu_lane1 <- CreateSeuratObject(counts = rna_lane1, assay = "RNA")
seu_lane1[["ATAC"]] <- CreateChromatinAssay(peakmat_lane1, fragments = frag_lane1)
seu_lane2 <- CreateSeuratObject(counts = rna_lane2, assay = "RNA")
seu_lane2[["ATAC"]] <- CreateChromatinAssay(peakmat_lane2, fragments = frag_lane2)
seu_merged <- merge(x = seu_lane1, y = list(seu_lane2), add.cell.ids = c("lane1", "lane2"))
# clear space: rm(seu_lane1, seu_lane2, peakmat_lane1, peakmat_lane2, rna_lane1, rna_lane2, temp, rna, meta_lane1, meta_lane2, frag_lane1, frag_lane2, peakwidths)

plan('default')
DefaultAssay(seu_merged) <- "RNA"
seu_merged <- NormalizeData(seu_merged)
seu_merged <- ScaleData(seu_merged)
seu_merged <- FindVariableFeatures(seu_merged)
seu_merged <- RunPCA(seu_merged)

DefaultAssay(seu_merged) <- "ATAC"
seu_merged <- FindTopFeatures(seu_merged, min.cutoff = 5)
seu_merged <- RunTFIDF(seu_merged)
seu_merged <- RunSVD(seu_merged)
seu_merged <- FindMultiModalNeighbors(seu_merged, reduction.list = list("pca", "lsi"), dims.list = list(1:30, 2:40), modality.weight.name = "RNA.weight", verbose = TRUE)
seu_merged <- RunUMAP(seu_merged, nn.name = "weighted.nn", assay = "RNA", verbose = TRUE)
seu_merged <- FindClusters(seu_merged, graph.name = "wsnn", algorithm = 3, resolution = 1, verbose = FALSE)
save(seu_merged, file='seu_merged.Robj')

seu_merged@meta.data[,'lane'] <- 1L
seu_merged@meta.data$lane[grep('lane2', colnames(seu_merged))] <- 2

## demultiplex using multi-seq fastqs -- replace 'XXX' with MULTI-seq BCs used in experiment (i.e., position in bar.ref)
cells <- rownames(seu_merged@meta.data)[which(seu_merged@meta.data$lane == 1)]
cells <- gsub(pattern = 'lane1_|-1', replacement = '', x = cells)
cells <- rna_wlist[match(cells, atac_wlist)]
read_table <- readTags(dir = ".", fastq.A = multi_r1_path[1], fastq.B = multi_r2_path[1], barcode.type = "MULTIseq", assay = "RNA", filter.cells = cells)
tag_mtx_lane1 <- alignTags(read_table, bar.ref)
rownames(tag_mtx_lane1) <- atac_wlist[match(rownames(tag_mtx_lane1), rna_wlist)]
rownames(tag_mtx_lane1) <- paste0('lane1_',rownames(tag_mtx_lane1),'-1'); save(tag_mtx_lane1, file='tag_mtx_lane1.Robj')
res_lane1 <- demultiplexTags(tag_mtx_lane1[,c('XXX')], plot.path = "./", plot.name = "lane1_", plot.diagnostics = F)

cells <- rownames(seu_merged@meta.data)[which(seu_merged@meta.data$lane == 2)]
cells <- gsub(pattern = 'lane2_|-1', replacement = '', x = cells)
cells <- rna_wlist[match(cells, atac_wlist)]
read_table <- readTags(dir = ".", fastq.A = multi_r1_path[2], fastq.B = multi_r2_path[2], barcode.type = "MULTIseq", assay = "RNA", filter.cells = cells)
tag_mtx_lane2 <- alignTags(read_table, bar.ref)
rownames(tag_mtx_lane2) <- atac_wlist[match(rownames(tag_mtx_lane2), rna_wlist)]
rownames(tag_mtx_lane2) <- paste0('lane2_',rownames(tag_mtx_lane2),'-1'); save(tag_mtx_lane2, file='tag_mtx_lane2.Robj')
res_lane2 <- demultiplexTags(tag_mtx_lane2[,c('XXX')], plot.path = "./", plot.name = "lane2_", plot.diagnostics = F); save(res_lane2, file='res_lane2.Robj')

## add classifications, parse, reprocess
seu_merged@meta.data[,'multi'] <- 0L
seu_merged@meta.data[names(res_lane1$final_assign),'multi'] <- res_lane1$final_assign
seu_merged@meta.data[names(res_lane2$final_assign),'multi'] <- res_lane2$final_assign

bad_clus <- c('XXX') # manually identify clusters with low RNA UMIs, ATAC UMIs, doublet/negative enrichment
bad_cells <- unique(c(rownames(seu_merged@meta.data)[which(seu_merged@meta.data$multi %in% c('multiplet','negative'))],
                      rownames(seu_merged@meta.data)[which(seu_merged@active.ident %in% bad_clus)]))
good_cells <- colnames(seu_merged)[which(colnames(seu_merged) %ni% bad_cells)]
seu_merged <- subset(seu_merged, cells=good_cells)

DefaultAssay(seu_merged) <- "RNA"
seu_merged <- NormalizeData(seu_merged)
seu_merged <- ScaleData(seu_merged); gc()
seu_merged <- FindVariableFeatures(seu_merged)
seu_merged <- RunPCA(seu_merged)

DefaultAssay(seu_merged) <- "ATAC"
seu_merged <- FindTopFeatures(seu_merged, min.cutoff = 5)
seu_merged <- RunTFIDF(seu_merged)
seu_merged <- RunSVD(seu_merged)
seu_merged <- FindMultiModalNeighbors(seu_merged, reduction.list = list("pca", "lsi"), dims.list = list(1:30, 2:40), modality.weight.name = "RNA.weight", verbose = TRUE)
seu_merged <- RunUMAP(seu_merged, nn.name = "weighted.nn", assay = "RNA", verbose = TRUE)
seu_merged <- FindClusters(seu_merged, graph.name = "wsnn", algorithm = 3, resolution = 0.25, verbose = FALSE)
save(seu_merged, file='seu_merged.Robj')