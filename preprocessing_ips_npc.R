
#Author: Franz AKE
#PLASS lab
#Barcelona

#imports Libs
library(data.table)
library(dplyr)
library(Seurat)
library(ggplot2)
library(ggthemes)


#import Datasets (A)
#===============
print("import...")
DGE.paths = c(
    "iPSC3" = "/home/fake/FAKE/PROJECTS/NEURODIFF/datas/processed_datas/iPSC3.counts.txt",
    "NPC3" = "/home/fake/FAKE/PROJECTS/NEURODIFF/datas/processed_datas/NPC3.counts.txt"
)
bamfiles.in = c(
    "iPSC3" = "/home/fake/FAKE/PROJECTS/NEURODIFF/datas/raw_datas/iPSC3_proc.bam",
    "NPC3" = "/home/fake/FAKE/PROJECTS/NEURODIFF/datas/raw_datas/NPC3_proc.bam"
)
samtools.exec = "samtools"

#reading...
myObjs = lapply(1:length(DGE.paths), function(x){
    path = DGE.paths[x]
    print(path)
    return(CreateSeuratObject(counts = read.table(path, header = T, sep="\t", row.names = 1), project = names(DGE.paths)[x],
           min.cells=3, min.features=1))
})
#merging...
myObj = merge(myObjs[[1]], myObjs[[2]])


#Quality control (B)
#===============
print("quality control...")
myObj[["percent.mt"]] <- PercentageFeatureSet(myObj, pattern = "^MT-")
VlnPlot(myObj, features = c("nCount_RNA","nFeature_RNA","percent.mt"))
VlnPlot(myObj %>% subset(nCount_RNA < 20000 & nCount_RNA > 300 & nFeature_RNA > 300 & nFeature_RNA < 7000 & percent.mt <5),
        features = c("nCount_RNA","nFeature_RNA","percent.mt"))
myObj.1 = subset(myObj, nCount_RNA < 20000 & nCount_RNA > 300 & nFeature_RNA > 300 & nFeature_RNA < 7000 & percent.mt <5)


#Dimensionality Reduction (C)
#========================
print("Dimensionality red...")
myObj.1 = NormalizeData(myObj.1, verbose = F, normalization.method = "LogNormalize")
myObj.1 = FindVariableFeatures(myObj.1, verbose = F)
myObj.1 = ScaleData(myObj.1, verbose = F,vars.to.regress = c("nCount_RNA","nFeature_RNA","percent.mt"))
myObj.1 = RunPCA(myObj.1, verbose = F)
ElbowPlot(myObj.1, ndims = 50)
nb.pcs = 9
myObj.1 = RunUMAP(myObj.1, dims = 1:nb.pcs, verbose = F)
DimPlot(myObj.1, cols = RColorBrewer::brewer.pal(name = "Set1", n=2), pt.size = 0.7)


#Generate BAM file And Barcodes associated (D)
#=========================================
print("BAM filtering..")
myObj.1$Barcodes = as.character(stringr::str_replace(colnames(myObj.1), pattern = "\\_", replacement = "-"))
Idents(myObj.1) = myObj.1$Barcodes

#create Barcode file(1)
lapply(1:length(bamfiles.in), function(x){
    sample = names(bamfiles.in)[x]
    print(sample)
    filter(myObj.1@meta.data, orig.ident == sample) %>%
        dplyr::select(Barcodes) %>% mutate(Barcodes = stringr::str_replace(Barcodes, "\\.","\\-")) %>%
        fwrite(paste0(sample, ".barcodes.txt"), sep = "\t", row.names = F, col.names = F)
    filter(myObj.1@meta.data, orig.ident == sample) %>%
        dplyr::select(Barcodes, orig.ident) %>%
        fwrite(paste0(sample, ".clusters.txt"), sep = "\t", row.names = F, col.names = F)
    FILTER.EXP = paste0(samtools.exec, " view -b -@ 10 -D XC:", sample, ".barcodes.txt ", bamfiles.in[x], " > ", sample, ".bam")
    print(FILTER.EXP)
    INDEX.EXP = paste0(samtools.exec, " index ", sample, ".bam")
    print(INDEX.EXP)
    print("filtering BAM...")
    #system(FILTER.EXP)
    print("indexing...")
    #system(INDEX.EXP)
    #write DGE
    print("writing DGE...")
    temp.obj = subset(myObj.1, orig.ident == sample)@assays$RNA@counts %>% data.frame()
    colnames(temp.obj) = stringr::str_replace(colnames(temp.obj), pattern="\\.", replacement="\\-")
    temp.obj %>% data.table(keep.rownames = T) %>%
        rename("GENE"="rn") %>%
        fwrite(file = paste0(sample, ".counts.txt"), sep = "\t", row.names = F, col.names = T)
})


#save processed Objects
saveRDS(myObj.1, file="processedObj.RDS")

#merge object
message("merging BAM files...")
system(paste0(samtools.exec, " merge -@ 4 -o ips_npc.bam iPSC3.bam NPC3.bam"))

message("DGE table...")
temp.obj = myObj.1@assays$RNA@counts %>% data.frame()
colnames(temp.obj) = stringr::str_replace(colnames(temp.obj), pattern="\\.", replacement="\\-")
temp.obj %>% data.table(keep.rownames = T) %>% rename("GENE"="rn") %>% 
    fwrite(file = paste0("ips_npc", ".counts.txt"), sep = "\t", row.names = F, col.names = T)

message("write clusters..")
myObj.1[[]] %>% mutate(Barcodes=stringr::str_replace(rownames(.),"\\.","\\-")) %>% distinct(Barcodes,orig.ident) %>% fwrite("ips_npc.clusters.txt",sep="\t",col.names=F)

message("write barcodes...")
myObj.1[[]] %>% mutate(Barcodes=stringr::str_replace(rownames(.),"\\.","\\-")) %>% distinct(Barcodes) %>% fwrite("ips_npc.barcodes.txt",sep="\t",col.names=F)









