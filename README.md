# Single-cell analysis of human Dropseq based scRNAseq and study of APA dynamic
Isoform quantification of scRNA-seq DropSeq neuronal differentiation  using SCALPEL


## A. Preprocessing of gene-level based scRNA-seq dataset

### 1. File importation
```
Author: Franz AKE
Barcelona, P-CRMC lab

#imports Libs
#------------
library(data.table)
library(dplyr)
library(Seurat)
library(ggplot2)
library(ggthemes)

#parameters
samtools.exec = "samtools"

#import files
#------------
#paths
DGE.paths = c(
    "iPSC3" = "/home/fake/FAKE/PROJECTS/NEURODIFF/datas/processed_datas/iPSC3.counts.txt",
    "NPC3" = "/home/fake/FAKE/PROJECTS/NEURODIFF/datas/processed_datas/NPC3.counts.txt"
)
bamfiles.in = c(
    "iPSC3" = "/home/fake/FAKE/PROJECTS/NEURODIFF/datas/raw_datas/iPSC3_proc.bam",
    "NPC3" = "/home/fake/FAKE/PROJECTS/NEURODIFF/datas/raw_datas/NPC3_proc.bam"
)

#reading
#-------
myObjs = lapply(1:length(DGE.paths), function(x){
    path = DGE.paths[x]
    
    #read counts (a)
   	#>>>>>>>>>>>
	counts = read.table(path, header = T, sep="\t", row.names = 1)
	
    #Create Seurat object (b) ~ minimum gene UMI = 4
	#>>>>>>>>>>>>>>>>>>>>>
    return(CreateSeuratObject(counts = read.table(path, header = T, sep="\t", row.names = 1), 
    project = names(DGE.paths)[x], min.cells=4, min.features=3))

})

#merging
#-------
myObj = merge(myObjs[[1]], myObjs[[2]])
```

### 2. Quality filtering

```
#Quality control
#---------------
print("quality control...")
myObj[["percent.mt"]] <- PercentageFeatureSet(myObj, pattern = "^MT-")
VlnPlot(myObj, features = c("nCount_RNA","nFeature_RNA","percent.mt"))
VlnPlot(myObj %>% subset(nCount_RNA < 20000 & nCount_RNA > 300 & nFeature_RNA > 300 & nFeature_RNA < 7000 & percent.mt <5),
        features = c("nCount_RNA","nFeature_RNA","percent.mt"))
myObj.1 = subset(myObj, nCount_RNA < 20000 & nCount_RNA > 300 & nFeature_RNA > 300 & nFeature_RNA < 7000 & percent.mt <5)

> VlnPlot(myObj, features = c("nCount_RNA","nFeature_RNA"))
```

***Figure1:** Quality filtering of samples*
![Quality filtering](https://data.cyverse.org/dav-anon/iplant/home/franzx5/scRNAseq_neuronal_differentiation/docs/Figure1.png)


### 3. Dimensionality Reduction

```
#Dimensionality Reduction
#------------------------
print("Dimensionality red...")
myObj.1 = NormalizeData(myObj.1, verbose = F, normalization.method = "LogNormalize")
myObj.1 = FindVariableFeatures(myObj.1, verbose = F)
myObj.1 = ScaleData(myObj.1, verbose = F,vars.to.regress = c("nCount_RNA","nFeature_RNA","percent.mt"))
myObj.1 = RunPCA(myObj.1, verbose = F)
ElbowPlot(myObj.1, ndims = 50)
nb.pcs = 9
myObj.1 = RunUMAP(myObj.1, dims = 1:nb.pcs, verbose = F)
> DimPlot(myObj.1, cols = RColorBrewer::brewer.pal(name = "Set1", n=2), pt.size = 0.7)
```
***Figure2:** UMAP plot of neuronal samples*
![UMAP plot of neuronal samples](https://data.cyverse.org/dav-anon/iplant/home/franzx5/scRNAseq_neuronal_differentiation/docs/Figure2.png)


### 4. Samples preprocessing

```
#Generate BAM file And Barcodes associated
#-----------------------------------------
print("BAM filtering..")
myObj.1$Barcodes = as.character(stringr::str_replace(colnames(myObj.1), pattern = "\\_", replacement = "-"))
Idents(myObj.1) = myObj.1$Barcodes

#create Barcode file(a)
lapply(1:length(bamfiles.in), function(x){
    sample = names(bamfiles.in)[x]
    print(sample)
    filter(myObj.1@meta.data, orig.ident == sample) %>%
        dplyr::select(Barcodes) %>% mutate(Barcodes = stringr::str_replace(Barcodes, "\\.","\\-")) %>%
        fwrite(paste0(sample, ".barcodes.txt"), sep = "\t", row.names = F, col.names = F)
    filter(myObj.1@meta.data, orig.ident == sample) %>%
        dplyr::select(Barcodes, orig.ident) %>%
        fwrite(paste0(sample, ".clusters.txt"), sep = "\t", row.names = F, col.names = F)

    #subprocessing(b)
    FILTER.EXP = paste0(samtools.exec, " view -b -@ 10 -D XC:", sample, ".barcodes.txt ", bamfiles.in[x], " > ", sample, ".bam")
    print(FILTER.EXP)
    INDEX.EXP = paste0(samtools.exec, " index ", sample, ".bam")
    print(INDEX.EXP)
    print("filtering BAM...")
    system(FILTER.EXP)
    print("indexing...")
    system(INDEX.EXP)
    
    #write DGE (c)
    print("writing DGE...")
    temp.obj = subset(myObj.1, orig.ident == sample)@assays$RNA@counts %>% data.frame()
    colnames(temp.obj) = stringr::str_replace(colnames(temp.obj), pattern="\\.", replacement="\\-")
    temp.obj %>% data.table(keep.rownames = T) %>%
        rename("GENE"="rn") %>%
        fwrite(file = paste0(sample, ".counts.txt"), sep = "\t", row.names = F, col.names = T)
})


#save processed seurat object (d)
saveRDS(myObj.1, file="processedObj.RDS")

#write meadata objects (e)
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

#stats about processed object
#----------------------------
An object of class Seurat 
19103 features across 2535 samples within 1 assay 
Active assay: RNA (19103 features, 2000 variable features)
 3 layers present: counts, data, scale.data
 2 dimensional reductions calculated: pca, umap

#number.cells = 2535
#nb.genes = 19.103 

```

Once processed, the current files are generated and can be used as input for SCALPEL execution

```
ips_npc.bam  ips_npc.bam.bai  ips_npc.barcodes.txt  ips_npc.clusters.txt  ips_npc.counts.txt  processedObj.RDS
```

## B. APA dynamic analysis using SCALPEL

All the information about SCALPEL installation and execution can be find on [SCALPEL](https://github.com/p-CMRC-LAB/SCALPEL) GitHub.

Following SCALPEL execution, the current folder **_results/final_results_** is generated

### 1. File importation

Using the iDGE counts matrix generated by SCALPEL, we can perform a downstream single-cell analysis by creating a Seurat object :

```
#sourcing to SCALPEL implemented functions scripts
#-------------------------------------------------
source("~/SCALPEL/src/scalpel_library.R")

#Opening (1)
#---------
scalpel.counts = fread("~/raw_neurofiles/results/final_results/ips_npc/ips_npc_APADGE.txt") %>%
  data.frame()
rownames(scalpel.counts) = scalpel.counts$gene_transcript
scalpel.counts = scalpel.counts[,2:ncol(scalpel.counts)]
colnames(scalpel.counts) = str_replace(colnames(scalpel.counts), pattern = "\\.", "\\-")

#get statistic about the number of genes and isoforms following SCALPEL
# execution
scalpel.counts.tab = data.table(rownames(scalpel.counts)) %>%
  tidyr::separate(V1, into=c("gene_name","transcript_name"), sep="\\*\\*\\*") %>%
  dplyr::filter(gene_name %in% rownames(in.obj)) %>%
  data.table()
summarise(scalpel.counts.tab, nb.genes=n_distinct(gene_name), nb.isofs = n_distinct(transcript_name))

#stats
#nb.cells = 2.535
#nb.genes = 18.797
#nb.isoforms = 78.268

#Create Seurat object (2)
#--------------------
scalpel.seurat = CreateSeuratObject(scalpel.counts, project = "ips_npc", min.cells=4, min.features=3)
scalpel.seurat$Barcodes = colnames(scalpel.seurat)
scalpel.seurat$Cell_Cluster = left_join(scalpel.seurat@meta.data, raw_cell.infos)$Cell_Cluster
scalpel.seurat 
```


### 2. Quality filtering
```
# ... (fixing issue arising from the fact that the samples are merged and min.cells=4 in CreateSeuratObject 
# fails to filter in each sample separately ...

#Discarding of isoforms expressed in less than 4 cells for each sample 
#---------------------------------------------------------------------
#get iPSC/NPCs (a)
scalpel.objs = lapply(c("iPSC3","NPC3"), function(x){
  print(x)
  sub.obj = subset(scalpel.seurat, Cell_Cluster==x)
  m <- sub.obj[["RNA"]]$counts
  t <- table(m@i)
  t.selected <- t[t >= 4]
  isoform.selected <- rownames(m)[as.numeric(names(t.selected)) + 1]
  myobj <- subset(scalpel.seurat, features = isoform.selected)
  return(data.table(gene_tr = rownames(myobj)))
})
scalpel.objs[[1]]$sample1 = "iPSC3"
scalpel.objs[[2]]$sample2 = "NPC3"

#get list of isoforms to keep (b)
data.table(gene_tr = rownames(scalpel.seurat)) %>%
  left_join(scalpel.objs[[1]]) %>%
  left_join(scalpel.objs[[2]]) %>%
  dplyr::filter((!(is.na(sample1)) & (!(is.na(sample2))))) -> isof.tokeep

#filtering (c)
scalpel.seurat = subset(scalpel.seurat, features =isof.tokeep$gene_tr)
scalpel.seurat

#stats following filtering
#-------------------------
> scalpel.seurat
An object of class Seurat 
49500 features across 2535 samples within 1 assay 
Active assay: RNA (49500 features, 4000 variable features)
 3 layers present: counts, data, scale.data
 2 dimensional reductions calculated: pca, umap

#nb.cells = 2535
#nb.genes = 13.431
#nb.isoforms = 49.500
```

### 3. Dimensionality reduction

The isoform quantification allows to perform dimensionality reduction analysis and visualize the sample using variable isoforms information :
```
#Dimensionality reduction analysis (1)
#---------------------------------
scalpel.seurat = NormalizeData(scalpel.seurat)
scalpel.seurat = FindVariableFeatures(scalpel.seurat, nfeatures = 4000)
scalpel.seurat = ScaleData(scalpel.seurat, verbose = F)
scalpel.seurat = RunPCA(scalpel.seurat, verbose = F)
ElbowPlot(scalpel.seurat, ndims = 50)
nb.pcs = 9
scalpel.seurat = RunUMAP(scalpel.seurat, dims = 1:nb.pcs, verbose = F)

#Visualization (2)
#-------------
Idents(scalpel.seurat) = scalpel.seurat$Cell_Cluster
DimPlot(in.obj, group.by = "Cell_Cluster", label = F, label.size = 0, pt.size = 1) + 
  ggtitle("gene based quantification") -> p1
DimPlot(scalpel.seurat, group.by = "Cell_Cluster", label = F, label.size = 10, pt.size = 1) +
  ggtitle("isoforms based quantification") -> p2
  
> p1 + p2
```

***Figure3:** DGE & iDGE based UMAP plot of neuronal samples*
![DGE & iDGE counts UMAP plots](https://data.cyverse.org/dav-anon/iplant/home/franzx5/scRNAseq_neuronal_differentiation/docs/Figure3.png)



### 3. Differential isoform usage in cell-types

Here, we perform a differential analysis to get differential isoform usage (DIU) genes in the iPS and NP cells. The analysis is only performed on genes with at least 2 isoforms :

```
#Differential isoform usage
#--------------------------
scalpel.dtu.sig = Find_isoforms(scalpel.seurat, pval_adjusted = 0.05, condition = "Cell_Cluster", threshold_tr_abundance = 0.1)
scalpel.dtu.sig$test = "iPSC3-NPC3"
scalpel.dtu.sig$iPSC3 = as.numeric(scalpel.dtu.sig$iPSC3) %>% round(2)
scalpel.dtu.sig$NPC3 = as.numeric(scalpel.dtu.sig$NPC3) %>% round(2)

#stats
> scalpel.dtu.sig %>%
  summarise(nb.genes = n_distinct(gene_name), nb.isoforms = n_distinct(transcript_name))

nb.DIUgenes = 2.095
nb.DIUisoforms = 6.208
```
Let's consider an example of DIU genes :

```
gene.in = "EIF1"
diu_trs = dplyr::filter(scalpel.dtu.sig, gene_name %in% gene.in)
> diu_trs
```
| iPS | NP | gene/isoform | gene_name | transcript_name | Pvalue | Pvalue.adjusted | test |
|--|--|--|--|--|--|--|--|
| 2566 | 1977 | EIF1***EIF1-201 | EIF1 | EIF1-201 |1.016049e-80 | 3.635588e-78 | iPS-NP |
| 5382 | 1950 | EIF1***EIF1-201 | EIF1 | EIF1-201 |1.016049e-80 | 3.635588e-78 | iPS-NP |

```
#Violin plot of the gene...
> VlnPlot(scalpel.seurat, features = trs.tab$gene)
```

***Figure4:** Violin plot of EIF1 isoforms expression in the samples*
![DIU genes](https://data.cyverse.org/dav-anon/iplant/home/franzx5/scRNAseq_neuronal_differentiation/docs/Figure4.png)


Let's visualize this gene reads coverage track on this genes :
```

#1 - first let's split the filtered BAM file by the defined samples
#------------------------------------------------------------------
filtered.BAM = "~/raw_neurofiles/final_results/ips_npc/ips_npc.filtered.bam"

# filtering (a)
bc.tags = dplyr::select(scalpel.seurat@meta.data, Barcodes, Cell_Cluster)
lapply(c("iPSC3","NPC3"), function(x){
  #current cluster
  print(x)
  #filter tab and get barcodes
  filter(bc.tags, Cell_Cluster==x) %>% dplyr::select(Barcodes) %>% fwrite(file = paste0(x,".txt"),sep = "\t")
  #samtools filtering
  filt.exp = paste0("/usr/local/Caskroom/mambaforge/base/envs/test/bin/samtools view -b -D XC:", 
                    x, ".txt ", filtered.BAM, " > ", x, ".bam")
  index.exp = paste0("/usr/local/Caskroom/mambaforge/base/envs/test/bin/samtools index ", x, ".bam")
  system(filt.exp)
  system(index.exp)
})

# Use the generated BAM file for each sample (b)
Cluster.bamfiles = c("iPS"="iPSC3.bam", "NPC"="NPC3.bam")

# import GTF reference file ...
GTF.path = "datas/gencode.v41.annotation.gtf"
gtf.gr = rtracklayer::import(GTF.path)


# Generate Coverage plot
# ----------------------
> CoverPlot(genome_gr = gtf.gr,
          gene_in = gene.in, genome_sp = "hg38",
          bamfiles = Cluster.bamfiles,
          samtools.bin = "/usr/local/Caskroom/mambaforge/base/envs/test/bin/samtools",
          transcripts_in = diu_trs$gene_tr, filter_trs = T)
```

***Figure5:** EIF1 coverage plot in iPS & NP samples*
![EIF1_coveragePlot](https://data.cyverse.org/dav-anon/iplant/home/franzx5/scRNAseq_neuronal_differentiation/docs/Figure5.png)


