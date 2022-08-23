DIALOGUE_SeuratExample<-function(){
  set.seed(123)
  # Simulating data based on Seurat object of PBMCs
  obj<-DIALOGUE_example.initialize(installation.flag = F)
  obj@meta.data$samples <- sample(c(paste0("sample",1:25)), size = ncol(pbmc3k), replace =TRUE)
  obj@meta.data$cell.subtypes<-obj@meta.data$seurat_annotations
  
  b1<-is.element(obj@meta.data$cell.subtypes,c("Naive CD4 T","CD14+ Mono"))
  b2<-is.element(obj@meta.data$cell.subtypes,c("Memory CD4 T","FCGR3A+ Mono"))
  b3<-is.element(obj@meta.data$cell.subtypes,c("CD8 T","DC"))
  obj@meta.data$samples[b1] <- sample(c(paste0("sample",1:10)), size = sum(b1), replace =TRUE)
  obj@meta.data$samples[b2] <- sample(c(paste0("sample",11:20)), size = sum(b2), replace =TRUE)
  obj@meta.data$samples[b3] <- sample(c(paste0("sample",21:25)), size = sum(b3), replace =TRUE)
  
  r1<- DIALOGUE_make.cell.type.seurat(obj, cell.subtypes = c("Naive CD4 T","Memory CD4 T","CD8 T"), name = "Tcell")
  r2<- DIALOGUE_make.cell.type.seurat(obj, cell.subtypes = c("CD14+ Mono","FCGR3A+ Mono","DC"), name = "Myeloid")
  
  rA<- list(Tcell = r1,Myeloid = r2)
  
  R <- DIALOGUE.run(rA = rA, # list of cell.type objects
                    main = "ToyExample",
                    k = 3, # number of MCPs to identify
                    results.dir = results.dir,
                    spatial.flag = F,
                    conf = "cellQ")# any potential confounders DIALOGUE needs to account for; default is "cellQ" 
  
  par(mfrow=c(1,2),oma = c(8, 1, 0, 5),xpd = T)
  # MCP1 marks CD8 T cells and DCs
  boxplot(R$scores$Tcell$MCP1~R$scores$Tcell$cell.subtypes,xlab = "",ylab = "MCP1",las=2)
  boxplot(R$scores$Myeloid$MCP1~R$scores$Myeloid$cell.subtypes,xlab = "",ylab = "MCP1",las=2)
  # MCP2 marks memory CD4 T cells and FCGR3A+ Mono
  boxplot(-R$scores$Tcell$MCP2~R$scores$Tcell$cell.subtypes,xlab = "",ylab = "MCP2",las=2)
  boxplot(-R$scores$Myeloid$MCP2~R$scores$Myeloid$cell.subtypes,xlab = "",ylab = "MCP2",las=2)
  # MCP3 marks naive CD4 T cells and CD14+ Mono
  boxplot(-R$scores$Tcell$MCP3~R$scores$Tcell$cell.subtypes,xlab = "",ylab = "MCP3",las=2)
  boxplot(-R$scores$Myeloid$MCP3~R$scores$Myeloid$cell.subtypes,xlab = "",ylab = "MCP3",las=2)
  
}

DIALOGUE_example.initialize<-function(installation.flag){
  if(installation.flag){
    devtools::install_github(repo = "https://github.com/livnatje/DIALOGUE")
    devtools::install_github('satijalab/seurat-data')
  }
  library(DIALOGUE)
  library(SeuratData)
  library(Seurat)
  InstallData("pbmc3k")
  data("pbmc3k")
  return(pbmc3k)
}

DIALOGUE_make.cell.type.seurat<- function(obj, cell.subtypes, name){
  # Given a Seurat object with "cell.subtypes" and "samples" information provided in the meta.data field
  
  if(is.null(obj@meta.data$cell.subtypes)){return("Cell subtypes information is missing")}
  if(is.null(obj@meta.data$samples)){return("Sample information is missing")}
  
  sub_obj<- subset(x = obj,cells = rownames(obj@meta.data)[is.element(obj@meta.data$cell.subtypes,cell.subtypes)])
  sub_obj<- NormalizeData(sub_obj,normalization.method = "LogNormalize", scale.factor = 100000)
  sub_obj<-FindVariableFeatures(sub_obj,selection.method = "vst", nfeatures = 2000)
  sub_obj<-ScaleData(sub_obj)
  sub_obj<-RunPCA(sub_obj,npcs = 30) 
  
  tpm<- as.matrix(sub_obj@assays$RNA@data)
  
  X<- sub_obj@reductions$pca@cell.embeddings
  assertthat::are_equal(rownames(X),colnames(tpm))
  n_cells<- ncol(tpm)
  samples<- sub_obj@meta.data$samples
  n_count<- as.vector(scale(log(sub_obj@meta.data$nCount_RNA)))
  metadata<- data.frame(nFeatures = sub_obj@meta.data$nFeature_RNA,
                        samples = samples,
                        cell.subtypes = as.character(sub_obj@meta.data$cell.subtypes))
  r<-make.cell.type(name = name ,tpm,samples,X,metadata, cellQ = n_count)
  return(r)
}



