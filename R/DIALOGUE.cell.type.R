#' DIALOGUE cell.type S4 class
#' @description An S4 class that represents a subset of cells (for example, a particular cell type).
#' 
#' @slot name cell type name;
#' @slot cells cell identifiers (n x 1);
#' @slot tpm gene expression matrix (m x n)
#' @slot genes genes (m x 1) represented in the [tpm] matrix;
#' @slot X features matrix (n x k), where n is the number of cells and k is the number of features, e.g., PCs, NMF components, tpm etc.
#' @slot samples the samples corresponding to the cells in [cells] (n x 1)
#' @seealso See \href{https://github.com/livnatje/DIALOGUE}{DIALOGUE GitHub page} for more details.
#' \code{\link{DIALOGUE.plot}}
#' @author Livnat Jerby-Arnon
#' @export

cell.type <- setClass(Class = "cell.type",
                      slots = c("name","cells","genes","cellQ",
                                "tpm","tpmAv","qcAv","zscores",
                                "X","samples",
                                "metadata",
                                "scores","scoresAv",
                                "tme","tme.qc","gene.pval",
                                "tme.OE",
                                "extra.scores",
                                "sig"))

#' make.cell.type
#'
#' This function generates a \linkS4class{cell.type} object for DIALOGUE.
#' @param name cell type name
#' @param tpm gene expression or any type of single-cell profiling (m x n)
#' @param X features matrix (n x k), e.g., PCs, NMF components, tpm etc.; these features will be used to identify the multicellular programs.
#' @param samples the sample of each cell (n x 1)
#' @param cellQ cell quality measures, e.g., number of reads/genes detected (n x 1)
#' @param metadata cell features (n x r)
#' @field cell.type a representation of a specific type of cells
#' @export

make.cell.type<-function(name,tpm,samples,X = NULL,metadata = NULL,
                         tpmAv = t(average.mat.rows(t(tpm),samples)),
                         cellQ){
  r<-cell.type(name = gsub("_","",name),
               cells = colnames(tpm),
               genes = rownames(tpm),
               cellQ = cellQ,
               tpm = tpm,
               tpmAv = tpmAv,
               qcAv = aggregate(x = cellQ,by = list(samples),FUN = mean),
               X = X,
               samples = samples,
               metadata = cbind.data.frame(cellQ = cellQ,metadata),
               extra.scores = list())
  print(paste("Cell type name:",r@name))
  if(!identical(r@cells,rownames(X))){
    print("Each row in X should include the original representation of the corresponding cell.")
    print("Error: Cell ids do not match the X input.")
    return("Error: Cell ids do not match the X input.")
  }
  rownames(r@qcAv)<-r@qcAv[,1]
  r@qcAv[,1]<-r@qcAv[,2]
  return(r)
}

cell.type.2.list<-function(r,idx){
  if(missing(idx)){
    idx<-slotNames(r)
  }
  r1<-lapply(idx,function(x) slot(r,name = x))
  names(r1)<-idx
  for(x in colnames(r@metadata)){r1[[x]]<-r@metadata[,x]}
  return(r1)
}
