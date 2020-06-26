#' DIALOGUE cell.type S4 class
#' @description An S4 class that represents a subset of cells (for example, a particular cell type).
#' 
#' @slot name cell type name;
#' @slot cells cell identifiers (n x 1);
#' @slot tpm gene expression matrix (m x n)
#' @slot genes genes (m x 1) represented in the [tpm] matrix;
#' @slot X features matrix (n x k), e.g., PCs, NMF components, tpm etc.
#' These features will be used to identify the multicellular programs.
#' @slot samples the samples corresponding to the cells in [cells] (n x 1)
#' @seealso See \href{https://github.com/livnatje/DIALOGUE/wiki}{DIALOGUE's wiki} for more details.
#' \code{\link{DIALOGUE.plot}}
#' @author Livnat Jerby-Arnon
cell.type <- setClass(Class = "cell.type",
                      slots = c("name","cells","genes",
                                "tpm","tpmAv","zscores",
                                "X","samples",
                                "metadata",
                                "scores","scoresAv",
                                "tme","gene.pval",
                                "tme.OE",
                                "extra.scores",
                                "sig"))

#' make.cell.type
#'
#' This function generates a \linkS4class{cell.type} object for DIALOGUE.
#' @param name cell type name
#' @param tpm gene expression or any type of single-cell profiling (mxn)
#' @param X features matrix (kxn), e.g., PCs, NMF components, tpm etc.; these features will be used to identify the multicellular programs.
#' @param samples the sample of each cell (1xn)
#' @param cellQ cell quality measures, e.g., number of reads/genes detected (1xn)
#' @field cell.type a representation of a specific type of cells
#' @export

make.cell.type<-function(name,tpm,samples,X = NULL,metadata = NULL,
                         tpmAv = t(average.mat.rows(t(tpm),samples))){
  r<-cell.type(name = name,
               cells = colnames(tpm),
               genes = rownames(tpm),
               tpm = tpm,
               tpmAv = tpmAv,
               X = X,
               samples = samples,
               metadata = metadata,
               extra.scores = list())
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