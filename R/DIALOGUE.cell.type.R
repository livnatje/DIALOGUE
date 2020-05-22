#' DIALOGUE cell.type S4 class
#' @description An S4 class that represents a subset of cells (for example, a particular cell type).
#' 
#' @slot name cell type name;
#' @slot cells cell identifiers (1xn);
#' @slot tpm gene expression matrix (mxn)
#' @slot genes genes (1xm) represented in the [tpm] matrix;
#' @slot X features matrix (kxn), e.g., PCs, NMF components, tpm etc.
#' These features will be used to identify the multicellular programs.
#' @slot samples the samples corresponding to the cells in [cells] (1xn)
#' @seealso See \href{https://github.com/livnatje/DIALOGUE}{DIALOGUE GitHub page} for more details.
#' \code{\link{DIALOGUE.plot}}
#' @author Livnat Jerby-Arnon
cell.type <- setClass(Class = "cell.type",
                      slots = c("name","cells","genes",
                                "tpm","tpmAv","zscores",
                                "X","samples","cellQ",
                                "scores","scoresAv",
                                "tme","gene.pval","conf",
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
#' @export
#' @examples
#' r<-make.cell.type(name = "CD8.T.cells",tpm = tpm,samples = c("sample1,"sample1","sample2","sample2"),
#' cellQ = colSumes(tpm>0),X = PCs)
#' @field cell.type a representation of a specific type of cells

make.cell.type<-function(name,tpm,samples,cellQ,X = NULL,conf = NULL,
                         tpmAv = t(average.mat.rows(t(tpm),samples))){
  r<-cell.type(name = name,
               cells = colnames(tpm),
               genes = rownames(tpm),
               tpm = tpm,
               tpmAv = tpmAv,
               X = X,
               samples = samples,
               cellQ = cellQ,
               conf = conf,
               extra.scores = list())
  return(r)
}

cell.type.2.list<-function(r,idx){
  if(missing(idx)){
    idx<-slotNames(r)
  }
  r1<-lapply(idx,function(x) slot(r,name = x))
  names(r1)<-idx
  if(is.null(r@conf)){return(r1)}
  for(x in colnames(r@conf)){r1[[x]]<-r@conf[,x]}
  return(r1)
}