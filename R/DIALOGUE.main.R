#' DIALOGUE.run
#'
#' DIALOGUE is a dimensionality reduction approach that uses cross-cell-type
#' associations to identify multicellular programs and map the cell transcriptome
#' as a function of its environment. 
#' 
#' @param rA list of \linkS4class{cell.type} objects.
#' @param k the number of multicellular programs to identify;
#' @param results.dir path to the results directory, where the output will be saved;
#' @param plot.flag if TRUE then \code{\link{DIALOGUE.plot}} will be called to plot the results; default is FALSE;
#' @param abn.c the minimal number of cells that a sample should have to be considered for MCP detection;
#' @param spatial.flag should be TRUE if working with spatial data with small niches, and TRUE if working with single cell data or larger tissue microenvironment niches. The default value is FALSE.
#' @param add.effects whether to add additional covariates to the multilevel model;
#' if TRUE then the covariates will extracted from the [conf] slot in the \linkS4class{cell.type} objects;
#' default is FALSE.
#' 
#' @return A list with the following components:
#' @return sig -  the multicellular programs, given as a list of signatures;
#' @return scores - the multicellular programs' scores in each cell;
#' @return gene.pval - the cross-cell-type p-values of each program;
#' @return pref - the correlation (R) and association (mixed-effects p-value) between the cell-type-sepcific
#' components of each multicellular programs.
#' @seealso See \href{https://github.com/livnatje/DIALOGUE}{DIALOGUE GitHub page} for more details.
#' \code{\link{DIALOGUE.plot}}
#' @author Livnat Jerby-Arnon
#' @export
#' 

DIALOGUE.run<-function(rA,main,k = 3,results.dir = getwd(),plot.flag = T,pheno = NULL,
                       PMD2 = F,conf = "cellQ",covar = c("cellQ","tme.qc"),n.genes = 200,
                       averaging.function = colMedians,p.anova = 0.05,specific.pair = NULL,
                       parallel.vs = F,center.flag = T,extra.sparse = F, bypass.emp = F, abn.c = 15, spatial.flag = F){
  full.version <- F
  names(rA)<-laply(rA,function(r) r@name)
  if(!grepl("\\/$", results.dir)){results.dir <- paste(results.dir, "/", sep = "")}
  R<-DIALOGUE1(rA = rA,k = k,main = main,
               results.dir = results.dir, PMD2 = PMD2,covar = covar,conf = conf,
               n.genes = n.genes,averaging.function = averaging.function,extra.sparse = extra.sparse,
               p.anova = p.anova,specific.pair = specific.pair,center.flag = center.flag,
               bypass.emp = bypass.emp,abn.c = abn.c,spatial.flag = spatial.flag)
  if(R$message=="No programs"){return(R)}
  if(!is.null(specific.pair)){
    main<-paste(main,paste(specific.pair,collapse = "."),sep = "_")
    rA<-rA[specific.pair]
  }
  R<-DIALOGUE2(rA = rA,main = main,results.dir = results.dir,parallel.vs = parallel.vs)
  R<-DIALOGUE3(rA = rA,main = main,results.dir = results.dir,full.version = full.version,pheno = pheno)
  if(plot.flag){
    DIALOGUE.plot(R,results.dir = results.dir,pheno = pheno)
  }
  return(R)
}

DIALOGUE.pheno<-function(R,pheno =  "clin.status",cca.flag = F,rA,frm,selected.samples = NULL){
  k<-R$k["DIALOGUE"]
  if(missing(frm)){
    frm<-gsub("+tme.qc","",R$frm,fixed = T)
    frm<-gsub("+ tme.qc","",frm,fixed = T)
    frm<-gsub(paste0("+",pheno),"",frm,fixed = T)
    frm<-gsub(paste0("+ ",pheno),"",frm,fixed = T)
  }
  
  f<-function(X){
    covar<-setdiff(R$covar,"tme.qc")
    r1<-lapply(covar, function(x) X[,x])
    names(r1)<-covar
    r1$pheno<-X[,pheno]
    if(!is.logical(r1$pheno)){r1$pheno<-r1$pheno==sort(unique(r1$pheno))[1]}
    r1<-c(r1,list(scores = as.matrix(X[,1:k]),samples = X$samples))
    if(any(is.na(r1$pheno))){
      print(paste("Identified",sum(is.na(r1$pheno)),"cells with no phenotype."))
      r1<-set.list(r1,!is.na(r1$pheno))
    }
    if(!is.null(selected.samples)){
      r1<-set.list(r1,is.element(r1$samples,selected.samples))
    }
    z<-apply.formula.HLM(r = r1,Y = r1$scores, X = r1$pheno,
                         MARGIN = 2,formula = frm)[,1]
    return(z)
  }
  if(cca.flag){
    if(is.null(R$metadata)){R<-DLG.add.metadata(R,rA = rA)}
    R$scores<-lapply(R$cell.types, function(x) cbind.data.frame(R$cca.scores[[x]],R$metadata[[x]]))
    names(R$scores)<-R$cell.types
  }
  Z<-laply(R$scores,f)
  X<-NULL;for(x in R$scores){X<-rbind(X,x)}
  R$covar<-c(R$covar,"cell.type")
  colnames(Z)<-colnames(R$scores[[1]])[grepl("MCP",colnames(R$scores[[1]]))]
  Z<-rbind(Z,f(X))
  rownames(Z)<-c(names(R$scores),"All")
  return(Z)
}

DIALOGUE1<-function(rA,k = 5,main,results.dir = "~/Desktop/DIALOGUE.results/",conf = "cellQ",
                    covar = c("cellQ","tme.qc","sex","pathology"),n.genes = 200,PMD2 = F,extra.sparse = F,
                    averaging.function = colMeans,p.anova = 0.05,specific.pair = NULL,center.flag = F,
                    seed1 = 1234,bypass.emp = F,abn.c = 15,spatial.flag = F){
  
  print("#************DIALOGUE Step I: PMD ************#")
  dir.create(results.dir, show.warings=FALSE, recursive=TRUE)
  if(!grepl("\\/$", results.dir)){results.dir <- paste(results.dir, "/", sep = "")}
  X<-lapply(rA, function(r){
    X1<-average.mat.rows(r@X,r@samples,f = averaging.function)
    if(spatial.flag){return(X1)}
    b<-get.abundant(r@samples,abn.c = abn.c,boolean.flag = T)
    p<-p.adjust(apply.anova(X = r@X[b,],y = r@samples[b],MARGIN = 2),method = "BH")
    print(paste0(r@name,": Removing ",sum(p>p.anova)," of ",length(p)," features."))
    X1<-X1[,names(p)[p<p.anova]]
    return(X1)
  })
  cell.types<-names(rA)
  names(X)<-cell.types
  n1<-length(cell.types)
  
  # Finding shared samples
  samples<-unlist(lapply(cell.types, function(x) rownames(X[[x]])))
  samplesU<-get.abundant(samples,n1)
  if(length(samplesU)<5){
    return("Error: Cannot run DIALOGUE with less than 5 samples.")
  }
  
  # Centering and scalling (optional)
  f<-function(X1){
    if(center.flag){
      X1<-center.matrix(X1,dim = 2,sd.flag = T)
      X1<-cap.mat(X1,cap = 0.01,MARGIN = 2)
    }
    X1<-X1[samplesU,]
    return(X1)
  }
  X<-lapply(X, f)
  k1<-ncol(X[[1]])
  
  if(is.null(specific.pair)){
    out<-DIALOGUE1.PMD(X = X,k = k,PMD2 = PMD2,extra.sparse = extra.sparse,seed1 = seed1)
    emp.p<-DIALOGUE1.PMD.empirical(X,k,n1 = 100,extra.sparse = extra.sparse)
    if(bypass.emp){
      emp.p1<-emp.p
      emp.p[]<-0.05
      print("Not using empirical p-values!")
      print(emp.p)
    }
  }else{
    out<-DIALOGUE1.PMD.pairwise(X,k,specific.pair)
    if(out$message=="No programs"){return(out)}
    rA<-rA[specific.pair]
    X<-X[specific.pair]
    emp.p<-DIALOGUE1.PMD.empirical(X,k,n1 = 20,extra.sparse = extra.sparse)
    emp.p<-subset.matrix(emp.p,is.element(rownames(emp.p),colnames(out$ws[[1]])))
    rownames(emp.p)<-paste0("MCP",1:nrow(emp.p))
    main<-paste(main,paste(specific.pair,collapse = "."),sep = "_")
    cell.types<-names(rA)
    n1<-length(cell.types)
  }
  
  Y<-lapply(names(X), function(i) X[[i]]%*%out$ws[[i]])
  names(Y)<-names(X)
  pairs1<-t(combn(names(X),2))
  cca.cor<-apply(pairs1,1,function(x) diag(cor(Y[[x[1]]],Y[[x[2]]])))
  cca.cor.p<-apply(pairs1,1,function(x) diag(spearman.cor(Y[[x[1]]],Y[[x[2]]],method = "pearson")$p))
  colnames(cca.cor)<-paste(pairs1[,1],pairs1[,2],sep = "_")
  colnames(cca.cor.p)<-colnames(cca.cor)
  y<-list()
  param<-list(conf = conf,covar = covar,n.genes = n.genes,PMD2 = PMD2,extra.sparse = extra.sparse,
              averaging.function = averaging.function,p.anova = p.anova,specific.pair = specific.pair,
              center.flag = center.flag,seed1 = seed1)
  R<-list(name = paste0("DIALOGUE1_",main),param = param,
          cell.types = cell.types,k = c(k,laply(X,ncol)),
          samples = samplesU,
          sample.PCs = X,
          samples.cells = list(),
          conf = conf,
          covar = covar,
          emp.p = emp.p,
          cca = out,cca.cor = list(R = cca.cor,P = cca.cor.p),
          cca.scores = list(),cca.gene.cor = list(),
          cca.sig = list(),cca.redun.cor = list())
  R$MCP.cell.types <- DIALOGUE.identify.cell.types(R)
  names(R$k)<-c("DIALOGUE",paste0("original.",cell.types))
  R$message<-paste("DIALOGUE1 found",nrow(emp.p),"programs.")
  
  for(x in cell.types){
    r<-rA[[x]]
    y[[x]]<-r@X[,rownames(out$ws[[x]])]%*%out$ws[[x]]
    scores0<-as.matrix(y[[x]])
    conf.m<-r@metadata[,is.element(colnames(r@metadata),conf)]
    r@scores<-t(get.residuals(t(scores0),conf.m))
    R$cca.scores[[x]]<-r@scores
    # R$cca.gene.cor[[x]]<-cor(t(r@tpm),r@scores)
    # R$cca.sig[[x]]<-get.top.cor(R$cca.gene.cor[[x]],q = n.genes,min.ci = 0.05)
    
    R$cca.gene.cor1[[x]]<-cor(t(r@tpm),r@scores)
    g1<-sort(unique(unlist(get.top.cor(R$cca.gene.cor1[[x]],q = n.genes,min.ci = 0.05))))
    R$cca.gene.cor[[x]]<-pcor.mat(t(r@tpm[g1,]),r@scores,r@cellQ)
    C1<-R$cca.gene.cor[[x]]$R
    P1<-R$cca.gene.cor[[x]]$P
    C1[P1>(0.05/nrow(r@tpm))]<-0
    R$cca.sig[[x]]<-get.top.cor(C1,q = n.genes,min.ci = 0.05)
    
    R$cca.redun.cor[[x]]<-cor(r@scores[,1:k])
    R$samples.cells[[x]]<-r@samples
  }
  if(bypass.emp){R$emp.p1<-emp.p1}
  saveRDS(R,file = paste0(results.dir,R$name,".rds"))
  dir.create(paste0(results.dir,"DIALOGUE2_",main))
  return(R)
}

DIALOGUE1.PMD<-function(X,k,PMD2 = F,extra.sparse = F,seed1 = 1234){
  set.seed(seed1)
  if(extra.sparse){
    perm.out <- MultiCCA.permute(X,type=rep("standard",length(X)),trace = F) 
  }else{
    perm.out <- MultiCCA.permute(X,type=rep("standard",length(X)),trace = F,penalties = sqrt(ncol(X[[1]]))/2)
  }
  out <- MultiCCA(X, type=rep("standard",length(X)),
                  penalty=perm.out$bestpenalties,niter = 100,
                  ncomponents=k, ws=perm.out$ws.init,trace = F)
  names(out$ws)<-names(X)
  m<-laply(out$ws,function(x) colSums(x!=0))
  colnames(m)<-paste0("MCP",1:ncol(m))
  rownames(m)<-names(X)
  for(i in names(X)){
    colnames(out$ws[[i]])<-paste0("MCP",1:k)
    rownames(out$ws[[i]])<-colnames(X[[i]])
  }
  if(!PMD2){return(out)}
  
  print("PMD #1");print("Number of features");print(m);
  set.seed(seed1)
  perm.out <- MultiCCA.permute(X,type=rep("standard",length(X)),trace = F,penalties = perm.out$penalties[,4:10])
  out <- MultiCCA(X, type=rep("standard",length(X)),
                  penalty=perm.out$bestpenalties,niter = 100,
                  ncomponents=k, ws=perm.out$ws.init,trace = F)
  m<-laply(out$ws,function(x) colSums(x!=0));rownames(m)<-names(X)
  print("PMD #2");print("Number of features");print(m)
  return(out)
}

DIALOGUE1.PMD.empirical<-function(X,k,seed = 1234,n1 = 20,extra.sparse = F,full.output = F){
  if(n1>1){
    cca.cor<-lapply(1:k, function(x) NULL)
    v<-DIALOGUE1.PMD.empirical(X,k,seed = 0,n1 = 1,extra.sparse = extra.sparse)
    cca.cor<-lapply(1:k, function(x) rbind(cca.cor[[x]],v[x,]))
    for(x in 2:n1){
      v<-DIALOGUE1.PMD.empirical(X,k,seed = sample(1:10000,1),n1 = 1,extra.sparse = extra.sparse)
      cca.cor<-lapply(1:k, function(x) rbind(cca.cor[[x]],v[x,]))
    }
    names(cca.cor)<-paste0("MCP",1:k)
    P<-laply(cca.cor,function(m){
      p<-ranksum.test.mat(as.matrix(t(m)),c(T,rep(F,nrow(m)-1)))[,"more"]
      return(p)})
    if(is.null(dim(P))){P<-as.matrix(P)}
    rownames(P)<-names(cca.cor)
    return(P)
  }
  if(seed>0){
    set.seed(seed)
    X<-lapply(X,function(X1){
      X2<-apply(X1,2,function(x) sample(x,length(x)))
      rownames(X2)<-rownames(X1)
      return(X2)})
  }
  if(extra.sparse){
    perm.out <- MultiCCA.permute(X,type=rep("standard",length(X)),trace = F) 
  }else{
    perm.out <- MultiCCA.permute(X,type=rep("standard",length(X)),trace = F,penalties = sqrt(ncol(X[[1]]))/2)
  }
  out <- MultiCCA(X, type=rep("standard",length(X)),
                  penalty=perm.out$bestpenalties,niter = 100,
                  ncomponents=k, ws=perm.out$ws.init,trace = F)
  names(out$ws)<-names(X)
  for(i in names(X)){
    colnames(out$ws[[i]])<-paste0("MCP",1:ncol(out$ws[[i]]))
    rownames(out$ws[[i]])<-colnames(X[[i]])
  }
  Y<-lapply(names(X), function(i) X[[i]]%*%out$ws[[i]])
  names(Y)<-names(X)
  pairs1<-t(combn(names(X),2))
  cca.cor<-apply(pairs1,1,function(x) diag(cor(Y[[x[1]]],Y[[x[2]]])))
  colnames(cca.cor)<-paste(pairs1[,1],pairs1[,2],sep = "_")
  if(full.output){
    out$Y<-Y
    out$cor<-cca.cor
    return(out)
  }
  return(cca.cor)
}

DIALOGUE1.PMD.pairwise<-function(X,k,specific.pair){
  out<-DIALOGUE1.PMD.empirical(X = X,k = k,n1 = 1,full.output = T,seed = 0)
  X1<-X[specific.pair]
  x1<-specific.pair[1]
  x2<-specific.pair[2]
  out1<-DIALOGUE1.PMD.empirical(X = X1,k = k,n1 = 1,full.output = T,seed = 0)
  c1<-spearman.cor(out$Y[[x1]],out1$Y[[x1]],method = "pearson")
  c2<-spearman.cor(out$Y[[x2]],out1$Y[[x2]],method = "pearson")
  c1$padj<-p.adjust.mat(c1$p)
  c2$padj<-p.adjust.mat(c2$p)
  B<-c1$padj<0.05&c2$padj<0.05&abs(c1$cor)>0.3&abs(c2$cor)>0.3
  m<-laply(1:ncol(c1$cor),function(i){
    x<-c1$cor[,i]
    idx<-which(abs(x)==max(abs(x)))
    return(c(c1$cor[idx,i],c2$cor[idx,i],c1$padj[idx,i],c2$padj[idx,i]))
  })
  rownames(m)<-colnames(c1$cor)
  colnames(m)<-c("R1","R2","P1","P2")
  b<-colSums(B)==0
  if(!any(b)){
    print(paste("No unique programs specific to",specific.pair[1],"and",specific.pair[2]))
    rslts<-list(cor = m,message = "No programs")
    return(rslts)
  }else{
    print(paste("Identified",sum(b),"unique programs for",specific.pair[1],"and",specific.pair[2]))
  }
  out1$mcp.comp<-m
  out1$ws[[x1]]<-out1$ws[[x1]][,b]
  out1$ws[[x2]]<-out1$ws[[x2]][,b]
  out1$message<-paste(sum(b),"programs.")
  return(out1)
}

DIALOGUE2<-function(rA,main,results.dir = "~/Desktop/DIALOGUE.results/",subsample.flag = T,parallel.vs = F){
  print("#************DIALOGUE Step II: HLM ************#")
  if(!grepl("\\/$", results.dir)){results.dir <- paste(results.dir, "/", sep = "")}
  cell.types<-names(rA)
  if(missing(main)){main<-paste0(cell.types,collapse = "_")}
  file1<-paste0(results.dir,"DIALOGUE1_",main,".rds")
  file2<-paste0(results.dir,"DIALOGUE2_",main,".rds")
  
  R<-readRDS(file1)
  R$frm<-paste0("y ~ (1 | samples) + x + ",paste(R$covar,collapse = " + "))
  
  k2<-ncol(R$cca$ws[[1]])
  pairs1<-t(combn(cell.types,2))
  sig<-R$cca.sig
  
  f<-function(i){
    x1<-pairs1[i,1];x2<-pairs1[i,2]
    print(paste("#************DIALOGUE Step II (multilevel modeling):",x1,"vs.",x2,"************#"))
    p<-paste0(x1,".vs.",x2)
    rslts<-DIALOGUE2.pair(R,rA[[x1]],rA[[x2]],cell.types,results.dir)
    return(rslts)
  }
  
  if(parallel.vs){
    R1<-mclapply(1:nrow(pairs1),f)
    names(R1)<-paste(pairs1[,1],pairs1[,2],sep = ".vs.")
    R<-c(R,R1)
  }else{
    for(i in 1:nrow(pairs1)){
      x1<-pairs1[i,1];x2<-pairs1[i,2]
      print(paste("#************DIALOGUE Step II (multilevel modeling):",x1,"vs.",x2,"************#"))
      R[[paste0(x1,".vs.",x2)]]<-DIALOGUE2.pair(R,rA[[x1]],rA[[x2]],cell.types,results.dir)
    }
  }
  
  R$name<-paste0("DIALOGUE2_",main)
  saveRDS(R,file = file2)
  return(R)
}

DIALOGUE2.pair<-function(R,r1,r2,cell.types,results.dir){
  x1<-r1@name;x2<-r2@name
  MCP.names<-names(R$MCP.cell.types)[laply(R$MCP.cell.types,function(x) sum(is.element(c(x1,x2),x))==2)]
  if(is.null(MCP.names)|length(MCP.names)==0){
    print("No MCPs identified for these cell types.")
    return(NULL)
  }
  print(paste(length(MCP.names),"MCPs identified for these cell types."))
  main<-gsub("DIALOGUE1_","",R$name)
  
  saveFile<-paste0(results.dir,"DIALOGUE2_",main,"/",x1,".vs.",x2,".rds")
  if(file.exists(saveFile)){
    return(readRDS(saveFile))
  }
  sig1<-R$cca.sig[[x1]]
  sig2<-R$cca.sig[[x2]]
  idx<-intersect(get.abundant(r1@samples),get.abundant(r2@samples))
  r1<-set.cell.type(r1,is.element(r1@samples,idx))
  r2<-set.cell.type(r2,is.element(r2@samples,idx))
  r1@scores<-R$cca.scores[[x1]][r1@cells,]
  r2@scores<-R$cca.scores[[x2]][r2@cells,]
  r<-DLG.get.OE(r1,r2,sig1,sig2,compute.scores = F);r1<-r$r1;r2<-r$r2
  r1@tme<-r2@tpmAv[,as.character(r1@samples)]
  r2@tme<-r1@tpmAv[,as.character(r2@samples)]
  
  #### NEW (added 06/01/21)
  r1@tme.qc<-as.matrix(r2@qcAv)[as.character(r1@samples),2]
  r2@tme.qc<-as.matrix(r1@qcAv)[as.character(r2@samples),2]
  #### NEW (ended)
  
  r1a<-cell.type.2.list(r1)
  r2a<-cell.type.2.list(r2)
  
  # r1a<-set.list(r1a,sample.per.label(r1a$samples,50),sampleName = paste0(r1a$name,"_A"))
  # r2a<-set.list(r2a,sample.per.label(r2a$samples,50),sampleName = paste0(r2a$name,"_A"))
  
  f1<-function(sig1,sig2,x){
    p1<-DIALOGUE2.mixed.effects(r2a,x,sig1,R$frm)
    p2<-DIALOGUE2.mixed.effects(r1a,x,sig2,R$frm)
    sig1f<-intersect.list1(get.top.cor(p1[!is.na(p1$Z),],q = 100,idx = "Z",min.ci = 1),r1@genes)
    sig2f<-intersect.list1(get.top.cor(p2[!is.na(p2$Z),],q = 100,idx = "Z",min.ci = 1),r2@genes)
    names(sig1f)<-gsub("Z.",paste0(x,"."),names(sig1f))
    names(sig2f)<-gsub("Z.",paste0(x,"."),names(sig2f))
    p1$program<-x;p2$program<-x
    p1$genes<-rownames(p1);p2$genes<-rownames(p2)
    results<-list(p1 = p1,p2 = p2,sig1f = sig1f,sig2f = sig2f)
    return(results)
  }
  
  # idx<-unique(get.strsplit(names(sig1),".",1))
  R1<-lapply(MCP.names,function(x){f1(sig1,sig2,x)})
  names(R1)<-MCP.names
  R1$p1<-NULL;R1$p2<-NULL
  for(x in MCP.names){
    R1$p1<-rbind(R1$p1,R1[[x]]$p1)
    R1$p2<-rbind(R1$p2,R1[[x]]$p2)
  }
  R1$sig1<-lapply(R[MCP.names], function(x) x$sig1f)
  R1$sig2<-lapply(R[MCP.names], function(x) x$sig2f)
  R1$name<-paste0(x1,".vs.",x2)
  saveRDS(R1,file = saveFile)
  return(R1)
}

DIALOGUE2.mixed.effects<-function(r1,x,sig2,frm = "y ~ (1 | samples) + x + cellQ"){
  # r1 was a cell.type S4 that was converted to a list.
  genes<-unlist(sig2[paste0(x,c(".up",".down"))])
  b<-is.element(genes,rownames(r1$tme))
  p<-apply.formula.HLM(r1,r1$scores[,x],
                       X = r1$tme[genes[b],],
                       MARGIN = 1,formula = frm)
  p$pval<-p.adjust(p$P,method = "BH")
  p$up<-is.element(rownames(p),sig2[[paste0(x,".up")]])
  if(all(b)){return(p)}
  P<-get.mat(genes,colnames(p))
  P[b,]<-as.matrix(p)
  rownames(P)<-gsub(".","-",rownames(P),fixed = T)
  P<-as.data.frame(P)
  return(P)
}

DIALOGUE3<-function(rA,main,results.dir = "~/Desktop/DIALOGUE.results/",full.version = F,pheno = NULL){
  print("#************Finalizing the scores************#")
  cell.types<-names(rA)
  if(missing(main)){main<-paste0(cell.types,collapse = "_")}
  print(paste0(results.dir,"/DIALOGUE2_",main,".rds"))
  R<-readRDS(paste0(results.dir,"/DIALOGUE2_",main,".rds"))
  R$gene.pval<-lapply(R$cell.types, function(x){DLG.multi.get.gene.pval(x,R)})
  names(R$gene.pval)<-R$cell.types
  
  rA<-lapply(rA,function(r){
    r<-DLG.find.scoring(r,R)
    return(r)
  })
  names(rA)<-cell.types
  R$pref<-list()
  idx<-unique(get.strsplit(names(R$sig[[1]]),".",1))
  pairs1<-t(combn(cell.types,2))
  
  for(i in 1:nrow(pairs1)){
    x1<-pairs1[i,1];x2<-pairs1[i,2]
    x<-paste0(x1,".vs.",x2)
    r1<-rA[[x1]];r2<-rA[[x2]]
    r<-DLG.get.OE(r1,r2,plot.flag = F,compute.scores = F)
    r1<-r$r1;r2<-r$r2
    idx<-intersect(get.abundant(r1@samples),get.abundant(r2@samples))
    R$pref[[x]]<-cbind.data.frame(R = diag(cor(r1@scoresAv[idx,],r2@scoresAv[idx,])),
                                  hlm = DLG.hlm.pval(r1,r2,formula = R$frm))
  }
  
  R$gene.pval<-lapply(rA,function(r1) r1@gene.pval)
  R$sig1<-lapply(rA,function(r1) r1@sig$sig1)
  R$sig2<-lapply(rA,function(r1) r1@sig$sig2)
  R$scores<-lapply(rA,function(r1){
    X<-cbind.data.frame(r1@scores,samples = r1@samples,
                        cells = r1@cells, cell.type = r1@name,
                        r1@metadata)
    return(X)})
  names(R$gene.pval)<-cell.types
  names(R$sig1)<-cell.types
  names(R$sig2)<-cell.types
  names(R$scores)<-cell.types
  
  R$name<-paste0("DLG.output_",main)
  R$MCPs.full<-sig2MCP(R$sig1)
  R$MCPs<-sig2MCP(R$sig2)
  
  R$cca.fit<-laply(R$cell.types,function(x) diag(cor(R$cca.scores[[x]],R$scores[[x]][,1:R$k["DIALOGUE"]])))
  rownames(R$cca.fit)<-R$cell.types
  
  fileName<-paste0(results.dir,"DLG.full.output_",main,".rds")
  # if(full.version){saveRDS(R,file = fileName)}
  
  if(!is.null(pheno)){R$phenoZ<-DIALOGUE.pheno(R,pheno = pheno)}
  if(full.version){saveRDS(R,file = fileName)}
  
  R1<-R[intersect(names(R),c("cell.types","scores","gene.pval","param","MCP.cell.types","MCPs",
                             "pref","k","name","phenoZ",results.dir))]
  fileName<-paste0(results.dir,"DLG.output_",main,".rds")
  saveRDS(R1,file = fileName)
  
  if(!full.version){
    file.remove(paste0(results.dir,"DIALOGUE1_",main,".rds"))
    file.remove(paste0(results.dir,"DIALOGUE2_",main,".rds"))
    unlink(paste0(results.dir,"DIALOGUE2_",main,"/"),recursive = T)
  }
  return(R1)
}

DLG.get.OE<-function(r1,r2,sig1,sig2,plot.flag = F,compute.scores = T){
  if(compute.scores){
    r1@scores<-get.OE(r1,sig1,semi.flag = T)
    r2@scores<-get.OE(r2,sig2,semi.flag = T)
  }
  r1@scoresAv<-average.mat.rows(r1@scores,r1@samples,f = colMedians)
  r2@scoresAv<-average.mat.rows(r2@scores,r2@samples,f = colMedians)
  r1@tme.OE<-r2@scoresAv[match(r1@samples,rownames(r2@scoresAv)),]
  r2@tme.OE<-r1@scoresAv[match(r2@samples,rownames(r1@scoresAv)),]
  r<-list(r1 = r1,r2= r2)
  if(!plot.flag){return(r)}
  # r$cor<-DLG.cor.plot(r1,r2,sd.flag = F,q1 = 1/3,q2 = 2/3)
  DLG.cor.plot(r1,r2,sd.flag = F,q1 = 1/3,q2 = 2/3)
  return(r)
}

DLG.fix.sig.names<-function(sig){
  b<-!get.abundant(get.strsplit(names(sig),".",1),boolean.flag = T)
  names(sig)[b]<-get.strsplit(names(sig[b]),".",1)
  return(sig)
}

DLG.multi.get.gene.pval<-function(cell.type,R){
  b<-grepl("vs.",names(R))
  pairs1<-get.strsplit(names(R),".vs.",1:2)
  b1<-b&is.element(pairs1[,1],cell.type)
  b2<-b&is.element(pairs1[,2],cell.type)
  if((sum(b1)+sum(b2))==0){return(NULL)}
  f1<-function(m1,pi = "p1"){
    x<-m1[[pi]]
    rownames(x)<-paste0(x$program,ifelse(x$up,".up_",".down_"),x$genes)
    return(x)
  }
  m<-c(lapply(R[b1],f1),lapply(R[b2],function(m1) f1(m1,"p2")))
  g<-unique(unlist(lapply(m,rownames)))
  p<-get.strsplit(g,"_",1:2)
  colnames(p)<-c("programF","genes")
  rownames(p)<-g
  p<-cbind.data.frame(p,program = get.strsplit(g,".",1),
                      up = grepl("up",g))
  names(m)<-gsub(paste0(cell.type,".vs."),"",names(m))
  names(m)<-gsub(paste0(".vs.",cell.type),"",names(m))
  for(i in names(m)){
    x<-m[[i]]
    g1<-paste0(x$program,ifelse(x$up,".up_",".down_"),x$genes)
    idx<-match(g,g1)
    p[,i]<-x$Z[idx]
  }
  
  # AD was without adjustments
  p.up<-p.adjust.mat.per.label(get.pval.from.zscores(p[,5:ncol(p)]),p$programF)
  p.down<-p.adjust.mat.per.label(get.pval.from.zscores(-p[,5:ncol(p)]),p$programF)
  
  p.ub<-0.1
  if(!is.matrix(p.up)){p.up<-as.matrix(p.up);p.down<-as.matrix(p.down)}
  m<-cbind.data.frame(p[,c(names(m),"programF","genes")],
                      p.up = fisher.combine(p.up),
                      p.down = fisher.combine(p.down),
                      n.up = rowSums(p.up<p.ub,na.rm = T),
                      nf.up = rowMeans(p.up<p.ub,na.rm = T),
                      n.down = rowSums(p.down<p.ub,na.rm = T),
                      nf.down = rowMeans(p.down<p.ub,na.rm = T),
                      p[,c("program","up")])
  m$N<-m$n.up
  m$N[!m$up]<-m$n.down[!m$up]
  m$Nf<-m$nf.up
  m$Nf[!m$up]<-m$nf.down[!m$up]
  m$p.up[!m$up]<-1
  m$p.down[m$up]<-1
  return(m)
  
}

DLG.find.scoring<-function(r1,R){
  gene.pval<-R$gene.pval[[r1@name]]
  if(is.null(gene.pval)){
    print("No MCPs identified.")
    r1<-DLG.initialize(r1,R)
    return(r1)
  }
  gene.pval<-gene.pval[is.element(gene.pval$genes,r1@genes)&!is.na(gene.pval[,1]),]
  g<-sort(unique(gene.pval$genes))
  
  # r1@cca.scores0<-r1@X%*%R$cca$ws[[r1@name]]
  WS<-R$cca$ws[[r1@name]]
  # r1@extra.scores$cca0<-r1@X%*%R$cca$ws[[r1@name]]
  r1@extra.scores$cca0<-r1@X[,rownames(WS)]%*%WS
  r1@zscores<-center.large.matrix(r1@tpm[g,],sd.flag = T)
  
  f<-function(x){
    y<-r1@extra.scores$cca0[,x]
    b<-gene.pval$program==x
    if(!any(b)){
      return(list(gene.pval = NULL,scores = y))
    }
    gene.pval<-gene.pval[b,]
    X<-t(r1@zscores[gene.pval$genes,])
    b<-is.element(colnames(X),gene.pval$genes[!gene.pval$up])
    X[,b]<-(-X[,b])
    gene.pval<-DLG.iterative.nnls(X,y,gene.pval)
    scores<-X%*%gene.pval$coef
    return(list(gene.pval = gene.pval,scores = scores))
  }
  m<-lapply(colnames(r1@extra.scores$cca0), f)
  # r1@scores0<-t(laply(m,function(x) x$scores))
  conf.m<-r1@metadata[,is.element(colnames(r1@metadata),R$conf)]
  r1@extra.scores$nnl0<-t(laply(m,function(x){return(as.matrix(x$scores))}))
  colnames(r1@extra.scores$nnl0)<-colnames(r1@extra.scores$cca0)
  r1@extra.scores$cca<-t(get.residuals(t(r1@extra.scores$cca0),conf.m))
  r1@scores<-t(get.residuals(t(r1@extra.scores$nnl0),conf.m))
  r1@scoresAv<-average.mat.rows(r1@scores,r1@samples)
  r1@gene.pval<-NULL
  for(x in m){
    r1@gene.pval<-rbind(r1@gene.pval,x$gene.pval)
  }
  m1<-r1@gene.pval
  m2<-r1@gene.pval
  
  idx<-laply(R$MCP.cell.types,length)
  names(idx)<-names(R$MCP.cell.types)
  m1$n.cells<-idx[m1$program]
  lb<-ceil(m1$n.cells/2)
  m1<-m1[m1$coef>0|(m1$n.up>=lb&m1$p.up<0.05)|(m1$n.down>=lb&m1$p.down<0.05),]
  m2<-m2[m2$Nf==1&(m2$p.up<0.05|m2$p.down<0.05),]
  r1@sig<-list(sig1 = split(m1$genes,m1$programF),sig2 = split(m2$genes,m2$programF))
  return(r1)
}

DLG.initialize<-function(r1,R){
  gene.pval<-R$gene.pval[[r1@name]]
  WS<-R$cca$ws[[r1@name]]
  r1@extra.scores$cca0<-r1@X[,rownames(WS)]%*%WS
  conf.m<-r1@metadata[,is.element(colnames(r1@metadata),R$conf)]
  r1@extra.scores$cca<-t(get.residuals(t(r1@extra.scores$cca0),conf.m))
  r1@scores<-r1@extra.scores$cca
  r1@scoresAv<-average.mat.rows(r1@scores,r1@samples)
  r1@gene.pval<-NULL
  r1@sig<-NULL
  return(r1)
}

DLG.iterative.nnls<-function(X,y,gene.pval){
  set.seed(1234)
  f.rank<-gene.pval$Nf
  y1<-y;y.fit<-rep(0,length(y))
  v<-list()
  gene.pval$coef<-0
  idx<-sort(unique(f.rank),decreasing = T)
  idx<-idx[idx>=(1/3)]
  for(n1 in idx){
    b1<-f.rank==n1
    if(sum(b1,na.rm = T)<5){next()}
    X1<-X[,b1]
    main<-paste0("N",n1)
    v[[main]]<-nnls::nnls(X1,y)
    y<-v[[main]]$residuals
    y.fit<-y.fit+v[[main]]$fitted
    gene.pval$coef[b1]<-v[[main]]$x
    if(length(unique(y.fit))>10 & cor(y1,y.fit)>0.95){
      # cor.plot(y.fit,y1,main = paste("NNLS fitting",n1))
      return(gene.pval)
    }
  }
  if(sum(f.rank<n1,na.rm = T)<5){return(gene.pval)}
  X1<-X[,f.rank<n1]
  main<-paste0("Ns",n1)
  v[[main]]<-nnls::nnls(X1,y)
  y<-v[[main]]$residuals
  y.fit<-y.fit+v[[main]]$fitted
  gene.pval$coef[f.rank<n1]<-v[[main]]$x
  cor.plot(y.fit,y1,main = paste("NNLS fitting -",n1))
  return(gene.pval)
}

DLG.cor.plot<-function(r1,r2,idx,q1 = 0.25,q2 = 0.75,sd.flag = F){
  if(missing(idx)){
    idx<-sort(unique(get.strsplit(colnames(r1@scores),".",1)))
  }
  if(length(idx)>1){
    par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
    lapply(idx, function(x){
      DLG.cor.plot(r1,r2,x,q1 = q1,q2 = q2,sd.flag = sd.flag)
      return(x)
    })
    return()
  }
  idx1<-intersect(rownames(r1@scoresAv),rownames(r2@scoresAv))
  x<-r1@scoresAv[idx1,idx]
  y<-r2@scoresAv[idx1,idx]
  # pheno<-is.element(rownames(r1@scoresAv),r1@samples[r1@pat1])
  x0<-laply(idx1,function(x) quantile(r1@scores[r1@samples==x,idx],q1))
  x1<-laply(idx1,function(x) quantile(r1@scores[r1@samples==x,idx],q2))
  y0<-laply(idx1,function(x) quantile(r2@scores[r2@samples==x,idx],q1))
  y1<-laply(idx1,function(x) quantile(r2@scores[r2@samples==x,idx],q2))
  if(sd.flag){
    xd<-laply(idx1,function(x) sd(r1@scores[r1@samples==x,idx]))/2
    yd<-laply(idx1,function(x) sd(r2@scores[r2@samples==x,idx]))/2
    x0<-x-(xd);x1<-x+(xd)
    y0<-y-(yd);y1<-y+(yd)
  }
  cor.plot(x,y,main = paste("Program",idx),cex = 1,
           xlim = c(min(c(x,x0)),max(c(x,x1))),
           ylim = c(min(c(y,y0)),max(c(y,y1))))
  for(i in 1:length(x)){
    # col<-ifelse(pheno[i],"red","grey")
    col<-"grey"
    arrows(x0=x0[i], y0=y[i], x1=x1[i], y1=y[i], code=3, col=col, lwd=2,angle=0,length = 0)
    arrows(x0=x[i], y0=y0[i], x1=x[i], y1=y1[i], code=3, col=col, lwd=2,angle=0,length = 0)
  }
  # points(x,y,col = ifelse(pheno,"red","black"),pch = 16,xlim = c(min(x0),max(x1)),ylim = c(min(y0),max(y1)))
  points(x,y,col = "black",pch = 16,xlim = c(min(x0),max(x1)),ylim = c(min(y0),max(y1)))
  return(spearman.cor(x,y))
  
}

DLG.hlm.pval<-function(r1,r2,formula = "y ~ (1 | samples) + x + cellQ"){
  idx<-c("samples","scores","tme.OE",intersect(slotNames(r1),setdiff(gsub(" ","",get.strsplit(formula,"+ ",1:10)),NA)))
  r1<-cell.type.2.list(r1,idx = idx)
  r2<-cell.type.2.list(r2,idx = idx)
  f<-function(x){
    p<-c(p1 = formula.HLM(r1,y = r1$scores[,x],x = r1$tme.OE[,x],formula = formula),
         p2 = formula.HLM(r2,y = r2$scores[,x],x = r2$tme.OE[,x],formula = formula))
    return(p)
  }
  idx<-unique(get.strsplit(colnames(r1$scores),".",1))
  m<-laply(idx,f)[,c(2,4)]
  rownames(m)<-idx
  return(m)
}

DLG.plot1<-function(R,i,mark.samples = NULL,d = 1){
  R$cell.types<-names(R$cca.scores)
  k<-ncol(R$cca.scores[[1]])
  R$scoresAv<-lapply(R$cell.types,function(x){
    m<-R$cca.scores[[x]]
    X<-average.mat.rows(as.matrix(m[,1:k]),R$samples.cells[[x]],f = colMedians)
    return(X)
  })
  idx<-unlist(lapply(R$scoresAv, function(m) rownames(m)))
  idx<-get.abundant(idx,length(R$cell.types))
  R$scoresAv<-lapply(R$scoresAv,function(m) m[idx,])
  col<-rep("black",length(idx))
  col[is.element(idx,mark.samples)]<-"red"
  pch<-ifelse(any(col=="red"),21,16)
  f<-function(i){
    m1<-t(laply(R$scoresAv,function(m) m[,i]))*d
    rownames(m1)<-idx;colnames(m1)<-R$cell.types
    pairs.panels(m1,hist.col = "grey",breaks = 50,bg = col,pch = pch,ellipses = F,smooth = T,lm = T,stars = T)
    title(i)
    return(m1)
  }
  idx1<-paste0("MCP",1:R$k["DIALOGUE"])
  if(!missing(i)){m1<-f(i);return(cor(m1))}
  m<-lapply(idx1, f)
}

DLG.add.metadata<-function(R,rA){
  if(missing(rA)){
    R$metadata<-lapply(R$cell.types,function(x) R$scores[[x]][,(R$k[1]+1):ncol(R$scores[[x]])])
    names(R$metadata)<-R$cell.types
    return(R)
  }
  
  R$metadata<-lapply(rA,function(r1){
    X<-cbind.data.frame(samples = r1@samples,
                        cells = r1@cells,
                        cell.type = r1@name,
                        r1@metadata)
    return(X)})
  names(R$metadata)<-R$cell.types
  return(R)
}

sig2MCP<-function(R.sig,k = 5){
  sig1<-unlist(R.sig,recursive = F)
  sig1<-c(sig1[grepl("up",names(sig1))],sig1[grepl("down",names(sig1))])
  MCPs<-lapply(1:k,function(x){
    idx<-paste0("MCP",x,".")
    sig1<-sig1[grepl(idx,names(sig1))]
    names(sig1)<-gsub(idx,"",names(sig1))
    return(sig1)
  })
  names(MCPs)<-paste0("MCP",1:k)
  return(MCPs)
}

DIALOGUE.identify.cell.types<-function(R){
  if(length(R$cell.types)==2){
    MCP.cell.types<-lapply(R$emp.p, function(x){
      if(x<0.1){return(R$cell.types)}
      return(NULL)
    })
    names(MCP.cell.types)<-rownames(R$emp.p)
    return(MCP.cell.types)
  }
  emp.p<-R$emp.p
  emp.p2<-R$emp.p
  colnames(emp.p2)<-paste(get.strsplit(colnames(emp.p),"_",2),
                          get.strsplit(colnames(emp.p),"_",1),sep = "_")
  emp.p<-cbind(emp.p,emp.p2)
  MCP.cell.types<-list()
  for(x in paste0("MCP",1:nrow(emp.p))){
    cell.types<-R$cell.types
    p1<-generic.vector2mat(emp.p[x,])
    diag(p1)<-0.05
    lb<-min(rowSums(p1<0.1))
    while(lb<ceil(length(cell.types)/2)){
      cell.types<-setdiff(cell.types,cell.types[rowSums(p1<0.1)<=lb][1])
      p1<-p1[cell.types,cell.types]
      lb<-min(rowSums(p1<0.1))
    }
    MCP.cell.types[[x]]<-cell.types
  }
  return(MCP.cell.types)
}

pcor.mat<-function(v1,v2,v3, method = 'spearman',use = "pairwise.complete.obs",alternative = "two.sided"){
  f<-function(x1,x2,x3){
    c.i<-tryCatch(pcor.test(x1,x2,x3,method = method),
                  error = function(err){return(NA)})
    if(is.list(c.i)){return(c(c.i$estimate,c.i$p.value))}
    return(c(NA,NA))
  }
  
  f<-function(x1,x2,x3){
    c.i<-pcor.test(x1,x2,x3,method = method)
    return(c(c.i$estimate,c.i$p.value))
  }
  
  P<-get.mat(colnames(v1),colnames(v2));R<-P
  for(x in 1:ncol(v2)){
    x2<-v2[,x]
    c1<-apply(v1,2,function(x1) f(x1,x2,v3))
    R[,x]<-c1[1,]
    P[,x]<-c1[2,]
  }
  padj<-p.adjust.mat(P,method = "BH")
  rslts<-list(R = R,P = P,padj = padj)
  return(rslts)
}

p.adjust.mat.per.label<-function(p,v){
  p1<-get.mat(rownames(p),colnames(p),data = NA)
  for(x in unique(v)){
    b<-is.element(v,x)
    if(is.null(ncol(p1))||ncol(p1)<2){
      p1[b]<-p.adjust(p[b])
    }else{
      p1[b,]<-p.adjust.mat(p[b,])
    }
    
  }
  return(p1)
}

generic.vector2mat<-function(v,rn = get.strsplit(names(v),"_",1),cn = get.strsplit(names(v),"_",2)){
  rnu <- sort(unique(rn))
  cnu <- sort(unique(cn))
  m<-get.mat(data = NA,m.rows = rnu,m.cols = cnu)
  for(x in rnu){
    v1<-v[rn==x]
    cn1 <- cn[rn==x]
    idx<-match(cnu,cn1)
    m[x,]<-v1[idx]
  }
  return(m)
}

