#' DIALOGUE.run
#'
#' DIALOGUE is a dimensionality reduction approach that uses cross-cell-type
#' associations to identify multicellular programs and map the cell transcriptome
#' as a function of its environment. 
#' 
#' @param rA list of \linkS4class{cell.type} objects.
#' @param k the number of multicellular programs to identify;
#' @param main results' name;
#' @param results.dir path to the results directory, where the output will be saved;
#' @param plot.flag if TRUE then \code{\link{DIALOGUE.plot}} will be called to plot the results; default is FALSE;
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
#' @example 
#' To run a toy example, first download
#' rA<-readRDS(system.file("extdata", "toy.example.rds", package = "DIALOGUE"))
#' summary(rA)
#' Length Class     Mode
#' TA2         1      cell.type S4  
#' Macrophages 1      cell.type S4
#' CD8         1      cell.type S4
#' Find multicellular programs:
#' R<-DIALOGUE.run(rA = rA,main = "toy.example",k = 2,results.dir = "~/Desktop/DIALOGUE.results/")
#' 
#' To regenerate the IBD results provided in our manuscript:
#' rA<-readRDS(system.file("extdata", "toy.example.rds", package = "DIALOGUE"))
#' R<-DIALOGUE.run(rA = rA,main = "IBD",k = 5,results.dir = "~/Desktop/DIALOGUE.results/")
#' 
#' @seealso See \href{https://github.com/livnatje/DIALOGUE}{DIALOGUE GitHub page} for more details.
#' \code{\link{DIALOGUE.plot}}
#' @author Livnat Jerby-Arnon
#' @references 
#' Jerby-Arnon and Regev, Mapping multicellular configurations using single-cell data
#' @export
#' 
DIALOGUE.run<-function(rA,main,k = 2,
                       results.dir = "~/Desktop/DIALOGUE.results/",
                       plot.flag = T,add.effects = F){
  full.version <- F
  R<-DIALOGUE1(rA,k = k,main = main,results.dir = results.dir)
  R<-DIALOGUE2(rA = rA,main = main,results.dir = results.dir,add.effects = add.effects)
  R<-DIALOGUE3(rA = rA,main = main,results.dir = results.dir,full.version = full.version)
  if(plot.flag){
    DIALOGUE.plot(R,results.dir = results.dir)
  }
  return(R)
}

DIALOGUE1<-function(rA,k2 = 5,main,results.dir = "~/Desktop/DIALOGUE.results/"){
  print("#************DIALOGUE Step I: Canonical Correlation Analysis (CCA)************#")
  X<-lapply(rA, function(r){
    X1<-average.mat.rows(r@X,r@samples,f = colMedians)
    return(X1)
  })
  cell.types<-names(rA)
  names(X)<-cell.types
  n1<-length(cell.types)

  # Finding shared samples
  samples<-unlist(lapply(cell.types, function(x) rownames(X[[x]])))
  samplesU<-get.abundant(samples,n1)

  # Centering and scalling
  f<-function(X1){
    X1<-center.matrix(X1,dim = 2,sd.flag = T)
    X1<-cap.mat(X1,cap = 0.01,MARGIN = 2)
    X1<-X1[samplesU,]
    return(X1)
  }
  X<-lapply(X, f)
  k1<-ncol(X[[1]])
  set.seed(1234)
  perm.out <- MultiCCA.permute(X,type=rep("standard",length(X)),trace = F)
  out <- MultiCCA(X, type=rep("standard",length(X)),
                  penalty=perm.out$bestpenalties,niter = 100,
                  ncomponents=k, ws=perm.out$ws.init,trace = F)
  names(out$ws)<-names(X)
  for(i in names(X)){
    colnames(out$ws[[i]])<-paste0("C",1:ncol(out$ws[[i]]))
    rownames(out$ws[[i]])<-colnames(X[[i]])
  }
  Y<-lapply(names(X), function(i) X[[i]]%*%out$ws[[i]])
  names(Y)<-names(X)
  pairs1<-t(combn(names(X),2))
  cca.cor<-apply(pairs1,1,function(x) diag(cor(Y[[x[1]]],Y[[x[2]]])))
  colnames(cca.cor)<-paste(pairs1[,1],pairs1[,2],sep = "_")
  y<-list()
  if(missing(main)){
    main<-paste0(names(out$ws),collapse = "_")
  }
  
  R<-list(name = paste0("DIALOGUE1_",main),
          cell.types = cell.types,k = c(k,laply(X,ncol)),
          samples = samplesU,sample.PCs = X,
          cca = out,cca.cor = cca.cor,
          cca.scores = list(),cca.gene.cor = list(),
          cca.sig = list(),cca.redun.cor = list())
  names(R$k)<-c("DIALOGUE",paste0("original.",cell.types))

  for(x in cell.types){
    r<-rA[[x]]
    y[[x]]<-r@X[,1:k1]%*%out$ws[[x]]
    scores0<-as.matrix(y[[x]])
    conf<-r@cellQ
    if(!is.null(r@conf)){conf<-cbind.data.frame(r@cellQ,r@conf)}
    r@scores<-t(get.residuals(t(scores0),conf))
    r@scoresAv<-average.mat.rows(r@scores,r@samples,f = colMedians)
    R$cca.scores[[x]]<-r@scores
    R$cca.gene.cor[[x]]<-cor(t(r@tpm),r@scores)
    R$cca.sig[[x]]<-get.top.cor(R$cca.gene.cor[[x]],100,min.ci = 0.05)
    R$cca.redun.cor[[x]]<-cor(r@scores[,1:k2])
  }

  saveRDS(R,file = paste0(results.dir,"/",R$name,".rds"))
  dir.create(paste0(results.dir,"/DIALOGUE2_",main))
  return(R)
}

DIALOGUE2<-function(rA,main,results.dir = "~/Desktop/DIALOGUE.results/",add.effects = T){
  cell.types<-names(rA)
  if(missing(main)){main<-paste0(cell.types,collapse = "_")}
  file1<-paste0(results.dir,"/DIALOGUE1_",main,".rds")
  file2<-paste0(results.dir,"/DIALOGUE2_",main,".rds")

  R<-readRDS(file1)
  if(add.effects & !is.null(rA[[1]]@conf)){
    R$frm<-paste("y ~ (1 | samples) + x + cellQ +",
                 paste(colnames(rA[[1]]@conf),collapse = " +"))
  }else{
    R$frm<-"y ~ (1 | samples) + x + cellQ"
  }

  k2<-ncol(R$cca$ws[[1]])
  pairs1<-t(combn(cell.types,2))
  sig<-R$cca.sig
  for(i in 1:nrow(pairs1)){
    x1<-pairs1[i,1];x2<-pairs1[i,2]
    print(paste("#************DIALOGUE Step II (multilevel modeling):",x1,"vs.",x2,"************#"))
    p<-paste0(x1,".vs.",x2)
    R[[p]]<-DIALOGUE2.pair(R,rA[[x1]],rA[[x2]],cell.types,results.dir)
  }
  R$name<-paste0("DIALOGUE2_",main)
  saveRDS(R,file = file2)
  return(R)
}

DIALOGUE2.pair<-function(R,r1,r2,cell.types,results.dir){
  main<-gsub("DIALOGUE1_","",R$name)
  x1<-r1@name;x2<-r2@name
  saveFile<-paste0(results.dir,"/DIALOGUE2_",main,"/",x1,".vs.",x2,".rds")
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

  r1a<-cell.type.2.list(r1)
  r2a<-cell.type.2.list(r2)
  r1a<-set.list(r1a,sample.per.label(r1a$samples,50),sampleName = paste0(r1a$name,"_A"))
  r2a<-set.list(r2a,sample.per.label(r2a$samples,50),sampleName = paste0(r2a$name,"_A"))

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

  idx<-unique(get.strsplit(names(sig1),".",1))
  R1<-lapply(idx,function(x){f1(sig1,sig2,x)})
  names(R1)<-idx
  R1$p1<-NULL;R1$p2<-NULL
  for(x in idx){
    R1$p1<-rbind(R1$p1,R1[[x]]$p1)
    R1$p2<-rbind(R1$p2,R1[[x]]$p2)
  }
  R1$sig1<-lapply(R[idx], function(x) x$sig1f)
  R1$sig2<-lapply(R[idx], function(x) x$sig2f)
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

DIALOGUE3<-function(rA,main,results.dir = "~/Desktop/DIALOGUE.results/",full.version = F){
  print("#************DIALOGUE Step III: Finalizing the scores************#")
  cell.types<-names(rA)
  if(missing(main)){main<-paste0(cell.types,collapse = "_")}
  R<-readRDS(paste0(results.dir,"/DIALOGUE2_",main,".rds"))
  R$gene.pval<-lapply(R$cell.types, function(x) DLG.multi.get.gene.pval(x,R))
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
    r<-DLG.get.OE(r1,r2,plot.flag = T,compute.scores = F)
    r1<-r$r1;r2<-r$r2
    idx<-intersect(get.abundant(r1@samples),get.abundant(r2@samples))
    R$pref[[x]]<-cbind.data.frame(R = diag(cor(r1@scoresAv[idx,],r2@scoresAv[idx,])),
                                  hlm = DLG.hlm.pval(r1,r2,formula = R$frm))
  }

  R$gene.pval<-lapply(rA,function(r1) r1@gene.pval)
  R$sig<-lapply(rA,function(r1) r1@sig)
  R$scores<-lapply(rA,function(r1){
    X<-cbind.data.frame(r1@scores,samples = r1@samples,cells = r1@cells, cell.type = r1@name)
    if(!is.null(r1@conf)){X<-cbind.data.frame(X,conf = r1@conf)}
    return(X)})
  names(R$gene.pval)<-cell.types
  names(R$sig)<-cell.types
  names(R$scores)<-cell.types
  R$name<-paste0("DLG.output_",main)

  if(!full.version){
    file.remove(paste0(results.dir,"DIALOGUE1_",main,".rds"))
    file.remove(paste0(results.dir,"DIALOGUE2_",main,".rds"))
    unlink(paste0(results.dir,"DIALOGUE2_",main,"/"),recursive = T)
  }
  fileName<-paste0(results.dir,"DLG.full.output_",main,".rds")
  saveRDS(R,file = fileName)
  R<-R[c("sig","scores","gene.pval","pref","k","cell.types","name")]
  fileName<-paste0(results.dir,"DLG.output_",main,".rds")
  saveRDS(R,file = fileName)
  return(R)
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

DLG.multi.get.gene.pval<-function(x,R){
  b<-grepl("vs.",names(R))
  pairs1<-get.strsplit(names(R),".vs.",1:2)
  b1<-b&is.element(pairs1[,1],x)
  b2<-b&is.element(pairs1[,2],x)
  m1<-NULL;m2<-NULL
  if(sum(b1)>1){m1<-t(laply(R[b1],function(m) m$p1[,"Z"]))}
  if(sum(b2)>1){m2<-t(laply(R[b2],function(m) m$p2[,"Z"]))}
  if(sum(b1)==1){m1<-R[[which(b1)]]$p1[,"Z"]}
  if(sum(b2)==1){m2<-R[[which(b2)]]$p2[,"Z"]}
  m<-cbind(m1,m2)
  colnames(m)<-c(pairs1[b1,2],pairs1[b2,1])
  if(sum(b1)>0){p<-R[[which(b1)[1]]]$p1}else{p<-R[[which(b2)[1]]]$p2}
  rownames(m)<-rownames(p)
  m<-cbind.data.frame(m,p[,c("up","program","genes")],
                      p.up = fisher.combine(get.pval.from.zscores(m)),
                      p.down = fisher.combine(get.pval.from.zscores(-m)),
                      n.up = rowSums(m>(-log10(0.1))),
                      n.down = rowSums(m<log10(0.1)))
  m$p.up[!m$up]<-1
  m$p.down[m$up]<-1
  return(m)

}

DLG.find.scoring<-function(r1,R){
  gene.pval<-R$gene.pval[[r1@name]]
  gene.pval<-gene.pval[rowSums(is.na(gene.pval))==0,]
  gene.pval<-gene.pval[is.element(gene.pval$genes,r1@genes),]
  g<-sort(unique(gene.pval$genes))

  # r1@cca.scores0<-r1@X%*%R$cca$ws[[r1@name]]
  r1@extra.scores$cca0<-r1@X%*%R$cca$ws[[r1@name]]
  r1@zscores<-center.large.matrix(r1@tpm[g,],sd.flag = T)

  f<-function(x){
    y<-r1@extra.scores$cca0[,x]
    gene.pval<-gene.pval[gene.pval$program==x,]
    X<-t(r1@zscores[gene.pval$genes,])
    b<-is.element(colnames(X),gene.pval$genes[!gene.pval$up])
    X[,b]<-(-X[,b])
    gene.pval<-DLG.iterative.nnls(X,y,gene.pval)
    scores<-X%*%gene.pval$coef
    return(list(gene.pval = gene.pval,scores = scores))
  }
  m<-lapply(colnames(r1@extra.scores$cca0), f)
  # r1@scores0<-t(laply(m,function(x) x$scores))
  conf<-r1@cellQ
  if(!is.null(r1@conf)){conf<-cbind.data.frame(r1@cellQ,r1@conf)}
  r1@extra.scores$nnl0<-t(laply(m,function(x) x$scores))
  colnames(r1@extra.scores$nnl0)<-colnames(r1@extra.scores$cca0)
  r1@extra.scores$cca<-t(get.residuals(t(r1@extra.scores$cca0),conf))
  r1@scores<-t(get.residuals(t(r1@extra.scores$nnl0),conf))
  r1@scoresAv<-average.mat.rows(r1@scores,r1@samples)
  r1@gene.pval<-NULL
  for(x in m){
    r1@gene.pval<-rbind(r1@gene.pval,x$gene.pval)
  }
  m<-r1@gene.pval
  m<-m[m$coef>0|m$p.down<1e-3|m$p.up<1e-3,]
  m$program<-paste0(m$program,ifelse(m$up,".up",".down"))
  r1@sig<-split(m$genes,m$program)
  return(r1)
}

DLG.iterative.nnls<-function(X,y,gene.pval){
  set.seed(1234)
  f.rank<-gene.pval$n.up
  f.rank[!gene.pval$up]<-gene.pval$n.down[!gene.pval$up]
  y1<-y;y.fit<-rep(0,length(y))
  v<-list()
  gene.pval$coef<-0
  idx<-sort(get.abundant(f.rank),decreasing = T)
  idx<-setdiff(idx,0)
  idx<-setdiff(idx,"0")
  for(n1 in idx){
    X1<-X[,f.rank==n1]
    main<-paste0("N",n1)
    v[[main]]<-nnls::nnls(X1,y)
    y<-v[[main]]$residuals
    y.fit<-y.fit+v[[main]]$fitted
    gene.pval$coef[f.rank==n1]<-v[[main]]$x
    if(cor(y1,y.fit)>0.95){
      cor.plot(y.fit,y1,main = paste("NNLS fitting",n1))
      return(gene.pval)
    }
  }
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
  idx<-c("samples","scores",intersect(slotNames(r1),setdiff(gsub(" ","",get.strsplit(formula,"+ ",1:10)),NA)))
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


