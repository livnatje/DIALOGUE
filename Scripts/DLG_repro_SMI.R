# Reproduce the results and plots for SMI lung cancer data.

DLG.get.file<-function(file1){
  workdir<-"change.to.your.work.directory"
  return(paste0(workdir,file1))
}

DLG.set.wd<-function(){
  dir.create(DLG.get.file("Results/"))
  dir.create(DLG.get.file("Tables/"))
  dir.create(DLG.get.file("Data/"))
  dir.create(DLG.get.file("Figures/"))
  return()
}

DLG_LungCancer.SMI_run<-function(){
  rA<-readRDS(file = DLG.get.file("Data/DLG.input_LungCancer.SMI.rds"))
  param<-DLG.get.param(k = 3,seed1 = 1234,
                       frm = "y ~ (1 | fov:slides) + (1 | slides) + x + cellQ + tme.qc",
                       averaging.function = colMeans,
                       center.flag = T,extra.sparse = F,
                       conf = c("cellQ"),covar = c("cellQ","tme.qc"),
                       results.dir = DLG.get.file("/Results/"),
                       plot.flag = F,n.genes = 50,
                       PMD2 = F,spatial.flag = T)
  R<-DIALOGUE.run(rA = rA,main = "LungCancer.SMI",param = param)
  DLG_LungCancer.SMI_TableS1C(R)
  DLG_LungCancer.SMI_Figure3(rA = rA)
  return(R)
}

DLG_LungCancer.SMI_TableS1C<-function(R){
  sig1<-R$MCPs$MCP1
  sig2<-R$MCPs$MCP2
  sig3<-R$MCPs$MCP3
  sig2<-setdiff.lists.by.idx(sig2[c(4:6,1:3)],sig1)[c(4:6,1:3)]
  sig3<-setdiff.lists.by.idx(sig3,sig1)
  X<-rbind(cbind(melt(sig1),MCP = "MCP1"),
           cbind(melt(sig2),MCP = "MCP2"),
           cbind(melt(sig3),MCP = "MCP3"))
  colnames(X)<-c("Genes","Compartment","MCP")
  write.csv(X,file = DLG.get.file("Tables/TableS1C.csv"),row.names = F)
}

DLG_LungCancer.SMI_Figure3<-function(rA,R){
  if(missing(rA)){
    rA<-readRDS(file = DLG.get.file("Data/DLG.input_LungCancer.SMI.rds"))
  }
  if(missing(R)){
    R<-readRDS(DLG.get.file("Results/DIALOGUE1_LungCancer.SMI.rds"))
  }
  
  R$scores<-lapply(R$cell.types,function(x) cbind.data.frame(R$cca.scores[[x]],rA[[x]]@metadata))
  names(R$scores)<-R$cell.types
  
  print("Reproducing Figure 3 from (Jerby and Regev, NBT 2022).")
  p1<-DLG_LungCancer.SMI.subplots(R,slide1 = 'Lung13',n.plots = 4,q1 = 0.9,MCPs = 1:2,both.sides  = T,
                                  fileName = DLG.get.file("/Figures/Figure3_Lung13.pdf"))
  p2<-DLG_LungCancer.SMI.subplots(R,slide1 = 'Lung9_Rep1',n.plots = 4,q1 = 0.9,MCPs = 1:2,both.sides  = F,
                                  fileName = DLG.get.file("/Figures/Figure3_Lung9_Rep1.pdf"))
  return(list(p1,p2))
}

DLG_LungCancer.SMI.subplots<-function(R,slide1 = "Lung5_Rep3",cex = 0.3,n.plots = 1,both.sides = F,
                                      fileName,q1 = 0.8,MCPs = seq(R$k[1],1,-1),r1){
  if(!missing(fileName)){pdf(fileName)}
  if(n.plots == 4){par(mfrow=c(2,2),oma = c(0, 0, 0, 0),xpd = F)}
  
  b<-lapply(R$scores,function(X) return(X$slides==slide1))
  scores<-lapply(names(R$scores), function(x) {return(R$scores[[x]][b[[x]],])})
  names(scores)<-names(R$scores)
  scoresB.up<-lapply(scores, function(X) apply(X[,1:R$k[1]],2,function(x){return(x>quantile(x,q1))}))
  scoresB.down<-lapply(scores, function(X) apply(X[,1:R$k[1]],2,function(x){return(x<quantile(x,1-q1))}))
  
  col1<-c("blue","red","green","orange")
  col2<-c("lightblue","bisque","grey","beige")
  names(col1)<-names(R$scores)
  names(col2)<-names(R$scores)
  coor<-NULL
  scoresBA.up<-NULL
  scoresCol.up<-NULL
  scoresBA.down<-NULL
  scoresCol.down<-NULL
  cell.types<-NULL
  scoresA<-NULL
  
  for(x in names(R$scores)){
    scoresA<-rbind(scoresA,center.matrix(scores[[x]][,MCPs],dim = 2,sd.flag = T))
    coor1<-scores[[x]][,c("coor.X","coor.Y")]
    coor<-rbind(coor,coor1)
    scoresBA.up<-rbind(scoresBA.up,scoresB.up[[x]])
    scoresBA.down<-rbind(scoresBA.down,scoresB.down[[x]])
    
    scoresCol.up<-rbind(scoresCol.up,apply(scoresB.up[[x]],2,function(v) ifelse(v,col1[x],col2[x])))
    scoresCol.down<-rbind(scoresCol.down,apply(scoresB.down[[x]],2,function(v) ifelse(v,col1[x],col2[x])))
    
    cell.types<-c(cell.types,rep(x,nrow(coor1)))
  }
  
  pch1<-c(16,21,22)[match(cell.types,unique(cell.types))]
  f1<-function(k1,scoresBA,scoresCol,str){
    call.plot.plus(coor,labels = paste0(cell.types,", ",ifelse(scoresBA[,k1],"High","Moderate/low")),
                   cex = cex+(0.05*scoresBA[,k1]),
                   b.top = scoresBA[,k1],
                   my.col = scoresCol[,k1],main = paste0(k1,", ",str),
                   xlab = "x-coor.",ylab = "y-coor.",add.N = T)
  }
  f2<-function(k1){
    par(mfrow=c(2,2),oma = c(0, 0, 0, 0),xpd = F)
    for(x in unique(cell.types)){
      b<-cell.types==x
      call.plot(coor[b,],labels = scoresAcap[b,k1],pch = 16,
                cex = 0.3,xlab = "x-coor.",ylab = "y-coor.",
                main = paste0("MCP",k1,": ",x))
    }
  }
  
  lapply(MCPs, function(x) f1(x,scoresBA = scoresBA.up,scoresCol = scoresCol.up,str = "UP"))
  if(both.sides){
    lapply(MCPs, function(x) f1(x,scoresBA = scoresBA.down,scoresCol = scoresCol.down,str = "DOWN"))
  }
  
  
  scoresAcap<-cap.mat(scoresA,cap = 0.01,MARGIN = 2)
  call.plot.multilabels(coor,labels = scoresAcap,cex = 0.3,main = "All cell types")
  lapply(MCPs, f2)
  
  
  if(!missing(fileName)){dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)}
  return(cbind.data.frame(coor = coor,scoresA = scoresA,up = scoresCol.up,down = scoresCol.down,cell.types = cell.types))
}


