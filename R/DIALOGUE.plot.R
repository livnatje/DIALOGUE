#' DIALOGUE.plot
#'
#' Plot DIALOGUE results.
#' The resulting plots will show the cell-type-specific components
#' of each multicellular program as a function of its other components;
#' and the composition of each multicellular program, depicting how many of its genes
#' are cell-type-specific, and how many are shared across multiple cell types.
#' @param rA list of cell type objects (see cell.type)
#' @examples
#' R<-DIALOGUE.run(rA,results.dir = "~/Desktop/Results/")
#' R<-DIALOGE.plot(R,results.dir = "~/Desktop/Figures/")
#' Alternatively
#' R<-DIALOGUE.run(rA,results.dir = "~/Desktop/Results/",plot.flag = T)
#' @author Livnat Jerby-Arnon
#' @export
#'

DIALOGUE.plot<-function(R,results.dir = "~/Desktop/DIALOGUE.results/",pheno = NULL){
  pdf(paste0(results.dir,"/",R$name,".pdf"))
  DIALOGUE.plot.av(R)
  DIALOGUE.plot.sig.comp(R)
  if(!is.null(pheno)){
    DIALOGUE.violin.pheno(R,pheno = pheno)
  }
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
}

DIALOGUE.plot.av<-function(R,i,mark.samples = NULL,d = 1){
  k<-length(R$MCPs)
  R$cell.types<-names(R$gene.pval)
  R$scoresAv<-lapply(R$scores,function(m) average.mat.rows(as.matrix(m[,1:k]),m$samples,f = colMedians))
  idx<-unlist(lapply(R$scoresAv, function(m) rownames(m)))
  idx<-get.abundant(idx,length(R$cell.types))
  R$scoresAv<-lapply(R$scoresAv,function(m) m[idx,])
  col<-rep("black",length(idx))
  col[is.element(idx,mark.samples)]<-"red"
  pch<-ifelse(any(col=="red"),21,16)
  f<-function(i){
    m1<-t(laply(R$scoresAv,function(m) m[,i]))*d
    rownames(m1)<-idx;colnames(m1)<-R$cell.types
    pairs.panels(m1,hist.col = "grey",breaks = 50,bg = col,pch = pch,ellipses = F,smooth = T,lm = T)
    title(i)
    return(m1)
  }
  idx1<-paste0("C",1:R$k["DIALOGUE"])
  if(!missing(i)){m1<-f(i);return(cor(m1))}
  m<-lapply(idx1, f)
}

DIALOGUE.plot.sig.comp<-function(R,main = ""){
  genes<-unique(sort(unlist(R$MCPs)))
  f<-function(sig1,d = 1){
    if(d==1){
      b<-!grepl(".down",names(sig1))
    }else{
      b<-grepl(".down",names(sig1))
    }
    if(!any(b)){return(rep(length(genes),0))}
    sig1<-sig1[b]
    names(sig1)<-get.strsplit(names(sig1),".",1)
    v<-list.2.ids(genes,sig1)
    return(c(v))
  }
  m1<-t(laply(R$MCPs,f))
  m2<-t(laply(R$MCPs,function(x) f(x,-1)))

  colnames(m1)<-paste0(names(R$MCPs),".up")
  colnames(m2)<-paste0(names(R$MCPs),".down")
  m<-cbind(m1,m2)
  idx<-R$cell.types
  idxA<-idx
  for(i in 2:length(idx)){
    idxA<-c(idxA,apply(combn(idx,i),2,function(x) paste(x,collapse = "&")))
  }
  idxA
  m1<-laply(idxA,function(x) colSums(m==x))
  rownames(m1)<-idxA
  b<-stri_count(str = rownames(m1),regex = "&")>1
  m2<-rbind(m1[!b,],colSums(subset.matrix(m1,b)))
  rownames(m2)[nrow(m2)]<-"> 2 cell types"
  if(all(m2["> 2 cell types",]==0)){m2<-m2[1:(nrow(m2)-1),]}
  m1<-melt(m2)
  colnames(m1)<-c("col","x","y")
  p<-ggplot(data=m1, aes(x=x, y=y, fill=col))+geom_bar(stat="identity")+
    labs(fill = "Cell type(s)", x = "Program", y = "No. of genes")
  p<-p+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black"))
  p<-p+theme(text = element_text(size=12),axis.text.x = element_text(angle=45, hjust=1))
  multiplot.util(list(NULL,p,NULL),cols = 1,nplots = 3)
  return(m2)
}

DIALOGUE.violin.pheno<-function(R,pheno = "clinical.status"){
  k<-R$k["DIALOGUE"]
  X<-NULL
  for(x in R$scores){
    x[,1:k]<-cap.mat(center.matrix(x[,1:k],dim = 2,sd.flag = T),cap = 0.01,MARGIN = 2)
    X<-rbind(X,x)
  }

  par(mfrow=c(2,1),oma = c(5, 0, 0, 7))
  laply(1:k,function(x){
    violin.split(scores = X[,x],treatment = X[,pheno],
                 conditions = X$cell.type,
                 main = paste0("MCP",x))
    return(x)
  })

  return()

}

multiplot.util<-function(plotlist,nplots = 4,cols = 2){
  flag<-F
  while(!is.null(plotlist)&!flag){
    print(multiplot(plotlist = plotlist[1:min(nplots,length(plotlist))],cols = cols))
    flag<-(min(nplots,length(plotlist))+1)>length(plotlist)
    plotlist<-plotlist[(min(nplots,length(plotlist))+1):length(plotlist)]
  }
}

list.2.ids<-function(ids,l,single.flag = F){
  if(length(l)==1){
    m<-ifelse(is.element(ids,l[[1]]),names(l),"")
    return(m)
  }
  B<-t(laply(l,function(x) is.element(ids,x)))
  colnames(B)<-names(l)
  m<-get.mat(ids,"Anno")
  m[]<-"ID"
  for(i in names(l)){
    m[B[,i]]<-paste(m[B[,i]],i,sep = "&")
  }
  m<-gsub("ID&","",m)
  m<-gsub("ID","",m)
  if(single.flag){
    m<-cbind(m,get.strsplit(m,"&",1))
    m[m[,1]=="",2]<-""
  }
  return(m)
}

get.mat<-function(m.rows,m.cols,data = NA){
  m<-matrix(data = data, nrow = length(m.rows),ncol = length(m.cols),
            dimnames = list(m.rows,m.cols))
  return(m)
}

multiplot<-function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

violin.split<-function(scores, treatment, conditions, main = "",xlab = "",ylab = "Scores"){
  # require(beanplot)
  treatment<-as.factor(treatment)
  beanplot(scores ~ treatment*conditions, ll = 0.0,las = 2,
           main = main, side = "both", xlab=xlab,ylab = ylab,
           col = list(c("lightblue", "black"),"gray"),
           axes=T,cex.main = 1)
  legend("bottomright", fill = c("lightblue","gray"),
         legend = levels(treatment), box.lty=0)
  return()
}

