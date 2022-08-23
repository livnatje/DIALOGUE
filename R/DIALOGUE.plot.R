#' DIALOGUE.plot
#'
#' Plot DIALOGUE results.
#' The resulting plots will show the cell-type-specific components
#' of each multicellular program as a function of its other components;
#' and the composition of each multicellular program, depicting how many of its genes
#' are cell-type-specific, and how many are shared across multiple cell types.
#'
#' In case there is a specific feature or phenotype of interest DIALOGUE will plot the MCP
#' expression when stratifying the cells according to that classification.
#'
#' @param R DIALOGUE output
#' @param results.dir the directory where the PDF figure file will be located
#' @param pheno (optional) the name of a binary feature of interest to visualize in relation to the MCPs.
#'
#' @examples
#' # Run DIALOGUE
#' R<-DIALOGUE.run(rA,results.dir = "~/Desktop/Results/",plot.flag = F)
#' # Plot the results
#' DIALOGE.plot(R,results.dir = "~/Desktop/Figures/",pheno = R$pheno)
#' # Alternatively
#' R<-DIALOGUE.run(rA,results.dir = "~/Desktop/Results/",plot.flag = T)
#'
#' @author Livnat Jerby-Arnon
#' @export
#'

DIALOGUE.plot<-function(R,results.dir = "~/Desktop/DIALOGUE.results/",
                        pheno = NULL,mark.samples = NULL,metadata = NULL,d = 1, MCPs = 1:R$k["DIALOGUE"]){
  
  pdf(paste0(results.dir,"/",R$name,".pdf"))
  DIALOGUE.plot.av(R,mark.samples = mark.samples,metadata = metadata,d = d,MCPs = MCPs)
  DIALOGUE.plot.sig.comp(R)
  if(!is.null(pheno)){
    DIALOGUE.violin.pheno(R,pheno = pheno,MCPs = MCPs,d = d)
  }
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
}

DIALOGUE.plot.av<-function(R,MCPs,mark.samples = NULL,d = 1,k = R$k["DIALOGUE"],
                           select.samples = NULL,averaging.function = colMeans,metadata = NULL){
  R$scoresAv<-lapply(R$scores,function(m) average.mat.rows(as.matrix(m[,1:k]),
                                                           m$samples,f = R$param$averaging.function))
  idx<-unlist(lapply(R$scoresAv, function(m) rownames(m)))
  idx<-get.abundant(idx,length(R$cell.types))
  if(!is.null(select.samples)){
    idx<-intersect(idx,select.samples)
  }
  R$scoresAv<-lapply(R$scoresAv,function(m) m[idx,])
  if(!is.null(metadata)){
    metadata<-metadata[idx,]
    col1<-metadata$col
    if(length(unique(metadata$size))){size1<-1}else{
      size1<-1.5*metadata$size/max(metadata$size)
      size1[size1<1]<-1
    }
  }else{
    col1<-rep("black",length(idx))
    col1[is.element(idx,mark.samples)]<-"red"
    size1<-1
  }
  pch<-ifelse(any(col1=="red"),21,16)
  f<-function(i){
    m1<-t(laply(R$scoresAv,function(m) m[,i]))*d
    # rownames(m1)<-idx
    colnames(m1)<-R$cell.types
    m1<-m1[,R$MCP.cell.types[[i]]]
    if(length(R$MCP.cell.types[[i]])<2){return()}
    pairs.panels(m1,hist.col = "grey",breaks = 50,bg = col1,pch = pch,ellipses = F,
                 smooth = F,lm = T,stars = T,method = "pearson",
                 cex = size1,cex.cor = 1)
    title(i)
    return(m1)
  }
  if(missing(MCPs)){
    MCPs<-paste0("MCP",1:R$k["DIALOGUE"])
  }
  m<-lapply(MCPs, f)
}

DIALOGUE.plot.sig.comp<-function(R,main = ""){
  R$MCPs<-R$MCPs[laply(R$MCPs,length)>1]
  genes<-unique(sort(unlist(R$MCPs)))
  f<-function(sig1,d = 1){
    if(d==1){
      b<-!grepl(".down",names(sig1))
    }else{
      b<-grepl(".down",names(sig1))
    }
    if(!any(b)){return(rep("",length(genes)))}
    sig1<-sig1[b]
    names(sig1)<-gsub(".up","",names(sig1))
    names(sig1)<-gsub(".down","",names(sig1))
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

DIALOGUE.violin.pheno<-function(R,pheno = "pathology",MCPs,selected.samples,d = 1){
  k<-R$k["DIALOGUE"]
  X<-NULL
  for(x in R$scores){
    x[,1:k]<-cap.mat(center.matrix(x[,1:k],dim = 2,sd.flag = T),cap = 0.01,MARGIN = 2)
    X<-rbind(X,x)
    X<-X[!is.na(X[,pheno]),]
  }
  if(!is.element("id",colnames(X))){X$id<-X$cell.type}
  # if(is.logical(X[,pheno])){X[,pheno]<-ifelse(X[,pheno],"Disease","Control")}
  if(!missing(selected.samples)){X<-X[is.element(X$samples,selected.samples),]}
  
  par(mfrow=c(1,1),oma = c(5, 0, 0, 7))
  f<-function(x){
    b<-is.element(X$cell.type,R$MCP.cell.types[[x]])
    violin.split(scores = d*X[b,x],treatment = X[b,pheno],
                 conditions = X$id[b],
                 main = x)
    return(x)
  }
  if(missing(MCPs)){MCPs<-paste0("MCP",1:k)}
  laply(MCPs,f)
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

violin.split<-function(scores, treatment, conditions, main = "",
                       xlab = "Sample",ylab = "Scores",legend.flag = T,show.pval = T){
  # require(beanplot)
  if(length(unique(conditions))==1){
    p<-t.test.mat(m = rbind(scores,scores),b = treatment == sort(treatment,decreasing = T)[1])[1,1]
  }else{
    p<-t.test.groups(x = rbind(scores,scores),b = treatment == sort(treatment,decreasing = T)[1],g = conditions)[1,]
    p<-p[sort(names(p))]
  }
  p[p<(-30)]<-(-30);p[p>30]<-30
  if(show.pval){
    conditions<-paste0(conditions,"\n",laply(10^-abs(p[conditions]),my.format.pval))
  }
  treatment<-as.factor(treatment)
  beanplot(scores ~ treatment*conditions, ll = 0.0,las = 2,
           main = main, side = "both", xlab=xlab,ylab = ylab,
           col = list(c("lightblue", "black"),"gray"),
           axes=T,cex.main = 1)
  if(legend.flag){
    legend("bottomright", fill = c("lightblue","gray"),
           legend = levels(treatment), box.lty=0)
  }
  return(p)
}

