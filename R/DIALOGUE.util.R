average.mat.rows<-function(m,ids,f = colMeans){
  ids.u<-sort(get.abundant(ids))
  m1<-laply(ids.u,function(x){return(f(m[is.element(ids,x),]))})
  rownames(m1)<-ids.u
  colnames(m1)<-colnames(m)
  
  ids.u1<-setdiff(unique(ids),ids.u)
  if(length(ids.u1)==0){return(m1)}
  b<-is.element(ids,ids.u1)
  m0<-m[b,]
  if(sum(b)==1){m0<-t(as.matrix(m0))}
  rownames(m0)<-ids[b]
  
  m2<-rbind(m1,m0)
  m2<-m2[sort(rownames(m2)),]
  return(m2)
}

prep4semiOE<-function(r,geneAv = NULL){
  paste("Preparing",r@name,"for OE computations.")
  r@zscores<-center.large.matrix(r@tpm,sd.flag = T,v = geneAv)
  r<-log.comp(r)
  return(r)
}

center.large.matrix<-function(m,sd.flag,v = NULL){
  if(is.null(v)){
    v<-rowMeans(m,na.rm = T)
  }
  if(!sd.flag){
    for(i in 1:nrow(m)){
      m[i,]<-m[i,]-v[i]
    }
  }else{
    for(i in 1:nrow(m)){
      x<-m[i,]
      m[i,]<-(x-v[i])/sd(x,na.rm = T)
    }
  }
  return(m)
}

get.abundant<-function(v,abn.c = 2,boolean.flag = F,top,decreasing = T){
  m<-as.matrix(table(v))
  m<-as.matrix(m[order(m,decreasing = decreasing),])
  if(!missing(top)){
    abn.c<-m[top]
  }
  m<-m[m>=abn.c,]
  abn.names<-names(m)
  if(boolean.flag){
    b<-is.element(v,abn.names)
    return(b)
  }
  return(abn.names)
}

center.matrix<-function(m,dim = 1,sd.flag = F){
  if(dim == 1){
    zscores<-sweep(m,1,rowMeans(m,na.rm = T),FUN = '-')
  }else{
    zscores<-sweep(m,2,colMeans(m,na.rm = T),FUN = '-')
  }
  if(sd.flag){
    zscores<-sweep(zscores,dim,apply(m,dim,function(x) (sd(x,na.rm = T))),FUN = '/')
  }
  return(zscores)
}

cap.mat<-function(M,cap = 0.01,MARGIN = 1){
  Z<-apply(M,MARGIN = MARGIN,function(x){
    q9<-quantile(x,1-cap)
    q1<-quantile(x,cap)
    x[x>q9]<-q9;x[x<q1]<-q1
    return(x)
  })
  if(MARGIN==1){Z<-t(Z)}
  return(Z)
}

get.residuals<-function(X,g,MARGIN = 1){
  if(MARGIN == 2){return(t(get.residuals(t(X),g,MARGIN = 1)))}
  f<-function(y){return(lm(y~.,data = as.data.frame(g))$residuals)}
  residuals<-t(apply(X,1,f))
  rownames(residuals)<-rownames(X)
  colnames(residuals)<-colnames(X)
  return(residuals)
}

get.top.cor<-function(m,q = 100,min.ci = 0,idx = NULL, add.prefix =""){
  m<-as.matrix(m)
  if(is.null(colnames(m))){colnames(m)<-1:ncol(m)}
  m.pos<-(-m);m.neg<-m
  
  colnames(m.pos)<-paste0(colnames(m.pos),".up")
  colnames(m.neg)<-paste0(colnames(m.neg),".down")
  v<-get.top.elements(cbind(m.pos,m.neg),q,min.ci = (-abs(min.ci)))
  names(v)<-c(colnames(m.pos),colnames(m.neg))
  if(!is.null(idx)){
    v<-v[paste(idx,c("up","down"),sep = ".")]
  }
  names(v)<-paste0(add.prefix,names(v))
  return(v)
}

get.top.elements<-function (m,q = 100,min.ci = NULL,main = ""){
  top.l<-list()
  v<-rownames(m)
  for (i in 1:ncol(m)){
    mi<-m[,i];mi<-mi[!is.na(mi)]
    idx<-order(mi,decreasing = F)
    ci <- mi[idx[min(q,length(mi))]]
    ci <- min(ci,min.ci)
    b <- m[,i]<=ci
    b[is.na(m[,i])]<-F
    top.l[[i]]<-sort(v[b])
  }
  if(main!=""){main<-paste0(main,".")}
  names(top.l)<-paste0(main,colnames(m))
  return(top.l)
}

set.list<-function (r,b,sampleName = NULL){
  rn<-lapply(r, set.field, b = b)
  if(!is.null(sampleName)){
    rn$sampleName<-sampleName
  }
  return(rn)
}

set.cell.type<-function(r,b,name = r@name){
  rn<-r
  for(x in slotNames(r)){
    slot(rn,x)<-set.field(v = slot(r,x),b = b)
  }
  return(rn)
}

set.field<-function (v,b){
  d <- dim(v)
  d.b<-length(b)
  if(!is.null(d)){
    if(d[1]==d.b){v <- subset(v,subset = b)}
    if(d[2]==d.b){v <- v[,b]}
  }else{if(length(v)==d.b){v <- v[b]}}
  return(v)
}

sample.per.label<-function(labels,size,boolean.flag = T,v,remove.flag = F){
  if(missing(v)){
    v<-1:length(labels)
  }
  ul<-unique(labels)
  vr<-NULL
  for(i in 1:length(ul)){
    b<-is.element(labels,ul[i])
    if(size<1){
      vr<-c(vr,sample(v[b],size = sum(b)*size))
    }else{
      vr<-c(vr,sample(v[b],size = min(size,sum(b))))
    }
  }
  b<-is.element(v,vr)
  if(remove.flag){
    b.small<-!is.element(labels,get.abundant(labels[b],size-1))
    b[b.small]<-F
    vr<-v[b]
  }
  if(boolean.flag){return(b)}
  return(vr)
}

get.strsplit<-function(v,sep,idx){
  v<-as.character(v)
  vi<-laply(strsplit(v,split = sep,fixed = T),function(x) x[idx])
  return(vi)
}

apply.formula.HLM<-function(r,X,Y,MARGIN = 1,formula = "y ~ (1 | samples) + x",ttest.flag = F){
  if(is.matrix(Y)){
    if(ttest.flag){
      m1<-t.test.mat(Y,X)
      b<-rowSums(p.adjust.mat(m1[,1:2])<0.1,na.rm = T)>0
      m1<-m1[b,];Y<-Y[b,]
    }
    m<-t(apply(Y,MARGIN = MARGIN,function(y){formula.HLM(y,X,r,formula = formula)}))
  }else{
    m<-t(apply(X,MARGIN = MARGIN,function(x){formula.HLM(Y,x,r,formula = formula)}))
  }
  colnames(m)<-c("Estimate","P")
  m<-cbind.data.frame(Z = get.cor.zscores(m[,"Estimate"],m[,"P"]),m)
  if(ttest.flag){
    m<-cbind.data.frame(m,ttest = m1)
  }
  return(m)
}

formula.HLM<-function(y,x,r0, formula = "y ~ (1 | samples) + x",val = ifelse(is.numeric(x),"","TRUE"),return.all = F){
  r0$x<-x;r0$y<-y
  f<-function(r0){
    M1 <- with(r0, lmer (formula = formula))
    if(return.all){
      c1<-summary(M1)$coef[,c("Estimate","Pr(>|t|)")]
    }else{
      c1<-summary(M1)$coef[paste0("x",val),]
      idx<-match(c("Estimate","Pr(>|t|)"),names(c1))
      c1<-c1[idx]
    }
    return(c1)
  }
  c1<-tryCatch({f(r0)},
               error = function(err){return(c(NA,NA))})
  return(c1)
}

get.cor.zscores<-function(c,p){
  v<-cbind(get.onesided.p.value(c,p),get.onesided.p.value(-c,p))
  z<-get.p.zscores(v)
  return(z)
}

get.onesided.p.value<-function(c,p){
  p[p==0] = min(p[p>0],na.rm = T)
  p.one.side <- p
  p.one.side[] <- NA
  b<-c>0&!is.na(c)
  p.one.side[b]=p[b]/2
  b<-c<=0&!is.na(c)
  p.one.side[b]=1-(p[b]/2)
  return(p.one.side)
}

get.p.zscores<-function(p){
  b<-p[,1]>0.5
  b[is.na(b)]<-F
  zscores<-(-log10(p[,1]))
  zscores[b]<-log10(p[b,2])
  # signficiant in p[,1] will be positive
  # signficiant in p[,2] will be negative
  return(zscores)
}

intersect.list1<-function(l,g,n1=0,HG.universe = NULL,prf = ""){
  l1<-lapply(l, function(x) intersect(x,g))
  l1<-l1[laply(l1,length)>n1]
  if(prf!=""){
    names(l1)<-paste(prf,names(l1),sep = ".")
  }
  if(!is.null(HG.universe)){
    p<-GO.enrichment.lapply(l[names(l1)],genes = HG.universe,list(g))
    names(l1)<-paste(names(l1),format(p,scientific = T,digits= 3),sep = "P = ")
  }
  return(l1)
}

fisher.combine<-function(p){
  p.f<-apply(p,1,get.fisher.p.value)
  return(p.f)
}

get.fisher.p.value<-function(p){
  p<-p[!is.na(p)]
  if(length(p)==1){
    p.fisher=p
  }else{
    p.fisher<- 1 - pchisq(-2*sum(log(p),na.rm = T), 2*sum(!is.na(p)))
  }
  return(p.fisher)
}

get.pval.from.zscores<-function(z){
  p<-10^(-abs(z))
  b<-z<0&!is.na(z)
  p[b]<-1-p[b]
  return(p)
}

cor.plot<-function(x,y = NULL,main = '',ylab = '', xlab = '',regression.flag = F,cex = 0.3,
                   xlim = NULL,ylim = NULL){
  if(is.null(y)){
    v<-colnames(x)
    xlab<-v[1];ylab<-v[2]
    y<-x[,2];x<-x[,1]
  }
  v<-spearman.cor(x,y)
  main <- paste(main,"\nR =",format(v[1],digits = 2),"P =",format(v[2],scientific = T,digits = 2))
  plot(x,y,main = main, xlab = xlab, ylab = ylab,cex = cex,pch=16,xlim = xlim,ylim = ylim)
  b<-!is.na(x)&!is.na(y)
  v<-lowess(x[b],y[b])
  lines(v,col = "red")
  if(!regression.flag){return()}
  y.d<-y-v$y[match(x,v$x)]
  y.sd<-sd(y.d,na.rm = T)
  y.av<-mean(y.d,na.rm = T)
  labels<-matrix(data = "Moderate",nrow = length(y))
  labels[y.d>(y.av+y.sd)]<-"High"
  labels[y.d<(y.av-y.sd)]<-"Low"
  my.plot(x,y,labels = labels,main = main,xlab = xlab,ylab = ylab)
  lines(v)
  return(y.d)
}

spearman.cor<-function(v1,v2 = NULL,method = 'spearman',use = "pairwise.complete.obs",
                       match.flag = F,alternative = "two.sided",upper.tri.flag = F){
  if(is.null(v2)){
    v2<-v1
  }
  if(!is.matrix(v1)){v1<-as.matrix(v1)}
  if(!is.matrix(v2)){v2<-as.matrix(v2)}
  if(match.flag){
    n=ncol(v1)
    if(is.null(colnames(v1))){colnames(v1)<-1:ncol(v1)}
    results<-get.mat(m.cols = c("R","P"),m.rows = colnames(v1))
    for(i in 1:ncol(v1)){
      c.i <- cor.test(v1[,i],v2[,i],method = method,use = use, alternative = alternative)
      results[i,1] <- c.i$estimate
      results[i,2] <- c.i$p.value
    }
  }else{
    n1=ncol(v1)
    m<-matrix(nrow = n1,ncol = ncol(v2))
    rownames(m)<-colnames(v1)
    colnames(m)<-colnames(v2)
    results<-list(cor = m, p = m)
    for(i in 1:n1){
      f<-function(x){
        c.i<-cor.test(v1[,i],x,method = method,use = use, alternative = alternative);
        c(c.i$estimate,c.i$p.value)}
      c.i <- apply(v2,2,f)
      results$cor[i,] <- c.i[1,]
      results$p[i,] <- c.i[2,]
    }
    if(ncol(v2)==1){
      results<-cbind(results$cor,results$p)
      colnames(results)<-c('R','P')
    }
  }
  if(upper.tri.flag){
    results$up <- cbind(results$cor[upper.tri(results$cor)],
                        results$p[upper.tri(results$p)])
  }
  return(results)
}