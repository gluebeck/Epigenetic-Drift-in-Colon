#### gene expression RSEM v2 data downloaded from GDC Portal 
#### analyses for static or drift related genes
#### drift: Genes.isl.drift5 (part2)
#### static: Genes.isl.static5

#### download expr file (rda on GitHub)
load('./data/GDC/integrative/expr_coadread.rda',verbose=T)
expr = as.data.frame(expr)
class(expr) = "numeric"

cond.normal = grepl("-11A-",colnames(expr))
cond.cancer = grepl("-01A-",colnames(expr))

expr.normal = expr[,cond.normal]
expr.cancer = expr[,cond.cancer]

tmp = rownames(expr.normal)
genes.expr = sapply(strsplit(tmp,"|",fixed=T), `[`, 1)  # strip of Entrez gene nom

bcr.expr.normal = substr(colnames(expr.normal),1,12)
bcr.expr.cancer = substr(colnames(expr.cancer),1,12)

## choose data
dat = cbind(out.tcga.left$Mvals[iset,])
geneList = out.tcga.left$genes[iset]  
islandList = out.tcga.left$islands[iset]

dat = cbind(out.tcga.right$Mvals[iset,])
geneList = out.tcga.right$genes[iset]  
islandList = out.tcga.right$islands[iset]

dat = cbind(static.tcga.left$Mvals)
geneList = static.tcga.left$genes  
islandList = static.tcga.left$islands

dat = cbind(static.tcga.right$Mvals)
geneList = static.tcga.right$genes  
islandList = static.tcga.right$islands

# dat = cbind(out.tcga.norm$Mvals)
# geneList = out.tcga.norm$genes  
# islandList = out.tcga.norm$islands

#### map arrayIDs to barcode
tmp=colnames(mswan)  #minfi mswan for TCGA
bcr = pdata.tcga$bcr_patient_barcode[match(tmp,pdata.tcga$arrayID)]
ids = na.omit(match(bcr.expr.normal,bcr))
bcr.tcga.norm = bcr[ids]

M = nrow(dat) #number of CGI to test for gex data
ids = colnames(dat)
bcr = pdata.tcga$bcr_patient_barcode[match(ids,pdata.tcga$arrayID)]

ids.gex.normal = match(bcr,bcr.expr.normal)
ids.gex.cancer = match(bcr,bcr.expr.cancer)
ids.mex.normal = match(bcr,bcr.expr.normal)

#### put it together, entry by entry
glab = rlab = dat.mex = dat.gex.normal = dat.gex.cancer = NULL
for ( m in 1:M) {
  mex = dat[m,]
  mil = islandList[m]
  genes  = geneList[[m]]; len = length(genes)
  if(len>0) {
    for(l in 1:len) {
      gene = genes[l]
      if(gene%in%genes.expr) {
        # find gene expression for gene
        gex.cancer  = expr.cancer[match(gene,genes.expr),ids.gex.cancer]
        gex.normal  = expr.normal[match(gene,genes.expr),ids.gex.normal]
        # pair them up
        dat.mex = rbind(dat.mex,mex)
        dat.gex.cancer = rbind(dat.gex.cancer,gex.cancer)
        dat.gex.normal = rbind(dat.gex.normal,gex.normal)
        glab = c(glab,gene)
        rlab = c(rlab,paste(mil,gene))
      }
    }
  }
}
rownames(dat.mex) = rownames(dat.gex.normal) = rownames(dat.gex.cancer) = rlab

# dat.gex = as.matrix(dat.gex)
G = nrow(dat.gex.cancer)

#### isl-gene correlations (using log(G+1) transform for RSEM)
corr = pval = ttpval = gex.normal = gex.cancer = numeric(G)
for (i in 1:G) {
  tmp = cor.test(dat.mex[i,],log(dat.gex.cancer[i,]+1))

  # dum = quantile(dat.gex.normal[i,],probs=c(.2,.8),na.rm=T)
  # gex.normal[i] = dum[2]/dum[1]
  gex.normal[i] = mean(dat.gex.normal[i,]+.1,na.rm=T)

  # dum = quantile(dat.gex.cancer[i,],probs=c(.2,.8),na.rm=T)
  # gex.cancer[i] = dum[2]/dum[1]
  gex.cancer[i] = mean(dat.gex.cancer[i,]+.1,na.rm=T)
  # tmp = cor.test(dat.mex[sample(1:G,1),],log(dat.gex[i,]+1))
  # tmp = cor.test(static.tcga.right$Mvals[sample(1:1016,1),],log(dat.gex[sample(1:G,1),]+1))
  # tmp = cor.test(static.tcga.left$Mvals[sample(1:1016,1),],log(dat.gex[sample(1:G,1),]+1))
  # tmp = cor.test(out.tcga.right$Mvals[sample(1:871,1),],log(dat.gex[i,]+1))
  
  corr[i] = tmp$estimate
  pval[i] = tmp$p.value
  ttpval[i] = t.test(x=dat.gex.normal[i,],y=dat.gex.cancer[i,])$p.value
}
iex = which(is.na(pval))
# iex = c(46,246,619,620) #left
# iex = 417 # right

tmp = qvalue(pval[-iex],fdr.level=.01); sum(tmp$significant)
summary(tmp)
sum(tmp$qvalues < 0.01 & corr[-iex] < 0)

tmp1 = gex.cancer[-iex]; tmp2 = gex.normal[-iex]
logfc = log2(tmp1/tmp2)
genes = glab[-iex]
cgi.genes = rlab[-iex]
cond0 = (tmp$qvalues < 0.01)

#### drift left colon (TCGA)
gexmex.left = data.frame(gene = cgi.genes,r=corr[-iex],
                          pval=pval[-iex],qval=tmp$qvalues,logFC = logfc)

## for Table (left colon)
dum = gexmex.left[cond0,] # top ranked by q-values < 0.01
io = order(dum$qval)
gexmex.left.drift = dum[io,]
write.csv(gexmex.left.drift,file='./GitHub/gexmex.left.drift.csv')

# > cor.test(gexmex.left.drift$logFC,gexmex.left.drift$r)
# 
# Pearson's product-moment correlation
# 
# data:  gexmex.left.drift$logFC and gexmex.left.drift$r
# t = 7.4068, df = 214, p-value = 2.942e-12
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
# 0.3386449 0.5519283
# sample estimates:
# cor 
# 0.4517172 

#### drift right colon (TCGA)
gexmex.right = data.frame(gene = cgi.genes,r=corr[-iex],
                         pval=pval[-iex],qval=tmp$qvalues,logFC = logfc)

## for Table (right colon)
dum = gexmex.right[cond0,] # top ranked by q-values < 0.01
io = order(dum$qval)
gexmex.right.drift = dum[io,]
write.csv(gexmex.right.drift,file='./GitHub/gexmex.right.drift.csv')

# > cor.test(gexmex.right.drift$logFC,gexmex.right.drift$r)
# 
# Pearson's product-moment correlation
# 
# data:  gexmex.right.drift$logFC and gexmex.right.drift$r
# t = 5.8779, df = 405, p-value = 8.691e-09
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
# 0.1882883 0.3675484
# sample estimates:
# cor 
# 0.2803609 

# > sum(abs(gexmex.right.drift$logFC)>1)/nrow(gexmex.right.drift)
# [1] 0.6289926
# > sum(abs(gexmex.left.drift$logFC)>1)/nrow(gexmex.left.drift)
# [1] 0.6203704


## static right colon (TCGA)
gexmex.Right = data.frame(gene = cgi.genes,r=corr[-iex],
                          pval=pval[-iex],qval=tmp$qvalues,logFC = logfc)

## for Table (right colon)
dum = gexmex.Right[cond0,] # top ranked by q-values < 0.01
io = order(dum$qval)
gexmex.right.static = dum[io,]
write.csv(gexmex.right.static,file='./GitHub/gexmex.right.static.csv')

## static left colon (TCGA)
gexmex.Left = data.frame(gene = cgi.genes,r=corr[-iex],
                          pval=pval[-iex],qval=tmp$qvalues,logFC = logfc)

## for Table (right colon)
dum = gexmex.Left[cond0,] # top ranked by q-values < 0.01
io = order(dum$qval)
gexmex.left.static = dum[io,]
write.csv(gexmex.left.static,file='./GitHub/gexmex.left.static.csv')

# > sum(abs(gexmex.left.static$logFC)>1)/nrow(gexmex.left.static)
# [1] 0.2315271
# > sum(abs(gexmex.right.static$logFC)>1)/nrow(gexmex.right.static)
# [1] 0.253112


length(intersect(gexmex.left.drift$gene,gexmex.right.drift$gene)) 

#### gene names
tmp = sapply(strsplit(as.character(gexmex.left$gene)," "), `[`, 2)
tmp = unique(tmp)
write(t(tmp),file='geneList.txt')

#### plot some data
gen = "ZNF304"
isl = "chr19:57862441-57863236"

cond = grepl(isl,rlab)
# cond = grepl(gen,glab)

plot(ilogit2(dat.mex[cond,]),(dat.gex[cond,]+1),pch=19,cex=0.6,
     xlab='CGI mean beta-value',ylab='(RSEM+1)')
legend('topleft',rlab[cond])

## add normal TCGA samples
ID0 = which(out.tcga.norm$islands==isl)
ID1 = which(out.tcga.left$islands==isl)

# ID = which(out.tcga.norm$genes=="PCDH8")
x0 = ilogit2(out.tcga.norm$Mvals[ID0,])
x1 = ilogit2(out.tcga.left$Mvals[ID1,])
y = expr.normal[grepl(gen,rownames(expr.normal)),ids.mex.normal]
z = ilogit2(out.sms.norm$Mvals[grepl(isl,out.sms.norm$islands),])

# points(x0,y,pch='+',col=2,cex=.5)

plot(density(x0),lwd=2,col=2,ylim=c(0,10),xlim=c(0,1),lty=3,main=rlab[cond],xlab='beta-value')
lines(density(x1),lwd=2,col=2)
lines(density(z),lwd=2,col=1)
legend('topright',c("normal healthy", "normal (adj cancer)","cancer"),
       lwd=c(2,2,2),col=c(1,2,2),lty=c(1,3,1))

#### change in gene expression for genes that are expressed RSEM2 > 1 in normal tissue
aux0 = log10(as.vector(dat.gex.normal[cond0,])+.1)
aux1 = log10(as.vector(dat.gex.cancer[cond0,])+.1)
plot(density(aux0,na.rm=T),col=1,lwd=2,ylim=c(0,1.))
lines(density(aux1,na.rm=T),col=2,lwd=2)

aux0 =  gex.normal[-iex] #[cond0]
cond =  aux0 > .5

tmp1 = gexmex.left$r[cond]
tmp2 = gexmex.left$logFC[cond]

cor.test(tmp1,tmp2)
plot(tmp1,tmp2,pch=19,cex=.6)
abline(h=0,v=0)

#### gex.normal vs gex.cancer
library(ggplot2)

labelSize = theme(axis.title.x = element_text(color="black", size=14, face="plain"),
                  axis.title.y = element_text(color="black", size=14, face="plain"),
                  title =        element_text(size=14, face='plain'))

## for left side
mex.left = apply(dat.mex,1,median)[-iex][cond0]
a = (gex.normal)[-iex][cond0]; b = (gex.cancer)[-iex][cond0]; r.left = gexmex.left$r[cond0]; cond2 = log10(a) > 0;
df.left = data.frame(x=log10(a),y=log2(b/a),site=rep("left",216),r=r.left,beta=mex.left)
df.left = df.left[cond2,]

# df.left = data.frame(x=log10(a),y=log10(b),site=rep("left",216),r=r.left,beta=mex.left)

## for right side
mex.right = apply(dat.mex,1,median)[-iex][cond0]
a = (gex.normal)[-iex][cond0]; b = (gex.cancer)[-iex][cond0]; r.right = gexmex.right$r[cond0]; cond2 = log10(a) > 0;
df.right = data.frame(x=log10(a),y=log2(b/a),site=rep("right",407),r=r.right,beta=mex.right)
df.right = df.right[cond2,]
# df.right = data.frame(x=log10(a),y=log10(b),site=rep("right",407),r=r.right,beta=mex.right)

df = rbind(df.left,df.right)
df$beta = ilogit2(df$beta)

# plot(log10(a),log2(b/a),pch=19,cex=0.6,xlab="log10 normal gene expression", ylab="log2 FC")
# abline(h=0)


ggplot(df, aes(x=x, y=y, color=beta)) + 
  geom_point(size=2.0, alpha=.8) +
  scale_color_gradient(low="darkblue", high="red") +
  geom_hline(yintercept=0, linetype=1, color='black', size=.5) +
  # geom_abline(slope=1,intercept=0, linetype=1, size=1.) +
  # geom_smooth(method="lm") +
  geom_density_2d(colour="black", bins=8) +
  geom_text(x=3,y=6,label="r = -0.31 (p = 4.e-12)",color='black',size = 6) +
  xlab("normal gene expression (log10)") +
  ylab("log2 FC") +
  xlim(-.5,4) +
  ylim(-6,6) +
  theme(legend.position="bottom") +
  # scale_color_manual(values=c("red", "green")) +
  labelSize

## other results (right)

## Null (random pairs n=10000)
# > sum(tmp.null.right$qvalues < 0.01 & corr.null.right[-iex.null] > 0)
# [1] 203
# > 203/9988
# [1] 0.02032439
# > sum(tmp.null.right$qvalues < 0.01 & corr.null.right[-iex.null] < 0)
# [1] 54
# > 54/9988
# [1] 0.005406488

## drift (right)
# > sum(tmp$qvalues < 0.01 & corr[-iex] < 0)
# [1] 223
# > sum(tmp$qvalues < 0.01 & corr[-iex] > 0)
# [1] 14
# > 223/643
# [1] 0.3468118  !!!!
# > 14/643
# [1] 0.02177294

## static
# > sum(tmp$qvalues < 0.01 & corr[-iex] < 0)
# [1] 191
# > 191/1070
# [1] 0.1785047
# > sum(tmp$qvalues < 0.01 & corr[-iex] > 0)
# [1] 12
# > 12/1070
# [1] 0.01121495

## results (left)
## Null (n=10000 random pairs)
# > sum(tmp$qvalues < 0.01 & corr[-iex] > 0)
# [1] 434
# > 434/9944
# [1] 0.04364441
# > sum(tmp$qvalues < 0.01 & corr[-iex] < 0)
# [1] 26
# > 26/9944
# [1] 0.002614642

## drift
# > sum(tmp$qvalues < 0.01 & corr[-iex] < 0)
# [1] 49
# > 49/640
# [1] 0.0765625
# > sum(tmp$qvalues < 0.01 & corr[-iex] > 0)
# [1] 6
# > 6/640
# [1] 0.009375

## static
# > sum(tmp$qvalues < 0.01 & corr[-iex] < 0)
# [1] 149
# > 149/1072
# [1] 0.1389925
# Error in args(obj) : attempt to use zero-length variable name
# > sum(tmp$qvalues < 0.01 & corr[-iex] > 0)
# [1] 23
# > 23/1072
# [1] 0.02145522

## test for island-by-island level drift with age
dat = cbind(out.tcga.left$Mvals[iset,],out.tcga.right$Mvals[iset,])
dat = static.tcga.left$Mvals
age = c(age.tcga.left,age.tcga.right)
M = nrow(dat)
pval=corr=numeric(M)
for (m in 1:M) {
  tmp = cor.test(dat[m,],age)
  corr[m] = tmp$estimate
  pval[m] = tmp$p.value
}

tmp = qvalue(pval,fdr.level=.01); sum(tmp$significant)
summary(tmp)

## compute z-scores relative to gene expression in methylomically young patients
dum = apply(dat,2,function(x) {mean(x[x>-20])})
plot(sort(dum),type='s')
io = order(dum)
io50 = io[1:50]
# mean gene expression and sd for io50
g.mu = apply(dat.gex[,io50],1,mean,na.rm=T)
g.sd = apply(dat.gex[,io50],1,sd,na.rm=T)
#
dat.z = dat.gex
for (i in 1:G) {
  dat.z[i,] = (dat.gex[i,]-g.mu[i])/g.sd[i] 
}

io = order(tmp$qvalues)
corr[-iex][io[1:10]]
pval[-iex][io[1:10]]
rlab[-iex][io[1:10]]

# results for right colon: 643 gene-island pairs
# > sum(corr[-iex][tmp$qvalues < 0.05] > 0)
# [1] 32
# > sum(corr[-iex][tmp$qvalues < 0.05] < 0)
# [1] 341
# >

# results for right colon: 640 gene-island pairs
# > sum(corr[-iex][tmp$qvalues < 0.05] > 0)
# [1] 18
# > sum(corr[-iex][tmp$qvalues < 0.05] < 0)
# [1] 112
# 


