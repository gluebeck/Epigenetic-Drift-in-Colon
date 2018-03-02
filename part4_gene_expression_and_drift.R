#### gene expression RSEM v2 data downloaded from CBio Portal 
#### analyses for static or drift related genes
#### drift: Genes.isl.drift5 (part2)
#### static: Genes.isl.static5

dum = as.data.frame(GEXstaticRaw) 
dum = as.data.frame(GEXdriftRaw)

genes = dum$COMMON
rownames(dum) = genes
dum = dum[,-(1:2)]
dum = dum[,-397] # remove "X399"
dum[dum=="NaN"]=NA
iex = apply(dum,1,function(x){all(is.na(x))}) # remove genes that lack data
dum = dum[!iex,]

GEXdata = dum
dum = colnames(GEXdata);
dum = substr(dum,1,12)
colnames(GEXdata)=dum
GEXdata=as.data.frame(GEXdata)

## choose
GEXdrift  = GEXdata # or
GEXstatic = GEXdata

## choose
dat = cbind(out.tcga.left$Mvals[iset,],out.tcga.right$Mvals[iset,])
dat = cbind(out.tcga.right$Mvals[iset,])
dat = cbind(out.tcga.left$Mvals[iset,])
dat = cbind(static.tcga.left$Mvals)
dat = cbind(static.tcga.right$Mvals)

ids = colnames(dat)
bcr = pdata.tcga$bcr_patient_barcode[match(ids,pdata.tcga$arrayID)]

ids.gex = (match(bcr,colnames(GEXdata)))

M = nrow(dat) #number of CGI to test for gex data

## choose
geneList = out.tcga.left$genes  
islandList = out.tcga.left$islands

geneList = out.tcga.right$genes  
islandList = out.tcga.right$islands

geneList = static.tcga.left$genes  
islandList = static.tcga.left$islands

geneList = static.tcga.right$genes  
islandList = static.tcga.right$islands


rlab = dat.mex = dat.gex = NULL
for ( m in 1:M) {
  mex = dat[m,]
  mil = islandList[m]
  genes  = geneList[[m]]; len = length(genes)
  if(len>0) {
    for(l in 1:len) {
      gene = genes[l]
      if(gene%in%rownames(GEXdata)) {
        # find gene expression for gene
        gex  = GEXdata[gene,ids.gex]
        # pair them up
        dat.mex = rbind(dat.mex,mex)
        dat.gex = rbind(dat.gex,gex)
        rlab = c(rlab,paste(mil,gene))
      }
    }
  }
}
rownames(dat.mex) = rownames(dat.gex) = rlab

dat.gex = as.matrix(dat.gex)
G = nrow(dat.gex)

# isl-gene correlations
corr = pval = numeric(G)
for (i in 1:G) {
  tmp = cor.test(dat.mex[i,],log(dat.gex[i,]+1))
  # tmp = cor.test(dat.mex[sample(1:G,1),],log(dat.gex[i,]+1))
  # tmp = cor.test(static.tcga.right$Mvals[sample(1:1016,1),],log(dat.gex[sample(1:G,1),]+1))
  # tmp = cor.test(static.tcga.left$Mvals[sample(1:1016,1),],log(dat.gex[sample(1:G,1),]+1))
  # tmp = cor.test(out.tcga.right$Mvals[sample(1:871,1),],log(dat.gex[i,]+1))
  
  corr[i] = tmp$estimate
  pval[i] = tmp$p.value
}
iex = which(is.na(pval))
# iex = c(46,246,619,620) #left
# iex = 417 # right

tmp = qvalue(pval[-iex],fdr.level=.01); sum(tmp$significant)
summary(tmp)

genes = rlab[-iex]
io = order(tmp$qvalues)
cond = (tmp$qvalues < 0.01)

## drift left colon (TCGA)
gexmex.left = data.frame(gene = genes[io],r=corr[-iex][io],
                          pval=pval[-iex][io],qval=tmp$qvalues[io])
gexmex.left = gexmex.left[cond[io],] # top ranked by q-values < 0.01
write.csv(gexmex.left,file='./GitHub/gexmex.left.csv')

## drift right colon (TCGA)
gexmex.right = data.frame(gene = genes[io],r=corr[-iex][io],
                          pval=pval[-iex][io],qval=tmp$qvalues[io])
gexmex.right = gexmex.right[cond[io],] # top ranked by q-values < 0.01
write.csv(gexmex.right,file='./GitHub/gexmex.right.csv')

intersect(gexmex.left$gene,gexmex.right$gene) 

## static right colon (TCGA)
gexmex.Right = data.frame(gene = genes[io],r=corr[-iex][io],
                          pval=pval[-iex][io],qval=tmp$qvalues[io])
gexmex.Right = gexmex.Right[cond[io],] # top ranked by q-values < 0.01
write.csv(gexmex.Right,file='./GitHub/gexmex.Right.csv')

## static left colon (TCGA)
gexmex.Left = data.frame(gene = genes[io],r=corr[-iex][io],
                          pval=pval[-iex][io],qval=tmp$qvalues[io])
gexmex.Left = gexmex.Left[cond[io],] # top ranked by q-values < 0.01
write.csv(gexmex.Left,file='./GitHub/gexmex.left.static.csv')

intersect(gexmex.Left$gene,gexmex.Right$gene) 

## plot some data
cond = grepl("chr3:37034228-37035356",rlab)
plot(ilogit2(dat.mex[cond,]),log10(dat.gex[cond,]+1),pch=19,cex=0.6,
     xlab='CGI mean beta-value',ylab='log10(RSEM+1)')
legend('topright',rlab[cond])

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


