library(gplots)
library(grid)
library(pheatmap)

## left colon (easily modified for right colon)
dat1 = cbind(out.gic.left$Mvals[iset,],out.tcga.left$Mvals[iset,ids.tcga.left])
## just take first 300 static islands in the list
dat2 = cbind(static.gic.left$Mvals[1:300,],static.tcga.left$Mvals[1:300,ids.tcga.left])
# dat2 = cbind(static.gic.left$Mvals,static.tcga.left$Mvals)
dat = rbind(dat2,dat1)
# dat = cbind(out.gic.left$Mvals[iset,],out.tcga.left$Mvals[iset,])
sex = factor(c(sex.gic.left,sex.tcga.left[ids.tcga.left]))
age = c(age.gic.left,age.tcga.left[ids.tcga.left])
age = factor(10*floor(age/10))
dum = 1:length(age.gic.left)

tmp = colnames(dat) # patients
# truncate long IDs
tmp = substr(tmp,1,7)

# aux1.norm.static = apply(dat2[,dum],2,function(x) {mean(x[x > -20])})
aux1.norm.drift  = apply(dat1[,dum],2,function(x) {mean(x[x > -20])})
aux1.canc.drift  = apply(dat1[,-dum],2,function(x) {mean(x[x > -20])})
aux2 = apply(dat[,dum],1,function(x) {mean(x[x > -20])})

aux2.norm.static = apply(dat2[,dum],1,function(x) {mean(x[x > -20])})
aux2.norm.drift  = apply(dat1[,dum],1,function(x) {mean(x[x > -20])})
aux2.canc.static = apply(dat2[,-dum],1,function(x) {mean(x[x > -20])})
aux2.canc.drift  = apply(dat1[,-dum],1,function(x) {mean(x[x > -20])})

# > cor.test(aux2.norm,aux2.canc)
# 
# Pearson's product-moment correlation
# 
# data:  aux2.norm and aux2.canc
# t = 37.996, df = 779, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
# 0.7798744 0.8292049
# sample estimates:
# cor 
# 0.8059342 

col.order.norm = order(aux1.norm.drift)
col.order.canc = order(aux1.canc.drift)
col.order = c(col.order.norm,max(dum)+col.order.canc)

row.order.static = order(aux2.norm.static)
row.order.drift  = order(aux2.norm.drift)
row.order = c(row.order.static,length(row.order.static)+row.order.drift)

mat = ilogit2(dat)[row.order,col.order] #c(1:10,10+order(pat))]

# labels = c(colnames(sqM)[1:10],tmp[c(id.MM.BE,id.BM.BE,id.HM.BE)])
# labels.o = labels[type.order]

annotation_col = data.frame(Age=age[col.order],Sex=sex[col.order])
rownames(annotation_col) = labels.o = as.character(col.order) #paste0("Subject", col.order)
tmp = RColorBrewer::brewer.pal(8, "Greys")

# Age=c("30"=tmp[1],"40"=tmp[2],"50"=tmp[3],
#      "60"=tmp[4],"70"=tmp[5],"80"=tmp[6],"90"=tmp[7])

Age=c("20"=tmp[1],"30"=tmp[2],"40"=tmp[3],
      "50"=tmp[4],"60"=tmp[5],"70"=tmp[6],"80"=tmp[7],"90"=tmp[8])

ann_colors =  list(Sex=c(F="pink",M="cyan"), Age=Age)
col.pal <- rev(RColorBrewer::brewer.pal(11, "RdYlBu"))
colnames(mat) = labels.o

# hmap.left.colon = 
pheatmap((mat), cluster_cols=F, cluster_rows=F, 
         # cutree_rows = 4,
         # kmeans_k = 4,
         # clustering_distance_rows = "correlation",#Pearson's
         # clustering_method = "ward.D2",
         main = " ",
         color=col.pal,
         show_colnames = T, show_rownames = F,
         annotation_names_col = F,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         labels_col = labels.o,
         gaps_col = max(dum),
         fontsize_col = 6,
         cellwidth = 1.8,cellheight=.27)

grid.text(x=unit(0.18,"npc"),y=unit(0.97,"npc"),"Normal Colon")
grid.text(x=unit(0.5,"npc"),y=unit(0.97,"npc"),"Carcinoma (left colon)")
grid.text(x=unit(0.02,"npc"),y=unit(0.75,"npc"),"static CGI",rot=90)
grid.text(x=unit(0.02,"npc"),y=unit(0.30,"npc"),"drift CGI",rot=90)

grid.text(x=unit(0.1,"npc"),y=unit(0.97,"npc"),"Normal Colon")
grid.text(x=unit(0.5,"npc"),y=unit(0.97,"npc"),"Carcinoma (right colon)")
grid.text(x=unit(0.02,"npc"),y=unit(0.75,"npc"),"static CGI",rot=90)
grid.text(x=unit(0.02,"npc"),y=unit(0.30,"npc"),"drift CGI",rot=90)

grid.text("beta-value",x=.845,y=.36,rot = 90)

