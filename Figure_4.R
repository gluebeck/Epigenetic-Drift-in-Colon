library(ggplot2)
library(plyr)

#### ggplot A: sojourn time predictions
#### use code part3_model_predictions_....

a = seq(20,90,1)

dat1 = data.frame(age=a, y=adT.left,  se=adT2.left, site=rep("Colon (dist.)",71))
dat2 = data.frame(age=a, y=adT.right, se=adT2.right, site=rep("Colon (prox.)",71))
dat3 = data.frame(age=a, y=adT.rect, se=adT2.rect, site=rep("Rectum",71))
df1 = rbind(dat1,dat2,dat3)

labelSize = theme(axis.title.x = element_text(color="black", size=14, face="plain"),
                  axis.title.y = element_text(color="black", size=14, face="plain"),
                  title =        element_text(size=14, face='plain'))

plot1 = ggplot(df1, aes(x=age)) +
  geom_line(aes(y=y,color=site),size=1.5) +
  # labs(title = "model predictions for sojourn time") +
  labs(title = "A") +
  geom_ribbon(aes(ymax=df1$y+df1$se, ymin=df1$y-df1$se,color=site), alpha=.15, size=0) +
  ylab("predicted sojourn time (years)") +
  xlab("age") +
  xlim(20,90) +
  ylim(0,85) +
  theme(legend.position="bottom") +
  labelSize

#### model fits to data assuming constant neoplastic mDNA drift rate

## non-linear regression: using adT on age - left colon excluding rectum
source('./Simulations/parameters_males_left_colon.R') 
dat = cbind(out.luo.left$Mvals[iset,],out.tcga.left$Mvals[iset,ids.tcga.left])
age = c(age.luo.left,age.tcga.left[ids.tcga.left])
adT = adT2 = numeric(length(age)); i=1

for (s in age) {
  u1 = seq(0,s-tlag,1); dens = f.adenoma(u1)
  u1.mean = sum(u1*dens)/sum(dens)
  u1.var =  sum(((u1-u1.mean)^2)*dens)/sum(dens)
  adT[i] = s - u1.mean 
  adT2[i] = 1.96*sqrt(u1.var)
  i=i+1
}

Mp = nrow(dat)
b.left = numeric(Mp)
for (m in 1:Mp) {
  # cat(m,'\n')
  cond = dat[m,] > -2
  dum = lm(dat[m,cond] ~ 0 + adT[cond], offset = rep(eps.left[m],length(adT[cond])))
  tmp = summary(dum)
  b.left[m] = tmp$coefficients[1,1]  # eps = tmp$coefficients[1,1]
}
tmp = quantile(b.left,probs=c(.025,.5,.975))

y.left = mean(eps.left)+tmp[2]*adT.left
y.left.up = mean(eps.left)+tmp[2]*(adT.left+adT2.left)
y.left.dn = mean(eps.left)+tmp[2]*(adT.left-adT2.left)

## non-linear regression: using adT on age - rectum only
source('./Simulations/parameters_males_rectum.R') 
dat = cbind(out.luo.rect$Mvals[iset,],out.tcga.left$Mvals[iset,ids.tcga.rect])
age = c(age.luo.rect,age.tcga.left[ids.tcga.rect])
adT = adT2 = numeric(length(age)); i=1

for (s in age) {
  u1 = seq(0,s-tlag,1); dens = f.adenoma(u1)
  u1.mean = sum(u1*dens)/sum(dens)
  u1.var =  sum(((u1-u1.mean)^2)*dens)/sum(dens)
  adT[i] = s - u1.mean 
  adT2[i] = 1.96*sqrt(u1.var)
  i=i+1
}

Mp = nrow(dat)
b.rect = numeric(Mp)
for (m in 1:Mp) {
  # cat(m,'\n')
  cond = dat[m,] > -2.1 # small adjustment
  dum = lm(dat[m,cond] ~ 0 + adT[cond], offset = rep(eps.left[m],length(adT[cond])))
  tmp = summary(dum)
  b.rect[m] = tmp$coefficients[1,1]  # eps = tmp$coefficients[1,1]
}
tmp = quantile(b.rect,probs=c(.025,.5,.975))

y.rect = mean(eps.left)+tmp[2]*adT.rect
y.rect.up = mean(eps.left)+tmp[2]*(adT.rect+adT2.rect)
y.rect.dn = mean(eps.left)+tmp[2]*(adT.rect-adT2.rect)

## non-linear regression: using adT on age - right colon
source('./Simulations/parameters_males_right_colon.R')
dat = cbind(out.luo.right$Mvals[iset,],out.tcga.right$Mvals[iset,])
age = c(age.luo.right,age.tcga.right)
adT = adT2 = numeric(length(age)); i=1

for (s in age) {
  u1 = seq(0,s-tlag,1); dens = f.adenoma(u1)
  u1.mean = sum(u1*dens)/sum(dens)
  u1.var =  sum(((u1-u1.mean)^2)*dens)/sum(dens)
  adT[i] = s - u1.mean 
  adT2[i] = 1.96*sqrt(u1.var)
  i=i+1
}

Mp = nrow(dat)
b.right = numeric(Mp)
for (m in 1:Mp) {
  # cat(m,'\n')
  cond = dat[m,] > -2
  dum = lm(dat[m,cond] ~ 0 + adT[cond], offset = rep(eps.right[m],length(adT[cond])))
  tmp = summary(dum)
  b.right[m] = tmp$coefficients[1,1]  # eps = tmp$coefficients[1,1]
}
tmp = quantile(b.right,probs=c(.025,.5,.975))

y.right = mean(eps.right)+tmp[2]*adT.right
y.right.up = mean(eps.right)+tmp[2]*(adT.right+adT2.right)
y.right.dn = mean(eps.right)+tmp[2]*(adT.right-adT2.right)

# ggplot
r0 = 1:length(age.luo.rect)
r1 = 1:length(age.luo.left)
r2 = 1:length(age.luo.right)
pred.on = c(y.rect[r0],y.left[r1],y.right[r2],
            y.rect[-r0],y.left[-r1],y.right[-r2])
pred.up = c(y.rect.up[r0],y.left.up[r1],y.right.up[r2],
            y.rect.up[-r0],y.left.up[-r1],y.right.up[-r2])
pred.dn = c(y.rect.dn[r0],y.left.dn[r1],y.right.dn[r2],
            y.rect.dn[-r0],y.left.dn[-r1],y.right.dn[-r2])

a = seq(20,90,1)
pred = data.frame(age=a, y.right=y.right,y.right.up=y.right.up,y.right.dn=y.right.dn,
                  y.left =y.left, y.left.up = y.left.up, y.left.dn = y.left.dn,
                  y.rect =y.rect, y.rect.up = y.rect.up, y.rect.dn = y.rect.dn)

dat = cbind(out.luo.rect$Mvals[iset,],
            out.luo.left$Mvals[iset,],
            out.luo.right$Mvals[iset,],
            out.tcga.left$Mvals[iset,ids.tcga.rect],
            out.tcga.left$Mvals[iset,ids.tcga.left],
            out.tcga.right$Mvals[iset,])

site = c(rep("Rectum",9),rep("Colon (dist.)",27),rep("Colon (prox.)",44),rep("Rectum",43),
         rep("Colon (dist.)",141),rep("Colon (prox.)",138))

age = c(age.luo.rect,age.luo.left,age.luo.right,age.tcga.left[ids.tcga.rect],
        age.tcga.left[ids.tcga.left],age.tcga.right)

aux1 = apply(dat,2,function(x) {mean(x[x > -2.])})
aux2 = apply(dat,2,function(x) {sd(x[x > -2.])})
mat = data.frame(age=age,Mval=aux1,InvSD=aux2,site=site)
# mat.left = mat[mat$site=="Left",]
# mat.right = mat[mat$site=="Right",]

#### optional: plot mmr status as defined in ... Island_level_drift.R
# r0 = 77:398; r0 = r0[c(!is.na(mmr.tcga.left),!is.na(mmr.tcga.right))] # or
# # r0 = 261:398; r0 = r0[c(!is.na(mmr.tcga.right))]
# mat.mmr0 = data.frame(age=age[r0],Mval=aux1[r0],MMR_IHC=rep("tested",length(r0)))
# 
# r1 = 77:398; r1 = r1[c(cond.mmr.left,cond.mmr.right)] # or
# # r1 = 261:398; r1 = r1[c(cond.mmr.right)]
# mat.mmr1 = data.frame(age=age[r1],Mval=aux1[r1],MMR_IHC=rep("negative",length(r1)))
# 
# r2 = setdiff(r0,r1)
# mat.mmr2 = data.frame(age=age[r2],Mval=aux1[r2],MMR_IHC=rep("positive",length(r2)))
# 
# mat.mmr = rbind(mat.mmr0,mat.mmr1)
# mat.mmr$MMR_IHC = as.character(mat.mmr$MMR_IHC)

## highlight some outliers
# dum = ellipseFun(c(48,0.5),c(25,1.2),npoints = 100)

plot2 =ggplot() +
  geom_point(data=mat,aes(x=age,y=Mval,color=site),show.legend = T) +
  # labs(title = "measured drift in CRC with scaled prediction") +
  labs(title = "B") +
  ylab("mean M-value") +
  xlim(23,90) +
  ylim(-2.8,1.8) +
  geom_line(data=pred,aes(x=age, y=y.rect, color="Rectum"),size=2,show.legend = F) +
  geom_line(data=pred,aes(x=age, y=y.left, color="Colon (dist.)"),size=2,show.legend = F) +
  geom_line(data=pred,aes(x=age, y=y.right, color="Colon (prox.)"),size=2,show.legend = F) +
  geom_ribbon(aes(x=age, ymax=y.rect.up, ymin=y.rect.dn), data = pred, 
              alpha=.15, size=0) +
  geom_ribbon(aes(x=age, ymax=y.left.up, ymin=y.left.dn), data = pred, 
              alpha=.15, size=0) +
  geom_ribbon(aes(x=age, ymax=y.right.up, ymin=y.right.dn), data = pred, 
              alpha=.15, size=0) +
  theme(legend.position="bottom") +
  # geom_path(aes(x,y),data=dum,color='blue',size=1.3) +
  # geom_rect(aes(xmin=35,xmax=60,ymin=.0,ymax=1), alpha=0.2,fill="cyan",color='blue') + 
  labelSize

# p1 +
#  geom_point(data=mat.mmr0,aes(x=age,y=Mval),color="green",size=2) +
#  geom_point(data=mat.mmr1,aes(x=age,y=Mval),color="red",size=2.5)

multiplot(plot1,plot2,cols=2)

#### plot estimated drift rate distributions  
  
df  = data.frame(b=c(b.left,b.right),site=
                     c(rep("Left",length(b.left)),rep("Right",length(b.right))))
ggplot(df,aes(df$b, fill=site)) + 
  geom_histogram(alpha = 0.4, position="identity",breaks=seq(0, .095, by = .005),color="black") + 
  # geom_density(aes(fill=site),alpha=.1) +
  xlim(c(0.02,.1)) + 
  labs(title="estimated drift rates") +
  labs(x="drift rate", y="count") +
  labelSize

