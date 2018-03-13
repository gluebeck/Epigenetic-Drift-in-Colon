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
  geom_line(aes(y=y,color=site),size=2) +
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
age = age.left = c(age.luo.left,age.tcga.left[ids.tcga.left])
adT = adT2 = numeric(length(age)); i=1

for (s in age) {
  u1 = seq(0,s-tlag,1); dens = f.adenoma(u1)
  u1.mean = sum(u1*dens)/sum(dens)
  u1.var =  sum(((u1-u1.mean)^2)*dens)/sum(dens)
  adT[i] = s - u1.mean 
  adT2[i] = 1.96*sqrt(u1.var)
  i=i+1
}
AdT.left = adT; AdT2.left = adT2

aux = apply(dat,2,function(x) {mean(x[x > -20.])})
dum = lm(aux ~ 0 + AdT.left, offset = rep(eps.left,length(AdT.left)))
tmp = summary(dum)
b.left = tmp$coefficients[1,1]
b.left.sd = tmp$coefficients[1,2] 

y.left = mean(eps.left)+b.left*adT.left
y.left.up = mean(eps.left)+b.left*(adT.left+adT2.left)
y.left.dn = mean(eps.left)+b.left*(adT.left-adT2.left)

## non-linear regression: using adT on age - rectum only
source('./Simulations/parameters_males_rectum.R') 
dat = cbind(out.luo.rect$Mvals[iset,],out.tcga.left$Mvals[iset,ids.tcga.rect])
age = age.rect = c(age.luo.rect,age.tcga.left[ids.tcga.rect])
adT = adT2 = numeric(length(age)); i=1

for (s in age) {
  u1 = seq(0,s-tlag,1); dens = f.adenoma(u1)
  u1.mean = sum(u1*dens)/sum(dens)
  u1.var =  sum(((u1-u1.mean)^2)*dens)/sum(dens)
  adT[i] = s - u1.mean 
  adT2[i] = 1.96*sqrt(u1.var)
  i=i+1
}
AdT.rect = adT; AdT2.rect = adT2

aux = apply(dat,2,function(x) {mean(x[x > -20.])})
dum = lm(aux ~ 0 + AdT.rect, offset = rep(eps.left,length(AdT.rect)))
tmp = summary(dum)
b.rect = tmp$coefficients[1,1]
b.rect.sd = tmp$coefficients[1,2] 

y.rect = mean(eps.left)+b.rect*adT.rect
y.rect.up = mean(eps.left)+b.rect*(adT.rect+adT2.rect)
y.rect.dn = mean(eps.left)+b.rect*(adT.rect-adT2.rect)

## non-linear regression: using adT on age - right colon
source('./Simulations/parameters_males_right_colon.R')
dat = cbind(out.luo.right$Mvals[iset,],out.tcga.right$Mvals[iset,])
age = age.right = c(age.luo.right,age.tcga.right)
adT = adT2 = numeric(length(age)); i=1

for (s in age) {
  u1 = seq(0,s-tlag,1); dens = f.adenoma(u1)
  u1.mean = sum(u1*dens)/sum(dens)
  u1.var =  sum(((u1-u1.mean)^2)*dens)/sum(dens)
  adT[i] = s - u1.mean 
  adT2[i] = 1.96*sqrt(u1.var)
  i=i+1
}
AdT.right = adT; AdT2.right = adT2

aux = apply(dat,2,function(x) {mean(x[x > -20.])})
dum = lm(aux ~ 0 + AdT.right, offset = rep(eps.right,length(AdT.right)))
tmp = summary(dum)
b.right = tmp$coefficients[1,1]
b.right.sd = tmp$coefficients[1,2] 

y.right = mean(eps.right)+b.right*adT.right
y.right.up = mean(eps.right)+b.right*(adT.right+adT2.right)
y.right.dn = mean(eps.right)+b.right*(adT.right-adT2.right)

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
pred.rect = data.frame(age=a, 
                       y.rect =y.rect, y.rect.up = y.rect.up, y.rect.dn = y.rect.dn)

pred.left = data.frame(age=a, 
                       y.left =y.left, y.left.up = y.left.up, y.left.dn = y.left.dn)

pred.right = data.frame(age=a, 
                        y.right = y.right,y.right.up = y.right.up,y.right.dn = y.right.dn)

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

aux1 = apply(dat,2,function(x) {mean(x[x > -20.])})
aux2 = apply(dat,2,function(x) {sd(x[x > -20.])})
mat = data.frame(age=age,Mval=aux1,InvSD=aux2,site=site)

plot2 =ggplot() +
  geom_point(data=mat,aes(x=age,y=Mval,color=site),show.legend = T) +
  # labs(title = "measured drift in CRC with scaled prediction") +
  labs(title = "B") +
  ylab("mean M-value") +
  xlim(20,90) +
  ylim(-3.,1.5) +
  geom_line(data=pred.rect,aes(x=age, y=y.rect, color="Rectum"),size=2,show.legend = F) +
  geom_line(data=pred.left,aes(x=age, y=y.left, color="Colon (dist.)"),size=2,show.legend = F) +
  geom_line(data=pred.right,aes(x=age, y=y.right, color="Colon (prox.)"),size=2,show.legend = F) +
  geom_ribbon(aes(x=age, ymax=y.rect.up, ymin=y.rect.dn), data = pred.rect, 
              alpha=.15, size=0) +
  geom_ribbon(aes(x=age, ymax=y.left.up, ymin=y.left.dn), data = pred.left, 
              alpha=.15, size=0) +
  geom_ribbon(aes(x=age, ymax=y.right.up, ymin=y.right.dn), data = pred.right, 
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

