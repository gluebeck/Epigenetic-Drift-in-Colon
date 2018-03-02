#### island-level drift in normal tissues 
#### determine intercepts
dat = cbind(out.gic.left$Mvals[iset,],
            out.sms.norm$Mvals[iset,],
            out.gic.right$Mvals[iset,])

age = c(age.gic.left,age.sms,age.gic.right)
Mp = nrow(dat)
slope.normal = pval.normal = intercept.normal = numeric(Mp)
for(m in 1:Mp) {
  aux = summary(lm(dat[m,] ~ age))
  pval.normal[m] = aux$coef[2,4]
  slope.normal[m] = aux$coef[2,1]
  intercept.normal[m] = aux$coef[1,1]
}
## assuming intercepts are similar between left and right colon for the drift islands
eps.left = eps.right = intercept.normal

#### requires functions to compute sojourn time distributions (adenoma-to-carcinoma) 
source('./MSCE_conditional_sojourn_time.R')

####  prepare left-right sojourn time predictions
source('./Simulations/parameters_males_right_colon.R')
a = seq(20,90,1)
adT = adT2 = numeric(length(a)); i=1

for (s in a) {
  u1 = seq(0,s-tlag,1); dens = f.adenoma(u1)
  u1.mean = sum(u1*dens)/sum(dens)
  u1.var =  sum(((u1-u1.mean)^2)*dens)/sum(dens)
  adT[i] = s - u1.mean 
  adT2[i] = 1.96*sqrt(u1.var)
  i=i+1
}
adT.right=adT; adT2.right=adT2

#### non-linear regression: using adT on age - left colon
source('./Simulations/parameters_males_left_colon.R')
a = seq(20,90,1)
adT = adT2 = numeric(length(a)); i=1

for (s in a) {
  u1 = seq(0,s-tlag,1); dens = f.adenoma(u1)
  u1.mean = sum(u1*dens)/sum(dens)
  u1.var =  sum(((u1-u1.mean)^2)*dens)/sum(dens)
  adT[i] = s - u1.mean 
  adT2[i] = 1.96*sqrt(u1.var)
  i=i+1
}
adT.left=adT; adT2.left=adT2

## non-linear regression: using adT on age - rectum
source('./Simulations/parameters_males_rectum.R')
a = seq(20,90,1)
adT = adT2 = numeric(length(a)); i=1

for (s in a) {
  u1 = seq(0,s-tlag,1); dens = f.adenoma(u1)
  u1.mean = sum(u1*dens)/sum(dens)
  u1.var =  sum(((u1-u1.mean)^2)*dens)/sum(dens)
  adT[i] = s - u1.mean 
  adT2[i] = 1.96*sqrt(u1.var)
  i=i+1
}
adT.rect=adT; adT2.rect=adT2

#### Figure_4.R here


