#### island-level drift in normal tissues 
#### determine intercepts
dat = cbind(out.gic.left$Mvals[iset,],
            out.sms.norm$Mvals[iset,],
            out.gic.right$Mvals[iset,])
age = c(age.gic.left,age.sms,age.gic.right)
sex = c(sex.gic.left,sex.rec,sex.gic.right)

## if right side only
## dat = out.gic.right$Mvals[iset,]; age = age.gic.right

## joint
a = age
aux1 = apply(dat,2,function(x) {mean(x[x > -20.])})
aux = summary(lm(aux1 ~ a))
intercept.normal = aux$coef[1,1]
slope.normal = aux$coef[2,1]

## females
a = age[sex=="F"]
aux1 = apply(dat[,sex=="F"],2,function(x) {mean(x[x > -20.])})
aux = summary(lm(aux1 ~ a))
intercept.normal.F = aux$coef[1,1]
slope.normal.F = aux$coef[2,1]

## males
a = age[sex=="M"]
aux1 = apply(dat[,sex=="M"],2,function(x) {mean(x[x > -20.])})
aux = summary(lm(aux1 ~ a))
intercept.normal.M = aux$coef[1,1]
slope.normal.M = aux$coef[2,1]

# Mp = nrow(dat)
# slope.normal = pval.normal = intercept.normal = numeric(Mp)
# for(m in 1:Mp) {
#   aux = summary(lm(dat[m,] ~ age))
#   pval.normal[m] = aux$coef[2,4]
#   slope.normal[m] = aux$coef[2,1]
#   intercept.normal[m] = aux$coef[1,1]
# }

## assuming intercepts are similar between left and right colon for the drift islands
eps.left = eps.right = intercept.normal

#### requires functions to compute sojourn time distributions (adenoma-to-carcinoma) 
# source('./GitHub/MSCE_conditional_sojourn_time.R')

####  prepare left-right sojourn time predictions
source('./GitHub/parameters_males_right_colon.R')
source('./GitHub/parameters_males_left_colon.R')
source('./GitHub/parameters_males_rectum.R')

source('./GitHub/parameters_females_right_colon.R')
source('./GitHub/parameters_females_left_colon.R')
source('./GitHub/parameters_females_rectum.R')

a = seq(20,90,1)
adT = adTe = numeric(length(a)); i=1

for (s in a) {
  u1 = seq(0,s-tlag,1); dens = f.adenoma(u1)
  u1.mean = sum(u1*dens)/sum(dens)
  u1.var =  sum(((u1-u1.mean)^2)*dens)/sum(dens)
  adT[i] = s - u1.mean 
  adTe[i] = sqrt(u1.var)
  i=i+1
}

adT.right.m=adT; adTe.right.m=adTe
adT.left.m=adT; adTe.left.m=adTe
adT.rect.m=adT; adTe.rect.m=adTe

adT.right.f=adT; adTe.right.f=adTe
adT.left.f=adT; adTe.left.f=adTe
adT.rect.f=adT; adTe.rect.f=adTe

#### Figure_4.R here


