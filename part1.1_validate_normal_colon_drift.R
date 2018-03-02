#### validation of island-level drift: drift5 islands
#### males SMS, vs males GICares (left) vs all GICares (right) 

## left colon (GICares)
age = c(age.gic.left)
sex = c(sex.gic.left)
dat = cbind(out.gic.left$Mvals[iset,]); Mp = nrow(dat)

inter.m = inter.f = pval.inter.m = numeric(Mp)
pval.slope.f = pval.slope.m = slope.f = slope.m = numeric(Mp)

for (m in 1:Mp) {
  aux = summary(lm(dat[m,] ~ age*sex))
  pval.slope.f[m] = aux$coef[2,4]; pval.slope.m[m] = aux$coef[4,4]
  slope.f[m] = aux$coef[2,1]; slope.m[m] = slope.f[m]+aux$coef[4,1]
  inter.f[m] = aux$coef[1,1]
  inter.m[m] = inter.f[m]+aux$coef[3,1]; pval.inter.m[m] = aux$coef[3,4]
}

# island-level slopes
GIC.Slope.f = slope.f
GIC.Slope.m = slope.m

## right colon (both sexes comb)
age = c(age.gic.right) 
dat = cbind(out.gic.right$Mvals[iset,]); Mp = nrow(dat)

slope.normal = pval.normal = intercept.normal = numeric(Mp)
for(i in 1:Mp) {
  aux = summary(lm(dat[i,] ~ age))
  pval.normal[i] = aux$coef[2,4]
  slope.normal[i] = aux$coef[2,1]
  intercept.normal[i] = aux$coef[1,1]
}

Slope.right = slope.normal

#### show data in boxplot: see Figure_3.R


# ## right colon. GIC right samples + Kaz Ascending Colon (n=11) 
# cond.kaz.right = loc.kaz=="A"
# 
# age = c(age.gic.right,age.kaz[cond.kaz.right])
# sex = c(sex.gic.right,sex.kaz[cond.kaz.right])
# 
# dat = cbind(out.gic.right$Mvals[iset,],out.kaz$Mvals[iset,cond.kaz.right]); Mp = nrow(dat)
# 
# inter.m = inter.f = pval.inter.m = numeric(Mp)
# pval.slope.f = pval.slope.m = slope.f = slope.m = numeric(Mp)
# 
# for (m in 1:Mp) {
#   aux = summary(lm(dat[m,] ~ age*sex))
#   pval.slope.f[m] = aux$coef[2,4]; pval.slope.m[m] = aux$coef[4,4]
#   slope.f[m] = aux$coef[2,1]; slope.m[m] = slope.f[m]+aux$coef[4,1]
#   inter.f[m] = aux$coef[1,1]
#   inter.m[m] = inter.f[m]+aux$coef[3,1]; pval.inter.m[m] = aux$coef[3,4]
# }
# 
# # GiCares + Kaz
# Kaz.Slope.f = slope.f
# Kaz.Slope.m = slope.m

