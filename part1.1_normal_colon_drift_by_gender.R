#### set-up ANCOVA regression with gender as categorical variable
age.rec = c(age.rec.f,age.rec.m)
sex.rec = c(rep("F",90),rep("M",60))

## SMS data for drift5 CpG probes
smsM = getM(mset.swan[cpgs.isl.drift5,c(ID.rectum.f,ID.rectum.m)])

## single probe level
dat = cbind(smsM[cpgs.isl.drift5,]); Mp = nrow(dat)
dat = cbind(smsM[cpgs.down,]); Mp = nrow(dat)
age = c(age.rec); sex=c(sex.rec)
## CpG-island level1 (SMS only)
dat = out.sms.norm$Mvals[iset,]; Mp = nrow(dat)
age = age.rec; sex=sex.rec
## CpG-island level2 (SMS + GICares left colon)
dat = cbind(out.gic.left$Mvals[iset,],
            out.sms.norm$Mvals[iset,]); Mp = nrow(dat)
age = c(age.gic.left,age.rec); sex = c(sex.gic.left,sex.rec)

## ANCOVA (sexes combined)
inter = pval  = slope = numeric(Mp)

for (m in 1:Mp) {
  aux = summary(lm(dat[m,] ~ age))
  pval[m] = aux$coef[2,4]
  slope[m] = aux$coef[2,1]
  inter[m] = aux$coef[1,1]
}

## ANCOVA (sex-adjusted)
inter.m = inter.f = pval.inter.m = numeric(Mp)
pval.slope.f = pval.slope.m = slope.f = slope.m = numeric(Mp)

for (m in 1:Mp) {
  aux = summary(lm(dat[m,] ~ age*sex))
  pval.slope.f[m] = aux$coef[2,4]; pval.slope.m[m] = aux$coef[4,4]
  slope.f[m] = aux$coef[2,1]; slope.m[m] = slope.f[m]+aux$coef[4,1]
  inter.f[m] = aux$coef[1,1]
  inter.m[m] = inter.f[m]+aux$coef[3,1]; pval.inter.m[m] = aux$coef[3,4]
}

## save outputs
# SMS cpg-level
SMS.slope.f = slope.f
SMS.slope.m = slope.m
# SMS island-level
SMS.Slope.f = slope.f
SMS.Slope.m = slope.m

## look for significant differences in slopes between males and females at q-value < 0.1
tmp = qvalue(pval.slope.m,fdr.level=.1); sum(tmp$significant)
summary(tmp)

## Figure.2.R

