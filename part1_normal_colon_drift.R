#### 1. identify CpG probes that undergo drift in SMS (rectal samples) jointly for males and females
#### 2. validate in GICares set for left-sided cancers

## only autosomes
cpgAUTO = manifestDATA$Name[!grepl("X|Y" ,manifestDATA$CHR)]
CpGs = intersect(rownames(mset.swan),cpgAUTO)

#### step1: filter out 'hypomethylated' probes with beta < 50% 
#### use M-vales
theta = logit2(0.5); rmax = apply(smsM,1,max)
cpgs = CpGs[rmax < theta]

#### regressing methylation on age to check for significant age effects in SMS for all cpgs
age.rec = c(age.rec.f,age.rec.m)
sex.rec = c(rep("F",90),rep("M",60)) # for later use

dat = getM(mset.swan[cpgs,c(ID.rectum.f,ID.rectum.m)])
age = age.rec
Mp = nrow(dat)

# initialize vectors
pval = intercept = slope = numeric(Mp)

for (m in 1:Mp) {
  aux = summary(lm(dat[m,] ~ age))
  pval[m] = aux$coef[2,4]; intercept[m] = aux$coef[1,1]; slope[m] = aux$coef[2,1]
}

# control FDR (needs qvalue libarry)
tmp = qvalue(pval.cpgs,fdr.level=.0001); sum(tmp$significant)
summary(tmp)

cpgs.drift = cpgs[tmp$significant]
cpgs.static= cpgs[!tmp$significant]

## plot some Luo and SMS data for cpgs.drift or cpgs.static. Eg. drift:
# for ( i in 1:100) {
#   cpg = cpgs.drift[0+i]
#   plot(age.rec,smsM[cpg,],pch=19,cex=.6,ylim=c(-10,2),xlab='age',ylab='M-value')
#   points(age.luo.normal,luoM.norm[cpg,],pch=19,cex=0.6,col=4)
#   points(age.luo.adenoma,luoM.aden[cpg,],pch=15,cex=1,col='green')
#   points(age.luo.cancer,luoM.carc[cpg,],pch=15,cex=1,col='red')
#   legend('bottomleft',paste('drift ',cpg))
#   Sys.sleep(1)
# }

#### step2: validation using GICares data gicM.left
## filter out hypomethylated probes with beta < 50% in SMS
theta = logit2(0.5); rmax = apply(gicM.left,1,max)
# cpgs.gic = CpGs[rmax < theta]
# cpgs.gic = cpgs

## regressing methyl on age to check for significant age effect
dat = gicM.left[cpgs,]
age = age.gic.left
Mp = nrow(dat)

## initialize vectors for regression
pval.gic = intercept.gic = slope.gic = numeric(Mp)

for (m in 1:Mp) {
  aux = summary(lm(dat[m,] ~ age))
  pval.gic[m] = aux$coef[2,4]; intercept.gic[m] = aux$coef[1,1]; slope.gic[m] = aux$coef[2,1]
}

## control FDR (needs qvalue libarry)
tmp = qvalue(pval.gic,fdr.level=.05); sum(tmp$significant)
summary(tmp)

cond = (tmp$pvalues < 0.05) # validate if significant at 0.05 p-value nominal
cpgs.drift.gic = cpgs[cond]
ids = na.omit(match(cpgs.drift,cpgs.drift.gic))
iex = which(slope.gic[cond][ids] < 0) # validate if in addition > slope

cpgs.Drift = cpgs.drift.gic[ids][-iex]  #12,700

# # confirm
# > length(cpgs.Drift)
# [1] 12700


