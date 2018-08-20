#### biological clocks applied to TCGA

## get horvath csv
horvath$coef = horvath$CoefficientTraining
cpgs.horvath = as.character(horvath$CpGmarker)[-1]
coef1 = horvath$coef[-1]

horvath$coef = horvath$CoefficientTrainingShrunk
coef2 = horvath$coef[-1]
coef2[is.na(coef2)] = 0

cpgs.hannum = Hannum71$Marker
coef3 = Hannum71$Coefficient

## normal tissues (Horvath)
age.normal = c(age.sms,age.gic.left,gicM.right)
dum1 = match(cpgs.horvath,rownames(smsM))
dum2 = match(cpgs.horvath,rownames(gicM.left))
dum3 = match(cpgs.horvath,rownames(gicM.right))
dat.normal.hor = ilogit2(cbind(smsM[dum1,],gicM.left[dum2,],gicM.right[dum3,]))

## normal tissues  (Hannum)
age.normal = c(age.sms,age.gic.left,age.gic.right)
dum1 = match(cpgs.hannum,rownames(smsM))
dum2 = match(cpgs.hannum,rownames(gicM.left))
dum3 = match(cpgs.hannum,rownames(gicM.right))
dat.normal.han = ilogit2(cbind(smsM[dum1,],gicM.left[dum2,],gicM.right[dum3,]))

## impute expected (mean) mAge.normal at age.tcga.left (Horvath)
dum = match(cpgs.horvath,rownames(tcgaM.left))
dat.left.cancer.hor = ilogit2(tcgaM.left[dum,])
dat.right.cancer.hor = ilogit2(tcgaM.right[dum,])
dat.cancer.hor = cbind(dat.left.cancer.hor,dat.right.cancer.hor)

## impute expected (mean) mAge.normal at age.tcga.left (Hannum)
dum = match(cpgs.hannum,rownames(tcgaM.left))
dat.left.cancer.han = ilogit2(tcgaM.left[dum,])
dat.right.cancer.han = ilogit2(tcgaM.right[dum,])
dat.cancer.han = cbind(dat.left.cancer.han,dat.right.cancer.han)

## correct dat.cancer for normal tissue fraction
ind = as.integer(age.normal/10)
p1 = ncol(dat.cancer.hor)
p2 = ncol(dat.cancer.han)
aux1 = matrix(0,ncol=p1,nrow=nrow(dat.cancer.hor))
aux2 = aux1
aux3 = matrix(0,ncol=p2,nrow=nrow(dat.cancer.han))

age = c(age.tcga.left,age.tcga.right)
fN = Ncount/100

for (i in 1:p1) {
  ind.cancer = as.integer(age[i]/10)
  if(ind.cancer==9) ind.cancer=8
  ids = which(ind==ind.cancer)
  val = apply(dat.normal.hor[,ids],1,mean)
  aux1[,i] = (dat.cancer.hor[,i]-Ncount[i]*val/100)/(1-Ncount[i]/100)
}

for (i in 1:p2) {
  ind.cancer = as.integer(age[i]/10)
  if(ind.cancer==9) ind.cancer=8
  ids = which(ind==ind.cancer)
  val = apply(dat.normal.han[,ids],1,mean)
  aux2[,i] = (dat.cancer.han[,i]-Ncount[i]*val/100)/(1-Ncount[i]/100)
}

## Horvath/Hannum
age.normal.horvath1 = 21*(0.695507258+(t(dat.normal.hor) %*% coef1)) + 20
age.normal.horvath2 = 21*(0.886919795+(t(dat.normal.hor) %*% coef2)) + 20
age.normal.hannum = t(dat.normal.han) %*% coef3

age.cancer.horvath1 = 21*(0.695507258+(t(dat.cancer.hor) %*% coef1)) + 20
#age.cancer.horvath1.adj = 21*(0.695507258+(t(aux1) %*% coef1)) + 20
## adjustment
bXn = ((age-20)/21 - 0.695507258)
bXt = ((t(dat.cancer.hor) %*% coef1) - fN*bXn)/(1-fN)
age.cancer.horvath1.adj = 21*(0.695507258+bXt) + 20

lm(age.normal ~ age.normal.horvath1)
age.cancer.horvath1.adj2 = -7.4 + 1.043*age.cancer.horvath2.adj

  
age.cancer.horvath2 = 21*(0.886919795+(t(dat.cancer.hor) %*% coef2)) + 20
#age.cancer.horvath2.adj = 21*(0.886919795+(t(aux1) %*% coef2)) + 20
## adjustment
bXn = ((age-20)/21 - 0.886919795)
bXt = ((t(dat.cancer.hor) %*% coef2) - fN*bXn)/(1-fN)
age.cancer.horvath2.adj = 21*(0.886919795+bXt) + 20

lm(age.normal ~ age.normal.horvath2)
age.cancer.horvath2.adj2 = -10.94 + 1.06*age.cancer.horvath2.adj

age.cancer.hannum =  t(dat.cancer.han) %*% coef3
# age.cancer.hannum.adj =  t(aux2) %*% coef3
## adjustment
bXn = age
bXt = ((t(dat.cancer.han) %*% coef3) - fN*bXn)/(1-fN)
age.cancer.hannum.adj = bXt

## applying calibration for normal tissue 
lm(age.normal ~ age.normal.hannum)
age.cancer.hannum.adj2 = -8.743 + 0.859 * age.cancer.hannum.adj
age.cancer.hannum.adj2 =  13.90 + 1.06 * age.cancer.hannum.adj


#### results in HannumHorvathDrift.txt



####  correct TCGA methylation levels for normal tissue contamination
####  get fractions of normal (tumor) cells for adjustment 
dat = cbind(out.tcga.left$Mvals[iset,],out.tcga.right$Mvals[iset,])
arrayID = colnames(dat)
# bcr = pdata.tcga$bcr_patient_uuid[match(arrayID,pdata.tcga$arrayID)]
tcgaID = pdata.tcga$sample_array_barcode[na.omit(match(arrayID,as.character(pdata.tcga$arrayID)))]
tcgaID = substr(tcgaID,1,16)

theta = 80
Ncount = aux1 = aux2 = rep(NA,length(tcgaID))
i=1
coadreadID = substr(as.character(slides_coadread$slide.bc),1,16)
tcf = slides_coadread$perc.tumor
ncf = slides_coadread$perc.norm
scf = slides_coadread$perc.strom
for (id in tcgaID) {
  dum0 = tcf[coadreadID==id]
  dum1 = ncf[coadreadID==id]
  dum2 = scf[coadreadID==id]
  
  dum = 100-dum0
  aux1[i] = length(dum)
  cond = dum < theta  # renove extremes
  aux2[i] = sum(!cond,na.rm=T)
  dum = dum[cond]
  Ncount[i]=mean(dum,na.rm=T); i=i+1
}

dum = mean(Ncount,na.rm=T)
Ncount[is.na(Ncount)]=dum

dum = apply(dat,2,mean)
cor.test(dum,Ncount)

#### estimated sojourn times (drift)
#### uses above tumor/normal fraction Ncount/100
#### uses slope.normal.M/F;  intercept
#### drift rate alpha_T: b.left.m/f
dat = cbind(out.tcga.left$Mvals[iset,],out.tcga.right$Mvals[iset,])
age = c(age.tcga.left,age.tcga.right)
sex = c(sex.tcga.left,sex.tcga.right)
aux = apply(dat,2,mean)

N.left.m = sum(sex.tcga.left=="M"); N.left.f = sum(sex.tcga.left=="F")
N.right.m = sum(sex.tcga.right=="M"); N.right.f = sum(sex.tcga.right=="F")
N.left = length(age.tcga.left); N.right = length(age.tcga.right)
N = length(age)
fN = Ncount/100

slope.normal = rep(slope.normal.M,N); slope.normal[sex=="F"]=slope.normal.F
slope.tumor.left =  rep(b.left.m,N.left); slope.tumor.left[sex.tcga.left=="F"]=b.left.f
slope.tumor.right =  rep(b.right.m,N.right); slope.tumor.right[sex.tcga.right=="F"]=b.right.f
slope.tumor = c(slope.tumor.left,slope.tumor.right)
intercept.normal = rep(intercept.normal.M,N); intercept.normal[sex=="F"]=intercept.normal.F

s = (aux - intercept.normal - fN*slope.normal*age)/((1-fN)*slope.tumor - fN*slope.normal)

