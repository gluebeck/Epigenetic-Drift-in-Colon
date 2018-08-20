#### drift and Relation to Island (uses manifestDATA)

dum = manifestDATA[CpGs,"Relation_to_UCSC_CpG_Island"]
dum = table(dum)

# Island N_Shelf N_Shore OpenSea S_Shelf S_Shore 
# 145733   24180   61148  172720   21692   47821 
# 
# OpenSea:
# #> 172720/473294
# #[1] 0.3649317

dum = manifestDATA[cpgs.Drift,"Relation_to_UCSC_CpG_Island"]
dum = table(dum)

# Island N_Shelf N_Shore OpenSea S_Shelf S_Shore 
# 8289     167    2095     771     109    1269 
# > sum(dum)
# [1] 12700
# 
# OpenSea:
#   > 771/12700
# [1] 0.06070866

#### identify island-associated drift (islands with minimum of 5 drift CpG probes) 
isl.drift = manifestDATA[cpgs.Drift,"UCSC_CpG_Islands_Name"]
dum = table(isl.drift)
Isl.drift5 = names(dum[dum>4])
Isl.drift5 = Isl.drift5[Isl.drift5!=""] #remove probes that are nore related to islands

cpgs.isl.drift5 = cpgs.Drift[isl.drift%in%Isl.drift5]
isl.drift5 = manifestDATA[cpgs.isl.drift5,"UCSC_CpG_Islands_Name"]

## determine genes among drift5 islands
dum = manifestDATA[cpgs.isl.drift5,"UCSC_RefGene_Name"]
Genes.isl.drift5 = unique(unlist(lapply(dum,strsplit,';')))
Genes.isl.drift5 = sort(Genes.isl.drift5)

# > length(cpgs.isl.drift5)
# [1] 7204
# > length(Isl.drift5)
# [1] 871

#### scatter plot of drift slopes (SMS) vs varRatio (Luo)
cpgs.static= setdiff(cpgs,cpgs.Drift)
cond0 = match(cpgs.static,cpgs)
## goes with Figure 1
cpgs.Static = cpgs.static[abs(slope.cpgs[cond0])<.002 & varRatio[cond0]>0] #35,751
cond1 = match(cpgs.Drift,cpgs)
cond2 = match(cpgs.Static,cpgs)

# ## Figure Slope vs VarRatio Plot.pdf
# plot(varRatio,slope,xlab='variance ratio',pch='.',col='grey',xlim=c(-4,6.5))
# abline(h=0,v=0,col=1)
# points(varRatio[cond1],slope[cond1],pch='.',col=2) 
# points(varRatio[cond2],slope[cond2],pch='.',col=4) 
# legend('bottomright',c('drift CpGs (n=8,837)','static CpGs (n=19,857)'),bty='n',pch=c(19,19),col=c(2,4))

#### identify static CGI (islands with a minimum of 5 static CpG probes on them)
isl.static = manifestDATA[cpgs.Static,"UCSC_CpG_Islands_Name"]
dum = table(isl.static)
Isl.static5 = names(dum[dum>4])
Isl.static5 = Isl.static5[Isl.static5!=""] #remove probes that are nore related to islands
cpgs.isl.static5 = cpgs.Static[isl.static%in%Isl.static5]
isl.static5 = manifestDATA[cpgs.isl.static5,"UCSC_CpG_Islands_Name"]
# gene list
dum = manifestDATA[cpgs.isl.static5,"UCSC_RefGene_Name"]
Genes.isl.static5 = unique(unlist(lapply(dum,strsplit,';')))
a = sort(Genes.isl.static5)


#### compute mean M-values for drift5 CpG-islands (requires CpGisl.level.info function)
#### normal, adenoma, cancer tissues

## normal tissues
out.sms.norm = CpGisl.level.info(Isl.drift5,isl.drift5,cpgs.isl.drift5,
                                 dat=smsM,full=T,prom=F)
out.luo.norm = CpGisl.level.info(Isl.drift5,isl.drift5,cpgs.isl.drift5,
                                 dat=luoM.norm,full=T,prom=F)
out.gic.left = CpGisl.level.info(Isl.drift5,isl.drift5,cpgs.isl.drift5,
                                 dat=gicM.left,full=T,prom=F)
out.gic.right = CpGisl.level.info(Isl.drift5,isl.drift5,cpgs.isl.drift5,
                                 dat=gicM.right,full=T,prom=F)

## tcga normals n=34
tmp=colnames(mswan)  #minfi mswan for TCGA 
bcr = pdata.tcga$bcr_patient_barcode[match(tmp,pdata.tcga$arrayID)]
ids = na.omit(match(bcr.expr.normal,bcr))
bcr.tcga.norm = bcr[ids]

dat = getM(mswan[,ids])
dat = dat[na.omit(match(cpgs.Drift,rownames(dat))),]

out.tcga.norm = CpGisl.level.info(Isl.drift5,isl.drift5,cpgs.isl.drift5,
                                 dat=dat,full=T,prom=F)

## adenoma
out.luo.aden = CpGisl.level.info(Isl.drift5,isl.drift5,cpgs.isl.drift5,
                                 dat=luoM.aden,full=T,prom=F)

## cancer (Luo)
out.luo.rect = CpGisl.level.info(Isl.drift5,isl.drift5,cpgs.isl.drift5,
                                 dat=luoM.rect,full=T,prom=F)
out.luo.left = CpGisl.level.info(Isl.drift5,isl.drift5,cpgs.isl.drift5,
                                 dat=luoM.left,full=T,prom=F)
out.luo.right = CpGisl.level.info(Isl.drift5,isl.drift5,cpgs.isl.drift5,
                                 dat=luoM.right,full=T,prom=F)

## cancer TCGA (left+rectum) n=184
cond.hist = pdata.tcga$hist == "01" & 
  pdata.tcga$history_of_neoadjuvant_treatment == "No" &
  pdata.tcga$microsatellite_instability != "YES"
cond.hist = cond.hist & grepl("Adenocarcinoma",pdata.tcga$histological_type) & 
  !grepl("Mucinous",pdata.tcga$histological_type)
 cond.hist[is.na(cond.hist)] = F

cond.rectum = pdata.tcga$anatomic_neoplasm_subdivision == "Rectum"
cond.rectum[is.na(cond.rectum)] = F

cond.left = 
  pdata.tcga$anatomic_neoplasm_subdivision == "Descending Colon" |
  pdata.tcga$anatomic_neoplasm_subdivision == "Sigmoid Colon" | 
  pdata.tcga$anatomic_neoplasm_subdivision == "Rectosigmoid Junction" |
  pdata.tcga$anatomic_neoplasm_subdivision =="Splenic Flexure"
cond.left[is.na(cond.left)] = F

## cancer TCGA (right) n=138
cond.right = pdata.tcga$anatomic_neoplasm_subdivision == "Right" |
  pdata.tcga$anatomic_neoplasm_subdivision == "Ascending Colon" |
  pdata.tcga$anatomic_neoplasm_subdivision == "Cecum" | 
  pdata.tcga$anatomic_neoplasm_subdivision == "Hepatic Flexure" |
  pdata.tcga$anatomic_neoplasm_subdivision =="Transverse Colon"
cond.right[is.na(cond.right)] = F

SID.tcga =paste(as.character(pdata.tcga$slide),'_',pdata.tcga$array,sep="")

## site-specific identfiers for TCGA data
cond = cond.hist & cond.rectum
SID.tcga.rect = SID.tcga[cond]
cond = cond.hist & cond.left
SID.tcga.left = SID.tcga[cond]
cond = cond.hist & cond.right
SID.tcga.right = SID.tcga[cond]

# cpgs.isl.drift5.tcga = intersect(cpgs.isl.drift5,rownames(mswan))
# iex = which(is.na(match(cpgs.isl.drift5,cpgs.isl.drift5.tcga)))
# isl.drift5.tcga = isl.drift5[-iex]

tcgaM.rect  = getM(mswan[,SID.tcga.rect])
tcgaM.left  = getM(mswan[,SID.tcga.left])
tcgaM.right = getM(mswan[,SID.tcga.right])

SID.tcga.leftrect = colnames(tcgaM.left)
## split out rectal and left-sided samples
ids.tcga.rect = match(SID.tcga.rect,SID.tcga.leftrect)
ids.tcga.left = match(SID.tcga.left,SID.tcga.leftrect)

## use appropriate selction 'cond' from above for rectum, left (distal), right (proximal)
age.tcga.rectum = as.numeric(pdata.tcga$age_at_initial_pathologic_diagnosis)[cond]
mmr.tcga.rectum = as.character(
  pdata.tcga$loss_expression_of_mismatch_repair_proteins_by_ihc_result)[cond]
cond.mmr.rectum = grepl("MLH1-Not Expressed|MSH2-Not Expressed|PMS2-Not Expressed",mmr.tcga.rectum)

age.tcga.left = as.numeric(pdata.tcga$age_at_initial_pathologic_diagnosis)[cond]
mmr.tcga.left = as.character(
  pdata.tcga$loss_expression_of_mismatch_repair_proteins_by_ihc_result)[cond]
cond.mmr.left = grepl("MLH1-Not Expressed|MSH2-Not Expressed|PMS2-Not Expressed",mmr.tcga.left)

age.tcga.right= as.numeric(pdata.tcga$age_at_initial_pathologic_diagnosis)[cond]
mmr.tcga.right = as.character(
  pdata.tcga$loss_expression_of_mismatch_repair_proteins_by_ihc_result)[cond]
cond.mmr.right = grepl("MLH1-Not Expressed|MSH2-Not Expressed|PMS2-Not Expressed",mmr.tcga.right)

out.tcga.left =  CpGisl.level.info(Isl.drift5,isl.drift5.tcga,cpgs.isl.drift5.tcga,
                                  dat=tcgaM.left,full=T,prom=F)
out.tcga.right = CpGisl.level.info(Isl.drift5,isl.drift5.tcga,cpgs.isl.drift5.tcga,
                                  dat=tcgaM.right,full=T,prom=F)

#### for static islands (GICares and TCGA )
static.sms.rect  = CpGisl.level.info(Isl.static5,isl.static5,cpgs.isl.static5,
                                     dat=smsM,full=T,prom=F)
static.gic.left  = CpGisl.level.info(Isl.static5,isl.static5,cpgs.isl.static5,
                                   dat=gicM.left,full=T,prom=F)
static.gic.right = CpGisl.level.info(Isl.static5,isl.static5,cpgs.isl.static5,
                                    dat=gicM.right,full=T,prom=F)

static.luo.rect = CpGisl.level.info(Isl.static5,isl.static5,cpgs.isl.static5,
                                     dat=luoM.rect,full=T,prom=F)
static.luo.left = CpGisl.level.info(Isl.static5,isl.static5,cpgs.isl.static5,
                                    dat=luoM.left,full=T,prom=F)
static.luo.right = CpGisl.level.info(Isl.static5,isl.static5,cpgs.isl.static5,
                                     dat=luoM.right,full=T,prom=F)

static.tcga.left = CpGisl.level.info(Isl.static5,isl.static5,cpgs.isl.static5,
                                     dat=tcgaM.left,full=T,prom=F)
static.tcga.right = CpGisl.level.info(Isl.static5,isl.static5,cpgs.isl.static5,
                                      dat=tcgaM.right,full=T,prom=F)

####  overall correlation enrichemnet by significance for TCGA left and right
iset.tcga.left  = CorPurify(dat=out.tcga.left$Mvals,alpha=c(seq(0.8,0.9,.1),.95))
iset.tcga.right = CorPurify(dat=out.tcga.right$Mvals,alpha=c(seq(0.8,0.9,.1),.95))
iset = iset.tcga = intersect(iset.tcga.left,iset.tcga.right)
genes = out.tcga.left$genes[iset]
a = sort(unique(unlist(genes)))

#### plot clock metric for normal and neoplastic data

## normal tissues (SMS, GIC, ...)
age = c(age.gic.left,age.rec)
sex = c(sex.gic.left,sex.rec)
dat = cbind(out.gic.left$Mvals[iset,],out.sms.norm$Mvals[iset,])

aux = apply(dat,2,function(x) {median(x[x > -2])})
# aux = apply(dat,2,quantile,probs = 0.5)
cor.test(age,aux)
# ANCOVA in ValidateNormalDrift.R
plot(age,aux,pch=19,cex=.7,xlim=c(20,90),ylim=c(-2.5,2.5),col='grey',xlab='age',ylab='M-value')
# points(age,aux,pch=19,cex=0.7,col='grey')
tmp = lm(aux~age); summary(tmp)
abline(tmp)

## Luo, Adenomas and Carcinomas 
age = c(age.luo.left,age.luo.right,age.luo.adenoma)
dat = cbind(out.luo.left$Mvals[iset,],out.luo.right$Mvals[iset,],out.luo.aden$Mvals[iset,])
aux = apply(dat,2,function(x) {mean(x[x > -2])})
cols = rep(4,ncol(dat)); cols[94:(94-18)]=3
# aux = apply(dat,2,quantile,probs = 0.5)
# plot(age,aux,pch=19,cex=.7,xlim=c(20,90),ylim=c(0,500),col = col)
points(age,aux,pch='+',col=cols)
cor.test(age,aux)

## Left Colorectal Adenoca (TCGA, n=184)
cond = cond.hist & cond.left
dum = SID.tcga[cond]

age.tcga.left = as.numeric(pdata.tcga$age_at_initial_pathologic_diagnosis)[cond]
sex.tcga.left = pdata.tcga$gender[cond]
sex.tcga.left[sex.tcga.left=="MALE"] = "M"
sex.tcga.left[sex.tcga.left=="FEMALE"] = "F"

dat = cbind(out.luo.left$Mvals[iset,],out.tcga.left$Mvals[iset,])
age = c(age.luo.left,age.tcga.left); sex = c(sex.luo.left,sex.tcga.left)
aux = apply(dat,2,function(x) {median(x[x > -2.5])})
cor.test(age,aux)

col=rep(1,length(aux));col[1:32] = 3
# aux = apply(dat,2,quantile,probs = 0.5)
# plot(age,aux,pch=19,cex=.7,xlim=c(20,90),ylim=c(-3.5,2),col = col)
points(age, aux, pch=19,cex=.5,col='orange')
lines(loess.smooth(age,aux),lwd=3,col='orange')
legend('topleft',c("Normal Colon","TCGA","Luo"),col=c('grey','orange','blue'),pch=c(19,19,19),
       bty='n',y.intersp = 1.3)

## Right Colon Adenoca (TCGA, n=138)
cond = cond.hist & cond.right
dum = SID.tcga[cond]

sex.tcga.right = pdata.tcga$gender[cond]
sex.tcga.right[sex.tcga.right=="MALE"] = "M"
sex.tcga.right[sex.tcga.right=="FEMALE"] = "F"

# iex = (1:138)[pdata.tcga$anatomic_neoplasm_subdivision[cond]=="Cecum"]

####  check significance of correlation with age for every island in dat
dat = cbind(out.luo.right$Mvals[iset,],out.tcga.right$Mvals[iset,])
age = c(age.luo.right,age.tcga.right); sex = c(sex.luo.right,sex.tcga.right)

dat = cbind(out.luo.left$Mvals[iset,],out.tcga.left$Mvals[iset,])
age = c(age.luo.left,age.tcga.left); sex = c(sex.luo.left,sex.tcga.left)

Mp = nrow(dat)
dum1 = dum2 = numeric(Mp)
for (m in 1:Mp) {
  tmp = cor.test(age,dat[m,])
  dum1[m] = tmp$est
  dum2[m] = tmp$p.value
}
dum = qvalue(dum2,fdr.level=.01); sum(dum$significant)
summary(dum)

iset.significant.right = iset[dum$qvalues<0.01]
iset.significant.left = iset[dum$qvalues<0.01]  #!!!! will not give q-values

#### optional: improve correlation of age with island-level methylation  
# sum(dum1.l>0.2 & dum1.r>0.2)
# iset.significant = iset[dum1.l>0.2 & dum1.r>0.2]
# length(iset.significant)
# 
# dat = out.tcga.right$Mvals[iset.significant,] #-iex]
# age = c(age.tcga.right); sex = c(sex.tcga.right)
# # age = c(age.tcga.right[-iex]); sex = c(sex.tcga.right[-iex])
# 
# dat = cbind(out.luo.right$Mvals[iset.significant,],out.tcga.right$Mvals[iset.significant,])
# age = c(age.luo.right,age.tcga.right); sex = c(sex.luo.right,sex.tcga.right)
# aux = apply(dat,2,function(x) {median(x[x > -2.5])})
# cor.test(age,aux)

