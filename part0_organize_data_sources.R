#### Luo data

pdata.luo = pData(mluoset$mset.swan)
mset.luo = mluoset$mset.swan
cond.luo.rect = "Rectum|rectum|Retum"
cond.luo.left = "oid|Left|Desc|distal"
cond.luo.right = "Right|right|cec|Cec|Trans|trans"

ID.luo.normal = (1:124)[pdata.luo$Histololgy=='normal']
ID.luo.adenoma = (1:124)[pdata.luo$Histololgy=='adenoma']
ID.luo.cancer = (1:124)[pdata.luo$Histololgy=='cancer']

ID.luo.rect = (1:124)[pdata.luo$Histololgy!='normal' 
                      & grepl(cond.luo.rect,pdata.luo$Site.size)]
ID.luo.left = (1:124)[pdata.luo$Histololgy!='normal' 
                      & grepl(cond.luo.left,pdata.luo$Site.size)]
ID.luo.right = (1:124)[pdata.luo$Histololgy!='normal'
                       & grepl(cond.luo.right,pdata.luo$Site.size)]

age.luo = as.numeric(pdata.luo$Age)
age.luo.adenoma = age.luo[ID.luo.adenoma]
age.luo.normal = age.luo[ID.luo.normal]
age.luo.cancer = age.luo[ID.luo.cancer]

age.luo.rect = age.luo[ID.luo.rect]
age.luo.left = age.luo[ID.luo.left]
age.luo.right = age.luo[ID.luo.right]

sex.luo = toupper(as.character(pdata.luo$Gender))
sex.luo.adenoma = sex.luo[ID.luo.adenoma]
sex.luo.normal = sex.luo[ID.luo.normal]
sex.luo.cancer = sex.luo[ID.luo.cancer]

sex.luo.rect = sex.luo[ID.luo.rect]
sex.luo.left = sex.luo[ID.luo.left]
sex.luo.right = sex.luo[ID.luo.right]

luoM.norm = getM(mset.luo[,ID.luo.normal])
luoM.aden = getM(mset.luo[,ID.luo.adenoma])
luoM.carc = getM(mset.luo[,ID.luo.cancer])

luoM.rect = getM(mset.luo[,ID.luo.rect])
luoM.left = getM(mset.luo[,ID.luo.left])
luoM.right= getM(mset.luo[,ID.luo.right])


#### SMS data (needs object pdata for SMS)
ID.rectum.f = rownames(pdata)[pdata$Sample_Group=="SMS" & pdata$Gender=="F"]
ID.rectum.m = rownames(pdata)[pdata$Sample_Group=="SMS" & pdata$Gender=="M"]

# ID.rectum.f.polyp = rownames(pdata)[pdata$Sample_Group=="SMS" & pdata$Gender=="F" 
#                                     & pdata$Disease_status=="polyps_present"]
# ID.rectum.f.nopol = rownames(pdata)[pdata$Sample_Group=="SMS" & pdata$Gender=="F"
#                                     & pdata$Disease_status=="no_polyps"]
# ID.rectum.m.polyp = rownames(pdata)[pdata$Sample_Group=="SMS" & pdata$Gender=="M" 
#                                     & pdata$Disease_status=="polyps_present"]
# ID.rectum.m.nopol = rownames(pdata)[pdata$Sample_Group=="SMS" & pdata$Gender=="M"
#                                     & pdata$Disease_status=="no_polyps"]

#### GICares data

#### ColoGICares data 


#### define probe set to study

## only autosomes
cpgAUTO = manifestDATA$Name[!grepl("X|Y" ,manifestDATA$CHR)]
CpGs = intersect(rownames(mset.swan),cpgAUTO)

## methylation data for SMS
smsM = getM(mset.swan[CpGs,c(ID.rectum.f,ID.rectum.m)])
age.rec.f = as.numeric(pdata[ID.rectum.f,"Age"]) 
age.rec.m = as.numeric(pdata[ID.rectum.m,"Age"]) 

## methylation data for GICares (left n=68, right=16)
SID.GICaRes = rownames(pdata)[pdata$Sample_Group=="GICaRes" & pdata$colon_side=="Left"]
gicM.left = getM(mset.swan[CpGs,SID.GICaRes])
age.gic.left = as.numeric(pdata$Age[match(SID.GICaRes,rownames(pdata))])
sex.gic.left = pdata$Gender[match(SID.GICaRes,rownames(pdata))]

SID.GICaRes = rownames(pdata)[pdata$Sample_Group=="GICaRes" & pdata$colon_side=="Right"]
gicM.right = getM(mset.swan[CpGs,SID.GICaRes])
age.gic.right = as.numeric(pdata$Age[match(SID.GICaRes,rownames(pdata))])
sex.gic.right = pdata$Gender[match(SID.GICaRes,rownames(pdata))]


