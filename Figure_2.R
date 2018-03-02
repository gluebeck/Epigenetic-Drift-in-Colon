#### Figure 2
#### ggplot2

library(ggplot2)
library(plyr)

len = length(cpgs.isl.drift5)
dat1  = data.frame(rate=c(SMS.slope.m,SMS.slope.f),sex=
                     c(rep("males",len),rep("females",len)),group=c(rep("probe"),len))
len = length(iset)
dat2  = data.frame(rate=c(SMS.Slope.m,SMS.Slope.f),sex=
                     c(rep("males",len),rep("females",len)),group=c(rep("island"),len))

# aux = ddply(dat2, "sex", summarise, grp.mean=mean(rate))
# head(aux)
labelSize = theme(axis.title.x = element_text(color="black", size=14, face="plain"),
                  axis.title.y = element_text(color="black", size=14, face="plain"),
                  title =        element_text(size=14, face='plain'))

ggplot(dat2, aes(x=rate, color=sex, fill=sex)) +
  # geom_density(aes(x=rate, fill=sex, color=sex),data=dat2,alpha=0.4,size=1.1) +
  # geom_histogram(aes(y=..density..), alpha=0.5, position="identity") +
  geom_density(alpha=.2,size=1.05) +
  geom_density(aes(x=rate, fill=sex, color=sex),data=dat1,alpha=0,size=1.05,linetype="dashed") +
  # geom_vline(data=aux, aes(xintercept=grp.mean, color=sex),size=1.1,linetype="solid") +
  theme(legend.position="bottom") +
  xlab("drift rate") +
  xlim(0,0.07) +
  labelSize

## save as Figure_2.pdf



