#### slope vs bval (given mval) for all probes tested
bval = ilogit2(mval)
plot(bval,slope,pch='.',col='grey',xlab='beta-value')
cond1 = match(cpgs.Drift,cpgs)
cond2 = match(cpgs.Static,cpgs)
points(bval[cond1],slope[cond1],pch='.',col=2)
points(bval[cond2],slope[cond2],pch='.',col=4)

#### get Variance Ratio using Luo data (cancer vs normal)
#### use set cpgs (all autosomal CpG probes < 50% methylation in SMS)
dat = luoM.norm[cpgs,]
sval1 = apply(dat,1,var)
dat = cbind(luoM.left[cpgs,],luoM.right[cpgs,])
sval2 = apply(dat,1,var)
varRatio = log2((sval2/sval1))

## slppes from M-value~age regressions in SMS data for set cpgs
slope = slope.cpgs

## get a control group (non-drift)
cpgs.static= setdiff(cpgs,cpgs.Drift)
cond0 = match(cpgs.static,cpgs)
cpgs.Static = cpgs.static[abs(slope[cond0])<.002 & varRatio[cond0]>0]

## Figure Slope vs VarRatio Plot.pdf
par(bg="grey")
plot(varRatio,slope,xlab='variance ratio',pch='.',col='grey',xlim=c(-3,6.5),cex.lab=1.3)
abline(h=0,v=0,col=1)
cond1 = match(cpgs.Drift,cpgs)
cond2 = match(cpgs.Static,cpgs)
points(varRatio[cond1],slope[cond1],pch='.',col=2) 
points(varRatio[cond2],slope[cond2],pch='.',col=4) 
legend('bottomright',c('drift CpGs (n=13,525)','static CpGs (n=35,751)'),
       bty='n',pch=c(19,19),col=c(2,4))

## in ggplot
dat = data.frame(varRatio=varRatio,rate=slope,group=rep("static CpGs",length(cpgs)))
dat$group = as.character(dat$group)
cond1 = match(cpgs.Drift,cpgs)
cond2 = match(cpgs.Static,cpgs)
dat$group[cond1]="drift CpGs"
dat$group[cond2]="controls"

labelSize = theme(axis.title.x = element_text(color="black", size=14, face="plain"),
                  axis.title.y = element_text(color="black", size=14, face="plain"),
                  title =        element_text(size=14, face='plain'))

ggplot(dat, aes(x=varRatio)) +
  geom_point(aes(y=rate,color=group), size=.1) +
  scale_color_manual(values=c("lightgreen", "red", "orange")) +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  theme(legend.position="bottom") +
  geom_hline(yintercept = 0) +
  # labs(title = "...") +
  xlim(-3,7.5) +
  ylim(-.05,.1) +
  xlab("log2 varRatio") +
  ylab("drift rate") +
  labelSize

