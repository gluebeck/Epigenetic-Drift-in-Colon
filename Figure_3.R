len = length(iset)
x = factor(c(rep("SMS males",len),rep("SMS females",len),rep("GIC males",len),
             rep("GIC females",len),rep("Right Colon",len)), 
           levels = c("SMS males","SMS females","GIC males","GIC females","Right Colon"), ordered = TRUE)

# x = factor(c(rep("SMS males",len),rep("SMS females",len),rep("GIC males",len),
#              rep("GIC females",len),rep("Right Colon males",len), rep("Right Colon females",len)), 
#            levels = c("SMS males","SMS females","GIC males","GIC females","Right Colon males",
#                       "Right Colon females"), ordered = TRUE)

y = c(SMS.Slope.m,SMS.Slope.f,GIC.Slope.m,GIC.Slope.f,Slope.right)
sex = c(rep("M",len),rep("F",len),rep("M",len),rep("F",len),rep("N",len)) #,rep("F",len))
labels = c("SMS males\n (n=60)","SMS females\n (n=90)","GIC males\n (n=31)","GIC females\n (n=37)",
           "Right Colon\n (n=14)") #,"Right Colon females\n (n=8)")

datf = data.frame(set=x,slopes=y,sex=sex)

labelSize = theme(axis.title.x = element_text(color="black", size=14, face="plain"),
                  axis.title.y = element_text(color="black", size=14, face="plain"),
                  title =        element_text(size=14, face='plain'))

  ggplot(datf, aes(x=set, y=slopes)) +
  # geom_boxplot(aes(colour=sex)) + geom_jitter(aes(colour=sex),width = 0.15, size=.3) + 
  geom_boxplot() + geom_jitter(aes(colour=sex),width = 0.15, size=.3) +
  # scale_fill_hue(l=40, c=75) +
  scale_color_manual(values=c("#E69F00", "#56B4E9","#999999")) +
  xlab("data set") +
  scale_x_discrete(labels=labels) +
  ylab("drift rates") +
  theme(legend.position="none") +
  labelSize



