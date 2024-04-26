---
Title: Dew treatment of anthers
Author: "Anna Kampová, Jan Petrášek & Stanislav Vosolsobě"
output:
  html_document:
  df_print: paged
word_document: default
pdf_document: default
---
  
# Statistical analysis of anther dehiscence after dew treatment
  
### By Stanislav Vosolsobě
  
  
This script works with input format of this type

```r
time	treatment	AZ	AI	AO
8	    DEW	      49	1	   0
8	    CTR	       0	1	  49
9	    DEW	      48	2	   0
9	    CTR	       0	0	  50
.     ...       ..  .   ..
```



## Read input data

```r
ath <- read.table("at_roseni",header=T)
```

## Compute GLM model

These functions create full model and optimize it by `step` procedure.

### For anther splitting

```r
mAI <- glm(cbind(AO+AI,AZ) ~ time * treatment, data = ath, family = "binomial")
step(mAI) # interaction removed by AIC
smAI <- step(mAI,trace=0)  # suppress print by trace
print(coef(smAI))
print(anova(smAI, test = 'Chisq')) 
# final model:
mAI <- glm(cbind(AO+AI,AZ) ~ time + treatment, data = ath, family = "binomial")
```

### For anther opening

```r
mAO <- glm(cbind(AO,AI+AZ) ~ time * treatment, data = ath, family = "binomial")
step(mAO) # interaction removed by AIC
smAO <- step(mAO,trace=0)
print(coef(smAO))
print(anova(smAO, test = 'Chisq'))
# final model
mAO <- glm(cbind(AO,AI+AZ) ~ time + treatment, data = ath, family = "binomial")
```

## Draw diagram

```r
pdf("GLM_ath.pdf",width=6,height = 5)
x <- c(8:11,14,20)
plot(1,type="n",xlim=c(8,20),ylim=c(-0.1,1.1),xaxt="n",xlab="time [hours]",
     ylab=expression("Probability of developmental transition"),main="Anther development")
abline(v=x, col="gray")
axis(side = 1, at = x, labels = x)
tp <- (80:200)/10

AIne <- predict(mAI,newdata = data.frame(time=tp,treatment=rep("CTR",length(tp))),type="response",se.fit = T)
AOne <- predict(mAO,newdata = data.frame(time=tp,treatment=rep("CTR",length(tp))),type="response",se.fit = T)
AIano <- predict(mAI,newdata = data.frame(time=tp,treatment=rep("DEW",length(tp))),type="response",se.fit = T)
AOano <- predict(mAO,newdata = data.frame(time=tp,treatment=rep("DEW",length(tp))),type="response",se.fit = T)


lines(AIne$fit~tp,lwd=2,col="orange")
lines(AOne$fit~tp,lwd=2,col="sienna")
lines(AIano$fit~tp,lwd=2,col="orange")
lines(AOano$fit~tp,lwd=2,col="sienna")

lines(AIne$fit+1.96*AIne$se.fit~tp,lwd=1,col="orange",lty=2)
lines(AOne$fit+1.96*AOne$se.fit~tp,lwd=1,col="sienna",lty=2)
lines(AIano$fit+1.96*AIano$se.fit~tp,lwd=1,col="orange",lty=2)
lines(AOano$fit+1.96*AOano$se.fit~tp,lwd=1,col="sienna",lty=2)

lines(AIne$fit-1.96*AIne$se.fit~tp,lwd=1,col="orange",lty=2)
lines(AOne$fit-1.96*AOne$se.fit~tp,lwd=1,col="sienna",lty=2)
lines(AIano$fit-1.96*AIano$se.fit~tp,lwd=1,col="orange",lty=2)
lines(AOano$fit-1.96*AOano$se.fit~tp,lwd=1,col="sienna",lty=2)

points(ath$AO/50~ath$time,col="sienna")
points((ath$AI+ath$AO)/50~ath$time,col="orange")

legend(x=15,y=0.7,col = c("orange","sienna"),lwd=2, legend = c("Anther splitting","Anther opening"))

dev.off()
```

