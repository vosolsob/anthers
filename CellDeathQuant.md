---
Title: Anther dehiscence
Author: "Anna Kampová, Jan Petrášek & Stanislav Vosolsobě"
output:
  html_document:
  df_print: paged
word_document: default
pdf_document: default
---
  
# Quantification of cell death during anther dehiscence
  
### By Stanislav Vosolsobě
  
## Libraries and data import

```r
require(betareg) # for beta regression
require(lmtest)  # for LR test after BetaReg
require(StepBeta) # for step after BetaReg
require(plotfunctions) # for gradient legend

setwd("/media/Home/home/standa/Plocha/Anička/konfo/FML_PCD/")
fml <- read.table("FML_PCD",header = T)
```

## Analysis of non-connective tissues

```r
subfml <- fml[fml$side!="connective",]
br <- betareg(data=subfml,mortality~time*region*side*treatment*tissue)
br_red <- StepBeta(br)
# reduced model computed again, because of error in 'predict'
br_red <-  betareg(formula =  mortality ~ time + region + tissue + treatment + side + time:region:side:treatment:tissue + treatment:tissue + time:side + time:treatment + time:treatment:tissue, data = subfml )

# Testing of the effect of treatment
br_notrt <- betareg(formula =  mortality ~ time + region + tissue + side + time:region:side:tissue + time:side + time:tissue, data = subfml )
lrtest(br_red,br_notrt)  # < 2.2e-16 ***

# Testing of the effect of tissue
br_notis <- betareg(formula =  mortality ~ time + region + treatment + side + time:region:side:treatment + time:side + time:treatment, data = subfml )
lrtest(br_red,br_notis)  # < 2.2e-16 ***

# Testing of the effect of side
br_noside <- betareg(formula =  mortality ~ time + region + tissue + treatment + time:region:treatment:tissue + treatment:tissue + time:treatment + time:treatment:tissue, data = subfml )
lrtest(br_red,br_noside)  # 6.345e-08 ***

# Testing of the effect of region
br_noreg <- betareg(formula =  mortality ~ time * tissue * treatment * side, data = subfml )
lrtest(br,br_noreg)  # < 2.2e-16 ***, unable to compute reduced model without 'region'
```

## Analysis of the connective

```r
confml <- fml[fml$side=="connective",]

conbr <- betareg(data=confml,mortality~time*treatment*tissue)
conbr_red <- StepBeta(conbr)

conbr_red <- betareg(formula =  mortality ~ treatment + time, data = confml )

# Testing of the effect of treatment
conbr_notrt <- betareg(formula =  mortality ~ time, data = confml )
lrtest(conbr_red,conbr_notrt) # 0.07804 .

# Testing of the effect of time
conbr_notime <- betareg(formula =  mortality ~ treatment, data = confml )
lrtest(conbr_red,conbr_notime) # 0.1536

# Comparison with null model
lrtest(conbr_red)  # 0.08991 .
```

## Creation of colour scale for diagrams

```r
pdf("FML_PCD-scale.pdf",height = 10,width = 15)

plot(type="n",1,xlim=c(0,10),ylim=c(-1,1),xlab="SHORT: 0 epi/endo - 20 epi/endo | LONG: 0 epi/endo - 20 epi/endo",xaxt="n",ylab="ABAXIAL  | ADAXIAL", yaxt="n")
rc <- colorRampPalette(colors = c("white","gold","red"), space = "rgb")(100)
cex=0   # suppression of text

pos <- c(0,0,1,1)
### ADAXIAL
## SHORT
# TIME 0
prd <- predict(br_red, newdata=data.frame(time=rep(0,2), region=c(1,3), tissue=rep("epidermis",2), treatment=rep("short",2), side=rep("adaxial",2)), type="response")
gradientLegend(valRange = c(0,1)/100,color = rc[round(100*prd[1]):round(100*prd[2])],inside = T,side = 2,n.seg = 1,dec = 0,tick.col = NA,border.col = NA,pos = pos,coords = T)
pos <- pos + c(1,0,1,0)

prd <- predict(br_red, newdata=data.frame(time=rep(0,2), region=c(1,3), tissue=rep("endothecium",2), treatment=rep("short",2), side=rep("adaxial",2)), type="response")
gradientLegend(valRange = c(0,1)/100,color = rc[round(100*prd[1]):round(100*prd[2])],inside = T,side = 2,n.seg = 1,dec = 0,tick.col = NA,border.col = NA,pos = pos,coords = T)
pos <- pos + c(1,0,1,0)

# TIME 20
prd <- predict(br_red, newdata=data.frame(time=rep(20,2), region=c(1,3), tissue=rep("epidermis",2), treatment=rep("short",2), side=rep("adaxial",2)), type="response")
gradientLegend(valRange = c(0,1)/100,color = rc[round(100*prd[1]):round(100*prd[2])],inside = T,side = 2,n.seg = 1,dec = 0,tick.col = NA,border.col = NA,pos = pos,coords = T)
pos <- pos + c(1,0,1,0)

prd <- predict(br_red, newdata=data.frame(time=rep(20,2), region=c(1,3), tissue=rep("endothecium",2), treatment=rep("short",2), side=rep("adaxial",2)), type="response")
gradientLegend(valRange = c(0,1)/100,color = rc[round(100*prd[1]):round(100*prd[2])],inside = T,side = 2,n.seg = 1,dec = 0,tick.col = NA,border.col = NA,pos = pos,coords = T)
pos <- pos + c(1,0,1,0)

## LONG
# TIME 0
prd <- predict(br_red, newdata=data.frame(time=rep(0,2), region=c(1,3), tissue=rep("epidermis",2), treatment=rep("long",2), side=rep("adaxial",2)), type="response")
gradientLegend(valRange = c(0,1)/100,color = rc[round(100*prd[1]):round(100*prd[2])],inside = T,side = 2,n.seg = 1,dec = 0,tick.col = NA,border.col = NA,pos = pos,coords = T)
pos <- pos + c(1,0,1,0)

prd <- predict(br_red, newdata=data.frame(time=rep(0,2), region=c(1,3), tissue=rep("endothecium",2), treatment=rep("long",2), side=rep("adaxial",2)), type="response")
gradientLegend(valRange = c(0,1)/100,color = rc[round(100*prd[1]):round(100*prd[2])],inside = T,side = 2,n.seg = 1,dec = 0,tick.col = NA,border.col = NA,pos = pos,coords = T)
pos <- pos + c(1,0,1,0)

# TIME 20
prd <- predict(br_red, newdata=data.frame(time=rep(20,2), region=c(1,3), tissue=rep("epidermis",2), treatment=rep("long",2), side=rep("adaxial",2)), type="response")
gradientLegend(valRange = c(0,1)/100,color = rc[round(100*prd[1]):round(100*prd[2])],inside = T,side = 2,n.seg = 1,dec = 0,tick.col = NA,border.col = NA,pos = pos,coords = T)
pos <- pos + c(1,0,1,0)

prd <- predict(br_red, newdata=data.frame(time=rep(20,2), region=c(1,3), tissue=rep("endothecium",2), treatment=rep("long",2), side=rep("adaxial",2)), type="response")
gradientLegend(valRange = c(0,1)/100,color = rc[round(100*prd[1]):round(100*prd[2])],inside = T,side = 2,n.seg = 1,dec = 0,tick.col = NA,border.col = NA,pos = pos,coords = T)
pos <- pos + c(1,0,1,0)

### ABAXIAL
pos <- c(0,-1,1,0)
## SHORT
# TIME 0
prd <- predict(br_red, newdata=data.frame(time=rep(0,2), region=c(1,3), tissue=rep("epidermis",2), treatment=rep("short",2), side=rep("abaxial",2)), type="response")
gradientLegend(valRange = c(0,1)/100,color = rc[round(100*prd[1]):round(100*prd[2])],inside = T,side = 2,n.seg = 1,dec = 0,tick.col = NA,border.col = NA,pos = pos,coords = T)
pos <- pos + c(1,0,1,0)

prd <- predict(br_red, newdata=data.frame(time=rep(0,2), region=c(1,3), tissue=rep("endothecium",2), treatment=rep("short",2), side=rep("abaxial",2)), type="response")
gradientLegend(valRange = c(0,1)/100,color = rc[round(100*prd[1]):round(100*prd[2])],inside = T,side = 2,n.seg = 1,dec = 0,tick.col = NA,border.col = NA,pos = pos,coords = T)
pos <- pos + c(1,0,1,0)

# TIME 20
prd <- predict(br_red, newdata=data.frame(time=rep(20,2), region=c(1,3), tissue=rep("epidermis",2), treatment=rep("short",2), side=rep("abaxial",2)), type="response")
gradientLegend(valRange = c(0,1)/100,color = rc[round(100*prd[1]):round(100*prd[2])],inside = T,side = 2,n.seg = 1,dec = 0,tick.col = NA,border.col = NA,pos = pos,coords = T)
pos <- pos + c(1,0,1,0)

prd <- predict(br_red, newdata=data.frame(time=rep(20,2), region=c(1,3), tissue=rep("endothecium",2), treatment=rep("short",2), side=rep("abaxial",2)), type="response")
gradientLegend(valRange = c(0,1)/100,color = rc[round(100*prd[1]):round(100*prd[2])],inside = T,side = 2,n.seg = 1,dec = 0,tick.col = NA,border.col = NA,pos = pos,coords = T)
pos <- pos + c(1,0,1,0)

## LONG
# TIME 0
prd <- predict(br_red, newdata=data.frame(time=rep(0,2), region=c(1,3), tissue=rep("epidermis",2), treatment=rep("long",2), side=rep("abaxial",2)), type="response")
gradientLegend(valRange = c(0,1)/100,color = rc[round(100*prd[1]):round(100*prd[2])],inside = T,side = 2,n.seg = 1,dec = 0,tick.col = NA,border.col = NA,pos = pos,coords = T)
pos <- pos + c(1,0,1,0)

prd <- predict(br_red, newdata=data.frame(time=rep(0,2), region=c(1,3), tissue=rep("endothecium",2), treatment=rep("long",2), side=rep("abaxial",2)), type="response")
gradientLegend(valRange = c(0,1)/100,color = rc[round(100*prd[1]):round(100*prd[2])],inside = T,side = 2,n.seg = 1,dec = 0,tick.col = NA,border.col = NA,pos = pos,coords = T)
pos <- pos + c(1,0,1,0)

# TIME 20
prd <- predict(br_red, newdata=data.frame(time=rep(20,2), region=c(1,3), tissue=rep("epidermis",2), treatment=rep("long",2), side=rep("abaxial",2)), type="response")
gradientLegend(valRange = c(0,1)/100,color = rc[round(100*prd[1]):round(100*prd[2])],inside = T,side = 2,n.seg = 1,dec = 0,tick.col = NA,border.col = NA,pos = pos,coords = T)
pos <- pos + c(1,0,1,0)

prd <- predict(br_red, newdata=data.frame(time=rep(20,2), region=c(1,3), tissue=rep("endothecium",2), treatment=rep("long",2), side=rep("abaxial",2)), type="response")
gradientLegend(valRange = c(0,1)/100,color = rc[round(100*prd[1]):round(100*prd[2])],inside = T,side = 2,n.seg = 1,dec = 0,tick.col = NA,border.col = NA,pos = pos,coords = T)
pos <- pos + c(1,0,1,0)

# General legend bar
gradientLegend(valRange = c(0,1),color = rc[1:100],inside = T,side = 2,n.seg = 3,dec = 1,tick.col = NA,border.col = NA,pos = c(9,0.2,9.5,0.8),coords = T)

# For connective
points(x=9,y=-0.5,col=rc[round(100*mean(confml$mortality))],cex=15,pch=16)

dev.off()
```
