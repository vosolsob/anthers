---
  Title: Anther dehiscence
Author: "Anna Kampová, Jan Petrášek & Stanislav Vosolsobě"
output:
  html_document:
  df_print: paged
word_document: default
pdf_document: default
---
  
# Statistical analysis of anther dehiscence
  
### By Stanislav Vosolsobě
  


      


```r
setwd("/media/Home/home/standa/Plocha/Anička/Mutanti_kveteni/")


dta <- import_horiz(filename = "mutanti_join",header=T, dt=3)
anthers(ath = dta,formula = c(AZ, AI, AO) ~ t * I(t^2) * genotyp * variant, emformula = ~ variant | genotyp, dt = 3,col = c("blue","red","coral","darkred","gold"),pdf = F,name = "mutanti_join")
        
dta_wt <- droplevels(dta[dta$genotyp=="WT",])
anthers(ath = dta_wt,formula = c(AZ, AI, AO) ~ t * I(t^2) * variant, emformula = ~ variant, dt = 3,col = c("darkgreen","sienna"),pdf = F,name = "attach")

dta_mut <- droplevels(dta[dta$variant=="intact",])
anthers(ath = dta_mut,formula = c(AZ, AI, AO) ~ t * I(t^2) * genotyp, emformula = ~ genotyp, dt = 3,col = c("blue","red","coral","darkred","gold"),pdf = F,name = "mutanti")
```


Function for import of dataset in long horizontal format, e.g.:
  
```r
genotyp variant flowers anthers AZ  AI  AO  AZ  AI  AO  AZ  AI  AO  AZ  AI  AO
WT      cut     30       120   120  0   0  118   2   0  114  4   2  10 	13   4
```


```r
import_horiz <- function(filename,   # input data
                         header=T,   # header in dataset
                         dt          # interval between observation
                         ){
  dta <- read.table(filename, header = header)
  az <- which(colnames(dta)=="AZ")
  nt <- (1 + ncol(dta) - az )/3
  dfl <- nrow(dta)*nt
  ath <- setNames(data.frame(matrix(ncol = 1+which(colnames(dta)=="AO"), nrow = dfl)), c("t",colnames(dta)[1:which(colnames(dta)=="AO")])   )
  k <- 1
  for (i in 1:nrow(dta)){
    for(j in 0:(nt-1)){
      ath$t[k] <- j*dt + dt/2
      for(v in 1:(az-1)){
        ath[k,v+1] <- dta[i,v]
      }
      ath[k,az+1] <- dta[i,3*j+az]
      ath[k,az+2] <- dta[i,3*j+az+1]
      ath[k,az+3] <- dta[i,3*j+az+2]
      k <- k + 1
    }
  }
  return(ath)
}
```

Function for the logistic regression

```r
anthers <- function(ath,   # input data, data frame with rows: time, variantID, # of closed , # of initiated and # of open anthers
                    formula,  # input model formula, time variable must be at the first place in right part, in left hand side must be in order 'closed', 'initiated', 'open'
                    emformula,    # formula for emmeans
                    dt = NULL,   # time step for horizontal table
                    vars=NULL,   # list of unique variantID for each factor in specified sorting, 'unique(...)' will be used if missing
                    col=NULL,    # colours in the order corresponding to 'expand.grid(vars)', if only fingle vector supplied, it must correspond to first element of 'vars'
                    lgd=NULL,    # legend caption corresponding to 'expand.grid(vars)'
                    pdf=F,name=NULL,            # generate PDF output? Respective file name
                    CI=T,splitting=T,opening=T, # show CI? Show anther spliting curve? Show anther opening curve?
                    xl=0         # X0 position for legend in 'Rate' plot
                    ){
  require(emmeans)
  require(shades)
  yAO <- cbind(ath[,which(colnames(ath)==all.vars(formula[-3])[3])],ath[,which(colnames(ath)==all.vars(formula[-3])[2])]+ath[,which(colnames(ath)==all.vars(formula[-3])[1])])
  colnames(yAO) <- c("yAO1","yAO2")
  yAI <- cbind(ath[,which(colnames(ath)==all.vars(formula[-3])[2])]+ath[,which(colnames(ath)==all.vars(formula[-3])[3])],ath[,which(colnames(ath)==all.vars(formula[-3])[1])])
  colnames(yAI) <- c("yAI1","yAI2")
  ath <- cbind(ath,yAI,yAO)
  x <- unique(ath[,which(colnames(ath)==all.vars(formula[-2])[1])])
  if(is.null(vars)){
    vars <- list()
    for(l in 2:length(all.vars(formula[-2]))){
      vars[[l-1]] <- unique(ath[,which(colnames(ath)==all.vars(formula[-2])[l])])
      names(vars)[l-1] <- all.vars(formula[-2])[l]
    }
  }
  evars <- expand.grid(vars)
  moAI <- yAI[,1]/(yAI[,2]+yAI[,1])
  moAO <- yAO[,1]/(yAO[,2]+yAO[,1])
  # Formulae must be reformulated because the environment issue with emmeans
  chfm <- as.character(formula)
  iformula <- as.formula(paste(chfm[2],chfm[1],chfm[3]))
  mAI <- glm(update(iformula,cbind(yAI1,yAI2)~.), data=ath, family = "binomial")
  smAI <- step(mAI,trace=0)  # suppress print by trace
  print(coef(smAI))
  print(anova(smAI, test = 'Chisq')) 
  mAO <- glm(update(iformula,cbind(yAO1,yAO2)~.), data=ath, family = "binomial")
  smAO <- step(mAO,trace=0)
  print(coef(smAO))
  print(anova(smAO, test = 'Chisq'))
  if(pdf) pdf(paste(name,"_curves.pdf",sep=""),width=8,height = 7)
  if(is.null(lgd)) lgd <- apply(as.matrix(evars),1,paste,collapse=" ")
  par(las=1)
  xa <- c(x-min(x),max(x)+min(x))
  plot(1,type="n",xlim=range(xa),ylim=c(0,1),xaxt="n", yaxt="n",xlab="Time afret dew treatment interruption [min]",
       ylab=expression("Probability of developmental transition"),main="Anther development")
  #abline(v=x, col="gray")
  abline(h=0.5, col="gray")
  axis(side = 1, at = xa, labels = xa)
  axis(side = 2, at = c(0,0.5,1), labels = c("0.0","0.5","1.0"))
  if(is.null(col)|(length(col)!=length(vars))){
    cole <- evars
    for(c in 1:ncol(cole)){
      cole[,c] <- as.numeric(cole[,c])
      pcole <- 0.6+(cole[,c]-mean(cole[,c]))/(max(cole[,c])-min(cole[,c]))*0.8
      if(c==1) if(length(col) != length(vars[[1]])){
        col <- palette.colors(n=length(unique(cole[,c])),palette = "Polychrome 36")[cole[,c]]
      }else{  # expansion of basic cole to whole grid
        cole[,1] <- col
        col <- cole[,1]
      }
      if(c==2) col <- saturation(col,values = recycle(pcole))
      #if(c==3) col <- rgb(r=cole[,1], g=cole[,2], b=cole[,3])
    }
  }
  tp <- (0:(10*max(xa)))/10
  print("Prediction")
  for(id in 1:nrow(evars)){
    # predict fitted values
    fac <- data.frame(t=tp,matrix(as.character(unlist(evars[id,])),nrow = length(tp),ncol = ncol(evars),byrow = T,dimnames = list(NULL,colnames(evars))))
    pmAI <- predict(smAI,newdata = fac, type="response")
    pmAO <- predict(smAO,newdata = fac, type="response")
    lines(pmAI~tp,lwd=2,col=col[id],lty=2)
    lines(pmAO~tp,lwd=2,col=col[id])
    # predict link & SE
    emAI <- predict(smAI,newdata = fac, type="link",se.fit = T)
    emAO <- predict(smAO,newdata = fac, type="link",se.fit = T)
    lines(1/(1+exp(-(emAI$fit+1.96*emAI$se.fit)))~tp,lwd=0.5,col=col[id],lty=2)
    lines(1/(1+exp(-(emAI$fit-1.96*emAI$se.fit)))~tp,lwd=0.5,col=col[id],lty=2)
    lines(1/(1+exp(-(emAO$fit+1.96*emAI$se.fit)))~tp,lwd=0.5,col=col[id],lty=1)
    lines(1/(1+exp(-(emAO$fit-1.96*emAI$se.fit)))~tp,lwd=0.5,col=col[id],lty=1)
  }
  for(id in 1:nrow(evars)){
    if(ncol(evars)>1){
      sel <- which(apply(ath[,all.vars(formula[-2])[-1]],1,paste,collapse=" ")==paste(as.matrix(evars[id,]),collapse = " ") )
    }else{
      sel <- which(ath[,all.vars(formula[-2])[-1]]==paste(as.matrix(evars[id,]),collapse = " ") )
    }
    points(moAO[sel]~ath[sel,which(colnames(ath)==all.vars(formula[-2])[1])],col=col[id],pch=19,lwd=2)
    points(moAI[sel]~ath[sel,which(colnames(ath)==all.vars(formula[-2])[1])],col=col[id],pch=1,lwd=2)
  }
  legend(x=0,y=1,border=NA, seg.len=4, col = c(col,"grey50","grey50","grey50"),pch=c(rep(15,nrow(evars)),1,19,NA), pt.cex=c(rep(2,nrow(evars)),1,1,1), lwd=c(rep(NA,nrow(evars)),2,2,0.5), lty=c(rep(NA,nrow(evars)),2,1,1), legend = c(lgd,"Splitting","Opening","95% CI"))
  if(pdf) dev.off()
  print("EMMEANS")
  chfm <- as.character(emformula)
  emmformula <- as.formula(paste(chfm[1],chfm[2]))
  emAI <- emmeans(smAI, emmformula, type = "response")
  emAO <- emmeans(smAO, emmformula, type = "response")
  print(emAI)
  print(pairs(emAI))
  print(emAO)
  print(pairs(emAO))
  cemAI <- confint(emAI)
  cemAO <- confint(emAO)
  emvar <- apply(cemAO[all.vars(formula[-2])[-1]],1,paste,collapse=" ")
  if(pdf) pdf(paste(name,"_rate.pdf",sep=""),width=2*(3+length(vars)),height = 7)
  if(length(vars)>1) par(las=2,mar=c(10,4,4,2))
  plot(1,type="n",xlim=c(0,length(emvar)),ylim=c(0,1),xaxt="n",xlab="",
       ylab=expression("Midtime probability of developmental transition"),main="Anther opening")
  xp <- (1:length(emvar))-0.5
  xo <- match(apply(as.matrix(evars),1,paste,collapse=" "),emvar)
  axis(side = 1, at = xp, labels = emvar[xo])
  segments(xp-0.2,cemAI$asymp.LCL[xo],xp-0.2,cemAI$asymp.UCL[xo], col=col,lwd=3)
  segments(xp+0.2,cemAO$asymp.LCL[xo],xp+0.2,cemAO$asymp.UCL[xo], col=col,lwd=3)
  points(cemAI$prob[xo]~I(xp-0.2),pch=23, cex=1.5, col=col, bg= "white",lwd=3)
  points(cemAO$prob[xo]~I(xp+0.2),pch=23, cex=1.5, col=col, bg= col,lwd=3)
  legend(x=xl,y=1,horiz=F,border=NA, seg.len=4, fill = col, legend = lgd)
  if(pdf) dev.off()
}
```

