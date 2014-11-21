library("ggplot2")
library("dplyr")
library("tidyr")

source("models.R")

coef.heatmap = function(m, e, y, top=T){
  
  co = coef(m)[-1, ] %>%
    as.matrix 
  co = co[order(rownames(co), decreasing=T), , drop=F]
  e = e[ , rownames(co)]
  
  if(top) par(mar=c(3.1,6.1,4.1,1)) else  par(mar=c(5.1,6.1,2.1,1)) 
  barplot(-co[,1],
          horiz=T, 
          col="black", 
          xlab="Coefficient",
          names.arg=gsub(".CT", "", rownames(co)),
          las=1,
          xlim=c(-0.6, 0.6))
  if(top) par(mar=c(3.1,1,4.1,3))  else par(mar=c(5.1,1,2.1,3)) 
  image(-scale(e), axes=F, col=redgreen(75),
        breaks=seq(-4, 4, length.out=76),
        main=if(top) "Health Donor Expression" else "Gene Heart Expression")
  
}

tiff("coefficient_heatmaps.tiff", res=300, width=8, height=8,
     units="in")

par(cex=1.3)
layout(matrix(c(1,2,2,3,4,4), nrow=2, byrow=T))

coef.heatmap(m.hd, x.hd.normed)
coef.heatmap(m.gh, x.gh.normed, top=F)

dev.off()

tiff("rocs.tiff", res=300, width=8, height=8, units="in")
par(cex=1.3)
par(mfrow=c(1,1))
par(mar=c(4.5,4.5,4.5,4.5))

plot(roc.hd_hd@x.values[[1]],
     roc.hd_hd@y.values[[1]],
     bty="n",
     col="black",
     type="l",
     lwd=2,
     xlab="False Positive Rate",
     ylab="True Positive Rate")

points(roc.hd_gh@x.values[[1]],
       roc.hd_gh@y.values[[1]],
       bty="n",
       col="black",
       lwd=2,
       type="l",
       lty=2)

points(roc.gh_gh@x.values[[1]],
       roc.gh_gh@y.values[[1]],
       bty="n",
       col="red",
       lwd=2,
       type="l",
       lty=1)

points(roc.gh_hd@x.values[[1]],
       roc.gh_hd@y.values[[1]],
       bty="n",
       col="red",
       lwd=2,
       type="l",
       lty=2)

lines(c(0, 1), c(0,1), col="grey")

dev.off()
