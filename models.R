library("dplyr")
library("tidyr")
library("glmnet")
library("pROC")

set.seed(123)

#### LOAD AND NORMALIZE DATA ####

# Load qPCR data
# hd: healthy donor
# gh: gene heart controls

hd = read.delim("./data/healthy_donor.txt",
                na.strings="Undetermined") 


gh = read.delim("./data/gene_heart.txt",
                na.strings="Undetermined") 


# Select qPCR data from data frame

x.hd = hd %>%
  select(ends_with("CT")) %>%
  select(-GAPDH.CT) %>%
  as.matrix

x.gh = gh %>%
  select(ends_with("CT")) %>%
  select(-GAPDH.CT) %>%
  as.matrix

# Replace missing with highest observed

x.hd[is.na(x.hd)] = max(x.hd, na.rm=T)
x.gh[is.na(x.gh)] = max(x.gh, na.rm=T)

# Normalize for GAPDH.CT

x.hd.normed = apply(x.hd, 2, function(z) z - hd$GAPDH.CT)
x.gh.normed = apply(x.gh, 2, function(z) z - gh$GAPDH.CT)

# Status 

y.hd = factor(hd$Sample.Type)
y.gh = factor(gh$Sample.Type)

#### TRAIN LINEAR MODELS ####

# Train linear models

m.hd = cv.glmnet(x.hd.normed,
                y.hd,
                family="binomial",
                nfolds=nrow(x.hd.normed),
                keep=T,
                grouped=F,
                alpha=0)

m.gh = cv.glmnet(x.gh.normed,
                 y.gh,
                 family="binomial",
                 nfolds=nrow(x.gh.normed),
                 keep=T,
                 grouped=F,
                 alpha=0)

# predict HD on GH

get.loocv = function(m) m$fit.preval[ ,m$lambda == m$lambda.1se]

p1 = get.loocv(m.hd)
p2 = get.loocv(m.gh)
p3 = predict(m.hd, newx=x.gh.normed)
p4 = predict(m.gh, newx=x.hd.normed)

r1 = roc(y.hd, p1)
r2 = roc(y.gh, p2)
r3 = roc(y.gh, p3)
r4 = roc(y.hd, p4)

col = brewer.pal(4, "Paired")[c(2,4,1,3)]

plot(r1, col=col[1])
plot(r2, add=T, col=col[2])
plot(r3, add=T, col=col[3])
plot(r4, add=T, col=col[4])



#### BOX PLOTS ####

make.bp = function(i){
  v = x[,i]
  
  stripchart(v ~ y + cohort,
             pch=20,
             frame=F,
             main=gsub(".CT", "", colnames(x)[i]),
             vertical=T,
             method="jitter",
             ylab="Normalized CT")
}

x = rbind(x.hd.normed, x.gh.normed)

y = c(as.character(y.hd), 
      as.character(y.gh))

cohort = c(rep("HD", nrow(x.hd)),
           rep("GH", nrow(x.gh)))

pdf("boxplots.pdf")
par(cex=1.3)
sapply(1:10, make.bp)
dev.off()


