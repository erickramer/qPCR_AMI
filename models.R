library("dplyr")
library("tidyr")
library("glmnet")
library("pROC")
library("ROCR")
library("RColorBrewer")

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

#### PREDICTIONS ####

# format for these variables
# <set used for training>_<set used for prediction>
# e.g. hd_gh is the the model trained on hd,
#      and predicted on gh

# the elastic net automatically does 
# LOOCV cross-validation
get.loocv = function(m) m$fit.preval[ ,m$lambda == m$lambda.1se]
hd_hd = get.loocv(m.hd)
gh_gh = get.loocv(m.gh)
hd_gh = predict(m.hd, newx=x.gh.normed)
gh_hd = predict(m.gh, newx=x.hd.normed)

# AUCs

# we have to use the weird notation because there
# are a lot of functions called "auc" and "roc"
# gotta make sure we're calling the right one

auc.hd_hd = pROC::auc(pROC::roc(y.hd, hd_hd))
auc.gh_hd = pROC::auc(pROC::roc(y.hd, gh_hd))
auc.gh_gh = pROC::auc(pROC::roc(y.gh, gh_gh))
auc.hd_gh = pROC::auc(pROC::roc(y.gh, hd_gh))

# ROC curves

roc.hd_hd = performance(prediction(hd_hd, y.hd), "tpr", "fpr")
roc.gh_hd = performance(prediction(gh_hd, y.hd), "tpr", "fpr")
roc.gh_gh = performance(prediction(gh_gh, y.gh), "tpr", "fpr")
roc.hd_gh = performance(prediction(hd_gh, y.gh), "tpr", "fpr")

col = brewer.pal(4, "Paired")[c(2,4,1,3)]

pdf("4_roc_curves.pdf")

plot(roc.hd_hd, lwd=2, col=col[1], main="qPCR ROC")
plot(roc.gh_gh, lwd=2, col=col[2], add=T)
plot(roc.gh_hd, lwd=2, col=col[3], add=T)
plot(roc.hd_gh, lwd=2, col=col[4], add=T)

legend(0.2, 0.3, 
       c("LOOCV with Health Donors", 
         "LOOCV with Gene Heart",
         "Gene Heart predictions with Healthy Donor Model", 
         "Heathy Donor predictions with Gene Heart Model"),
       lwd=2,
       col=col)

dev.off()