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
                na.strings="Undetermined") %>%
  filter(ID != 132) # Patient 132 is missing values for 6 genes
                    # Discarding as low quality

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

ci.auc.hd_hd = ci.auc(pROC::roc(y.hd, hd_hd))
ci.auc.gh_hd = ci.auc(pROC::roc(y.hd, gh_hd))
ci.auc.gh_gh = ci.auc(pROC::roc(y.gh, gh_gh))
ci.auc.hd_gh = ci.auc(pROC::roc(y.gh, hd_gh))

lb.aucs = sapply(list(ci.auc.hd_hd,
                      ci.auc.hd_gh,
                      ci.auc.gh_gh,
                      ci.auc.gh_hd), function(x) x[1])

ub.aucs = sapply(list(ci.auc.hd_hd,
                      ci.auc.hd_gh,
                      ci.auc.gh_gh,
                      ci.auc.gh_hd), function(x) x[3])

aucs = data.frame("Training Set"=c("HD", "HD", "GH", "GH"),
                  "Testing Set"=c("HD", "GH", "GH", "HD"),
                  AUC=c(auc.hd_hd, auc.hd_gh, auc.gh_gh, auc.gh_hd),
                  AUC.Lower=lb.aucs,
                  AUC.Upper=ub.aucs)

write.table(aucs, file="aucs.csv", sep=",", col.names=T,
            row.names=F)

# ROC curves

roc.hd_hd = performance(prediction(hd_hd, y.hd), "tpr", "fpr")
roc.gh_hd = performance(prediction(gh_hd, y.hd), "tpr", "fpr")
roc.gh_gh = performance(prediction(gh_gh, y.gh), "tpr", "fpr")
roc.hd_gh = performance(prediction(hd_gh, y.gh), "tpr", "fpr")