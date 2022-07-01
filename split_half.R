library(haven)
library(lme4)
library(psych)
library(MLmetrics)
library(GPArotation)
library(mirt)
library(directlabels)
library(ggplot2)
library(ggeffects)
library(car)

### Notes for myself
# linear regression and pairwise correlations
# APA report of what you did
# rmarkdown
# simulation

# remember to set wd; why won't setwd work???
setwd("~/Documents/split_half")

# load data
study_1 <- read_sav("study1_separated.sav")
study_2 <- read.csv("study_2_dichotamized.csv")

# to delete rows
# study_1 <- study_1[-c()]

# save data to SPSS file
# save(study_1,file="study_1_edited.sav")


####### STUDY 1

# 5 number summary
summary(study_1)
names(study_1)


## split data into variables/scales

## STUDY 1 ##
# pre-test
pre_anxiety <- study_1[1]
pre_item_ME <- study_1[3:7]
pre_total_ME <- study_1[47]

# post-test
post_anxiety <- study_1[2]
post_item_ME <- study_1[18:22]
post_total_ME <- study_1[48]


# correlations
ME_corr_matrix <- cor(pre_item_ME,post_item_ME)
anxiety_corr <- cor(pre_anxiety,post_anxiety)
# r = 0.61 for anxiety pre and post

## STUDY 2 ##
# pre_test
pre_anxiety_2 <- study_2[21]
pre_item_attention <- study_2[1:10]
pre_item_attention_copy <- pre_item_attention

# post_test
post_anxiety_2 <- study_2[22]
post_item_attention <- study_2[11:20]
post_item_attention_copy <- post_item_attention

# combine
pre_item_attention_copy[is.na(pre_item_attention)] <- 0
pre_total_attention <- rowSums(pre_item_attention_copy)

post_item_attention_copy[is.na(post_item_attention)] <- 0
post_total_attention <- rowSums(post_item_attention_copy)

# correlations
attention_corr_matrix <- cor(pre_item_attention_copy,post_item_attention_copy)
anxiety_corr_2 <- cor(pre_anxiety_2,post_anxiety_2)
# r = 0.86 for anxiety pre and post


### BEGIN ANALYSES ###

## reliability

#### STUDY 1 ####
# ALPHA
# overall alpha
alpha(ME_corr_matrix,check.keys = TRUE)
# alpha = 0.86

# pre alpha
alpha(pre_item_ME)
# alpha = 0.4

# post alpha
alpha(post_item_ME)
# alpha = 0.45

# OMEGA
# overall omega
omega(ME_corr_matrix,nfactors = 1)
# omega = 0.87

# pre omega
omega(pre_item_ME,nfactors = 1)
# omega = 0.4

# post omega
omega(pre_item_ME,nfactors = 1)
# omega = 0.4


#### STUDY 2 ####
# ALPHA
# overall alpha
# alpha(attention_corr_matrix,check.keys = TRUE)
# alpha = 0.6

# pre alpha
# alpha(pre_item_attention_copy)
# alpha = 0.11
theta_se <- fscores(mirt_pre, full.scores.SE = TRUE) 
empirical_rxx(theta_se)
empirical_plot(pre_item_attention)
# 0.53
theta_se <- fscores(mirt_pre, full.scores.SE = TRUE, method = 'ML')
empirical_rxx(theta_se)

# post alpha
# alpha(post_item_attention_copy)
# alpha = 0.21

# OMEGA
# overall omega
omega(attention_corr_matrix,nfactors = 1)
# omega = 0.63

# pre omega
omega(pre_item_attention,nfactors = 1)
# omega = 0.68

# post omega
omega(post_item_attention,nfactors = 1)
# omega = 0.73

# IRT
# item difficulties/item response curves
# 1 factor solution
# oblimin rotation

#### STUDY 1 ####
# pre_ME
object <- irt.fa(pre_item_ME)
# items 1,2,3 are not highly discriminating; not good items
# centered around theta = 0; -1 to 1

# post_ME
post_fa_study_1 <- irt.fa(post_item_ME, cut = 0)
# item 4 is not highly disciminating; where is item 2?

#### STUDY 2 ####

# pre attention
# irt.fa(pre_item_attention_copy)
# items 2 and 8
# ICC
mirt_pre <- mirt(data = pre_item_attention, model = 1)
plot(mirt_pre,type = 'infotrace')
# trace lines
plot_pre <- plot(mirt_pre, type = 'trace')
direct.label(plot_pre, 'top.points')

# post attention
##
# irt.fa(post_item_attention_copy)
# only item 7 is highly discriminating
mirt_post <- mirt(data = post_item_attention, model = 1)
plot(mirt_post,type = 'infotrace')
# only 4 is bad item
# itemplot will plot individual items 

# descriptive
# study 1
describe(pre_item_ME,skew=FALSE)
describe(post_item_ME,skew=FALSE)

describe(pre_anxiety, skew = FALSE)
describe(post_anxiety, skew = FALSE)

# study 2
describe(pre_item_attention_copy,skew=FALSE)
describe(post_item_attention_copy,skew=FALSE)

describe(pre_anxiety_2, skew = FALSE)
describe(post_anxiety_2, skew = FALSE)


# Multiple Regression
# study 1
# pre anxiety, pre ME, and post ME predict post anxiety
lm_1 <- lm(unlist(post_anxiety) ~ unlist(pre_anxiety) + unlist(pre_total_ME) + unlist(post_total_ME))
plot(lm_1)

# study 2
# pre anxiety, pre ME, and post ME predict post anxiety
lm_2 <- lm(unlist(post_anxiety_2) ~ unlist(pre_anxiety_2) + unlist(pre_total_attention) + unlist(post_total_attention))
plot(lm_2)
scatterplot(unlist(post_anxiety_2)~ (unlist(pre_anxiety_2)))