## Code to analyze LEAP 2016 MS data
## Vincent Fugère 2019-2020
## This code fits linear models and exports a caterpillar plot

rm(list=ls())

library(tidyverse)
library(magrittr)
library(scales)
library(shape)
library(RColorBrewer)
library(viridis)
library(plotrix)

library(lme4)
library(performance)
library(sjPlot)
library(sjlabelled)
library(sjmisc)

devtools::source_url('https://raw.githubusercontent.com/VFugere/Rfuncs/master/vif.R')
devtools::source_url('https://raw.githubusercontent.com/VFugere/Rfuncs/master/utils.R')

#### load data ####

load('~/Google Drive/Recherche/LEAP Postdoc/2016/MSdata.RData')

#### colour palette & plotting parameters ####

pchs <- c(1,0)

gly.cols <- c('gray90',brewer.pal(9, 'Reds')[2:8])
imi.cols <- c('gray90',brewer.pal(9, 'Blues')[2:8])
both.cols <- c('gray90',brewer.pal(9, 'Greens')[2:8])

glycolfunc <- colorRampPalette(gly.cols)
imicolfunc <- colorRampPalette(imi.cols)
bothcolfunc <- colorRampPalette(both.cols)

allcols <- c(gly.cols,gly.cols,imi.cols,imi.cols,both.cols,both.cols)

##### format treatment variables for models ####

#adding ordered factor for GAMs
merged.data$o.nut <- as.ordered(merged.data$nut.fac)
#rescaling pesticides gradients from 0 to 1 to compare effect with nutrient factor
merged.data$sc.gly <- rescale(merged.data$gly, c(0,1))
merged.data$sc.imi <- rescale(merged.data$imi, c(0,1))
#adding site factor
merged.data$site.f <- as.factor(merged.data$site)
#adding date factor
merged.data$date.f <- as.factor(merged.data$date)
#reordering
merged.data <- select(merged.data, date:pond.id,o.nut:date.f,everything())

#adding RUE
merged.data$RUE <- with(merged.data, total_zoo/total)
merged.data$RUE[55:56] <- 0

##### Linear models to quantify effect sizes#####

merged.data$nut.num <- merged.data$nut.num-1

merged.data$sc.gly.std <- scale(merged.data$sc.gly)
merged.data$sc.imi.std <- scale(merged.data$sc.imi)

merged.data$nut.num.std2 <- arm::rescale(merged.data$nut.num, 'full')
merged.data$sc.gly.std2 <- arm::rescale(merged.data$sc.gly)
merged.data$sc.imi.std2 <- arm::rescale(merged.data$sc.imi)

merged.data$date.f.sc <- as.factor(rescale(merged.data$date, c(0,1)))

#all 3 stressors ranging from 0 to 1
m1 <- lmer(log10(total) ~ date.f:nut.num + date.f:sc.gly + date.f:sc.imi + date.f:sc.gly:sc.imi + date.f:sc.gly:nut.num +
                      date.f:sc.imi:nut.num + (1|site.f),merged.data)
#standardizing by 1 sd: wrong! (see Gelman 2008)
m2 <- lmer(log10(total) ~ date.f:nut.num + date.f:sc.gly.std + date.f:sc.imi.std + date.f:sc.gly.std:sc.imi.std + date.f:sc.gly.std:nut.num +
             date.f:sc.imi.std:nut.num + (1|site.f),merged.data)
#standarding by 2 sds instead, so that all range between -0.5 and 0.5
m3 <- lmer(log10(total) ~ date.f:nut.num.std2 + date.f:sc.gly.std2 + date.f:sc.imi.std2 + date.f:sc.gly.std2:sc.imi.std2 + date.f:sc.gly.std2:nut.num.std2 +
             date.f:sc.imi.std2:nut.num.std2 + (1|site.f),merged.data)
#gly and imi from 1:8 - definitely wrong!
m4 <- lmer(log10(total) ~ date.f:nut.num + date.f:gly + date.f:imi + date.f:gly:imi + date.f:gly:nut.num +
             date.f:imi:nut.num + (1|site.f),merged.data)

performance(m1)
performance(m3)

plot_model(m1)
plot_model(m3)

#does not change results very much to use std2, but R2 is way worse

#adding main effects to see if makes a diff
m1.m <- lmer(log10(total) ~ date.f + nut.num + sc.gly + sc.imi + date.f:nut.num + date.f:sc.gly + date.f:sc.imi + date.f:sc.gly:sc.imi + date.f:sc.gly:nut.num +
             date.f:sc.imi:nut.num + (1|site.f),merged.data)
#do I need to add a scaled version of date.f?
m1.m.sc <- lmer(log10(total) ~ date.f.sc + nut.num + sc.gly + sc.imi + date.f.sc:nut.num + date.f.sc:sc.gly + date.f.sc:sc.imi + date.f.sc:sc.gly:sc.imi + date.f.sc:sc.gly:nut.num +
               date.f.sc:sc.imi:nut.num + (1|site.f),merged.data)
#no definitely not. Does not change anything.

#standarding by 2 sds instead, so that all range between -0.5 and 0.5
m3.m <- lmer(log10(total) ~ date.f + nut.num.std2 + sc.gly.std2 + sc.imi.std2 + date.f:nut.num.std2 + date.f:sc.gly.std2 + date.f:sc.imi.std2 + date.f:sc.gly.std2:sc.imi.std2 + date.f:sc.gly.std2:nut.num.std2 +
             date.f:sc.imi.std2:nut.num.std2 + (1|site.f),merged.data)

performance(m1.m)
plot_model(m1.m)

performance(m1.m.sc)
plot_model(m1.m.sc)

performance(m3.m)
performance(m3)
plot_model(m3.m)
plot_model(m1)

#removing time:nut interaction
m4.m <- lmer(log10(total) ~ date.f + nut.num.std2 + sc.gly.std2 + sc.imi.std2 + date.f:sc.gly.std2 + date.f:sc.imi.std2 + date.f:sc.gly.std2:sc.imi.std2 + date.f:sc.gly.std2:nut.num.std2 +
               date.f:sc.imi.std2:nut.num.std2 + (1|site.f),merged.data)
anova(m4.m)
plot_model(m4.m, type='diag')
performance(m4.m)
plot_model(m4.m)

m5 <- lmer(log10(total) ~ date.f*nut.num*sc.gly*sc.imi + (1|site.f),merged.data)
anova(m5)
plot_model(m5, type='diag')
performance(m5)
plot_model(m5)

#how would things look without temporal patterns and interactions
sub <- filter(merged.data, date == 7)
sub.m <- lm(log10(total) ~ nut.num+sc.gly+sc.imi+nut.num:sc.gly+nut.num:sc.imi+sc.gly:sc.imi,sub)
sub.m2 <- lm(log10(total) ~ nut.num.std2+sc.gly.std2+sc.imi.std2+nut.num.std2:sc.gly.std2+nut.num.std2:sc.imi.std2+sc.gly.std2:sc.imi.std2,sub)
plot_model(sub.m)
plot_model(sub.m2)
car::vif(sub.m)
car::vif(sub.m2)
summary(sub.m2)
#sub.m2 is a perfectly adequate regression model for this single time point. How do the coefs compare with the global model of choice?
coef(sub.m) #+0.54 nut, -0.53 gly
coef(sub.m2) #+0.44 (0.09) nut , -0.48 (0.09) gly
compare_performance(sub.m,sub.m2) #the fit is absolutely identical.
fixef(m1) #+0.23 nut, -0.98 gly
fixef(m3) #+0.44 nut, -0.48 gly. Highly similar, without main effects in the model
fixef(m1.m) #if I add main + interaction, I get +0.54 & -0.53, the same coefs as in sub.m
fixef(m3.m) #if I add main + interaction, I get +0.44 & -0.48, the same coefs as in sub.m2
compare_performance(m1,m3,m1.m,m3.m) #so m3.m is the only one that recovers the right coefficient + has a good R2
plot_model(m3.m)
get_model_data(m3, type = 'est') #0.44(0.18) & -0.47(0.18) #standard errors are twice as large as lm
library(brms)
m3.b <- brm(log10(total) ~ date.f:nut.num.std2 + date.f:sc.gly.std2 + date.f:sc.imi.std2 + date.f:sc.gly.std2:sc.imi.std2 + date.f:sc.gly.std2:nut.num.std2 +
             date.f:sc.imi.std2:nut.num.std2 + (1|site.f),merged.data)
summary(m3.b) #0.44(0.18) & -0.47(0.18) : exactly the same coefs and standard errors

#manually checking unexpected imi effect at last time points
sub <- filter(merged.data, date == 43)
sub.m <- lm(log10(total) ~ nut.num+sc.gly+sc.imi+nut.num:sc.gly+nut.num:sc.imi+sc.gly:sc.imi,sub)
sub.m2 <- lm(log10(total) ~ nut.num.std2+sc.gly.std2+sc.imi.std2+nut.num.std2:sc.gly.std2+nut.num.std2:sc.imi.std2+sc.gly.std2:sc.imi.std2,sub)
plot_model(sub.m)
plot_model(sub.m2)
#no imi effect...
# a list of significant effect I should recover:
sub <- filter(merged.data, date == 30)
sub.m2 <- lm(log10(total) ~ nut.num.std2+sc.gly.std2+sc.imi.std2+nut.num.std2:sc.gly.std2+nut.num.std2:sc.imi.std2+sc.gly.std2:sc.imi.std2,sub)
summary(sub.m2)
#day 1: nut+
#day 7: nut+, gly-
#day 15: nut+
#day 30: nut+ m.s, gly+
#day 35: nut+ 
#day 43: nut+ m.s, gly+
get_model_data(m1, type='est')
get_model_data(m3, type='est')
get_model_data(m3.m, type='est')

#another model, with main effect of date but not others
m5.m <- lmer(log10(total) ~ date.f + date.f:nut.num.std2 + date.f:sc.gly.std2 + date.f:sc.imi.std2 + date.f:sc.gly.std2:sc.imi.std2 + date.f:sc.gly.std2:nut.num.std2 +
               date.f:sc.imi.std2:nut.num.std2 + (1|site.f),merged.data)
performance(m5.m) #good
plot_model(m5.m) #good
#day 7 linear model with scaled variables: #+0.44 (0.09) nut, -0.48 (0.09) gly
get_model_data(m5.m, type='est') #+0.44(0.13) nut, -0.48 (0.14) gly. Perfect but SEs somewhat larger
plot_model(m5.m, type='diag')

m6.m <- lm(log10(total) ~ date.f + date.f:nut.num.std2 + date.f:sc.gly.std2 + date.f:sc.imi.std2 + date.f:sc.gly.std2:sc.imi.std2 + date.f:sc.gly.std2:nut.num.std2 +
               date.f:sc.imi.std2:nut.num.std2,merged.data)
summary(m6.m) #removing random effect does not change fixed effects

#ok, lets go with m5.m then.

##### Fitting and validating all m5.m-type models #####

ba.mod.lin <- lmer(log10(BA) ~ date.f + date.f:nut.num.std2 + date.f:sc.gly.std2 + date.f:sc.imi.std2 + date.f:sc.gly.std2:sc.imi.std2 + date.f:sc.gly.std2:nut.num.std2 +
                     date.f:sc.imi.std2:nut.num.std2 + (1|site.f),merged.data)
plot_model(ba.mod.lin,type='est')
plot_model(ba.mod.lin,type='diag')
performance(ba.mod.lin)
plot(fitted(ba.mod.lin)~log10(merged.data$BA))

chla.mod.lin <- lmer(log10(total) ~ date.f + date.f:nut.num.std2 + date.f:sc.gly.std2 + date.f:sc.imi.std2 + date.f:sc.gly.std2:sc.imi.std2 + date.f:sc.gly.std2:nut.num.std2 +
                       date.f:sc.imi.std2:nut.num.std2 + (1|site.f),merged.data)
plot_model(chla.mod.lin,type='est')
plot_model(chla.mod.lin,type='diag')
performance(chla.mod.lin)
plot(fitted(chla.mod.lin)~log10(merged.data$total))

zoo.mod.lin <- lmer(log10p(total_zoo) ~ date.f + date.f:nut.num.std2 + date.f:sc.gly.std2 + date.f:sc.imi.std2 + date.f:sc.gly.std2:sc.imi.std2 + date.f:sc.gly.std2:nut.num.std2 +
                      date.f:sc.imi.std2:nut.num.std2 + (1|site.f),merged.data)
plot_model(zoo.mod.lin,type='est')
plot_model(zoo.mod.lin,type='diag')
performance(zoo.mod.lin)
plot(fitted(zoo.mod.lin)~log10p(merged.data$total_zoo_adult))

EP.mod.lin <- lmer(use ~ date.f + date.f:nut.num.std2 + date.f:sc.gly.std2 + date.f:sc.imi.std2 + date.f:sc.gly.std2:sc.imi.std2 + date.f:sc.gly.std2:nut.num.std2 +
                     date.f:sc.imi.std2:nut.num.std2 + (1|site.f),merged.data)
plot_model(EP.mod.lin,type='est')
plot_model(EP.mod.lin,type='diag')
performance(EP.mod.lin)
plot(fitted(EP.mod.lin)~merged.data$use)

nep.mod.lin <- lmer(log10p(NEP) ~ date.f + date.f:nut.num.std2 + date.f:sc.gly.std2 + date.f:sc.imi.std2 + date.f:sc.gly.std2:sc.imi.std2 + date.f:sc.gly.std2:nut.num.std2 +
                      date.f:sc.imi.std2:nut.num.std2 + (1|site.f),merged.data)
summary(nep.mod.lin)
plot_model(nep.mod.lin,type='est')
plot_model(nep.mod.lin,type='diag')
performance(nep.mod.lin) #no random effect variance...
plot(fitted(nep.mod.lin)~log10p(merged.data$NEP))
#can I trust these parameter estimates? Trying other packages
nep.mod.lin <- glmmTMB::glmmTMB(log10p(NEP) ~ date.f + date.f:nut.num.std2 + date.f:sc.gly.std2 + date.f:sc.imi.std2 + date.f:sc.gly.std2:sc.imi.std2 + date.f:sc.gly.std2:nut.num.std2 +
                      date.f:sc.imi.std2:nut.num.std2 + (1|site.f),merged.data)
nep.mod.lin.b <- brm(log10p(NEP) ~ date.f + date.f:nut.num.std2 + date.f:sc.gly.std2 + date.f:sc.imi.std2 + date.f:sc.gly.std2:sc.imi.std2 + date.f:sc.gly.std2:nut.num.std2 +
                                  date.f:sc.imi.std2:nut.num.std2 + (1|site.f),merged.data)
fixef(nep.mod.lin)
fixef(nep.mod.lin.b)
#almost identical parameter estimates so I'll trust the frequentist LMM even if sigma[pond] == 0

rue.mod.lin <- lmer(log10(RUE) ~ date.f + date.f:nut.num.std2 + date.f:sc.gly.std2 + date.f:sc.imi.std2 + date.f:sc.gly.std2:sc.imi.std2 + date.f:sc.gly.std2:nut.num.std2 +
                      date.f:sc.imi.std2:nut.num.std2 + (1|site.f),merged.data,subset=is.finite(log(RUE)))
summary(rue.mod.lin)
plot_model(rue.mod.lin,type='est')
plot_model(rue.mod.lin,type='diag')
performance(rue.mod.lin)
     