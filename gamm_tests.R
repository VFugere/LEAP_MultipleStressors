## Code to analyze LEAP 2016 MS data
## Vincent Fugère 2019-2020
## This code tries a number of gamms to find the optimal model structure
## also trying a number of linear models

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
library(mgcv)
library(itsadug)

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

#### finding the optimal gamm ####

#all smooths + interactions
ba.model <- gam(log10(BA) ~ o.nut + s(date,k=4) + s(sc.gly, k = 4) + s(sc.imi, k = 4) + ti(date,sc.gly, k=4) + ti(date,sc.gly, by = o.nut, k=4) + ti(date,sc.imi, k=4) + ti(date,sc.imi, by = o.nut, k=4) + ti(date,sc.gly,sc.imi, k=4) + ti(date,sc.gly,sc.imi, by = o.nut, k=4) + s(date, site.f, bs='fs',k=3, m=2), data=merged.data, method = 'REML')
fit <- fitted(ba.model)
summary(ba.model)
gam.check(ba.model)

#no 3 way interaction
ba.model2 <- gam(log10(BA) ~ o.nut + s(date,k=4) + s(sc.gly, k = 4) + s(sc.imi, k = 4) + ti(date,sc.gly, k=4) + ti(date,sc.gly, by = o.nut, k=4) + ti(date,sc.imi, k=4) + ti(date,sc.imi, by = o.nut, k=4) + ti(date,sc.gly,sc.imi, k=4) + s(date, site.f, bs='fs',k=3, m=2), data=merged.data, method = 'REML')
fit2 <- fitted(ba.model2)
summary(ba.model2)
gam.check(ba.model2)
plot(fit2~fit)

#not eactly sure what the te() in this one does...
ba.model3 <- gam(log10(BA) ~ nut.fac + te(date,sc.gly,sc.imi, by = nut.fac, k=4) + s(date, site.f, bs='fs',k=3, m=2), data=merged.data, method = 'REML')
summary(ba.model3)
gam.check(ba.model3)
fit3 <- fitted(ba.model3)
plot(fit3~fit)

#no smooth for overall effect of pesticides, only interactions with time. Not sure this is appropriate
ba.model4 <- gam(log10(BA) ~ o.nut + s(date,k=4) + ti(date,sc.gly, k=4) + ti(date,sc.gly, by = o.nut, k=4) + ti(date,sc.imi, k=4) + ti(date,sc.imi, by = o.nut, k=4) + ti(date,sc.gly,sc.imi, k=4) + ti(date,sc.gly,sc.imi, by = o.nut, k=4) + s(date, site.f, bs='fs',k=3, m=2), data=merged.data, method = 'REML')
fit4 <- fitted(ba.model4)
summary(ba.model4)
gam.check(ba.model4)
plot(fit4~fit)

#te approach instead
ba.model5 <- gam(log10(BA) ~ nut.fac + te(date,sc.gly, by = nut.fac, k=4) + te(date,sc.imi, by = nut.fac, k=4) + ti(date,sc.gly,sc.imi, by = nut.fac, k=4) + s(date, site.f, bs='fs',k=3, m=2), data=merged.data, method = 'REML')
fit5 <- fitted(ba.model5)
summary(ba.model5)
gam.check(ba.model5)
plot(fit5~fit)

#changing s to ti
ba.model6 <- gam(log10(BA) ~ o.nut + ti(date,k=4) + ti(sc.gly, k = 4) + ti(sc.imi, k = 4) + ti(date,sc.gly, k=4) + ti(date,sc.gly, by = o.nut, k=4) + ti(date,sc.imi, k=4) + ti(date,sc.imi, by = o.nut, k=4) + ti(date,sc.gly,sc.imi, k=4) + ti(date,sc.gly,sc.imi, by = o.nut, k=4) + s(date, site.f, bs='fs',k=3, m=2), data=merged.data, method = 'REML')
summary(ba.model6)
gam.check(ba.model6)
plot(predict(ba.model6)~predict(ba.model))

AIC(ba.model,ba.model2,ba.model3,ba.model4,ba.model5,ba.model6)