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

#adding RUE
merged.data$RUE <- with(merged.data, total_zoo/total)
merged.data$RUE[55:56] <- 0

#### finding the optimal gamm ####

#' chla.model <- gam(log10(total) ~ o.nut + ti(date,k=4) + ti(sc.gly, k = 4) + ti(sc.imi, k = 4) + ti(date,sc.gly, k=6) + ti(date,sc.gly, by = o.nut, k=3) + ti(date,sc.imi, k=6) + ti(date,sc.imi, by = o.nut, k=3) + ti(date,sc.gly,sc.imi, k=6) + ti(date,sc.gly,sc.imi, by = o.nut, k=3) + s(site.f, bs='re'), data=merged.data, method = 'REML')
#' gam.check(chla.model)
#' summary(chla.model)
#' plot(chla.model)
#' 
#' ba.model <- gam(log10(BA) ~ o.nut + ti(date,k=6) + ti(sc.gly, k = 3) + ti(sc.imi, k = 3) + ti(date,sc.gly, k=6) + ti(date,sc.gly, by = o.nut, k=3) + ti(date,sc.imi, k=6) + ti(date,sc.imi, by = o.nut, k=3) + ti(date,sc.gly,sc.imi, k=6) + ti(date,sc.gly,sc.imi, by = o.nut, k=3) + s(site.f, bs='re'), data=merged.data, method = 'REML')
#' gam.check(ba.model)
#' summary(ba.model)
#' plot(ba.model)
#' 
#' m1 <- gam(log10(BA) ~ o.nut + ti(date,k=6) + ti(sc.gly, k = 3) + ti(sc.imi, k = 3) + ti(date,sc.gly, k=6) + ti(date,sc.gly, by = o.nut, k=3) + ti(date,sc.imi, k=6) + ti(date,sc.imi, by = o.nut, k=3) + ti(date,sc.gly,sc.imi, k=6) + ti(date,sc.gly,sc.imi, by = o.nut, k=3) + s(site.f, bs='re'), data=merged.data, method = 'REML')
#' m2 <- gam(log10(BA) ~ o.nut + ti(date,k=4) + ti(sc.gly, k = 3) + ti(sc.imi, k = 3) + ti(date,sc.gly, k=4) + ti(date,sc.gly, by = o.nut, k=3) + ti(date,sc.imi, k=4) + ti(date,sc.imi, by = o.nut, k=3) + ti(date,sc.gly,sc.imi, k=4) + ti(date,sc.gly,sc.imi, by = o.nut, k=3) + s(date, site.f, bs='fs',k=4), data=merged.data, method = 'REML')
#' AIC(m1,m2)
#' plot(m2)
#' 
#' #' model with random smooths is doing a lot better (R2 and AIC). However, the basis dimension
#' #' is so low that random smooths are in fact straight lines. Moreover, the need to reduce k
#' #' considerably for all smooths leads to a poorer representation of pesticide effects.
#' #' If I include random smooths, need to find a way to increase k by reducing number of effects
#' 
#' m3 <- gam(log10(BA) ~ o.nut + ti(date,k=5) + ti(date,sc.gly, k=6) + ti(date,sc.gly, by = o.nut, k=3) + ti(date,sc.imi, k=6) + ti(date,sc.imi, by = o.nut, k=3) + ti(date,sc.gly,sc.imi, k=5) + ti(date,sc.gly,sc.imi, by = o.nut, k=3) + s(date, site.f, bs='fs',k=3), data=merged.data, method = 'REML')
#' summary(m3)
#' gam.check(m3)
#' plot(m3)
#' 
#' AIC(m3,m2)
#' 
#' m4 <- gam(log10(total) ~ o.nut + ti(date,k=5) + ti(date,sc.gly, k=6) + ti(date,sc.gly, by = o.nut, k=3) + ti(date,sc.imi, k=6) + ti(date,sc.imi, by = o.nut, k=3) + ti(date,sc.gly,sc.imi, k=5) + ti(date,sc.gly,sc.imi, by = o.nut, k=3) + s(date, site.f, bs='fs',k=3), data=merged.data, method = 'REML')
#' summary(m4)
#' gam.check(m4)
#' plot(m4)
#' 
#' AIC(m4,chla.model) # model with random smooths is definitey better. All gam.check() calls suggest basis dimensions are ok

ba.m <- gam(log10(BA) ~ o.nut + ti(date,k=5) + ti(date,sc.gly, k=6) + ti(date,sc.gly, by = o.nut, k=3) + ti(date,sc.imi, k=6) + ti(date,sc.imi, by = o.nut, k=3) + ti(date,sc.gly,sc.imi, k=5) + ti(date,sc.gly,sc.imi, by = o.nut, k=3) + s(date, site.f, bs='fs',k=3), data=merged.data, method = 'REML')
chla.m <- gam(log10(total) ~ o.nut + ti(date,k=5) + ti(date,sc.gly, k=6) + ti(date,sc.gly, by = o.nut, k=3) + ti(date,sc.imi, k=6) + ti(date,sc.imi, by = o.nut, k=3) + ti(date,sc.gly,sc.imi, k=5) + ti(date,sc.gly,sc.imi, by = o.nut, k=3) + s(date, site.f, bs='fs',k=3), data=merged.data, method = 'REML')
zoo.m <- gam(log10p(total_zoo) ~ o.nut + ti(date,k=5) + ti(date,sc.gly, k=6) + ti(date,sc.gly, by = o.nut, k=3) + ti(date,sc.imi, k=6) + ti(date,sc.imi, by = o.nut, k=3) + ti(date,sc.gly,sc.imi, k=5) + ti(date,sc.gly,sc.imi, by = o.nut, k=3) + s(date, site.f, bs='fs',k=3), data=merged.data, method = 'REML')

## finding a good distribution for 'use' variable
# 
# use.model <- gam(use ~ o.nut + ti(date,k=5) + ti(date,sc.gly, k=6) + ti(date,sc.gly, by = o.nut, k=3) + ti(date,sc.imi, k=6) + ti(date,sc.imi, by = o.nut, k=3) + ti(date,sc.gly,sc.imi, k=5) + ti(date,sc.gly,sc.imi, by = o.nut, k=3) + s(date, site.f, bs='fs',k=3), data=merged.data, method = 'REML')
# poisson.mod <- gam(use ~ o.nut + ti(date,k=5) + ti(date,sc.gly, k=6) + ti(date,sc.gly, by = o.nut, k=3) + ti(date,sc.imi, k=6) + ti(date,sc.imi, by = o.nut, k=3) + ti(date,sc.gly,sc.imi, k=5) + ti(date,sc.gly,sc.imi, by = o.nut, k=3) + s(date, site.f, bs='fs',k=3), data=merged.data, method = 'REML',family=poisson)
# negbin.mod <- gam(use ~ o.nut + ti(date,k=5) + ti(date,sc.gly, k=6) + ti(date,sc.gly, by = o.nut, k=3) + ti(date,sc.imi, k=6) + ti(date,sc.imi, by = o.nut, k=3) + ti(date,sc.gly,sc.imi, k=5) + ti(date,sc.gly,sc.imi, by = o.nut, k=3) + s(date, site.f, bs='fs',k=3), data=merged.data, method = 'REML',family=nb)
# compare_performance(use.model,poisson.mod,negbin.mod)
# gam.check(use.model)
# gam.check(poisson.mod) #worse
# gam.check(negbin.mod)
# #better to stick to gaussian model after all

use.m <- gam(use ~ o.nut + ti(date,k=5) + ti(date,sc.gly, k=6) + ti(date,sc.gly, by = o.nut, k=3) + ti(date,sc.imi, k=6) + ti(date,sc.imi, by = o.nut, k=3) + ti(date,sc.gly,sc.imi, k=5) + ti(date,sc.gly,sc.imi, by = o.nut, k=3) + s(date, site.f, bs='fs',k=3), data=merged.data, method = 'REML')
nep.m <- gam(log10p(NEP) ~ o.nut + ti(date,k=5) + ti(date,sc.gly, k=6) + ti(date,sc.gly, by = o.nut, k=3) + ti(date,sc.imi, k=6) + ti(date,sc.imi, by = o.nut, k=3) + ti(date,sc.gly,sc.imi, k=5) + ti(date,sc.gly,sc.imi, by = o.nut, k=3) + s(date, site.f, bs='fs',k=3), data=merged.data, method = 'REML')
rue.m <- gam(log10(RUE) ~ o.nut + ti(date,k=5) + ti(date,sc.gly, k=6) + ti(date,sc.gly, by = o.nut, k=3) + ti(date,sc.imi, k=6) + ti(date,sc.imi, by = o.nut, k=3) + ti(date,sc.gly,sc.imi, k=5) + ti(date,sc.gly,sc.imi, by = o.nut, k=3) + s(date, site.f, bs='fs',k=3), data=merged.data, subset=is.finite(log(RUE)), method = 'REML')
diatoms.m <- gam(log10p(diatoms) ~ o.nut + ti(date,k=5) + ti(date,sc.gly, k=6) + ti(date,sc.gly, by = o.nut, k=3) + ti(date,sc.imi, k=6) + ti(date,sc.imi, by = o.nut, k=3) + ti(date,sc.gly,sc.imi, k=5) + ti(date,sc.gly,sc.imi, by = o.nut, k=3) + s(date, site.f, bs='fs',k=3), data=merged.data, method = 'REML')
prop.diatoms.m <- gam((diatoms/total) ~ o.nut + ti(date,k=5) + ti(date,sc.gly, k=6) + ti(date,sc.gly, by = o.nut, k=3) + ti(date,sc.imi, k=6) + ti(date,sc.imi, by = o.nut, k=3) + ti(date,sc.gly,sc.imi, k=5) + ti(date,sc.gly,sc.imi, by = o.nut, k=3) + s(date, site.f, bs='fs',k=3), data=merged.data, method = 'REML', family=betar)
 
save(ba.m,chla.m,zoo.m,use.m,nep.m,rue.m,diatoms.m,prop.diatoms.m, file='~/Google Drive/Recherche/LEAP Postdoc/2016/GAMMs.RData')

##### tweaking k's to better capture data, as shows on scattergams ####

ba.m <- gam(log10(BA) ~ o.nut + ti(date,k=5) + ti(date,sc.gly, k=6) + ti(date,sc.gly, by = o.nut, k=3) + ti(date,sc.imi, k=6) + ti(date,sc.imi, by = o.nut, k=3) + ti(date,sc.gly,sc.imi, k=5) + ti(date,sc.gly,sc.imi, by = o.nut, k=3) + s(date, site.f, bs='fs',k=3), data=merged.data, method = 'REML')
ba.m2 <- gam(log10(BA) ~ o.nut + ti(date,k=4) + ti(date,sc.gly, k=6) + ti(date,sc.gly, by = o.nut, k=3) + ti(date,sc.imi, k=6) + ti(date,sc.imi, by = o.nut, k=4) + ti(date,sc.gly,sc.imi, k=5) + ti(date,sc.gly,sc.imi, by = o.nut, k=3) + s(date, site.f, bs='fs',k=3), data=merged.data, method = 'REML')
ba.m3 <- gam(log10(BA) ~ o.nut + ti(date,k=3) + ti(date,sc.gly, k=6) + ti(date,sc.gly, by = o.nut, k=3) + ti(date,sc.imi, k=6) + ti(date,sc.imi, by = o.nut, k=5) + ti(date,sc.gly,sc.imi, k=4) + ti(date,sc.gly,sc.imi, by = o.nut, k=3) + s(date, site.f, bs='fs',k=3), data=merged.data, method = 'REML')
ba.m4 <- gam(log10(BA) ~ o.nut + ti(date,k=3) + ti(date,sc.gly, k=6) + ti(date,sc.gly, by = o.nut, k=3) + ti(date,sc.imi, k=6) + ti(date,sc.imi, by = o.nut, k=6) + ti(date,sc.gly,sc.imi, k=4) + ti(date,sc.gly,sc.imi, by = o.nut, k=4) + s(date, site.f, bs='fs',k=3), data=merged.data, method = 'REML')
ba.m5 <- gam(log10(BA) ~ o.nut + ti(date,k=5) + ti(date,sc.gly, k=5) + ti(date,sc.gly, by = o.nut, k=3) + ti(date,sc.imi, k=6) + ti(date,sc.imi, by = o.nut, k=4) + ti(date,sc.gly,sc.imi, k=5) + ti(date,sc.gly,sc.imi, by = o.nut, k=3) + s(date, site.f, bs='fs',k=3), data=merged.data, method = 'REML')
ba.m6 <- gam(log10(BA) ~ o.nut + ti(date,k=6) + ti(date,sc.gly, k=6) + ti(date,sc.gly, by = o.nut, k=4) + ti(date,sc.imi, k=6) + ti(date,sc.imi, by = o.nut, k=6) + ti(date,sc.gly,sc.imi, k=4) + ti(date,sc.gly,sc.imi, by = o.nut, k=3) + s(date, site.f, bs='fs',k=3), data=merged.data, method = 'REML')

compare_performance(ba.m,ba.m2,ba.m3,ba.m4,ba.m5,ba.m6)
scattergam(var=vars[1],model.name=ba.m, varname=var.names[1])
scattergam(var=vars[1],model.name=ba.m2, varname=var.names[1])
scattergam(var=vars[1],model.name=ba.m3, varname=var.names[1])
scattergam(var=vars[1],model.name=ba.m4, varname=var.names[1])
scattergam(var=vars[1],model.name=ba.m5, varname=var.names[1])
scattergam(var=vars[1],model.name=ba.m6, varname=var.names[1])

ba.m6 -> ba.m

chla.m <- gam(log10(total) ~ o.nut + ti(date,k=5) + ti(date,sc.gly, k=6) + ti(date,sc.gly, by = o.nut, k=3) + ti(date,sc.imi, k=6) + ti(date,sc.imi, by = o.nut, k=3) + ti(date,sc.gly,sc.imi, k=5) + ti(date,sc.gly,sc.imi, by = o.nut, k=3) + s(date, site.f, bs='fs',k=3), data=merged.data, method = 'REML')
chla.m2 <- gam(log10(total) ~ o.nut + ti(date,k=6) + ti(date,sc.gly, k=6) + ti(date,sc.gly, by = o.nut, k=6) + ti(date,sc.imi, k=4) + ti(date,sc.imi, by = o.nut, k=3) + ti(date,sc.gly,sc.imi, k=5) + ti(date,sc.gly,sc.imi, by = o.nut, k=3) + s(date, site.f, bs='fs',k=3), data=merged.data, method = 'REML')
chla.m3 <- gam(log10(total) ~ o.nut + ti(date,k=6) + ti(date,sc.gly, k=6) + ti(date,sc.gly, by = o.nut, k=6) + ti(date,sc.imi, k=6) + ti(date,sc.imi, by = o.nut, k=3) + ti(date,sc.gly,sc.imi, k=4) + ti(date,sc.gly,sc.imi, by = o.nut, k=4) + s(date, site.f, bs='fs',k=3), data=merged.data, method = 'REML')
compare_performance(chla.m,chla.m2,chla.m3)
scattergam(var=vars[2],model.name=chla.m, varname=var.names[2])
scattergam(var=vars[2],model.name=chla.m2, varname=var.names[2])
scattergam(var=vars[2],model.name=chla.m3, varname=var.names[2])

chla.m3 -> chla.m

zoo.m <- gam(log10p(total_zoo) ~ o.nut + ti(date,k=5) + ti(date,sc.gly, k=6) + ti(date,sc.gly, by = o.nut, k=3) + ti(date,sc.imi, k=6) + ti(date,sc.imi, by = o.nut, k=3) + ti(date,sc.gly,sc.imi, k=5) + ti(date,sc.gly,sc.imi, by = o.nut, k=3) + s(date, site.f, bs='fs',k=3), data=merged.data, method = 'REML')
zoo.m2 <- gam(log10p(total_zoo) ~ o.nut + ti(date,k=6) + ti(date,sc.gly, k=6) + ti(date,sc.gly, by = o.nut, k=6) + ti(date,sc.imi, k=6) + ti(date,sc.imi, by = o.nut, k=3) + ti(date,sc.gly,sc.imi, k=4) + ti(date,sc.gly,sc.imi, by = o.nut, k=3) + s(date, site.f, bs='fs',k=3), data=merged.data, method = 'REML')
zoo.m3 <- gam(log10p(total_zoo) ~ o.nut + ti(date,k=6) + ti(date,sc.gly, k=6) + ti(date,sc.gly, by = o.nut, k=4) + ti(date,sc.imi, k=6) + ti(date,sc.imi, by = o.nut, k=3) + ti(date,sc.gly,sc.imi, k=5) + ti(date,sc.gly,sc.imi, by = o.nut, k=3) + s(date, site.f, bs='fs',k=3), data=merged.data, method = 'REML')
zoo.m4 <- gam(log10p(total_zoo) ~ o.nut + ti(date,k=3) + ti(date,sc.gly, k=3) + ti(date,sc.imi, k=3) + ti(date,sc.imi, by = o.nut, k=3) + ti(date,sc.gly,sc.imi, k=6) + s(date, site.f, bs='fs',k=3), data=merged.data, method = 'REML')
compare_performance(zoo.m,zoo.m2,zoo.m3,zoo.m4)
scattergam(var=vars[3],model.name=zoo.m, varname=var.names[3])
scattergam(var=vars[3],model.name=zoo.m2, varname=var.names[3])
scattergam(var=vars[3],model.name=zoo.m3, varname=var.names[3])
scattergam(var=vars[3],model.name=zoo.m4, varname=var.names[3])

zoo.m3 -> zoo.m

use.m <- gam(use ~ o.nut + ti(date,k=5) + ti(date,sc.gly, k=6) + ti(date,sc.gly, by = o.nut, k=3) + ti(date,sc.imi, k=6) + ti(date,sc.imi, by = o.nut, k=3) + ti(date,sc.gly,sc.imi, k=5) + ti(date,sc.gly,sc.imi, by = o.nut, k=3) + s(date, site.f, bs='fs',k=3), data=merged.data, method = 'REML')
use.m2 <- gam(use ~ o.nut + ti(date,k=6) + ti(date,sc.gly, k=6) + ti(date,sc.gly, by = o.nut, k=4) + ti(date,sc.imi, k=6) + ti(date,sc.imi, by = o.nut, k=3) + ti(date,sc.gly,sc.imi, k=5) + ti(date,sc.gly,sc.imi, by = o.nut, k=3) + s(date, site.f, bs='fs',k=3), data=merged.data, method = 'REML')
use.m3 <- gam(use ~ o.nut + ti(date,k=6) + ti(date,sc.gly, k=5) + ti(date,sc.gly, by = o.nut, k=3) + ti(date,sc.imi, k=5) + ti(date,sc.imi, by = o.nut, k=4) + ti(date,sc.gly,sc.imi, k=5) + ti(date,sc.gly,sc.imi, by = o.nut, k=4) + s(date, site.f, bs='fs',k=3), data=merged.data, method = 'REML')
compare_performance(use.m,use.m2,use.m3)
scattergam(var=vars[4],model.name=use.m, varname=var.names[4])
scattergam(var=vars[4],model.name=use.m2, varname=var.names[4])
scattergam(var=vars[4],model.name=use.m3, varname=var.names[4])

use.m3 -> use.m

nep.m <- gam(log10p(NEP) ~ o.nut + ti(date,k=5) + ti(date,sc.gly, k=6) + ti(date,sc.gly, by = o.nut, k=3) + ti(date,sc.imi, k=6) + ti(date,sc.imi, by = o.nut, k=3) + ti(date,sc.gly,sc.imi, k=5) + ti(date,sc.gly,sc.imi, by = o.nut, k=3) + s(date, site.f, bs='fs',k=3), data=merged.data, method = 'REML')
nep.m2 <- gam(log10p(NEP) ~ o.nut + ti(date,k=6) + ti(date,sc.gly, k=6) + ti(date,sc.gly, by = o.nut, k=4) + ti(date,sc.imi, k=6) + ti(date,sc.imi, by = o.nut, k=3) + ti(date,sc.gly,sc.imi, k=5) + ti(date,sc.gly,sc.imi, by = o.nut, k=3) + s(date, site.f, bs='fs',k=3), data=merged.data, method = 'REML')
nep.m3 <- gam(log10p(NEP) ~ o.nut + ti(date,k=6) + ti(date,sc.gly, k=6) + ti(date,sc.gly, by = o.nut, k=6) + ti(date,sc.imi, k=3) + ti(date,sc.imi, by = o.nut, k=3) + ti(date,sc.gly,sc.imi, k=3) + ti(date,sc.gly,sc.imi, by = o.nut, k=5) + s(date, site.f, bs='fs',k=3), data=merged.data, method = 'REML')
compare_performance(nep.m,nep.m2,nep.m3)
scattergam(var=vars[5],model.name=nep.m, varname=var.names[5])
scattergam(var=vars[5],model.name=nep.m2, varname=var.names[5])
scattergam(var=vars[5],model.name=nep.m3, varname=var.names[5])

nep.m3 -> nep.m

rue.m <- gam(log10(RUE) ~ o.nut + ti(date,k=5) + ti(date,sc.gly, k=6) + ti(date,sc.gly, by = o.nut, k=3) + ti(date,sc.imi, k=6) + ti(date,sc.imi, by = o.nut, k=3) + ti(date,sc.gly,sc.imi, k=5) + ti(date,sc.gly,sc.imi, by = o.nut, k=3) + s(date, site.f, bs='fs',k=3), data=merged.data, subset=is.finite(log(RUE)), method = 'REML')
rue.m2 <- gam(log10(RUE) ~ o.nut + ti(date,k=6) + ti(date,sc.gly, k=6) + ti(date,sc.gly, by = o.nut, k=4) + ti(date,sc.imi, k=6) + ti(date,sc.imi, by = o.nut, k=3) + ti(date,sc.gly,sc.imi, k=5) + ti(date,sc.gly,sc.imi, by = o.nut, k=3) + s(date, site.f, bs='fs',k=3), data=merged.data, subset=is.finite(log(RUE)), method = 'REML')
rue.m3 <- gam(log10(RUE) ~ o.nut + ti(date,k=6) + ti(date,sc.gly, k=6) + ti(date,sc.gly, by = o.nut, k=6) + ti(date,sc.imi, k=6) + ti(date,sc.imi, by = o.nut, k=5) + ti(date,sc.gly,sc.imi, k=6) + ti(date,sc.gly,sc.imi, by = o.nut, k=3) + s(date, site.f, bs='re',k=2), data=merged.data, subset=is.finite(log(RUE)), method = 'REML')
compare_performance(rue.m,rue.m2,rue.m3)
scattergam(var=vars[6],model.name=rue.m, varname=var.names[6])
scattergam(var=vars[6],model.name=rue.m2, varname=var.names[6])
scattergam(var=vars[6],model.name=rue.m3, varname=var.names[6])

rue.m2 -> rue.m

save(ba.m,chla.m,zoo.m,use.m,nep.m,rue.m,file='~/Google Drive/Recherche/LEAP Postdoc/2016/GAMMs.RData')

#### scattergams 3 rows ####

vars <- c('BA','total','total_zoo','use','NEP','RUE')

var.names <- c(expression(log[10]~bact.~(cells/mu*L)),
               expression(log[10]~chl.~italic(a)~(mu*g/L)),
               expression(log[10](1+zoo)~(mu*g/L)),
               expression(subtrates~used~by~bact.),
               expression(log[10]~(1+Delta*DO)~(mu*g/L)),
               expression(log[10]~(Zoo:Phyto)~(mu*g/mu*g)))

pchs.2 <- c(21,22)

plot.data <- select(merged.data, date:date.f, vars)
plot.data <- plot.data %>% mutate_at(vars(BA,total,RUE), log10)
plot.data <- plot.data %>% mutate_at(vars(total_zoo,NEP), log10p)

scattergam <- function(var, varname, model.name){
tmp <- plot.data[,c(colnames(plot.data)[1:13],var)]
tmp <- drop_na(tmp)
tmp <- filter(tmp, is.finite(tmp[,14]))
ylims <- range(tmp[,14])
for(letter in c('C|D','E|H','J|K')){
  for(date.x in Sampling.dates){
    sub <- tmp %>% filter(date == date.x, str_detect(site, letter))
    sub$pesticide <- rescale(as.numeric(str_remove(sub$site, letter)),c(0,1))
    if(date.x == 1 & letter == 'J|K'){
      plot(y=sub[,14],x=sub[,15],xlab=NULL,ylab=NULL,bty='o',type='p',ylim=ylims,pch=pchs.2[sub$nut.num],col=1,bg=allcols[sub$pond.id],xaxt='n')
      axis(1, lwd=0, lwd.ticks = 1, at = seq(0,1,length.out = 8), labels=1:8)
    }else if(date.x == 1 & letter != 'J|K'){
      plot(y=sub[,14],x=sub[,15],xlab=NULL,ylab=NULL,bty='o',type='p',ylim=ylims,pch=pchs.2[sub$nut.num],col=1,bg=allcols[sub$pond.id],xaxt='n')
    }else if(date.x != 1 & letter == 'J|K'){
      plot(y=sub[,14],x=sub[,15],xlab=NULL,ylab=NULL,bty='o',type='p',ylim=ylims,pch=pchs.2[sub$nut.num],col=1,bg=allcols[sub$pond.id],yaxt='n',xaxt='n')
      axis(1, lwd=0, lwd.ticks = 1, at = seq(0,1,length.out = 8), labels=1:8)
    }else{
      plot(y=sub[,14],x=sub[,15],xlab=NULL,ylab=NULL,bty='o',type='p',ylim=ylims,pch=pchs.2[sub$nut.num],col=1,bg=allcols[sub$pond.id],yaxt='n',xaxt='n')
    }
    if(letter == 'C|D'){
      plot_smooth(model.name, view="sc.gly", cond=list('date'=date.x,'sc.imi'=0,'o.nut'='low'), rm.ranef=TRUE, rug=F,col=alpha(gly.cols[8],0.7),lty=1,add=T,print.summary = F)
      plot_smooth(model.name, view="sc.gly", cond=list('date'=date.x,'sc.imi'=0,'o.nut'='high'), rm.ranef=TRUE, rug=F,col=alpha(gly.cols[8],0.7),lty=2,lwd=1.5,add=T,print.summary = F)
    }else if(letter == 'E|H'){
      plot_smooth(model.name, view="sc.imi", cond=list('date'=date.x,'sc.gly'=0,'o.nut'='low'), rm.ranef=TRUE, rug=F,col=alpha(imi.cols[8],0.7),lty=1,add=T,print.summary = F)
      plot_smooth(model.name, view="sc.imi", cond=list('date'=date.x,'sc.gly'=0,'o.nut'='high'), rm.ranef=TRUE, rug=F,col=alpha(imi.cols[8],0.7),lty=2,lwd=1.5,add=T,print.summary = F)
    }else{
      fitted <- as.data.frame(predict.gam(model.name, newdata = list('date' = rep(date.x,200), 'sc.gly' = rep(seq(0,1,length.out=100),2), 'sc.imi' = rep(seq(0,1,length.out=100),2), 'o.nut' = c(rep('low',100),rep('high',100)), 'site.f' = rep('C1',200)), exclude = s(date,site.f), se.fit = T))
      fitted$lwr <- fitted$fit - 1.96*fitted$se.fit
      fitted$upr <- fitted$fit + 1.96*fitted$se.fit
      poly(x=seq(0,1,length.out=100),upper=fitted$upr[1:100],lower=fitted$lwr[1:100],fill=alpha(both.cols[8],0.1))
      points(fitted$fit[1:100]~seq(0,1,length.out=100),type='l',col=alpha(both.cols[8],0.7),lty=1)
      poly(x=seq(0,1,length.out=100),upper=fitted$upr[101:200],lower=fitted$lwr[101:200],fill=alpha(both.cols[8],0.1))
      points(fitted$fit[101:200]~seq(0,1,length.out=100),type='l',col=alpha(both.cols[8],0.7),lty=2,lwd=1.5)
    }
  }
}
mtext(varname,side=2,outer=T,line=2.5,cex=1.2)
mtext('pesticide nominal concentration (dose 1 to 8)',side=1,outer=T,line=2.5,cex=1.2)
mtext(paste('day',Sampling.dates,' '),side=3,outer=T,line=0.1,at=seq(0.1,0.93,length.out = 6),adj=0.5)
}

pdf('~/Desktop/FigS5-10_scattergams.pdf',width=5.5,height = 3,pointsize = 8,onefile = T)
par(mfrow=c(3,6),mar=c(0.1,0.1,0.1,0.1),oma=c(4,4,2,0.5),cex=1,xpd=T)
scattergam(var=vars[1],model.name=ba.m, varname=var.names[1])
scattergam(var=vars[2],model.name=chla.m, varname=var.names[2])
scattergam(var=vars[3],model.name=zoo.m, varname=var.names[3])
scattergam(var=vars[4],model.name=use.m, varname=var.names[4])
scattergam(var=vars[5],model.name=nep.m, varname=var.names[5])
scattergam(var=vars[6],model.name=rue.m, varname=var.names[6])
dev.off()


###### ugly visualization of gamms ####

VF.gam.plot <- function(model.name, varname){
  #fitted <- as.data.frame(predict(model.name))
  ylims <- range(fitted(model.name))
  vals <- seq(0,1,length.out=30) #how many lines to draw?
  gly.cols.plot <- glycolfunc(length(vals))
  imi.cols.plot <- imicolfunc(length(vals))
  both.cols.plot <- bothcolfunc(length(vals))
  plot_smooth(model.name, view="date", cond=list('sc.gly'=0,'sc.imi'=0,'o.nut'='low'), rm.ranef=TRUE, rug=F,col=gly.cols.plot[1],print.summary = F,hide.label = T,cex.axis=1,ann=F,bty='l',legend_plot_all = F, h0=NA,ylim = ylims,se=0)
  for(i in 2:length(vals)){plot_smooth(model.name, view="date", cond=list('sc.gly'=vals[i],'sc.imi'=0,'o.nut'='low'), rm.ranef=TRUE, rug=F,col=gly.cols.plot[i],add=T,print.summary = F,se=0)}
  legend('topleft',legend='low',bty='n')
  mtext(varname,2,line=2.5)
  plot_smooth(model.name, view="date", cond=list('sc.gly'=0,'sc.imi'=0,'o.nut'='high'), rm.ranef=TRUE, rug=F,col=gly.cols.plot[1],print.summary = F,hide.label = T,cex.axis=1,ann=F,bty='l',legend_plot_all = F, h0=NA,ylim = ylims,se=0)
  for(i in 2:length(vals)){plot_smooth(model.name, view="date", cond=list('sc.gly'=vals[i],'sc.imi'=0,'o.nut'='high'), rm.ranef=TRUE, rug=F,col=gly.cols.plot[i],add=T,print.summary = F,se=0)}
  legend('topleft',legend='high',bty='n')
  plot_smooth(model.name, view="date", cond=list('sc.gly'=0,'sc.imi'=0,'o.nut'='low'), rm.ranef=TRUE, rug=F,col=imi.cols.plot[1],print.summary = F,hide.label = T,cex.axis=1,ann=F,bty='l',legend_plot_all = F, h0=NA,ylim = ylims,se=0)
  for(i in 2:length(vals)){plot_smooth(model.name, view="date", cond=list('sc.gly'=0,'sc.imi'=vals[i],'o.nut'='low'), rm.ranef=TRUE, rug=F,col=imi.cols.plot[i],add=T,print.summary = F,se=0)}
  legend('topleft',legend='low',bty='n')
  plot_smooth(model.name, view="date", cond=list('sc.gly'=0,'sc.imi'=0,'o.nut'='high'), rm.ranef=TRUE, rug=F,col=imi.cols.plot[1],print.summary = F,hide.label = T,cex.axis=1,ann=F,bty='l',legend_plot_all = F, h0=NA,ylim = ylims,se=0)
  for(i in 2:length(vals)){plot_smooth(model.name, view="date", cond=list('sc.gly'=0,'sc.imi'=vals[i],'o.nut'='high'), rm.ranef=TRUE, rug=F,col=imi.cols.plot[i],add=T,print.summary = F,se=0)}
  legend('topleft',legend='high',bty='n')
  plot_smooth(model.name, view="date", cond=list('sc.gly'=0,'sc.imi'=0,'o.nut'='low'), rm.ranef=TRUE, rug=F,col=both.cols.plot[1],print.summary = F,hide.label = T,cex.axis=1,ann=F,bty='l',legend_plot_all = F, h0=NA,ylim = ylims,se=0)
  for(i in 2:length(vals)){plot_smooth(model.name, view="date", cond=list('sc.gly'=vals[i],'sc.imi'=vals[i],'o.nut'='low'), rm.ranef=TRUE, rug=F,col=both.cols.plot[i],add=T,print.summary = F,se=0)}
  legend('topleft',legend='low',bty='n')
  plot_smooth(model.name, view="date", cond=list('sc.gly'=0,'sc.imi'=0,'o.nut'='high'), rm.ranef=TRUE, rug=F,col=both.cols.plot[1],print.summary = F,hide.label = T,cex.axis=1,ann=F,bty='l',legend_plot_all = F, h0=NA,ylim = ylims,se=0)
  for(i in 2:length(vals)){plot_smooth(model.name, view="date", cond=list('sc.gly'=vals[i],'sc.imi'=vals[i],'o.nut'='high'), rm.ranef=TRUE, rug=F,col=both.cols.plot[i],add=T,print.summary = F,se=0)}
  legend('topleft',legend='high',bty='n')
}

pdf('~/Desktop/gamms.pdf',width=15,height=10,pointsize = 12)
par(mfrow=c(3,6),mar=c(2,2,2,2),oma=c(2,2,0,0),cex=1)
VF.gam.plot(model.name=ba.m, varname=expression(log[10]~bacterial~abundance~(cells/mu*L)))
VF.gam.plot(model.name=chla.m, varname=expression(log[10]~chlorophyll~italic(a)~(mu*g/L)))
VF.gam.plot(model.name=zoo.m, varname=expression(log[10](1+zooplankton~biomass)~(mu*g/L)))
mtext('date',1,outer=T,line=0.5)
dev.off()

pdf('~/Desktop/gamms_ef.pdf',width=15,height=10,pointsize = 12)
par(mfrow=c(3,6),mar=c(2,2,2,2),oma=c(2,2,0,0),cex=1)
VF.gam.plot(model.name=use.m, varname=expression(carbon~use~potential))
VF.gam.plot(model.name=nep.m, varname=expression(log[10]~daytime~Delta*DO~(mu*g/L)))
VF.gam.plot(model.name=rue.m, varname=expression(log[10]~RUE~(mu*g/mu*g)))
mtext('date',1,outer=T,line=0.5)
dev.off()

#### only max vs. min: see main code for final version of this ####

VF.gam.plot.simpler <- function(model.name, varname){
  fitted <- as.data.frame(predict(model.name, se.fit = T,exclude='s(date,site.f)'))
  #fitted <- as.data.frame(predict(model.name, se.fit = T,exclude='s(site.f)'))
  fitted$lwr <- fitted$fit - 1.96*fitted$se.fit
  fitted$upr <- fitted$fit + 1.96*fitted$se.fit
  ylims <- range(c(fitted$lwr,fitted$upr))
  plot_smooth(model.name, view="date", cond=list('sc.gly'=0,'sc.imi'=0,'o.nut'='low'), rm.ranef=TRUE, rug=F,col='gray50',print.summary = F,hide.label = T,cex.axis=1,ann=F,bty='l',legend_plot_all = F, h0=NA,ylim = ylims,lwd=2)
  plot_smooth(model.name, view="date", cond=list('sc.gly'=0.5,'sc.imi'=0,'o.nut'='low'), rm.ranef=TRUE, rug=F,col=gly.cols[5],add=T,print.summary = F,lwd=2)
  plot_smooth(model.name, view="date", cond=list('sc.gly'=1,'sc.imi'=0,'o.nut'='low'), rm.ranef=TRUE, rug=F,col=gly.cols[8],add=T,print.summary = F,lwd=2)
  mtext(varname,2,line=2.5)
  plot_smooth(model.name, view="date", cond=list('sc.gly'=0,'sc.imi'=0,'o.nut'='high'), rm.ranef=TRUE, rug=F,col='gray50',print.summary = F,hide.label = T,cex.axis=1,ann=F,bty='l',legend_plot_all = F, h0=NA,ylim = ylims,lwd=2)
  plot_smooth(model.name, view="date", cond=list('sc.gly'=0.5,'sc.imi'=0,'o.nut'='high'), rm.ranef=TRUE, rug=F,col=gly.cols[5],add=T,print.summary = F,lwd=2)
  plot_smooth(model.name, view="date", cond=list('sc.gly'=1,'sc.imi'=0,'o.nut'='high'), rm.ranef=TRUE, rug=F,col=gly.cols[8],add=T,print.summary = F,lwd=2)
  plot_smooth(model.name, view="date", cond=list('sc.gly'=0,'sc.imi'=0,'o.nut'='low'), rm.ranef=TRUE, rug=F,col='gray50',print.summary = F,hide.label = T,cex.axis=1,ann=F,bty='l',legend_plot_all = F, h0=NA,ylim = ylims,lwd=2)
  plot_smooth(model.name, view="date", cond=list('sc.gly'=0,'sc.imi'=0.5,'o.nut'='low'), rm.ranef=TRUE, rug=F,col=imi.cols[5],add=T,print.summary = F,lwd=2)
  plot_smooth(model.name, view="date", cond=list('sc.gly'=0,'sc.imi'=1,'o.nut'='low'), rm.ranef=TRUE, rug=F,col=imi.cols[8],add=T,print.summary = F,lwd=2)
  plot_smooth(model.name, view="date", cond=list('sc.gly'=0,'sc.imi'=0,'o.nut'='high'), rm.ranef=TRUE, rug=F,col='gray50',print.summary = F,hide.label = T,cex.axis=1,ann=F,bty='l',legend_plot_all = F, h0=NA,ylim = ylims,lwd=2)
  plot_smooth(model.name, view="date", cond=list('sc.gly'=0,'sc.imi'=0.5,'o.nut'='high'), rm.ranef=TRUE, rug=F,col=imi.cols[5],add=T,print.summary = F,lwd=2)
  plot_smooth(model.name, view="date", cond=list('sc.gly'=0,'sc.imi'=1,'o.nut'='high'), rm.ranef=TRUE, rug=F,col=imi.cols[8],add=T,print.summary = F,lwd=2)
  plot_smooth(model.name, view="date", cond=list('sc.gly'=0,'sc.imi'=0,'o.nut'='low'), rm.ranef=TRUE, rug=F,col='gray50',print.summary = F,hide.label = T,cex.axis=1,ann=F,bty='l',legend_plot_all = F, h0=NA,ylim = ylims,lwd=2)
  plot_smooth(model.name, view="date", cond=list('sc.gly'=0.5,'sc.imi'=0.5,'o.nut'='low'), rm.ranef=TRUE, rug=F,col=both.cols[5],add=T,print.summary = F,lwd=2)
  plot_smooth(model.name, view="date", cond=list('sc.gly'=1,'sc.imi'=1,'o.nut'='low'), rm.ranef=TRUE, rug=F,col=both.cols[8],add=T,print.summary = F,lwd=2)
  plot_smooth(model.name, view="date", cond=list('sc.gly'=0,'sc.imi'=0,'o.nut'='high'), rm.ranef=TRUE, rug=F,col='gray50',print.summary = F,hide.label = T,cex.axis=1,ann=F,bty='l',legend_plot_all = F, h0=NA,ylim = ylims,lwd=2)
  plot_smooth(model.name, view="date", cond=list('sc.gly'=0.5,'sc.imi'=0.5,'o.nut'='high'), rm.ranef=TRUE, rug=F,col=both.cols[5],add=T,print.summary = F,lwd=2)
  plot_smooth(model.name, view="date", cond=list('sc.gly'=1,'sc.imi'=1,'o.nut'='high'), rm.ranef=TRUE, rug=F,col=both.cols[8],add=T,print.summary = F,lwd=2)
}

pdf('~/Desktop/gamms_simpler.pdf',width=15,height=10,pointsize = 12)
par(mfrow=c(3,6),mar=c(2,2,2,2),oma=c(2,2,0,0),cex=1)
VF.gam.plot.simpler(model.name=ba.m, varname=expression(log[10]~bacterial~abundance~(cells/mu*L)))
VF.gam.plot.simpler(model.name=chla.m, varname=expression(log[10]~chlorophyll~italic(a)~(mu*g/L)))
VF.gam.plot.simpler(model.name=zoo.m, varname=expression(log[10](1+zooplankton~biomass)~(mu*g/L)))
mtext('date',1,outer=T,line=0.5)
dev.off()

pdf('~/Desktop/gamms_simpler_ef.pdf',width=15,height=10,pointsize = 12)
par(mfrow=c(3,6),mar=c(2,2,2,2),oma=c(2,2,0,0),cex=1)
VF.gam.plot.simpler(model.name=use.m, varname=expression(C~sources~used))
VF.gam.plot.simpler(model.name=nep.m, varname=expression(log[10]~daytime~Delta*DO~(mu*g/L)))
VF.gam.plot.simpler(model.name=rue.m, varname=expression(log[10]~RUE~(mu*g/mu*g)))
mtext('date',1,outer=T,line=0.5)
dev.off()

#### contour plots ####

# heatcolfunc <- colorRampPalette(brewer.pal(11, 'RdYlBu')[11:1])
# heat.cols <- heatcolfunc(50)
#heat.cols <- wes_palette('Zissou1', 50, 'continuous')
heat.cols <- viridis(50)

VF.cont.plot <- function(model.name,varname){
  zlims <- range( fitted <- as.data.frame(predict(model.name, se.fit = F,exclude='s(date,site.f)')))
  for(nut.lvl in c('low','high')){
    for(date.x in Sampling.dates){
      if(date.x == 43 & nut.lvl == 'low'){
        fvisgam(model.name, view = c('sc.gly','sc.imi'), too.far=0.3,cond = list('date' = date.x, 'o.nut' = nut.lvl), zlim=zlims, add.color.legend=T,dec=1,hide.label=T,plot.type = 'contour', lwd=1.5,color = heat.cols, main = NULL,rm.ranef = T,print.summary = F,yaxt='n',xaxt='n',xlab=NULL,ylab=NULL)
        }else{
        fvisgam(model.name, view = c('sc.gly','sc.imi'), too.far=0.3,cond = list('date' = date.x, 'o.nut' = nut.lvl), zlim=zlims, add.color.legend=F,hide.label=T,plot.type = 'contour', lwd=1.5,color = heat.cols, main = NULL,rm.ranef = T,print.summary = F,yaxt='n',xaxt='n',xlab=NULL,ylab=NULL)
      }
    }
  }
  mtext('glyphosate',1,outer=T,line=0.5)
  mtext('imidacloprid',2,outer=T,line=0.5)
  mtext(paste('day',Sampling.dates,' '),side=3,outer=T,line=0.5,at=seq(0.09,0.92,length.out = 6),adj=0.5)
  mtext(c('low nut.','high nut.'),side=4,outer=T,line=0.5,at=c(0.75,0.25))
  mtext(varname,side=3,outer=T,line=2,adj=0.5,cex=1.2)
}

pdf('~/Desktop/contourplots.pdf',width=7.5,height=3.5,pointsize = 10,onefile = T)
par(mfrow=c(2,6),mar=c(0.1,0.1,0.1,0.1),oma=c(2,2,3.5,2),cex=1,xpd=T)
VF.cont.plot(model.name=ba.m, varname=expression(log[10]~bacterial~abundance~(cells/mu*L)))
VF.cont.plot(model.name=chla.m, varname=expression(log[10]~chlorophyll~italic(a)~(mu*g/L)))
VF.cont.plot(model.name=zoo.m, varname=expression(log[10](1+zooplankton~biomass)~(mu*g/L)))
VF.cont.plot(model.name=use.m, varname=expression(carbon~use~potential))
VF.cont.plot(model.name=nep.m, varname=expression(log[10]~daytime~Delta*DO~(mu*g/L)))
VF.cont.plot(model.name=rue.m, varname=expression(log[10]~zooplankton:phytoplankton~(mu*g/mu*g)))
dev.off()

#### other gamms:ef, diatoms ####

par(mfrow=c(2,6),mar=c(0.1,0.1,0.1,0.1),oma=c(2,2,2,2),cex=1,xpd=T)
VF.cont.plot(use.m)
summary(use.m)
VF.cont.plot(nep.m)
summary(nep.m)
VF.cont.plot(diatoms.m)
summary(diatoms.m)
VF.cont.plot(prop.diatoms.m)
summary(prop.diatoms.m)

tmp <- merged.data[,c(colnames(merged.data)[1:13])]
tmp$prop.diatoms <- with(merged.data, diatoms/total)
tmp <- drop_na(tmp)
ylims <- range(tmp[,14])
par(mfrow=c(3,6),mar=c(0.1,0.1,0.1,0.1),oma=c(4,4,2,0.5),cex=1,xpd=T)
pchs.2 <- c(16,15)

for(letter in c('C|D','E|H','J|K')){
  for(date.x in Sampling.dates){
    sub <- tmp %>% filter(date == date.x, str_detect(site, letter))
    sub$pesticide <- as.numeric(str_remove(sub$site, letter))
    if(date.x == 1 & letter == 'J|K'){
      plot(y=sub[,14],x=sub[,15],xlab=NULL,ylab=NULL,bty='o',type='p',ylim=ylims,pch=pchs.2[sub$nut.num],col=allcols[sub$pond.id])
    }else if(date.x == 1 & letter != 'J|K'){
      plot(y=sub[,14],x=sub[,15],xlab=NULL,ylab=NULL,bty='o',type='p',ylim=ylims,pch=pchs.2[sub$nut.num],col=allcols[sub$pond.id],xaxt='n')
    }else if(date.x != 1 & letter == 'J|K'){
      plot(y=sub[,14],x=sub[,15],xlab=NULL,ylab=NULL,bty='o',type='p',ylim=ylims,pch=pchs.2[sub$nut.num],col=allcols[sub$pond.id],yaxt='n')
    }else{
      plot(y=sub[,14],x=sub[,15],xlab=NULL,ylab=NULL,bty='o',type='p',ylim=ylims,pch=pchs.2[sub$nut.num],col=allcols[sub$pond.id],yaxt='n',xaxt='n')
    }
    #   sub.low <- filter(sub, nut==1)
    #   lines(predict(loess(sub.low[,14]~sub.low[,15])), col=alpha(1,0.5), lwd=1.5,lty=1)
    #   sub.high <- filter(sub, nut==2)
    #   lines(predict(loess(sub.high[,14]~sub.high[,15])), col=alpha(1,0.5), lwd=1.5,lty=2)
  }
}
mtext('proportion diatoms',side=2,outer=T,line=2.5,cex=1.2)
mtext('pesticide concentration',side=1,outer=T,line=2.5,cex=1.2)
mtext(paste('day',Sampling.dates,' '),side=3,outer=T,line=0.1,at=seq(0.1,0.93,length.out = 6),adj=0.5)


##### 3D plots ####

VF.3D.plot <- function(model.name,label){
  fitted <- as.data.frame(predict(model.name))
  zlims <- range(fitted)
  for(nut.lvl in c('low','high')){
    for(date.x in Sampling.dates){
      vis.gam(model.name, view = c('sc.gly','sc.imi'), cond = list('date' = date.x, 'o.nut' = nut.lvl), zlim=zlims, plot.type = 'persp', color = 'topo', main = NULL,xlab='glyphosate',ylab='imidacloprid',zlab=label,theta=45)
    }
  }
}

pdf('~/Desktop/3dplots.pdf',width=14,height=6,pointsize = 8,onefile = T)
par(mfrow=c(2,6),mar=c(0,0,0,0),oma=c(2,2,2,2),cex=1,xpd=T)
VF.3D.plot(ba.m,'BA')
VF.3D.plot(chla.m,'chlorophyll')
VF.3D.plot(zoo.m,'zoo biomass')
dev.off()

pdf('~/Desktop/3dplots_ef.pdf',width=14,height=6,pointsize = 8,onefile = T)
par(mfrow=c(2,6),mar=c(0,0,0,0),oma=c(2,2,2,2),cex=1,xpd=T)
VF.3D.plot(use.m,'CUP')
VF.3D.plot(nep.m,'NEP')
VF.3D.plot(rue.m,'RUE')
dev.off()

