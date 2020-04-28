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

ba.model <- gam(log10(BA) ~ o.nut + ti(date,k=5) + ti(date,sc.gly, k=6) + ti(date,sc.gly, by = o.nut, k=3) + ti(date,sc.imi, k=6) + ti(date,sc.imi, by = o.nut, k=3) + ti(date,sc.gly,sc.imi, k=5) + ti(date,sc.gly,sc.imi, by = o.nut, k=3) + s(date, site.f, bs='fs',k=3), data=merged.data, method = 'REML')
chla.model <- gam(log10(total) ~ o.nut + ti(date,k=5) + ti(date,sc.gly, k=6) + ti(date,sc.gly, by = o.nut, k=3) + ti(date,sc.imi, k=6) + ti(date,sc.imi, by = o.nut, k=3) + ti(date,sc.gly,sc.imi, k=5) + ti(date,sc.gly,sc.imi, by = o.nut, k=3) + s(date, site.f, bs='fs',k=3), data=merged.data, method = 'REML')
zoo.model <- gam(log10p(total_zoo) ~ o.nut + ti(date,k=5) + ti(date,sc.gly, k=6) + ti(date,sc.gly, by = o.nut, k=3) + ti(date,sc.imi, k=6) + ti(date,sc.imi, by = o.nut, k=3) + ti(date,sc.gly,sc.imi, k=5) + ti(date,sc.gly,sc.imi, by = o.nut, k=3) + s(date, site.f, bs='fs',k=3), data=merged.data, method = 'REML')

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

###### a simple visualization of gamms ####

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
VF.gam.plot(model.name=ba.model, varname=expression(log[10]~bacterial~abundance~(cells/mu*L)))
VF.gam.plot(model.name=chla.model, varname=expression(log[10]~chlorophyll~italic(a)~(mu*g/L)))
VF.gam.plot(model.name=zoo.model, varname=expression(log[10](1+zooplankton~biomass)~(mu*g/L)))
mtext('date',1,outer=T,line=0.5)
dev.off()

#### only max vs. min ####

VF.gam.plot.simpler <- function(model.name, varname){
  fitted <- as.data.frame(predict(model.name, se.fit = T,exclude='s(date,site.f)'))
  #fitted <- as.data.frame(predict(model.name, se.fit = T,exclude='s(site.f)'))
  fitted$lwr <- fitted$fit - 1.96*fitted$se.fit
  fitted$upr <- fitted$fit + 1.96*fitted$se.fit
  ylims <- range(c(fitted$lwr,fitted$upr))
  plot_smooth(model.name, view="date", cond=list('sc.gly'=0,'sc.imi'=0,'o.nut'='low'), rm.ranef=TRUE, rug=F,col='gray50',print.summary = F,hide.label = T,cex.axis=1,ann=F,bty='l',legend_plot_all = F, h0=NA,ylim = ylims)
  plot_smooth(model.name, view="date", cond=list('sc.gly'=0.5,'sc.imi'=0,'o.nut'='low'), rm.ranef=TRUE, rug=F,col=gly.cols[4],add=T,print.summary = F)
  plot_smooth(model.name, view="date", cond=list('sc.gly'=1,'sc.imi'=0,'o.nut'='low'), rm.ranef=TRUE, rug=F,col=gly.cols[8],add=T,print.summary = F)
  legend('topleft',legend='low',bty='n')
  mtext(varname,2,line=2.5)
  plot_smooth(model.name, view="date", cond=list('sc.gly'=0,'sc.imi'=0,'o.nut'='high'), rm.ranef=TRUE, rug=F,col='gray50',print.summary = F,hide.label = T,cex.axis=1,ann=F,bty='l',legend_plot_all = F, h0=NA,ylim = ylims)
  plot_smooth(model.name, view="date", cond=list('sc.gly'=0.5,'sc.imi'=0,'o.nut'='high'), rm.ranef=TRUE, rug=F,col=gly.cols[4],add=T,print.summary = F)
  plot_smooth(model.name, view="date", cond=list('sc.gly'=1,'sc.imi'=0,'o.nut'='high'), rm.ranef=TRUE, rug=F,col=gly.cols[8],add=T,print.summary = F)
  legend('topleft',legend='high',bty='n')
  plot_smooth(model.name, view="date", cond=list('sc.gly'=0,'sc.imi'=0,'o.nut'='low'), rm.ranef=TRUE, rug=F,col='gray50',print.summary = F,hide.label = T,cex.axis=1,ann=F,bty='l',legend_plot_all = F, h0=NA,ylim = ylims)
  plot_smooth(model.name, view="date", cond=list('sc.gly'=0,'sc.imi'=0.5,'o.nut'='low'), rm.ranef=TRUE, rug=F,col=imi.cols[4],add=T,print.summary = F)
  plot_smooth(model.name, view="date", cond=list('sc.gly'=0,'sc.imi'=1,'o.nut'='low'), rm.ranef=TRUE, rug=F,col=imi.cols[8],add=T,print.summary = F)
  legend('topleft',legend='low',bty='n')
  plot_smooth(model.name, view="date", cond=list('sc.gly'=0,'sc.imi'=0,'o.nut'='high'), rm.ranef=TRUE, rug=F,col='gray50',print.summary = F,hide.label = T,cex.axis=1,ann=F,bty='l',legend_plot_all = F, h0=NA,ylim = ylims)
  plot_smooth(model.name, view="date", cond=list('sc.gly'=0,'sc.imi'=0.5,'o.nut'='high'), rm.ranef=TRUE, rug=F,col=imi.cols[4],add=T,print.summary = F)
  plot_smooth(model.name, view="date", cond=list('sc.gly'=0,'sc.imi'=1,'o.nut'='high'), rm.ranef=TRUE, rug=F,col=imi.cols[8],add=T,print.summary = F)
  legend('topleft',legend='high',bty='n')
  plot_smooth(model.name, view="date", cond=list('sc.gly'=0,'sc.imi'=0,'o.nut'='low'), rm.ranef=TRUE, rug=F,col='gray50',print.summary = F,hide.label = T,cex.axis=1,ann=F,bty='l',legend_plot_all = F, h0=NA,ylim = ylims)
  plot_smooth(model.name, view="date", cond=list('sc.gly'=0.5,'sc.imi'=0.5,'o.nut'='low'), rm.ranef=TRUE, rug=F,col=both.cols[4],add=T,print.summary = F)
  plot_smooth(model.name, view="date", cond=list('sc.gly'=1,'sc.imi'=1,'o.nut'='low'), rm.ranef=TRUE, rug=F,col=both.cols[8],add=T,print.summary = F)
  legend('topleft',legend='low',bty='n')
  plot_smooth(model.name, view="date", cond=list('sc.gly'=0,'sc.imi'=0,'o.nut'='high'), rm.ranef=TRUE, rug=F,col='gray50',print.summary = F,hide.label = T,cex.axis=1,ann=F,bty='l',legend_plot_all = F, h0=NA,ylim = ylims)
  plot_smooth(model.name, view="date", cond=list('sc.gly'=0.5,'sc.imi'=0.5,'o.nut'='high'), rm.ranef=TRUE, rug=F,col=both.cols[4],add=T,print.summary = F)
  plot_smooth(model.name, view="date", cond=list('sc.gly'=1,'sc.imi'=1,'o.nut'='high'), rm.ranef=TRUE, rug=F,col=both.cols[8],add=T,print.summary = F)
  legend('topleft',legend='high',bty='n')
}

pdf('~/Desktop/gamms_simpler.pdf',width=15,height=10,pointsize = 12)
par(mfrow=c(3,6),mar=c(2,2,2,2),oma=c(2,2,0,0),cex=1)
VF.gam.plot.simpler(model.name=ba.model, varname=expression(log[10]~bacterial~abundance~(cells/mu*L)))
VF.gam.plot.simpler(model.name=chla.model, varname=expression(log[10]~chlorophyll~italic(a)~(mu*g/L)))
VF.gam.plot.simpler(model.name=zoo.model, varname=expression(log[10](1+zooplankton~biomass)~(mu*g/L)))
mtext('date',1,outer=T,line=0.5)
dev.off()

#### contour plots ####

heatcolfunc <- colorRampPalette(brewer.pal(11, 'RdYlBu')[11:1])
heat.cols <- heatcolfunc(50)
#heat.cols <- wes_palette('Zissou1', 50, 'continuous')
#heat.cols <- viridis(50)

VF.cont.plot <- function(model.name){
  zlims <- range( fitted <- as.data.frame(predict(model.name, se.fit = F,exclude='s(date,site.f)')))
  for(nut.lvl in c('low','high')){
    for(date.x in Sampling.dates){
      fvisgam(model.name, view = c('sc.gly','sc.imi'), too.far=0.3,cond = list('date' = date.x, 'o.nut' = nut.lvl), zlim=zlims, add.color.legend=F,hide.label=T,plot.type = 'contour', lwd=1.5,color = heat.cols, main = NULL,rm.ranef = T,dec=1,print.summary = F,yaxt='n',xaxt='n',xlab=NULL,ylab=NULL)
    }
  }
  mtext('glyphosate',1,outer=T,line=0.5)
  mtext('imidacloprid',2,outer=T,line=0.5)
  mtext(paste('day',Sampling.dates,' '),side=3,outer=T,line=0.5,at=seq(0.09,0.92,length.out = 6),adj=0.5)
  mtext(c('low nut','high nut'),side=4,outer=T,line=0.5,at=c(0.75,0.25))
}

pdf('~/Desktop/contourplots.pdf',width=7.5,height=3,pointsize = 12,onefile = T)
par(mfrow=c(2,6),mar=c(0.1,0.1,0.1,0.1),oma=c(2,2,2,2),cex=1,xpd=T)
VF.cont.plot(ba.model)
VF.cont.plot(chla.model)
VF.cont.plot(zoo.model)
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
VF.3D.plot(ba.model,'BA')
VF.3D.plot(chla.model,'chlorophyll')
VF.3D.plot(zoo.model,'zoo biomass')
dev.off()

##### Linear models to quantify effect sizes#####

merged.data$nut.num <- merged.data$nut.num-1

ba.mod.lin <- lmer(log10(BA) ~ date.f:nut.num + date.f:sc.gly + date.f:sc.imi + date.f:sc.gly:sc.imi + date.f:sc.gly:nut.num +
                     date.f:sc.imi:nut.num + (1|site.f),merged.data)

plot_model(ba.mod.lin,type='est')
plot_model(ba.mod.lin,type='diag')
performance(ba.mod.lin)
plot(fitted(ba.mod.lin)~log10(merged.data$BA))
ba.coefs <- get_model_data(ba.mod.lin, type = 'est')

chla.mod.lin <- lmer(log10(total) ~ date.f:nut.num + date.f:sc.gly + date.f:sc.imi + date.f:sc.gly:sc.imi + date.f:sc.gly:nut.num +
                       date.f:sc.imi:nut.num + (1|site.f),merged.data)
plot_model(chla.mod.lin,type='est')
plot_model(chla.mod.lin,type='diag')
performance(chla.mod.lin)
plot(fitted(chla.mod.lin)~log10(merged.data$BA))
chla.coefs <- get_model_data(chla.mod.lin, type = 'est')

zoo.mod.lin <- lmer(log10p(total_zoo_adult) ~ date.f:nut.num + date.f:sc.gly + date.f:sc.imi + date.f:sc.gly:sc.imi + date.f:sc.gly:nut.num +
                       date.f:sc.imi:nut.num + (1|site.f),merged.data)
plot_model(zoo.mod.lin,type='est')
plot_model(zoo.mod.lin,type='diag')
performance(zoo.mod.lin)
plot(fitted(zoo.mod.lin)~log10p(merged.data$total_zoo_adult))
zoo.coefs <- get_model_data(zoo.mod.lin, type = 'est')

EP.mod.lin <- lmer(use ~ date.f:nut.num + date.f:sc.gly + date.f:sc.imi + date.f:sc.gly:sc.imi + date.f:sc.gly:nut.num +
                     date.f:sc.imi:nut.num + (1|site.f),merged.data)
plot_model(EP.mod.lin,type='est')
plot_model(EP.mod.lin,type='diag')
performance(EP.mod.lin)
plot(fitted(EP.mod.lin)~merged.data$use)

nep.mod.lin <- lmer(log10p(NEP) ~ date.f:nut.num + date.f:sc.gly + date.f:sc.imi + date.f:sc.gly:sc.imi + date.f:sc.gly:nut.num +
                      date.f:sc.imi:nut.num + (1|site.f),merged.data)
summary(nep.mod.lin)
plot_model(nep.mod.lin,type='est')
plot_model(nep.mod.lin,type='diag')
performance(nep.mod.lin)
plot(fitted(nep.mod.lin)~log10p(merged.data$NEP))

# coef plot

plot_models(ba.mod.lin,chla.mod.lin,std.est = "std2")
tab_model(ba.mod.lin,chla.mod.lin)

par(mfrow=c(7,1),mar=c(2,4,1,1),oma=c(2,0,0,0),cex=1)

emptyplot(xlim=c())
