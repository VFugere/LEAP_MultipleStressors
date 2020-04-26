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

chla.model <- gam(log10(total) ~ o.nut + ti(date,k=4) + ti(sc.gly, k = 4) + ti(sc.imi, k = 4) + ti(date,sc.gly, k=6) + ti(date,sc.gly, by = o.nut, k=3) + ti(date,sc.imi, k=6) + ti(date,sc.imi, by = o.nut, k=3) + ti(date,sc.gly,sc.imi, k=6) + ti(date,sc.gly,sc.imi, by = o.nut, k=3) + s(site.f, bs='re'), data=merged.data, method = 'REML')
gam.check(chla.model)
summary(chla.model)
plot(chla.model)

ba.model <- gam(log10(BA) ~ o.nut + ti(date,k=6) + ti(sc.gly, k = 3) + ti(sc.imi, k = 3) + ti(date,sc.gly, k=6) + ti(date,sc.gly, by = o.nut, k=3) + ti(date,sc.imi, k=6) + ti(date,sc.imi, by = o.nut, k=3) + ti(date,sc.gly,sc.imi, k=6) + ti(date,sc.gly,sc.imi, by = o.nut, k=3) + s(site.f, bs='re'), data=merged.data, method = 'REML')
gam.check(ba.model)
summary(ba.model)
plot(ba.model)

m1 <- gam(log10(BA) ~ o.nut + ti(date,k=6) + ti(sc.gly, k = 3) + ti(sc.imi, k = 3) + ti(date,sc.gly, k=6) + ti(date,sc.gly, by = o.nut, k=3) + ti(date,sc.imi, k=6) + ti(date,sc.imi, by = o.nut, k=3) + ti(date,sc.gly,sc.imi, k=6) + ti(date,sc.gly,sc.imi, by = o.nut, k=3) + s(site.f, bs='re'), data=merged.data, method = 'REML')
m2 <- gam(log10(BA) ~ o.nut + ti(date,k=4) + ti(sc.gly, k = 3) + ti(sc.imi, k = 3) + ti(date,sc.gly, k=4) + ti(date,sc.gly, by = o.nut, k=3) + ti(date,sc.imi, k=4) + ti(date,sc.imi, by = o.nut, k=3) + ti(date,sc.gly,sc.imi, k=4) + ti(date,sc.gly,sc.imi, by = o.nut, k=3) + s(date, site.f, bs='fs',k=4), data=merged.data, method = 'REML')
AIC(m1,m2)
plot(m2)

#' model with random smooths is doing a lot better (R2 and AIC). However, the basis dimension
#' is so low that random smooths are in fact straight lines. Moreover, the need to reduce k
#' considerably for all smooths leads to a poorer representation of pesticide effects.
#' If I include random smooths, need to find a way to increase k by reducing number of effects

m2 <- gam(log10(BA) ~ o.nut + ti(date,k=4) + ti(sc.gly, k = 3) + ti(sc.imi, k = 3) + ti(date,sc.gly, k=4) + ti(date,sc.gly, by = o.nut, k=3) + ti(date,sc.imi, k=4) + ti(date,sc.imi, by = o.nut, k=3) + ti(date,sc.gly,sc.imi, k=4) + ti(date,sc.gly,sc.imi, by = o.nut, k=3) + s(date, site.f, bs='fs',k=4), data=merged.data, method = 'REML')


###### a simple visualization of gamms ####

VF.gam.plot <- function(model.name, varname){
  #fitted <- as.data.frame(predict(model.name, se.fit = T,exclude='s(date,site)'))
  fitted <- as.data.frame(predict(model.name, se.fit = T,exclude='s(site.f)'))
  fitted$lwr <- fitted$fit - 1.96*fitted$se.fit
  fitted$upr <- fitted$fit + 1.96*fitted$se.fit
  #ylims <- range(c(fitted$lwr,fitted$upr))
  ylims <- range(fitted$fit)
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
mtext('date',1,outer=T,line=0.5)
dev.off()

VF.gam.plot.simpler <- function(model.name, varname){
  #fitted <- as.data.frame(predict(model.name, se.fit = T,exclude='s(date,site)'))
  fitted <- as.data.frame(predict(model.name, se.fit = T,exclude='s(site.f)'))
  fitted$lwr <- fitted$fit - 1.96*fitted$se.fit
  fitted$upr <- fitted$fit + 1.96*fitted$se.fit
  ylims <- range(c(fitted$lwr,fitted$upr))
  plot_smooth(model.name, view="date", cond=list('sc.gly'=0,'sc.imi'=0,'o.nut'='low'), rm.ranef=TRUE, rug=F,col=gly.cols[1],print.summary = F,hide.label = T,cex.axis=1,ann=F,bty='l',legend_plot_all = F, h0=NA,ylim = ylims)
  plot_smooth(model.name, view="date", cond=list('sc.gly'=1,'sc.imi'=0,'o.nut'='low'), rm.ranef=TRUE, rug=F,col=gly.cols[8],add=T,print.summary = F)
  legend('topleft',legend='low',bty='n')
  mtext(varname,2,line=2.5)
  plot_smooth(model.name, view="date", cond=list('sc.gly'=0,'sc.imi'=0,'o.nut'='high'), rm.ranef=TRUE, rug=F,col=gly.cols[1],print.summary = F,hide.label = T,cex.axis=1,ann=F,bty='l',legend_plot_all = F, h0=NA,ylim = ylims)
  plot_smooth(model.name, view="date", cond=list('sc.gly'=1,'sc.imi'=0,'o.nut'='high'), rm.ranef=TRUE, rug=F,col=gly.cols[8],add=T,print.summary = F)
  legend('topleft',legend='high',bty='n')
  plot_smooth(model.name, view="date", cond=list('sc.gly'=0,'sc.imi'=0,'o.nut'='low'), rm.ranef=TRUE, rug=F,col=imi.cols[1],print.summary = F,hide.label = T,cex.axis=1,ann=F,bty='l',legend_plot_all = F, h0=NA,ylim = ylims)
  plot_smooth(model.name, view="date", cond=list('sc.gly'=0,'sc.imi'=1,'o.nut'='low'), rm.ranef=TRUE, rug=F,col=imi.cols[8],add=T,print.summary = F)
  legend('topleft',legend='low',bty='n')
  plot_smooth(model.name, view="date", cond=list('sc.gly'=0,'sc.imi'=0,'o.nut'='high'), rm.ranef=TRUE, rug=F,col=imi.cols[1],print.summary = F,hide.label = T,cex.axis=1,ann=F,bty='l',legend_plot_all = F, h0=NA,ylim = ylims)
  plot_smooth(model.name, view="date", cond=list('sc.gly'=0,'sc.imi'=1,'o.nut'='high'), rm.ranef=TRUE, rug=F,col=imi.cols[8],add=T,print.summary = F)
  legend('topleft',legend='high',bty='n')
  plot_smooth(model.name, view="date", cond=list('sc.gly'=0,'sc.imi'=0,'o.nut'='low'), rm.ranef=TRUE, rug=F,col=both.cols[1],print.summary = F,hide.label = T,cex.axis=1,ann=F,bty='l',legend_plot_all = F, h0=NA,ylim = ylims)
  plot_smooth(model.name, view="date", cond=list('sc.gly'=1,'sc.imi'=1,'o.nut'='low'), rm.ranef=TRUE, rug=F,col=both.cols[8],add=T,print.summary = F)
  legend('topleft',legend='low',bty='n')
  plot_smooth(model.name, view="date", cond=list('sc.gly'=0,'sc.imi'=0,'o.nut'='high'), rm.ranef=TRUE, rug=F,col=both.cols[1],print.summary = F,hide.label = T,cex.axis=1,ann=F,bty='l',legend_plot_all = F, h0=NA,ylim = ylims)
  plot_smooth(model.name, view="date", cond=list('sc.gly'=1,'sc.imi'=1,'o.nut'='high'), rm.ranef=TRUE, rug=F,col=both.cols[8],add=T,print.summary = F)
  legend('topleft',legend='high',bty='n')
}

pdf('~/Desktop/gamms_simpler.pdf',width=15,height=10,pointsize = 12)
par(mfrow=c(3,6),mar=c(2,2,2,2),oma=c(2,2,0,0),cex=1)
VF.gam.plot.simpler(model.name=ba.model, varname=expression(log[10]~bacterial~abundance~(cells/mu*L)))
VF.gam.plot.simpler(model.name=chla.model, varname=expression(log[10]~chlorophyll~italic(a)~(mu*g/L)))
mtext('date',1,outer=T,line=0.5)
dev.off()

#### contour plots ####

heatcolfunc <- colorRampPalette(brewer.pal(11, 'RdYlBu')[11:1])
heat.cols <- heatcolfunc(50)
#heat.cols <- wes_palette('Zissou1', 50, 'continuous')
#heat.cols <- viridis(50)

VF.cont.plot <- function(model.name){
  zlims <- range(fitted(model.name))
  for(nut.lvl in c('low','high')){
    for(date.x in Sampling.dates){
      fvisgam(model.name, view = c('sc.gly','sc.imi'), too.far=0.25,cond = list('date' = date.x, 'o.nut' = nut.lvl), zlim=zlims, add.color.legend=F,hide.label=T,plot.type = 'contour', lwd=1.5,color = heat.cols, main = NULL,rm.ranef = T,dec=1,print.summary = F,yaxt='n',xaxt='n',xlab=NULL,ylab=NULL)
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
dev.off()

#####

vis.gam(ba.model)

#chla

zlims <- range(fitted(chla.model))

fvisgam(chla.model, view = c('sc.gly','sc.imi'), cond = list('date' = 1, 'o.nut' = 'low'), zlim=zlims, add.color.legend=F,hide.label=T,plot.type = 'contour', lwd=1.5,color = magma(100), main = NULL,rm.ranef = T,dec=1,print.summary = F)
mtext('LN day 1',3,cex=1)
fvisgam(chla.model, view = c('sc.gly','sc.imi'), cond = list('date' = 7, 'o.nut' = 'low'), zlim=zlims, add.color.legend=F,hide.label=T,plot.type = 'contour', lwd=1.5,color = magma(100), main = NULL,rm.ranef = T,dec=1,print.summary = F)
mtext('LN day 7',3,cex=1)
fvisgam(chla.model, view = c('sc.gly','sc.imi'), cond = list('date' = 15, 'o.nut' = 'low'), zlim=zlims, add.color.legend=F,hide.label=T,plot.type = 'contour', lwd=1.5,color = magma(100), main = NULL,rm.ranef = T,dec=1,print.summary = F)
mtext('LN day 15',3,cex=1)
fvisgam(chla.model, view = c('sc.gly','sc.imi'), cond = list('date' = 30, 'o.nut' = 'low'), zlim=zlims, add.color.legend=F,hide.label=T,plot.type = 'contour', lwd=1.5,color = magma(100), main = NULL,rm.ranef = T,dec=1,print.summary = F)
mtext('LN day 30',3,cex=1)
fvisgam(chla.model, view = c('sc.gly','sc.imi'), cond = list('date' = 35, 'o.nut' = 'low'), zlim=zlims, add.color.legend=F,hide.label=T,plot.type = 'contour', lwd=1.5,color = magma(100), main = NULL,rm.ranef = T,dec=1,print.summary = F)
mtext('LN day 35',3,cex=1)
fvisgam(chla.model, view = c('sc.gly','sc.imi'), cond = list('date' = 43, 'o.nut' = 'low'), zlim=zlims, add.color.legend=T,hide.label=T,plot.type = 'contour', lwd=1.5,color = magma(100), main = NULL,rm.ranef = T,dec=1,print.summary = F)
mtext('LN day 43',3,cex=1)
fvisgam(chla.model, view = c('sc.gly','sc.imi'), cond = list('date' = 1, 'o.nut' = 'high'), zlim=zlims, add.color.legend=F,hide.label=T,plot.type = 'contour', lwd=1.5,color = magma(100), main = NULL,rm.ranef = T,dec=1,print.summary = F)
mtext('HN day 1',3,cex=1)
fvisgam(chla.model, view = c('sc.gly','sc.imi'), cond = list('date' = 7, 'o.nut' = 'high'), zlim=zlims, add.color.legend=F,hide.label=T,plot.type = 'contour', lwd=1.5,color = magma(100), main = NULL,rm.ranef = T,dec=1,print.summary = F)
mtext('HN day 7',3,cex=1)
fvisgam(chla.model, view = c('sc.gly','sc.imi'), cond = list('date' = 15, 'o.nut' = 'high'), zlim=zlims, add.color.legend=F,hide.label=T,plot.type = 'contour', lwd=1.5,color = magma(100), main = NULL,rm.ranef = T,dec=1,print.summary = F)
mtext('HN day 15',3,cex=1)
fvisgam(chla.model, view = c('sc.gly','sc.imi'), cond = list('date' = 30, 'o.nut' = 'high'), zlim=zlims, add.color.legend=F,hide.label=T,plot.type = 'contour', lwd=1.5,color = magma(100), main = NULL,rm.ranef = T,dec=1,print.summary = F)
mtext('HN day 30',3,cex=1)
fvisgam(chla.model, view = c('sc.gly','sc.imi'), cond = list('date' = 35, 'o.nut' = 'high'), zlim=zlims, add.color.legend=F,hide.label=T,plot.type = 'contour', lwd=1.5,color = magma(100), main = NULL,rm.ranef = T,dec=1,print.summary = F)
mtext('HN day 35',3,cex=1)
fvisgam(chla.model, view = c('sc.gly','sc.imi'), cond = list('date' = 43, 'o.nut' = 'high'), zlim=zlims, add.color.legend=T,hide.label=T,plot.type = 'contour', lwd=1.5,color = magma(100), main = NULL,rm.ranef = T,dec=1,print.summary = F)
mtext('HN day 43',3,cex=1)

mtext('glyphosate',1,outer=T,line=0.5)
mtext('imidacloprid',2,outer=T,line=0.5)

dev.off()

##### trying out mgcViz package ####



#### linear model ####

ba.mod.lin <- lmer(log10(BA) ~ date.f + nut.fac:date.f + date.f:sc.gly + date.f:sc.imi + date.f:sc.gly:sc.imi + date.f:sc.gly:nut.fac +
                     date.f:sc.imi:nut.fac + date.f:sc.gly:sc.imi:nut.fac + (1|site.f),merged.data)

plot_model(ba.mod.lin,type='est')
plot_model(ba.mod.lin,type='re')
plot_model(ba.mod.lin,type='diag')
performance(ba.mod.lin)
plot(fitted(ba.mod.lin)~log10(merged.data$BA))
ba.coefs <- get_model_data(ba.mod.lin, type = 'est')

chla.mod.lin <- lmer(log10(total) ~ date.f + nut.fac:date.f + date.f:sc.gly + date.f:sc.imi + date.f:sc.gly:sc.imi + date.f:sc.gly:nut.fac +
                       date.f:sc.imi:nut.fac + date.f:sc.gly:sc.imi:nut.fac + (1|site.f),merged.data)

plot_model(chla.mod.lin,type='est')
plot_model(chla.mod.lin,type='re')
plot_model(chla.mod.lin,type='diag')
performance(chla.mod.lin)
plot(fitted(chla.mod.lin)~log10(merged.data$BA))
chla.coefs <- get_model_data(chla.mod.lin, type = 'est')

EP.mod.lin <- lmer(use ~ date.f + nut.fac + date.f:sc.gly + date.f:sc.imi + date.f:sc.gly:sc.imi + date.f:sc.gly:nut.fac +
                     date.f:sc.imi:nut.fac + date.f:sc.gly:sc.imi:nut.fac + (1|site.f),merged.data)

plot_model(EP.mod.lin,type='est')
plot_model(EP.mod.lin,type='re')
plot_model(EP.mod.lin,type='diag')
performance(EP.mod.lin)
plot(fitted(EP.mod.lin)~merged.data$use);abline(a=0,b=1,lty=3)

# coef plot

plot_models(ba.mod.lin,chla.mod.lin,std.est = "std2")
tab_model(ba.mod.lin,chla.mod.lin)

par(mfrow=c(7,1),mar=c(2,4,1,1),oma=c(2,0,0,0),cex=1)

emptyplot(xlim=c())
