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

ba.mod.lin <- lmer(log10(BA) ~ date.f:nut.num + date.f:sc.gly + date.f:sc.imi + date.f:sc.gly:sc.imi + date.f:sc.gly:nut.num +
                     date.f:sc.imi:nut.num + (1|site.f),merged.data)
# plot_model(ba.mod.lin,type='est')
# plot_model(ba.mod.lin,type='diag')
# performance(ba.mod.lin)
# plot(fitted(ba.mod.lin)~log10(merged.data$BA))
ba.coefs <- get_model_data(ba.mod.lin, type = 'est')

chla.mod.lin <- lmer(log10(total) ~ date.f:nut.num + date.f:sc.gly + date.f:sc.imi + date.f:sc.gly:sc.imi + date.f:sc.gly:nut.num +
                       date.f:sc.imi:nut.num + (1|site.f),merged.data)
# plot_model(chla.mod.lin,type='est')
# plot_model(chla.mod.lin,type='diag')
# performance(chla.mod.lin)
# plot(fitted(chla.mod.lin)~log10(merged.data$BA))
chla.coefs <- get_model_data(chla.mod.lin, type = 'est')

zoo.mod.lin <- lmer(log10p(total_zoo) ~ date.f:nut.num + date.f:sc.gly + date.f:sc.imi + date.f:sc.gly:sc.imi + date.f:sc.gly:nut.num +
                      date.f:sc.imi:nut.num + (1|site.f),merged.data)
# plot_model(zoo.mod.lin,type='est')
# plot_model(zoo.mod.lin,type='diag')
# performance(zoo.mod.lin)
# plot(fitted(zoo.mod.lin)~log10p(merged.data$total_zoo_adult))
zoo.coefs <- get_model_data(zoo.mod.lin, type = 'est')

EP.mod.lin <- lmer(use ~ date.f:nut.num + date.f:sc.gly + date.f:sc.imi + date.f:sc.gly:sc.imi + date.f:sc.gly:nut.num +
                     date.f:sc.imi:nut.num + (1|site.f),merged.data)
# plot_model(EP.mod.lin,type='est')
# plot_model(EP.mod.lin,type='diag')
# performance(EP.mod.lin)
# plot(fitted(EP.mod.lin)~merged.data$use)
use.coefs <- get_model_data(EP.mod.lin, type = 'est')


nep.mod.lin <- lmer(log10p(NEP) ~ date.f:nut.num + date.f:sc.gly + date.f:sc.imi + date.f:sc.gly:sc.imi + date.f:sc.gly:nut.num +
                      date.f:sc.imi:nut.num + (1|site.f),merged.data)
# summary(nep.mod.lin)
# plot_model(nep.mod.lin,type='est')
# plot_model(nep.mod.lin,type='diag')
# performance(nep.mod.lin)
# plot(fitted(nep.mod.lin)~log10p(merged.data$NEP))
nep.coefs <- get_model_data(nep.mod.lin, type = 'est')

rue.mod.lin <- lmer(log10(RUE) ~ date.f:nut.num + date.f:sc.gly + date.f:sc.imi + date.f:sc.gly:sc.imi + date.f:sc.gly:nut.num +
                      date.f:sc.imi:nut.num + (1|site.f),merged.data,subset=is.finite(log(RUE)))
# summary(rue.mod.lin)
# plot_model(rue.mod.lin,type='est')
# plot_model(rue.mod.lin,type='diag')
# performance(rue.mod.lin)
rue.coefs <- get_model_data(rue.mod.lin, type = 'est')

#tab_model(ba.mod.lin,chla.mod.lin)

#### coef plot ####

pdf('~/Desktop/coef_plot.pdf',width=6,height = 6,pointsize = 12)

# var.names <- c(expression(log[10]~BA~(cells/mu*L)),
#                expression(log[10]~chl.~italic(a)~(mu*g/L)),
#                expression(log[10](1+zoo)~(mu*g/L)),
#                expression(C~sources~used),
#                expression(log[10]~(1+Delta*DO)~(mu*g/L)),
#                expression(log[10]~RUE~(mu*g/mu*g)))

var.names <- c('Bact.','Chl. a','Zoo.','C use',expression(Delta*DO),'RUE')

par(mfrow=c(6,6),mar=c(0.1,0.1,0.1,0.1),oma=c(4,4.5,3,0.5),cex=1)

#BA#

ba.coefs$alpha <- 0.2
ba.coefs$alpha[str_detect(ba.coefs$p.stars, '\\*')] <- 1
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(ba.coefs$conf.low,ba.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
abline(h=0,col='gray80');abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=ba.coefs$estimate[1:6],pch=15,col=scales::alpha('black',ba.coefs$alpha[1:6]))
arrows(x0=1:6,y0=ba.coefs$conf.low[1:6],y1=ba.coefs$conf.high[1:6],length=0,col=scales::alpha('black',ba.coefs$alpha[1:6]))
mtext(var.names[1],side=2,outer=F,line=1,cex=1)
axis(2,lwd=0,lwd.ticks = 1,at=0,labels='')
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(ba.coefs$conf.low,ba.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
abline(h=0,col='gray80');abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=ba.coefs$estimate[7:12],pch=16,col=scales::alpha(gly.cols[8],ba.coefs$alpha[7:12]))
arrows(x0=1:6,y0=ba.coefs$conf.low[7:12],y1=ba.coefs$conf.high[7:12],length=0,col=scales::alpha(gly.cols[8],ba.coefs$alpha[7:12]))
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(ba.coefs$conf.low,ba.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
abline(h=0,col='gray80');abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=ba.coefs$estimate[25:30],pch=15,col=scales::alpha(gly.cols[8],ba.coefs$alpha[25:30]))
arrows(x0=1:6,y0=ba.coefs$conf.low[25:30],y1=ba.coefs$conf.high[25:30],length=0,col=scales::alpha(gly.cols[8],ba.coefs$alpha[25:30]))
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(ba.coefs$conf.low,ba.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
abline(h=0,col='gray80');abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=ba.coefs$estimate[13:18],pch=16,col=scales::alpha(imi.cols[8],ba.coefs$alpha[13:18]))
arrows(x0=1:6,y0=ba.coefs$conf.low[13:18],y1=ba.coefs$conf.high[13:18],length=0,col=scales::alpha(imi.cols[8],ba.coefs$alpha[13:18]))
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(ba.coefs$conf.low,ba.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
abline(h=0,col='gray80');abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=ba.coefs$estimate[31:36],pch=15,col=scales::alpha(imi.cols[8],ba.coefs$alpha[31:36]))
arrows(x0=1:6,y0=ba.coefs$conf.low[31:36],y1=ba.coefs$conf.high[31:36],length=0,col=scales::alpha(imi.cols[8],ba.coefs$alpha[31:36]))
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(ba.coefs$conf.low,ba.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
abline(h=0,col='gray80');abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=ba.coefs$estimate[19:24],pch=16,col=scales::alpha(both.cols[8],ba.coefs$alpha[19:24]))
arrows(x0=1:6,y0=ba.coefs$conf.low[19:24],y1=ba.coefs$conf.high[19:24],length=0,col=scales::alpha(both.cols[8],ba.coefs$alpha[19:24]))

mtext('independent variable',side=3,outer=T,line=1.7,cex=1.3)
mtext(c('nut','gly','gly:nut','imi','imi:nut','gly:imi'),side=3,outer=T,line=0.2,at=seq(0.09,0.91,length.out = 6),adj=0.5,cex=1)
mtext('sampling occasion',side=1,outer=T,line=2.5,cex=1.3)
mtext('response variable',side=2,outer=T,line=3,cex=1.3)

#Chla#

chla.coefs$alpha <- 0.2
chla.coefs$alpha[str_detect(chla.coefs$p.stars, '\\*')] <- 1
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(chla.coefs$conf.low,chla.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
abline(h=0,col='gray80');abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=chla.coefs$estimate[1:6],pch=15,col=scales::alpha('black',chla.coefs$alpha[1:6]))
arrows(x0=1:6,y0=chla.coefs$conf.low[1:6],y1=chla.coefs$conf.high[1:6],length=0,col=scales::alpha('black',chla.coefs$alpha[1:6]))
mtext(var.names[2],side=2,outer=F,line=1,cex=1)
axis(2,lwd=0,lwd.ticks = 1,at=0,labels='')
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(chla.coefs$conf.low,chla.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
abline(h=0,col='gray80');abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=chla.coefs$estimate[7:12],pch=16,col=scales::alpha(gly.cols[8],chla.coefs$alpha[7:12]))
arrows(x0=1:6,y0=chla.coefs$conf.low[7:12],y1=chla.coefs$conf.high[7:12],length=0,col=scales::alpha(gly.cols[8],chla.coefs$alpha[7:12]))
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(chla.coefs$conf.low,chla.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
abline(h=0,col='gray80');abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=chla.coefs$estimate[25:30],pch=15,col=scales::alpha(gly.cols[8],chla.coefs$alpha[25:30]))
arrows(x0=1:6,y0=chla.coefs$conf.low[25:30],y1=chla.coefs$conf.high[25:30],length=0,col=scales::alpha(gly.cols[8],chla.coefs$alpha[25:30]))
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(chla.coefs$conf.low,chla.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
abline(h=0,col='gray80');abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=chla.coefs$estimate[13:18],pch=16,col=scales::alpha(imi.cols[8],chla.coefs$alpha[13:18]))
arrows(x0=1:6,y0=chla.coefs$conf.low[13:18],y1=chla.coefs$conf.high[13:18],length=0,col=scales::alpha(imi.cols[8],chla.coefs$alpha[13:18]))
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(chla.coefs$conf.low,chla.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
abline(h=0,col='gray80');abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=chla.coefs$estimate[31:36],pch=15,col=scales::alpha(imi.cols[8],chla.coefs$alpha[31:36]))
arrows(x0=1:6,y0=chla.coefs$conf.low[31:36],y1=chla.coefs$conf.high[31:36],length=0,col=scales::alpha(imi.cols[8],chla.coefs$alpha[31:36]))
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(chla.coefs$conf.low,chla.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
abline(h=0,col='gray80');abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=chla.coefs$estimate[19:24],pch=16,col=scales::alpha(both.cols[8],chla.coefs$alpha[19:24]))
arrows(x0=1:6,y0=chla.coefs$conf.low[19:24],y1=chla.coefs$conf.high[19:24],length=0,col=scales::alpha(both.cols[8],chla.coefs$alpha[19:24]))

#zoo#

zoo.coefs$alpha <- 0.2
zoo.coefs$alpha[str_detect(zoo.coefs$p.stars, '\\*')] <- 1
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(zoo.coefs$conf.low,zoo.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
abline(h=0,col='gray80');abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=zoo.coefs$estimate[1:6],pch=15,col=scales::alpha('black',zoo.coefs$alpha[1:6]))
arrows(x0=1:6,y0=zoo.coefs$conf.low[1:6],y1=zoo.coefs$conf.high[1:6],length=0,col=scales::alpha('black',zoo.coefs$alpha[1:6]))
mtext(var.names[3],side=2,outer=F,line=1,cex=1)
axis(2,lwd=0,lwd.ticks = 1,at=0,labels='')
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(zoo.coefs$conf.low,zoo.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
abline(h=0,col='gray80');abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=zoo.coefs$estimate[7:12],pch=16,col=scales::alpha(gly.cols[8],zoo.coefs$alpha[7:12]))
arrows(x0=1:6,y0=zoo.coefs$conf.low[7:12],y1=zoo.coefs$conf.high[7:12],length=0,col=scales::alpha(gly.cols[8],zoo.coefs$alpha[7:12]))
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(zoo.coefs$conf.low,zoo.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
abline(h=0,col='gray80');abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=zoo.coefs$estimate[25:30],pch=15,col=scales::alpha(gly.cols[8],zoo.coefs$alpha[25:30]))
arrows(x0=1:6,y0=zoo.coefs$conf.low[25:30],y1=zoo.coefs$conf.high[25:30],length=0,col=scales::alpha(gly.cols[8],zoo.coefs$alpha[25:30]))
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(zoo.coefs$conf.low,zoo.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
abline(h=0,col='gray80');abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=zoo.coefs$estimate[13:18],pch=16,col=scales::alpha(imi.cols[8],zoo.coefs$alpha[13:18]))
arrows(x0=1:6,y0=zoo.coefs$conf.low[13:18],y1=zoo.coefs$conf.high[13:18],length=0,col=scales::alpha(imi.cols[8],zoo.coefs$alpha[13:18]))
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(zoo.coefs$conf.low,zoo.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
abline(h=0,col='gray80');abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=zoo.coefs$estimate[31:36],pch=15,col=scales::alpha(imi.cols[8],zoo.coefs$alpha[31:36]))
arrows(x0=1:6,y0=zoo.coefs$conf.low[31:36],y1=zoo.coefs$conf.high[31:36],length=0,col=scales::alpha(imi.cols[8],zoo.coefs$alpha[31:36]))
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(zoo.coefs$conf.low,zoo.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
abline(h=0,col='gray80');abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=zoo.coefs$estimate[19:24],pch=16,col=scales::alpha(both.cols[8],zoo.coefs$alpha[19:24]))
arrows(x0=1:6,y0=zoo.coefs$conf.low[19:24],y1=zoo.coefs$conf.high[19:24],length=0,col=scales::alpha(both.cols[8],zoo.coefs$alpha[19:24]))

#use#

use.coefs$alpha <- 0.2
use.coefs$alpha[str_detect(use.coefs$p.stars, '\\*')] <- 1
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(use.coefs$conf.low,use.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
abline(h=0,col='gray80');abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=use.coefs$estimate[1:6],pch=15,col=scales::alpha('black',use.coefs$alpha[1:6]))
arrows(x0=1:6,y0=use.coefs$conf.low[1:6],y1=use.coefs$conf.high[1:6],length=0,col=scales::alpha('black',use.coefs$alpha[1:6]))
mtext(var.names[4],side=2,outer=F,line=1,cex=1)
axis(2,lwd=0,lwd.ticks = 1,at=0,labels='')
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(use.coefs$conf.low,use.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
abline(h=0,col='gray80');abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=use.coefs$estimate[7:12],pch=16,col=scales::alpha(gly.cols[8],use.coefs$alpha[7:12]))
arrows(x0=1:6,y0=use.coefs$conf.low[7:12],y1=use.coefs$conf.high[7:12],length=0,col=scales::alpha(gly.cols[8],use.coefs$alpha[7:12]))
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(use.coefs$conf.low,use.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
abline(h=0,col='gray80');abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=use.coefs$estimate[25:30],pch=15,col=scales::alpha(gly.cols[8],use.coefs$alpha[25:30]))
arrows(x0=1:6,y0=use.coefs$conf.low[25:30],y1=use.coefs$conf.high[25:30],length=0,col=scales::alpha(gly.cols[8],use.coefs$alpha[25:30]))
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(use.coefs$conf.low,use.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
abline(h=0,col='gray80');abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=use.coefs$estimate[13:18],pch=16,col=scales::alpha(imi.cols[8],use.coefs$alpha[13:18]))
arrows(x0=1:6,y0=use.coefs$conf.low[13:18],y1=use.coefs$conf.high[13:18],length=0,col=scales::alpha(imi.cols[8],use.coefs$alpha[13:18]))
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(use.coefs$conf.low,use.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
abline(h=0,col='gray80');abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=use.coefs$estimate[31:36],pch=15,col=scales::alpha(imi.cols[8],use.coefs$alpha[31:36]))
arrows(x0=1:6,y0=use.coefs$conf.low[31:36],y1=use.coefs$conf.high[31:36],length=0,col=scales::alpha(imi.cols[8],use.coefs$alpha[31:36]))
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(use.coefs$conf.low,use.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
abline(h=0,col='gray80');abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=use.coefs$estimate[19:24],pch=16,col=scales::alpha(both.cols[8],use.coefs$alpha[19:24]))
arrows(x0=1:6,y0=use.coefs$conf.low[19:24],y1=use.coefs$conf.high[19:24],length=0,col=scales::alpha(both.cols[8],use.coefs$alpha[19:24]))

#nep#

nep.coefs$alpha <- 0.2
nep.coefs$alpha[str_detect(nep.coefs$p.stars, '\\*')] <- 1
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(nep.coefs$conf.low,nep.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
abline(h=0,col='gray80');abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=nep.coefs$estimate[1:6],pch=15,col=scales::alpha('black',nep.coefs$alpha[1:6]))
arrows(x0=1:6,y0=nep.coefs$conf.low[1:6],y1=nep.coefs$conf.high[1:6],length=0,col=scales::alpha('black',nep.coefs$alpha[1:6]))
mtext(var.names[5],side=2,outer=F,line=1,cex=1)
axis(2,lwd=0,lwd.ticks = 1,at=0,labels='')
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(nep.coefs$conf.low,nep.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
abline(h=0,col='gray80');abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=nep.coefs$estimate[7:12],pch=16,col=scales::alpha(gly.cols[8],nep.coefs$alpha[7:12]))
arrows(x0=1:6,y0=nep.coefs$conf.low[7:12],y1=nep.coefs$conf.high[7:12],length=0,col=scales::alpha(gly.cols[8],nep.coefs$alpha[7:12]))
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(nep.coefs$conf.low,nep.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
abline(h=0,col='gray80');abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=nep.coefs$estimate[25:30],pch=15,col=scales::alpha(gly.cols[8],nep.coefs$alpha[25:30]))
arrows(x0=1:6,y0=nep.coefs$conf.low[25:30],y1=nep.coefs$conf.high[25:30],length=0,col=scales::alpha(gly.cols[8],nep.coefs$alpha[25:30]))
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(nep.coefs$conf.low,nep.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
abline(h=0,col='gray80');abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=nep.coefs$estimate[13:18],pch=16,col=scales::alpha(imi.cols[8],nep.coefs$alpha[13:18]))
arrows(x0=1:6,y0=nep.coefs$conf.low[13:18],y1=nep.coefs$conf.high[13:18],length=0,col=scales::alpha(imi.cols[8],nep.coefs$alpha[13:18]))
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(nep.coefs$conf.low,nep.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
abline(h=0,col='gray80');abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=nep.coefs$estimate[31:36],pch=15,col=scales::alpha(imi.cols[8],nep.coefs$alpha[31:36]))
arrows(x0=1:6,y0=nep.coefs$conf.low[31:36],y1=nep.coefs$conf.high[31:36],length=0,col=scales::alpha(imi.cols[8],nep.coefs$alpha[31:36]))
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(nep.coefs$conf.low,nep.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
abline(h=0,col='gray80');abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=nep.coefs$estimate[19:24],pch=16,col=scales::alpha(both.cols[8],nep.coefs$alpha[19:24]))
arrows(x0=1:6,y0=nep.coefs$conf.low[19:24],y1=nep.coefs$conf.high[19:24],length=0,col=scales::alpha(both.cols[8],nep.coefs$alpha[19:24]))

#rue#

rue.coefs$alpha <- 0.2
rue.coefs$alpha[str_detect(rue.coefs$p.stars, '\\*')] <- 1
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(rue.coefs$conf.low,rue.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
abline(h=0,col='gray80');abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=rue.coefs$estimate[1:6],pch=15,col=scales::alpha('black',rue.coefs$alpha[1:6]))
arrows(x0=1:6,y0=rue.coefs$conf.low[1:6],y1=rue.coefs$conf.high[1:6],length=0,col=scales::alpha('black',rue.coefs$alpha[1:6]))
mtext(var.names[6],side=2,outer=F,line=1,cex=1)
axis(2,lwd=0,lwd.ticks = 1,at=0,labels='')
axis(1,lwd=0,lwd.ticks = 1)
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(rue.coefs$conf.low,rue.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
abline(h=0,col='gray80');abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=rue.coefs$estimate[7:12],pch=16,col=scales::alpha(gly.cols[8],rue.coefs$alpha[7:12]))
arrows(x0=1:6,y0=rue.coefs$conf.low[7:12],y1=rue.coefs$conf.high[7:12],length=0,col=scales::alpha(gly.cols[8],rue.coefs$alpha[7:12]))
axis(1,lwd=0,lwd.ticks = 1)
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(rue.coefs$conf.low,rue.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
axis(1,lwd=0,lwd.ticks = 1)
abline(h=0,col='gray80');abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=rue.coefs$estimate[25:30],pch=15,col=scales::alpha(gly.cols[8],rue.coefs$alpha[25:30]))
arrows(x0=1:6,y0=rue.coefs$conf.low[25:30],y1=rue.coefs$conf.high[25:30],length=0,col=scales::alpha(gly.cols[8],rue.coefs$alpha[25:30]))
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(rue.coefs$conf.low,rue.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
axis(1,lwd=0,lwd.ticks = 1)
abline(h=0,col='gray80');abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=rue.coefs$estimate[13:18],pch=16,col=scales::alpha(imi.cols[8],rue.coefs$alpha[13:18]))
arrows(x0=1:6,y0=rue.coefs$conf.low[13:18],y1=rue.coefs$conf.high[13:18],length=0,col=scales::alpha(imi.cols[8],rue.coefs$alpha[13:18]))
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(rue.coefs$conf.low,rue.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
axis(1,lwd=0,lwd.ticks = 1)
abline(h=0,col='gray80');abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=rue.coefs$estimate[31:36],pch=15,col=scales::alpha(imi.cols[8],rue.coefs$alpha[31:36]))
arrows(x0=1:6,y0=rue.coefs$conf.low[31:36],y1=rue.coefs$conf.high[31:36],length=0,col=scales::alpha(imi.cols[8],rue.coefs$alpha[31:36]))
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(rue.coefs$conf.low,rue.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
axis(1,lwd=0,lwd.ticks = 1)
abline(h=0,col='gray80');abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=rue.coefs$estimate[19:24],pch=16,col=scales::alpha(both.cols[8],rue.coefs$alpha[19:24]))
arrows(x0=1:6,y0=rue.coefs$conf.low[19:24],y1=rue.coefs$conf.high[19:24],length=0,col=scales::alpha(both.cols[8],rue.coefs$alpha[19:24]))

dev.off()