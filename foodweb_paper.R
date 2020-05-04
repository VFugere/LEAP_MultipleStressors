## Code to analyze LEAP 2016 MS data
## Vincent Fugère 2019-2020

rm(list=ls())

library(tidyverse)
library(skimr)
library(magrittr)

library(scales)
library(shape)
library(RColorBrewer)
library(viridis)
library(wesanderson)
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

#### Figure panel: experimental design ####

gly.doses <- c(0,0.04,0.1,0.3,0.7,2,5.5,15) #ppm
imi.doses <- c(0,0.15,0.4,1,3,8,22,60) #ppb

#pdf('~/Desktop/figure1b.pdf',width=5.5,height=6,pointsize = 12)
pdf('~/Desktop/figure1.pdf',width=5.5,height=9.5,pointsize = 12)
layout(rbind(1,2),heights = c(6,3.5))
par(cex=1)

plot(x=0,xlim=c(0,9),ylim=c(0,9),type='n',ann=F,yaxt='n',xaxt='n',bty='n')
points(x=1:8,y=rep(8,8),pch=16,cex=4,col=gly.cols)
points(x=1:8,y=rep(7,8),pch=15,cex=4,col=gly.cols)
points(x=1:8,y=rep(5,8),pch=16,cex=4,col=imi.cols)
points(x=1:8,y=rep(4,8),pch=15,cex=4,col=imi.cols)
points(x=1:8,y=rep(2,8),pch=16,cex=4,col=both.cols)
points(x=1:8,y=rep(1,8),pch=15,cex=4,col=both.cols)
title(ylab='nutrient treatment',line=0.3,cex.lab=1.1)
text(x=rep(-0.4,6),y=c(8,7,5,4,2,1),labels=rep(c('low','high'),3),pos=4,cex=0.8)
text(x=1:8,y=rep(8,8),labels=as.character(gly.doses),cex=0.7)
text(x=1:8,y=rep(7,8),labels=as.character(gly.doses),cex=0.7)
text(x=1:8,y=rep(5,8),labels=as.character(imi.doses),cex=0.7)
text(x=1:8,y=rep(4,8),labels=as.character(imi.doses),cex=0.7)
text(x=1:8,y=rep(2.15,8),labels=as.character(gly.doses),cex=0.7)
text(x=1:8,y=rep(1.85,8),labels=as.character(imi.doses),cex=0.7)
text(x=1:8,y=rep(1.15,8),labels=as.character(gly.doses),cex=0.7)
text(x=1:8,y=rep(0.85,8),labels=as.character(imi.doses),cex=0.7)
text(x=rep(4.5,3),y=c(8.4,5.4,2.4),pos=3,labels=c('glyphosate only','imidacloprid only','both pesticides'),cex=1.1)
text(x=0.5,y=6.3,pos=4,labels=make.italic(expression(numbers~indicate~target~concentrations~'in'~mg/L)),cex=0.7)
text(x=0.5,y=3.3,pos=4,labels=make.italic(expression(numbers~indicate~target~concentrations~'in'~mu*g/L)),cex=0.7)
text(x=0.5,y=0.3,pos=4,labels=make.italic('top number = glyphosate (ppm), bottom number = imidacloprid (ppb)'),cex=0.7)
#dev.off()

#pdf('~/Desktop/figure1c.pdf',width=5.5,height=3.5,pointsize = 12)
plot(x=0,xlim=c(0,44),ylim=c(-1,1),type='n',ann=F,yaxt='n',xaxt='n',bty='n')
abline(h=0)
segments(x0=pulse.dates[1:2],y0=rep(-0.12,2),y1=rep(0.5,2),lty=3)
segments(x0=Sampling.dates,y0=rep(0,6),y1=rep(-0.12,6))
points(x=1:43,y=rep(0,43),pch=21,col=1,bg='white',cex=0.5)
points(x=Sampling.dates,y=rep(0,6),pch=21,col=1,bg=1)
text(x=Sampling.dates,y=rep(-0.25,6),labels=as.character(Sampling.dates),cex=0.8)
text(x=c(2,42),y=rep(-0.5,2),labels=c('(17-Aug)','(28-Sep)'),cex=0.8)
text(x=pulse.dates[1],y=0.5,label=expression(italic(dose~1)),pos=4)
text(x=pulse.dates[2],y=0.5,label=expression(italic(dose~2)),pos=4)
text(x=22,y=-0.7,label='Day of experiment')

dev.off()

#### format nutrient and pesticide data ####

nut <- left_join(nut, treat)

imi.lake <- filter(imi, site == 'LAKE')
gly.lake <- filter(gly, site == 'LAKE')

gly <- gly %>% left_join(treat, by = 'site') %>% filter(site %!in% c('E1','H1','LAKE'), date < 45)
gly$gly.measured.ppm <- gly$gly.measured.ppb/1000
gly$log.gly.measured <- log10p(gly$gly.measured.ppb)
gly$log.gly.target <- log10p(gly$gly.target.ppb)

imi <- imi %>% left_join(treat, by = 'site') %>% filter(site %!in% c('C5','D1','LAKE'), date < 45)
imi$log.imi.measured <- log10p(imi$imi.measured.ppb)
imi$log.imi.target <- log10p(imi$imi.target.ppb)

#### Figure: did treatments work? ####

pdf('~/Desktop/FigS1.pdf',width=11.5,height = 8.5,pointsize = 12)
layout(rbind(c(1,1,2,2),c(3,3,4,5),c(6,6,7,8)))
par(cex=1,mar=c(4,4,1,1))

## nutrients

tmp <- filter(nut, gly == 1, imi == 1, site != 'LAKE')
lake <- filter(nut, site == 'LAKE')

#TN

plotfunctions::emptyPlot(xlim=c(5,43),ylim=range(tmp$TN)*c(0.9,1.1),yaxt='n',xaxt='n',ann=F, bty='l', log = 'y')
axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
title(ylab=expression(paste('TN (',mu,g~L^-1,')')),line=2.8)
title(xlab="day", cex.lab=1,line=2.8)
abline(v=pulse.dates[1:2],lty=3)
for(i in 1:6){
  pond <- unique(tmp$site)[i]
  sub.sub <- filter(tmp, site == pond)
  points(x=sub.sub$date,y=sub.sub$TN,type='o',lwd=1,pch=pchs[sub.sub$nut.num[1]],col=gly.cols[1])
}
points(x=lake$date,y=lake$TN,type='o',lwd=1,pch=16,col=1)
legend(x=7,y=max(tmp$TN),bty='n',legend=c('low nutrient','high nutrient','source lake'),pch=c(1,0,16),lty=1,col=c('gray80','gray80','black'))
text(x=pulse.dates[1],y=max(tmp$TN)*1.05,label=expression(italic(dose~1)),pos=4)
text(x=pulse.dates[2],y=max(tmp$TN)*1.05,label=expression(italic(dose~2)),pos=4)

#TP

plotfunctions::emptyPlot(xlim=c(5,43),ylim=range(tmp$TP)*c(0.9,1.1),yaxt='n',xaxt='n',ann=F, bty='l', log = 'y')
axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
title(ylab=expression(paste('TP (',mu,g~L^-1,')')),line=2.8)
title(xlab="day", cex.lab=1,line=2.8)
abline(v=pulse.dates[1:2],lty=3)
for(i in 1:6){
  pond <- unique(tmp$site)[i]
  sub.sub <- filter(tmp, site == pond)
  points(x=sub.sub$date,y=sub.sub$TP,type='o',lwd=1,pch=pchs[sub.sub$nut.num[1]],col=gly.cols[1])
}
points(x=lake$date,y=lake$TP,type='o',lwd=1,pch=16,col=1)
#legend(x=7,y=max(tmp$TP),bty='n',legend=c('low nutrient','high nutrient','source lake'),pch=c(1,0,16),lty=1,col=c('gray80','gray80','black'))
text(x=pulse.dates[1],y=max(tmp$TP)*1.05,label=expression(italic(dose~1)),pos=4)
text(x=pulse.dates[2],y=max(tmp$TP)*1.05,label=expression(italic(dose~2)),pos=4)

# #molar ratio
#
# # N 14.0067 g/moles
# # P 30.97376 g/moles
# nut$NP <- (nut$TN/14.0067)/(nut$TP/30.97376)
# tmp <- filter(nut, gly == 1, imi == 1, site != 'LAKE')
# lake <- filter(nut, site == 'LAKE')
# plotfunctions::emptyPlot(xlim=c(5,43),ylim=range(tmp$NP),yaxt='n',xaxt='n',ann=F, bty='l')
# axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
# axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
# title(ylab='N:P molar ratio',line=2.8)
# title(xlab="day", cex.lab=1,line=2.8)
# abline(v=pulse.dates[1:2],lty=3)
# for(i in 1:6){
#   pond <- unique(tmp$site)[i]
#   sub.sub <- filter(tmp, site == pond)
#   points(x=sub.sub$date,y=sub.sub$NP,type='o',lwd=1,pch=pchs[sub.sub$nut.num[1]],col=gly.cols[1])
# }
# points(x=lake$date,y=lake$NP,type='o',lwd=1,pch=16,col=1)
# legend(x=8,y=max(tmp$NP),bty='n',legend=c('low nutrient','high nutrient','source lake'),pch=c(1,0,16),lty=1,col=c('gray80','gray80','black'))

## pesticides

ts.ponds.gly <- c('C1','C4','C8','D1','D4','D8','J4','J8','K4','K8') #ponds with 4-5 gly measurements over time 
ts.ponds.imi <- c('E1','E4','E8','H1','H4','H8','J4','J8','K4','K8') #ponds with 4-5 imi measurements over time

glyp1 <- filter(gly, date == 6)
glyp2 <- filter(gly, date == 34)
glyts <- filter(gly, site %in% ts.ponds.gly)
imip1 <- filter(imi, date == 6)
imip2 <- filter(imi, date == 34)
imits <- filter(imi, site %in% ts.ponds.imi)

# gly

plotfunctions::emptyPlot(xlim=c(5,43),ylim=c(-0.5,max(glyts$log.gly.measured)*1.1),yaxt='n',xaxt='n',ann=F, bty='l')
title(ylab=log[10](1+glyphosate)~(mu*g~L^-1),line=2.5)
axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
title(xlab="day", cex.lab=1,line=2.8)
abline(v=pulse.dates[1:2],lty=3)
for(i in 1:10){
  pond <- unique(glyts$site)[i]
  sub.sub <- filter(glyts, site == pond)
  points(x=sub.sub$date,y=sub.sub$log.gly.measured,type='o',lwd=1,pch=pchs[sub.sub$nut.num[1]],col=allcols[sub.sub$pond.id])
}
points(x=gly.lake$date,y=log10p(gly.lake$gly.measured.ppb),pch=16,col=1)
text(x=pulse.dates[1],y=max(glyts$log.gly.measured)*1.05,label=expression(italic(dose~1)),pos=4)
text(x=pulse.dates[2],y=max(glyts$log.gly.measured)*1.05,label=expression(italic(dose~2)),pos=4)
xseqs <- seq(37,43,length.out = 100)
yrg <- c(-0.5,max(glyts$log.gly.measured)*1.1)
ypts <- seq(from=yrg[1],to=yrg[2],length.out = 10)
points(rep(ypts[3],100)~xseqs, col=bothcolfunc(100), pch=16)
points(rep(ypts[6],100)~xseqs, col=glycolfunc(100), pch=16)
text(x=c(40,40),y=c(ypts[3],ypts[6]),labels = make.italic(c('both pesticides','glyphosate')), pos=3,cex=0.9)

plotfunctions::emptyPlot(xlim=range(glyp1$log.gly.target),ylim=range(glyp1$log.gly.measured),yaxt='n',xaxt='n',ann=F, bty='o')
title(ylab=measured~(log[10]~1+ppb),line=2.5)
title(xlab=target~(log[10]~1+ppb), cex.lab=1,line=2.8)
axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
abline(a=0,b=1,lty=2)
points(log.gly.measured~log.gly.target,glyp1,pch=pchs[glyp1$nut.num],col=allcols[glyp1$pond.id])
legend('topleft',bty='n',legend='dose 1',cex=0.95)

plotfunctions::emptyPlot(xlim=range(glyp2$log.gly.target),ylim=range(glyp2$log.gly.measured),yaxt='n',xaxt='n',ann=F, bty='o')
title(ylab=measured~(log[10]~1+ppb),line=2.5)
title(xlab=target~(log[10]~1+ppb), cex.lab=1,line=2.8)
axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
abline(a=0,b=1,lty=2)
points(log.gly.measured~log.gly.target,glyp2,pch=pchs[glyp2$nut.num],col=allcols[glyp2$pond.id])
legend('topleft',bty='n',legend='dose 2',cex=0.95)

#imi

plotfunctions::emptyPlot(xlim=c(5,43),ylim=c(-0.1,max(imits$log.imi.measured)*1.1),yaxt='n',xaxt='n',ann=F, bty='l')
title(ylab=log[10](1+imidacloprid)~(mu*g~L^-1),line=2.5)
axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
title(xlab="day", cex.lab=1,line=2.8)
abline(v=pulse.dates[1:2],lty=3)
for(i in 1:10){
  pond <- unique(imits$site)[i]
  sub.sub <- filter(imits, site == pond)
  points(x=sub.sub$date,y=sub.sub$log.imi.measured,type='o',lwd=1,pch=pchs[sub.sub$nut.num[1]],col=allcols[sub.sub$pond.id])
}
points(x=imi.lake$date,y=log10p(imi.lake$imi.measured.ppb),pch=16,col=1)
text(x=pulse.dates[1],y=max(imits$log.imi.measured)*1.05,label=expression(italic(dose~1)),pos=4)
text(x=pulse.dates[2],y=max(imits$log.imi.measured)*1.05,label=expression(italic(dose~2)),pos=4)
xseqs <- seq(37,43,length.out = 100)
yrg <- c(-0.1,max(imits$log.imi.measured)*1.1)
ypts <- seq(from=yrg[1],to=yrg[2],length.out = 10)
points(rep(ypts[3],100)~xseqs, col=bothcolfunc(100), pch=16)
points(rep(ypts[6],100)~xseqs, col=imicolfunc(100), pch=16)
text(x=c(40,40),y=c(ypts[3],ypts[6]),labels = make.italic(c('both pesticides','imidacloprid')), pos=3,cex=0.9)

plotfunctions::emptyPlot(xlim=range(imip1$log.imi.target),ylim=range(imip1$log.imi.measured),yaxt='n',xaxt='n',ann=F, bty='o')
title(ylab=measured~(log[10]~1+ppb),line=2.5)
title(xlab=target~(log[10]~1+ppb), cex.lab=1,line=2.8)
axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
abline(a=0,b=1,lty=2)
points(log.imi.measured~log.imi.target,imip1,pch=pchs[imip1$nut.num],col=allcols[imip1$pond.id])
legend('topleft',bty='n',legend='dose 1',cex=0.95)

plotfunctions::emptyPlot(xlim=range(imip2$log.imi.target),ylim=range(imip2$log.imi.measured),yaxt='n',xaxt='n',ann=F, bty='o')
title(ylab=measured~(log[10]~1+ppb),line=2.5)
title(xlab=target~(log[10]~1+ppb), cex.lab=1,line=2.8)
axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
abline(a=0,b=1,lty=2)
points(log.imi.measured~log.imi.target,imip2,pch=pchs[imip2$nut.num],col=allcols[imip2$pond.id])
legend('topleft',bty='n',legend='dose 2',cex=0.95)

dev.off()

#### supplementary figure: other physico-chem variables ####

pdf('~/Desktop/FigS2.pdf',width=7,height = 5,pointsize = 8)
par(mfrow=c(4,3),mar=c(2,2,1,1),oma=c(2.5,2.5,1,1),cex=1)
pchs.2 <- c(21,22)
vars <- c('pH.mean','DO.mean','SPC.mean')
var.names <- c(expression(pH),expression(DO~'(mg/L)'),expression(SPC~(mu*g/cm)))

for(v in 1:3){
  tmp <- merged.data[,c(colnames(merged.data)[1:13],vars[v])]
  tmp$date.idx <- tmp$date - 0.7
  tmp$date.idx[tmp$nut.num == 2] <- tmp$date[tmp$nut.num == 2] + 0.7
  for(letter in c('C|D','E|H','J|K')){
    plotfunctions::emptyPlot(xlim=c(1,43),ylim=range(tmp[,14]),yaxt='n',xaxt='n',ann=F, bty='l')
    axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
    axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
    abline(v=pulse.dates[1:2],lty=3)
    tmp.sub <- filter(tmp, str_detect(site, letter))
    for(i in 1:n_distinct(tmp.sub$site)){
      pond <- unique(tmp.sub$site)[i]
      sub.sub <- filter(tmp.sub, site == pond)
      points(x=sub.sub$date.idx,y=sub.sub[,14],type='l',lwd=1.2, pch=pchs.2[sub.sub$nut.num], col=alpha(allcols[sub.sub$pond.id],1))
    }
    for(i in 1:n_distinct(tmp.sub$site)){
      pond <- unique(tmp.sub$site)[i]
      sub.sub <- filter(tmp.sub, site == pond)
      points(x=sub.sub$date.idx,y=sub.sub[,14],type='p',pch=pchs.2[sub.sub$nut.num], col=1, bg=alpha(allcols[sub.sub$pond.id],1))
    }
    if(letter == 'C|D'){mtext(var.names[v],side=2,line=3)}
  }
}

xmax <- 2064
hobodat <- hobodat %>% select(-date, -time) %>% gather('site', 'temp', -time.num)
hobodat <- left_join(hobodat, treat, by = 'site')
for(letter in c('C|D','E|H','J|K')){
  plotfunctions::emptyPlot(xlim=c(1,xmax),ylim=range(hobodat$temp),yaxt='n',xaxt='n',ann=F, bty='l')
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
  hobodat.sub <- filter(hobodat, str_detect(site, letter))
  for(i in 1:n_distinct(hobodat.sub$site)){
    pond <- unique(hobodat.sub$site)[i]
    sub.sub <- filter(hobodat.sub, site == pond)
    points(x=sub.sub$time.num,y=sub.sub$temp,type='l',lwd=1, col=alpha(allcols[sub.sub$pond.id],0.9))
  }
  abline(v=(pulse.dates[1:2])*48,lty=3)
  if(letter == 'C|D'){mtext("temperature (°C)",side=2,line=3)}
}

mtext('day of experiment',side=1,outer=T,line=1)

dev.off()

##### format variables for models ####

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

#### Fig. 2 & 3: time series of response variables ####

vars <- c('BA','total','total_zoo','use','NEP','RUE')

var.names <- c(expression(log[10]~bact.~(cells/mu*L)),
               expression(log[10]~chl.~italic(a)~(mu*g/L)),
               expression(log[10](1+zoo)~(mu*g/L)),
               expression(C~sources~used),
               expression(log[10]~(1+Delta*DO)~(mu*g/L)),
               expression(log[10]~(Z:P)~(mu*g/mu*g)))

plot.data <- select(merged.data, date:date.f, vars)
plot.data <- plot.data %>% mutate_at(vars(BA,total,RUE), log10)
plot.data <- plot.data %>% mutate_at(vars(total_zoo,NEP), log10p)

pdf('~/Desktop/Fig2.pdf',width=9,height = 6,pointsize = 10)
par(mfrow=c(3,3),mar=c(2,2,1,1),oma=c(2.5,2.5,1,1),cex=1)

for(v in 1:3){
  
  tmp <- plot.data[,c(colnames(plot.data)[1:13],vars[v])]
  tmp <- drop_na(tmp)
  tmp <- filter(tmp, is.finite(tmp[,14]))
  
  tmp$date.idx <- tmp$date - 0.7
  tmp$date.idx[tmp$nut.num == 2] <- tmp$date[tmp$nut.num == 2] + 0.7
  for(letter in c('C|D','E|H','J|K')){
    plotfunctions::emptyPlot(xlim=c(1,43),ylim=range(tmp[,14]),yaxt='n',xaxt='n',ann=F, bty='l')
    axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
    axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
    tmp.sub <- filter(tmp, str_detect(site, letter))
    abline(v=pulse.dates[1:2],lty=3)
    for(i in 1:n_distinct(tmp.sub$site)){
      pond <- unique(tmp.sub$site)[i]
      sub.sub <- filter(tmp.sub, site == pond)
      points(x=sub.sub$date.idx,y=sub.sub[,14],type='l',lwd=1.2, pch=pchs.2[sub.sub$nut.num], col=alpha(allcols[sub.sub$pond.id],1))
    }
    for(i in 1:n_distinct(tmp.sub$site)){
      pond <- unique(tmp.sub$site)[i]
      sub.sub <- filter(tmp.sub, site == pond)
      points(x=sub.sub$date.idx,y=sub.sub[,14],type='p',pch=pchs.2[sub.sub$nut.num], col=1, bg=alpha(allcols[sub.sub$pond.id],1))
    }
    if(letter == 'C|D'){mtext(var.names[v],side=2,line=3,cex=1.1)}
  }
  if(v == 3){mtext('day of experiment',side=1,outer=T,line=1,cex=1.1)}
}

dev.off()

pdf('~/Desktop/Fig3.pdf',width=9,height = 6,pointsize = 10)
par(mfrow=c(3,3),mar=c(2,2,1,1),oma=c(2.5,2.5,1,1),cex=1)

for(v in 4:6){
  
  tmp <- plot.data[,c(colnames(plot.data)[1:13],vars[v])]
  tmp <- drop_na(tmp)
  tmp <- filter(tmp, is.finite(tmp[,14]))
  
  tmp$date.idx <- tmp$date - 0.7
  tmp$date.idx[tmp$nut.num == 2] <- tmp$date[tmp$nut.num == 2] + 0.7
  for(letter in c('C|D','E|H','J|K')){
    plotfunctions::emptyPlot(xlim=c(1,43),ylim=range(tmp[,14]),yaxt='n',xaxt='n',ann=F, bty='l')
    axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
    axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
    tmp.sub <- filter(tmp, str_detect(site, letter))
    abline(v=pulse.dates[1:2],lty=3)
    for(i in 1:n_distinct(tmp.sub$site)){
      pond <- unique(tmp.sub$site)[i]
      sub.sub <- filter(tmp.sub, site == pond)
      points(x=sub.sub$date.idx,y=sub.sub[,14],type='l',lwd=1.2, pch=pchs.2[sub.sub$nut.num], col=alpha(allcols[sub.sub$pond.id],1))
    }
    for(i in 1:n_distinct(tmp.sub$site)){
      pond <- unique(tmp.sub$site)[i]
      sub.sub <- filter(tmp.sub, site == pond)
      points(x=sub.sub$date.idx,y=sub.sub[,14],type='p',pch=pchs.2[sub.sub$nut.num], col=1, bg=alpha(allcols[sub.sub$pond.id],1))
    }
    if(letter == 'C|D'){mtext(var.names[v],side=2,line=3,cex=1.1)}
  }
  if(v == 6){mtext('day of experiment',side=1,outer=T,line=1,cex=1.1)}
}

dev.off()

#### Linear models ####

merged.data$nut <- merged.data$nut-1

ba.mod.lin <- lmer(log10(BA) ~ date.f:nut + date.f:sc.gly + date.f:sc.imi + date.f:sc.gly:sc.imi + date.f:sc.gly:nut +
                     date.f:sc.imi:nut + (1|site.f),merged.data)
ba.coefs <- get_model_data(ba.mod.lin, type = 'est')

chla.mod.lin <- lmer(log10(total) ~ date.f:nut + date.f:sc.gly + date.f:sc.imi + date.f:sc.gly:sc.imi + date.f:sc.gly:nut +
                       date.f:sc.imi:nut + (1|site.f),merged.data)
chla.coefs <- get_model_data(chla.mod.lin, type = 'est')

zoo.mod.lin <- lmer(log10p(total_zoo) ~ date.f:nut + date.f:sc.gly + date.f:sc.imi + date.f:sc.gly:sc.imi + date.f:sc.gly:nut +
                      date.f:sc.imi:nut + (1|site.f),merged.data)
zoo.coefs <- get_model_data(zoo.mod.lin, type = 'est')

EP.mod.lin <- lmer(use ~ date.f:nut + date.f:sc.gly + date.f:sc.imi + date.f:sc.gly:sc.imi + date.f:sc.gly:nut +
                     date.f:sc.imi:nut + (1|site.f),merged.data)
use.coefs <- get_model_data(EP.mod.lin, type = 'est')

nep.mod.lin <- lmer(log10p(NEP) ~ date.f:nut + date.f:sc.gly + date.f:sc.imi + date.f:sc.gly:sc.imi + date.f:sc.gly:nut +
                      date.f:sc.imi:nut + (1|site.f),merged.data)
nep.coefs <- get_model_data(nep.mod.lin, type = 'est')

rue.mod.lin <- lmer(log10(RUE) ~ date.f:nut + date.f:sc.gly + date.f:sc.imi + date.f:sc.gly:sc.imi + date.f:sc.gly:nut +
                      date.f:sc.imi:nut + (1|site.f),merged.data,subset=is.finite(log(RUE)))
rue.coefs <- get_model_data(rue.mod.lin, type = 'est')

tab_model(ba.mod.lin,chla.mod.lin,zoo.mod.lin, dv.labels = c('bacterial abundance','chlorophyll a','zooplankton biomass'), file='~/Desktop/TableS1.html')
tab_model(EP.mod.lin,nep.mod.lin,rue.mod.lin, dv.labels = c('C substrates used','oxygen production','zooplankton:phytoplankton'), file='~/Desktop/TableS2.html')

#### Figure 4: coef plot ####

pdf('~/Desktop/Fig4.pdf',width=6,height = 6,pointsize = 12)

var.names <- c('Bact.','Chl. a','Zoo.','C use',expression(Delta*DO),'Z:P')

par(mfrow=c(6,6),mar=c(0.1,0.1,0.1,0.1),oma=c(4,4.5,3,0.5),cex=1)

#BA#

ba.coefs$alpha <- 0.2
ba.coefs$alpha[str_detect(ba.coefs$p.stars, '\\*')] <- 1
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(ba.coefs$conf.low,ba.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col = alpha(gly.cols[1],0.4), border = NA)
abline(h=0,col=1);abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=ba.coefs$estimate[1:6],pch=16,col=scales::alpha('black',ba.coefs$alpha[1:6]))
arrows(x0=1:6,y0=ba.coefs$conf.low[1:6],y1=ba.coefs$conf.high[1:6],length=0,col=scales::alpha('black',ba.coefs$alpha[1:6]))
mtext(var.names[1],side=2,outer=F,line=1,cex=1)
axis(2,lwd=0,lwd.ticks = 1,at=0,labels='')
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(ba.coefs$conf.low,ba.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col = alpha(gly.cols[4],0.2), border = NA)
abline(h=0,col=1);abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=ba.coefs$estimate[7:12],pch=16,col=scales::alpha('black',ba.coefs$alpha[7:12]))
arrows(x0=1:6,y0=ba.coefs$conf.low[7:12],y1=ba.coefs$conf.high[7:12],length=0,col=scales::alpha('black',ba.coefs$alpha[7:12]))
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(ba.coefs$conf.low,ba.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col = alpha(gly.cols[4],0.2), border = NA)
abline(h=0,col=1);abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=ba.coefs$estimate[25:30],pch=16,col=scales::alpha('black',ba.coefs$alpha[25:30]))
arrows(x0=1:6,y0=ba.coefs$conf.low[25:30],y1=ba.coefs$conf.high[25:30],length=0,col=scales::alpha('black',ba.coefs$alpha[25:30]))
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(ba.coefs$conf.low,ba.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col = alpha(imi.cols[6],0.2), border = NA)
abline(h=0,col=1);abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=ba.coefs$estimate[13:18],pch=16,col=scales::alpha('black',ba.coefs$alpha[13:18]))
arrows(x0=1:6,y0=ba.coefs$conf.low[13:18],y1=ba.coefs$conf.high[13:18],length=0,col=scales::alpha('black',ba.coefs$alpha[13:18]))
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(ba.coefs$conf.low,ba.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col = alpha(imi.cols[6],0.2), border = NA)
abline(h=0,col=1);abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=ba.coefs$estimate[31:36],pch=16,col=scales::alpha('black',ba.coefs$alpha[31:36]))
arrows(x0=1:6,y0=ba.coefs$conf.low[31:36],y1=ba.coefs$conf.high[31:36],length=0,col=scales::alpha('black',ba.coefs$alpha[31:36]))
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(ba.coefs$conf.low,ba.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col = alpha(both.cols[6],0.2), border = NA)
abline(h=0,col=1);abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=ba.coefs$estimate[19:24],pch=16,col=scales::alpha('black',ba.coefs$alpha[19:24]))
arrows(x0=1:6,y0=ba.coefs$conf.low[19:24],y1=ba.coefs$conf.high[19:24],length=0,col=scales::alpha('black',ba.coefs$alpha[19:24]))

mtext('independent variable',side=3,outer=T,line=1.7,cex=1.3)
mtext(c('nut','gly','gly:nut','imi','imi:nut','gly:imi'),side=3,outer=T,line=0.2,at=seq(0.09,0.91,length.out = 6),adj=0.5,cex=1)
#mtext(c('nut.','gly.','gly. × nut.','imi.','imi. × nut.','gly. × imi.'),side=3,outer=T,line=0.2,at=seq(0.09,0.91,length.out = 6),adj=0.5,cex=1)
mtext('time point',side=1,outer=T,line=2.5,cex=1.3)
mtext('response variable',side=2,outer=T,line=3,cex=1.3)

#Chla#

chla.coefs$alpha <- 0.2
chla.coefs$alpha[str_detect(chla.coefs$p.stars, '\\*')] <- 1
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(chla.coefs$conf.low,chla.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col = alpha(gly.cols[1],0.4), border = NA)
abline(h=0,col=1);abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=chla.coefs$estimate[1:6],pch=16,col=scales::alpha('black',chla.coefs$alpha[1:6]))
arrows(x0=1:6,y0=chla.coefs$conf.low[1:6],y1=chla.coefs$conf.high[1:6],length=0,col=scales::alpha('black',chla.coefs$alpha[1:6]))
mtext(var.names[2],side=2,outer=F,line=1,cex=1)
axis(2,lwd=0,lwd.ticks = 1,at=0,labels='')
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(chla.coefs$conf.low,chla.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col = alpha(gly.cols[4],0.2), border = NA)
abline(h=0,col=1);abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=chla.coefs$estimate[7:12],pch=16,col=scales::alpha('black',chla.coefs$alpha[7:12]))
arrows(x0=1:6,y0=chla.coefs$conf.low[7:12],y1=chla.coefs$conf.high[7:12],length=0,col=scales::alpha('black',chla.coefs$alpha[7:12]))
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(chla.coefs$conf.low,chla.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col = alpha(gly.cols[4],0.2), border = NA)
abline(h=0,col=1);abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=chla.coefs$estimate[25:30],pch=16,col=scales::alpha('black',chla.coefs$alpha[25:30]))
arrows(x0=1:6,y0=chla.coefs$conf.low[25:30],y1=chla.coefs$conf.high[25:30],length=0,col=scales::alpha('black',chla.coefs$alpha[25:30]))
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(chla.coefs$conf.low,chla.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col = alpha(imi.cols[6],0.2), border = NA)
abline(h=0,col=1);abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=chla.coefs$estimate[13:18],pch=16,col=scales::alpha('black',chla.coefs$alpha[13:18]))
arrows(x0=1:6,y0=chla.coefs$conf.low[13:18],y1=chla.coefs$conf.high[13:18],length=0,col=scales::alpha('black',chla.coefs$alpha[13:18]))
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(chla.coefs$conf.low,chla.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col = alpha(imi.cols[6],0.2), border = NA)
abline(h=0,col=1);abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=chla.coefs$estimate[31:36],pch=16,col=scales::alpha('black',chla.coefs$alpha[31:36]))
arrows(x0=1:6,y0=chla.coefs$conf.low[31:36],y1=chla.coefs$conf.high[31:36],length=0,col=scales::alpha('black',chla.coefs$alpha[31:36]))
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(chla.coefs$conf.low,chla.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col = alpha(both.cols[6],0.2), border = NA)
abline(h=0,col=1);abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=chla.coefs$estimate[19:24],pch=16,col=scales::alpha('black',chla.coefs$alpha[19:24]))
arrows(x0=1:6,y0=chla.coefs$conf.low[19:24],y1=chla.coefs$conf.high[19:24],length=0,col=scales::alpha('black',chla.coefs$alpha[19:24]))

#zoo#

zoo.coefs$alpha <- 0.2
zoo.coefs$alpha[str_detect(zoo.coefs$p.stars, '\\*')] <- 1
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(zoo.coefs$conf.low,zoo.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col = alpha(gly.cols[1],0.4), border = NA)
abline(h=0,col=1);abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=zoo.coefs$estimate[1:6],pch=16,col=scales::alpha('black',zoo.coefs$alpha[1:6]))
arrows(x0=1:6,y0=zoo.coefs$conf.low[1:6],y1=zoo.coefs$conf.high[1:6],length=0,col=scales::alpha('black',zoo.coefs$alpha[1:6]))
mtext(var.names[3],side=2,outer=F,line=1,cex=1)
axis(2,lwd=0,lwd.ticks = 1,at=0,labels='')
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(zoo.coefs$conf.low,zoo.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col = alpha(gly.cols[4],0.2), border = NA)
abline(h=0,col=1);abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=zoo.coefs$estimate[7:12],pch=16,col=scales::alpha('black',zoo.coefs$alpha[7:12]))
arrows(x0=1:6,y0=zoo.coefs$conf.low[7:12],y1=zoo.coefs$conf.high[7:12],length=0,col=scales::alpha('black',zoo.coefs$alpha[7:12]))
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(zoo.coefs$conf.low,zoo.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col = alpha(gly.cols[4],0.2), border = NA)
abline(h=0,col=1);abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=zoo.coefs$estimate[25:30],pch=16,col=scales::alpha('black',zoo.coefs$alpha[25:30]))
arrows(x0=1:6,y0=zoo.coefs$conf.low[25:30],y1=zoo.coefs$conf.high[25:30],length=0,col=scales::alpha('black',zoo.coefs$alpha[25:30]))
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(zoo.coefs$conf.low,zoo.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col = alpha(imi.cols[6],0.2), border = NA)
abline(h=0,col=1);abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=zoo.coefs$estimate[13:18],pch=16,col=scales::alpha('black',zoo.coefs$alpha[13:18]))
arrows(x0=1:6,y0=zoo.coefs$conf.low[13:18],y1=zoo.coefs$conf.high[13:18],length=0,col=scales::alpha('black',zoo.coefs$alpha[13:18]))
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(zoo.coefs$conf.low,zoo.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col = alpha(imi.cols[6],0.2), border = NA)
abline(h=0,col=1);abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=zoo.coefs$estimate[31:36],pch=16,col=scales::alpha('black',zoo.coefs$alpha[31:36]))
arrows(x0=1:6,y0=zoo.coefs$conf.low[31:36],y1=zoo.coefs$conf.high[31:36],length=0,col=scales::alpha('black',zoo.coefs$alpha[31:36]))
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(zoo.coefs$conf.low,zoo.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col = alpha(both.cols[6],0.2), border = NA)
abline(h=0,col=1);abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=zoo.coefs$estimate[19:24],pch=16,col=scales::alpha('black',zoo.coefs$alpha[19:24]))
arrows(x0=1:6,y0=zoo.coefs$conf.low[19:24],y1=zoo.coefs$conf.high[19:24],length=0,col=scales::alpha('black',zoo.coefs$alpha[19:24]))

#use#

use.coefs$alpha <- 0.2
use.coefs$alpha[str_detect(use.coefs$p.stars, '\\*')] <- 1
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(use.coefs$conf.low,use.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col = alpha(gly.cols[1],0.4), border = NA)
abline(h=0,col=1);abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=use.coefs$estimate[1:6],pch=16,col=scales::alpha('black',use.coefs$alpha[1:6]))
arrows(x0=1:6,y0=use.coefs$conf.low[1:6],y1=use.coefs$conf.high[1:6],length=0,col=scales::alpha('black',use.coefs$alpha[1:6]))
mtext(var.names[4],side=2,outer=F,line=1,cex=1)
axis(2,lwd=0,lwd.ticks = 1,at=0,labels='')
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(use.coefs$conf.low,use.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col = alpha(gly.cols[4],0.2), border = NA)
abline(h=0,col=1);abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=use.coefs$estimate[7:12],pch=16,col=scales::alpha('black',use.coefs$alpha[7:12]))
arrows(x0=1:6,y0=use.coefs$conf.low[7:12],y1=use.coefs$conf.high[7:12],length=0,col=scales::alpha('black',use.coefs$alpha[7:12]))
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(use.coefs$conf.low,use.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col = alpha(gly.cols[4],0.2), border = NA)
abline(h=0,col=1);abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=use.coefs$estimate[25:30],pch=16,col=scales::alpha('black',use.coefs$alpha[25:30]))
arrows(x0=1:6,y0=use.coefs$conf.low[25:30],y1=use.coefs$conf.high[25:30],length=0,col=scales::alpha('black',use.coefs$alpha[25:30]))
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(use.coefs$conf.low,use.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col = alpha(imi.cols[6],0.2), border = NA)
abline(h=0,col=1);abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=use.coefs$estimate[13:18],pch=16,col=scales::alpha('black',use.coefs$alpha[13:18]))
arrows(x0=1:6,y0=use.coefs$conf.low[13:18],y1=use.coefs$conf.high[13:18],length=0,col=scales::alpha('black',use.coefs$alpha[13:18]))
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(use.coefs$conf.low,use.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col = alpha(imi.cols[6],0.2), border = NA)
abline(h=0,col=1);abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=use.coefs$estimate[31:36],pch=16,col=scales::alpha('black',use.coefs$alpha[31:36]))
arrows(x0=1:6,y0=use.coefs$conf.low[31:36],y1=use.coefs$conf.high[31:36],length=0,col=scales::alpha('black',use.coefs$alpha[31:36]))
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(use.coefs$conf.low,use.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col = alpha(both.cols[6],0.2), border = NA)
abline(h=0,col=1);abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=use.coefs$estimate[19:24],pch=16,col=scales::alpha('black',use.coefs$alpha[19:24]))
arrows(x0=1:6,y0=use.coefs$conf.low[19:24],y1=use.coefs$conf.high[19:24],length=0,col=scales::alpha('black',use.coefs$alpha[19:24]))

#nep#

nep.coefs$alpha <- 0.2
nep.coefs$alpha[str_detect(nep.coefs$p.stars, '\\*')] <- 1
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(nep.coefs$conf.low,nep.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col = alpha(gly.cols[1],0.4), border = NA)
abline(h=0,col=1);abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=nep.coefs$estimate[1:6],pch=16,col=scales::alpha('black',nep.coefs$alpha[1:6]))
arrows(x0=1:6,y0=nep.coefs$conf.low[1:6],y1=nep.coefs$conf.high[1:6],length=0,col=scales::alpha('black',nep.coefs$alpha[1:6]))
mtext(var.names[5],side=2,outer=F,line=1,cex=1)
axis(2,lwd=0,lwd.ticks = 1,at=0,labels='')
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(nep.coefs$conf.low,nep.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col = alpha(gly.cols[4],0.2), border = NA)
abline(h=0,col=1);abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=nep.coefs$estimate[7:12],pch=16,col=scales::alpha('black',nep.coefs$alpha[7:12]))
arrows(x0=1:6,y0=nep.coefs$conf.low[7:12],y1=nep.coefs$conf.high[7:12],length=0,col=scales::alpha('black',nep.coefs$alpha[7:12]))
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(nep.coefs$conf.low,nep.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col = alpha(gly.cols[4],0.2), border = NA)
abline(h=0,col=1);abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=nep.coefs$estimate[25:30],pch=16,col=scales::alpha('black',nep.coefs$alpha[25:30]))
arrows(x0=1:6,y0=nep.coefs$conf.low[25:30],y1=nep.coefs$conf.high[25:30],length=0,col=scales::alpha('black',nep.coefs$alpha[25:30]))
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(nep.coefs$conf.low,nep.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col = alpha(imi.cols[6],0.2), border = NA)
abline(h=0,col=1);abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=nep.coefs$estimate[13:18],pch=16,col=scales::alpha('black',nep.coefs$alpha[13:18]))
arrows(x0=1:6,y0=nep.coefs$conf.low[13:18],y1=nep.coefs$conf.high[13:18],length=0,col=scales::alpha('black',nep.coefs$alpha[13:18]))
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(nep.coefs$conf.low,nep.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col = alpha(imi.cols[6],0.2), border = NA)
abline(h=0,col=1);abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=nep.coefs$estimate[31:36],pch=16,col=scales::alpha('black',nep.coefs$alpha[31:36]))
arrows(x0=1:6,y0=nep.coefs$conf.low[31:36],y1=nep.coefs$conf.high[31:36],length=0,col=scales::alpha('black',nep.coefs$alpha[31:36]))
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(nep.coefs$conf.low,nep.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col = alpha(both.cols[6],0.2), border = NA)
abline(h=0,col=1);abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=nep.coefs$estimate[19:24],pch=16,col=scales::alpha('black',nep.coefs$alpha[19:24]))
arrows(x0=1:6,y0=nep.coefs$conf.low[19:24],y1=nep.coefs$conf.high[19:24],length=0,col=scales::alpha('black',nep.coefs$alpha[19:24]))

#rue#

rue.coefs$alpha <- 0.2
rue.coefs$alpha[str_detect(rue.coefs$p.stars, '\\*')] <- 1
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(rue.coefs$conf.low,rue.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
axis(1,lwd=0,lwd.ticks = 1)
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col = alpha(gly.cols[1],0.4), border = NA)
abline(h=0,col=1);abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=rue.coefs$estimate[1:6],pch=16,col=scales::alpha('black',rue.coefs$alpha[1:6]))
arrows(x0=1:6,y0=rue.coefs$conf.low[1:6],y1=rue.coefs$conf.high[1:6],length=0,col=scales::alpha('black',rue.coefs$alpha[1:6]))
mtext(var.names[6],side=2,outer=F,line=1,cex=1)
axis(2,lwd=0,lwd.ticks = 1,at=0,labels='')
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(rue.coefs$conf.low,rue.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
axis(1,lwd=0,lwd.ticks = 1)
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col = alpha(gly.cols[4],0.2), border = NA)
abline(h=0,col=1);abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=rue.coefs$estimate[7:12],pch=16,col=scales::alpha('black',rue.coefs$alpha[7:12]))
arrows(x0=1:6,y0=rue.coefs$conf.low[7:12],y1=rue.coefs$conf.high[7:12],length=0,col=scales::alpha('black',rue.coefs$alpha[7:12]))
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(rue.coefs$conf.low,rue.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
axis(1,lwd=0,lwd.ticks = 1)
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col = alpha(gly.cols[4],0.2), border = NA)
abline(h=0,col=1);abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=rue.coefs$estimate[25:30],pch=16,col=scales::alpha('black',rue.coefs$alpha[25:30]))
arrows(x0=1:6,y0=rue.coefs$conf.low[25:30],y1=rue.coefs$conf.high[25:30],length=0,col=scales::alpha('black',rue.coefs$alpha[25:30]))
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(rue.coefs$conf.low,rue.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
axis(1,lwd=0,lwd.ticks = 1)
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col = alpha(imi.cols[6],0.2), border = NA)
abline(h=0,col=1);abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=rue.coefs$estimate[13:18],pch=16,col=scales::alpha('black',rue.coefs$alpha[13:18]))
arrows(x0=1:6,y0=rue.coefs$conf.low[13:18],y1=rue.coefs$conf.high[13:18],length=0,col=scales::alpha('black',rue.coefs$alpha[13:18]))
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(rue.coefs$conf.low,rue.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
axis(1,lwd=0,lwd.ticks = 1)
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col = alpha(imi.cols[6],0.2), border = NA)
abline(h=0,col=1);abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=rue.coefs$estimate[31:36],pch=16,col=scales::alpha('black',rue.coefs$alpha[31:36]))
arrows(x0=1:6,y0=rue.coefs$conf.low[31:36],y1=rue.coefs$conf.high[31:36],length=0,col=scales::alpha('black',rue.coefs$alpha[31:36]))
plot(x=NULL,xlim=c(0.5,6.5),ylim=range(c(rue.coefs$conf.low,rue.coefs$conf.high)),bty='o',yaxt='n',xaxt='n')
axis(1,lwd=0,lwd.ticks = 1)
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col = alpha(both.cols[6],0.2), border = NA)
abline(h=0,col=1);abline(v=c(1.5,4.5),lty=3)
points(x=1:6,y=rue.coefs$estimate[19:24],pch=16,col=scales::alpha('black',rue.coefs$alpha[19:24]))
arrows(x0=1:6,y0=rue.coefs$conf.low[19:24],y1=rue.coefs$conf.high[19:24],length=0,col=scales::alpha('black',rue.coefs$alpha[19:24]))

dev.off()

#### fitting GAMMs ####

# ba.mo <- gam(log10(BA) ~ o.nut + ti(date,k=5) + ti(date,sc.gly, k=6) + ti(date,sc.gly, by = o.nut, k=3) + ti(date,sc.imi, k=6) + ti(date,sc.imi, by = o.nut, k=3) + ti(date,sc.gly,sc.imi, k=5) + ti(date,sc.gly,sc.imi, by = o.nut, k=3) + s(date, site.f, bs='fs',k=3), data=merged.data, method = 'REML')
# chla.m <- gam(log10(total) ~ o.nut + ti(date,k=5) + ti(date,sc.gly, k=6) + ti(date,sc.gly, by = o.nut, k=3) + ti(date,sc.imi, k=6) + ti(date,sc.imi, by = o.nut, k=3) + ti(date,sc.gly,sc.imi, k=5) + ti(date,sc.gly,sc.imi, by = o.nut, k=3) + s(date, site.f, bs='fs',k=3), data=merged.data, method = 'REML')
# zoo.m <- gam(log10p(total_zoo) ~ o.nut + ti(date,k=5) + ti(date,sc.gly, k=6) + ti(date,sc.gly, by = o.nut, k=3) + ti(date,sc.imi, k=6) + ti(date,sc.imi, by = o.nut, k=3) + ti(date,sc.gly,sc.imi, k=5) + ti(date,sc.gly,sc.imi, by = o.nut, k=3) + s(date, site.f, bs='fs',k=3), data=merged.data, method = 'REML')
# 
# use.m <- gam(use ~ o.nut + ti(date,k=5) + ti(date,sc.gly, k=6) + ti(date,sc.gly, by = o.nut, k=3) + ti(date,sc.imi, k=6) + ti(date,sc.imi, by = o.nut, k=3) + ti(date,sc.gly,sc.imi, k=5) + ti(date,sc.gly,sc.imi, by = o.nut, k=3) + s(date, site.f, bs='fs',k=3), data=merged.data, method = 'REML')
# nep.m <- gam(log10p(NEP) ~ o.nut + ti(date,k=5) + ti(date,sc.gly, k=6) + ti(date,sc.gly, by = o.nut, k=3) + ti(date,sc.imi, k=6) + ti(date,sc.imi, by = o.nut, k=3) + ti(date,sc.gly,sc.imi, k=5) + ti(date,sc.gly,sc.imi, by = o.nut, k=3) + s(date, site.f, bs='fs',k=3), data=merged.data, method = 'REML')
# rue.m <- gam(log10(RUE) ~ o.nut + ti(date,k=5) + ti(date,sc.gly, k=6) + ti(date,sc.gly, by = o.nut, k=3) + ti(date,sc.imi, k=6) + ti(date,sc.imi, by = o.nut, k=3) + ti(date,sc.gly,sc.imi, k=5) + ti(date,sc.gly,sc.imi, by = o.nut, k=3) + s(date, site.f, bs='fs',k=3), data=merged.data, method = 'REML')

load('~/Google Drive/Recherche/LEAP Postdoc/2016/GAMMs.RData')

#### Figure 5: GAMM summary ####

sdcdcs

gamtabs(ba.m)


#### Figure S3-S4: contour plots ###

fvdfv

#### Figures S5-S10: scatterplots with GAMM fit ####

pdf('~/Desktop/scatterplots.pdf',width=5.5,height = 3,pointsize = 8,onefile = T)
par(mfrow=c(3,6),mar=c(0.1,0.1,0.1,0.1),oma=c(4,4,2,0.5),cex=1,xpd=T)
pchs.2 <- c(16,15)

for(v in 1:length(vars)){
  tmp <- plot.data[,c(colnames(plot.data)[1:13],vars[v])]
  tmp <- drop_na(tmp)
  tmp <- filter(tmp, is.finite(tmp[,14]))
  ylims <- range(tmp[,14])
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
  mtext(var.names[v],side=2,outer=T,line=2.5,cex=1.2)
  mtext('nominal pesticide concentration (dose 1 to 8)',side=1,outer=T,line=2.5,cex=1.2)
  mtext(paste('day',Sampling.dates,' '),side=3,outer=T,line=0.1,at=seq(0.1,0.93,length.out = 6),adj=0.5)
}

dev.off()


scattergam <- function(var, varname, model.name, xticks=F){
  tmp <- plot.data[,c(colnames(plot.data)[1:13],var)]
  tmp <- drop_na(tmp)
  tmp <- filter(tmp, is.finite(tmp[,14]))
  ylims <- range(tmp[,14])
  for(date.x in Sampling.dates){
    sub <- tmp %>% filter(date == date.x)
    sub$pesticide <- rescale(as.numeric(str_remove(sub$site, 'C|D|E|H|J|K')),c(0,1))
    if(date.x == 1 & xticks==T){
      plot(y=sub[,14],x=sub[,15],bty='o',xlim=c(-0.05,1.05),type='p',ylim=ylims,pch=pchs[sub$nut.num],col=allcols[sub$pond.id],xlab=NULL,ylab=NULL, xaxt='n')
      mtext(varname,side=2,outer=F,line=2.5,cex=1.1)
      axis(1,at=rescale(1:8,c(0,1)),labels=1:8,lwd=0,lwd.ticks = 1)
    }else if(date.x == 1 & xticks==F){
      plot(y=sub[,14],x=sub[,15],bty='o',xlim=c(-0.05,1.05),type='p',ylim=ylims,pch=pchs[sub$nut.num],col=allcols[sub$pond.id],xlab=NULL,ylab=NULL, xaxt='n')
      mtext(varname,side=2,outer=F,line=2.5,cex=1.1)
    }else if(date.x != 1 & xticks==T){
      plot(y=sub[,14],x=sub[,15],bty='o',xlim=c(-0.05,1.05),type='p',ylim=ylims,pch=pchs[sub$nut.num],col=allcols[sub$pond.id],xlab=NULL,ylab=NULL,yaxt='n', xaxt='n')
      axis(1,at=rescale(1:8,c(0,1)),labels=1:8,lwd=0,lwd.ticks = 1)
    }else{
      plot(y=sub[,14],x=sub[,15],bty='o',xlim=c(-0.05,1.05),type='p',ylim=ylims,pch=pchs[sub$nut.num],col=allcols[sub$pond.id],xlab=NULL,ylab=NULL,yaxt='n',xaxt='n')
    }
    plot_smooth(model.name, view="sc.gly", cond=list('date'=date.x,'sc.imi'=0,'o.nut'='low'), rm.ranef=TRUE, rug=F,col=alpha(gly.cols[8],0.7),lty=1,add=T,print.summary = F,se=0)
    plot_smooth(model.name, view="sc.gly", cond=list('date'=date.x,'sc.imi'=0,'o.nut'='high'), rm.ranef=TRUE, rug=F,col=alpha(gly.cols[8],0.7),lty=2,lwd=1.5,add=T,print.summary = F,se=0)
    plot_smooth(model.name, view="sc.imi", cond=list('date'=date.x,'sc.gly'=0,'o.nut'='low'), rm.ranef=TRUE, rug=F,col=alpha(imi.cols[8],0.7),lty=1,add=T,print.summary = F,se=0)
    plot_smooth(model.name, view="sc.imi", cond=list('date'=date.x,'sc.gly'=0,'o.nut'='high'), rm.ranef=TRUE, rug=F,col=alpha(imi.cols[8],0.7),lty=2,lwd=1.5,add=T,print.summary = F,se=0)
    fitvals <- predict.gam(model.name, newdata = list('date' = rep(date.x,200), 'sc.gly' = rep(seq(0,1,length.out=100),2), 'sc.imi' = rep(seq(0,1,length.out=100),2), 'o.nut' = c(rep('low',100),rep('high',100)), 'site.f' = rep('C1',200)), exclude = s(date,site.f))
    points(fitvals[1:100]~seq(0,1,length.out=100),type='l',col=alpha(both.cols[8],0.7),lty=1)  
    points(fitvals[101:200]~seq(0,1,length.out=100),type='l',col=alpha(both.cols[8],0.7),lty=2,lwd=1.5)  
  }
}

