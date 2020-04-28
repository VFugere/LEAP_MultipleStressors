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

pdf('~/Desktop/figure1.pdf',width=5.5,height=6,pointsize = 12)
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

pdf('~/Desktop/figure2.pdf',width=11.5,height = 8.5,pointsize = 12)
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

pdf('~/Desktop/FigS1.pdf',width=7,height = 5,pointsize = 8)
par(mfrow=c(4,3),mar=c(2,2,1,1),oma=c(2.5,2.5,1,1),cex=1)
pchs.2 <- c(16,15)
vars <- c('pH.mean','DO.mean','SPC.mean')
var.names <- c(expression(pH),expression(DO~'(mg/L)'),expression(SPC~(mu*g/cm)))

for(v in 1:3){
  tmp <- merged.data[,c(colnames(merged.data)[1:13],vars[v])]
  tmp$date.idx <- tmp$date - 0.6
  tmp$date.idx[tmp$nut.num == 2] <- tmp$date[tmp$nut.num == 2] + 0.6
  for(letter in c('C|D','E|H','J|K')){
    plotfunctions::emptyPlot(xlim=c(1,43),ylim=range(tmp[,14]),yaxt='n',xaxt='n',ann=F, bty='l')
    axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
    axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
    tmp.sub <- filter(tmp, str_detect(site, letter))
    for(i in 1:n_distinct(tmp.sub$site)){
      pond <- unique(tmp.sub$site)[i]
      sub.sub <- filter(tmp.sub, site == pond)
      points(x=sub.sub$date.idx,y=sub.sub[,14],type='o',lwd=1.2, pch=pchs.2[sub.sub$nut.num], col=alpha(allcols[sub.sub$pond.id],1))
    }
    abline(v=pulse.dates[1:2],lty=3)
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
merged.data$RUE <- with(merged.data, total/total_zoo)
merged.data$RUE[is.na(merged.data$RUE)] <- 0

#### fitting GAMMs ####

ba.mo <- gam(log10(BA) ~ o.nut + ti(date,k=5) + ti(date,sc.gly, k=6) + ti(date,sc.gly, by = o.nut, k=3) + ti(date,sc.imi, k=6) + ti(date,sc.imi, by = o.nut, k=3) + ti(date,sc.gly,sc.imi, k=5) + ti(date,sc.gly,sc.imi, by = o.nut, k=3) + s(date, site.f, bs='fs',k=3), data=merged.data, method = 'REML')
chla.m <- gam(log10(total) ~ o.nut + ti(date,k=5) + ti(date,sc.gly, k=6) + ti(date,sc.gly, by = o.nut, k=3) + ti(date,sc.imi, k=6) + ti(date,sc.imi, by = o.nut, k=3) + ti(date,sc.gly,sc.imi, k=5) + ti(date,sc.gly,sc.imi, by = o.nut, k=3) + s(date, site.f, bs='fs',k=3), data=merged.data, method = 'REML')
zoo.m <- gam(log10p(total_zoo) ~ o.nut + ti(date,k=5) + ti(date,sc.gly, k=6) + ti(date,sc.gly, by = o.nut, k=3) + ti(date,sc.imi, k=6) + ti(date,sc.imi, by = o.nut, k=3) + ti(date,sc.gly,sc.imi, k=5) + ti(date,sc.gly,sc.imi, by = o.nut, k=3) + s(date, site.f, bs='fs',k=3), data=merged.data, method = 'REML')

use.m <- gam(use ~ o.nut + ti(date,k=5) + ti(date,sc.gly, k=6) + ti(date,sc.gly, by = o.nut, k=3) + ti(date,sc.imi, k=6) + ti(date,sc.imi, by = o.nut, k=3) + ti(date,sc.gly,sc.imi, k=5) + ti(date,sc.gly,sc.imi, by = o.nut, k=3) + s(date, site.f, bs='fs',k=3), data=merged.data, method = 'REML')
nep.m <- gam(log10p(NEP) ~ o.nut + ti(date,k=5) + ti(date,sc.gly, k=6) + ti(date,sc.gly, by = o.nut, k=3) + ti(date,sc.imi, k=6) + ti(date,sc.imi, by = o.nut, k=3) + ti(date,sc.gly,sc.imi, k=5) + ti(date,sc.gly,sc.imi, by = o.nut, k=3) + s(date, site.f, bs='fs',k=3), data=merged.data, method = 'REML')
rue.m <- gam(log10p(RUE) ~ o.nut + ti(date,k=5) + ti(date,sc.gly, k=6) + ti(date,sc.gly, by = o.nut, k=3) + ti(date,sc.imi, k=6) + ti(date,sc.imi, by = o.nut, k=3) + ti(date,sc.gly,sc.imi, k=5) + ti(date,sc.gly,sc.imi, by = o.nut, k=3) + s(date, site.f, bs='fs',k=3), data=merged.data, method = 'REML')

