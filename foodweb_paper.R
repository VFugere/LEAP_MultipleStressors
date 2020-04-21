## Code to analyze LEAP 2016 MS data
## Vincent Fugère 2019-2020

rm(list=ls())

library(tidyverse)
library(skimr)
library(magrittr)
library(readxl)

library(scales)
library(shape)
library(RColorBrewer)
library(viridis)
library(plotrix)

library(corrgram)
library(lme4)
library(performance)
library(sjPlot)
library(sjlabelled)
library(sjmisc)
library(mgcv)
library(itsadug)

devtools::source_url('https://raw.githubusercontent.com/VFugere/Rfuncs/master/vif.R')
devtools::source_url('https://raw.githubusercontent.com/VFugere/Rfuncs/master/utils.R')

#### basic data

Sampling.dates <-c('17/08/16','23/08/16','31/08/16','15/09/16','20/09/16','28/09/16')
Sampling.dates <- as.Date(Sampling.dates, format = '%d/%m/%y')
Sampling.dates <- format(Sampling.dates, '%j')
Sampling.dates <- as.numeric(Sampling.dates)
Sampling.dates <- Sampling.dates - 229 #day 1 of exp is Julian day 230

pulse.dates <- as.Date(c('22/08/2016','19/09/2016','30/09/2016'), format = '%d/%m/%Y')
pulse.dates <- format(pulse.dates, '%j')
pulse.dates <- as.numeric(pulse.dates) - 229

treat <- read.csv('~/Google Drive/Recherche/LEAP Postdoc/2016/raw data/LEAP2016treatments.csv', stringsAsFactors = F) %>%
  rename('site' = pond,'nut' = nut.f, 'gly' = gly.lvl, 'imi' = imi.lvl, 'gly.target.ppb' = gly.conc, 'imi.target.ppb' = imi.conc) %>%
  select(-nut.ug.P,-nut.rel,-array,-nb)
treat$nut.fac <- factor(treat$nut, levels = c('low','high'))
treat$nut.num <- as.numeric(treat$nut.fac)
treat$pond.id <- 1:48

#three.cols <- c('#058A51','#130CA6','#EB0C00')
#library(wesanderson)
#three.cols <- wes_palette('Zissou1',5)[c(3,1,5)]
three.cols <- brewer.pal(3,'Dark2')[c(2,3,1)]
glycolfunc <- colorRampPalette(c("gray80", three.cols[1]))
gly.cols <- glycolfunc(8)
imicolfunc <- colorRampPalette(c("gray80", three.cols[2]))
imi.cols <- imicolfunc(8)
bothcolfunc <- colorRampPalette(c("gray80", three.cols[3]))
both.cols <- bothcolfunc(8)
allcols <- c(gly.cols,gly.cols,imi.cols,imi.cols,both.cols,both.cols)
pchs <- c(16,15)

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

#### nutrient and pesticide data ####

nut <- read_xlsx('~/Google Drive/Recherche/LEAP Postdoc/2016/raw data/2016nutrients.xlsx') %>%
  group_by(site,time.point,var) %>% summarize(day = mean(day), conc = mean(conc.ug.per.L, na.rm=T)) %>%
  ungroup %>% filter(var != 'SRP', time.point != 4) %>% select(-time.point) %>% spread(var,conc) %>% rename(date = day)

gly <- read_xlsx('~/Google Drive/Recherche/LEAP Postdoc/2016/raw data/Contaminants/gly_clean.xlsx') %>% select(-time.point, -gly.expected.ppb, -gly.max.ppb)
gly$date <- gly$date %>% as.Date(format = '%d.%m.%Y') %>% format('%j') %>% as.numeric
gly$date <- gly$date - 229

imi <- read_xlsx('~/Google Drive/Recherche/LEAP Postdoc/2016/raw data/Contaminants/imi_clean.xlsx') %>% select(-time.point, -notes)
imi$date <- imi$date %>% as.Date(format = '%d.%m.%Y') %>% format('%j') %>% as.numeric
imi$date <- imi$date - 229

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

#### YSI & metab data ####

YSI <- read.csv('~/Google Drive/Recherche/LEAP Postdoc/2016/raw data/YSI_exp.csv',header=T,stringsAsFactors=F)
YSI$date <- sapply(strsplit(YSI$time, '[ ]'), function(x) x[1])
YSI$time <- sapply(strsplit(YSI$time, '[ ]'), function(x) x[2])
YSI[YSI$pond == 'F1','pond'] <- 'H3'
names(YSI)[7] <- 'site' #to be consistent with other data frames
YSI$date <- as.Date(YSI$date, format = '%d/%m/%Y')
YSI$date <- format(YSI$date, '%j')
YSI$date <- as.numeric(YSI$date)
YSI$date <- YSI$date - 229 #day 1 of exp is Julian day 230
YSI <- YSI[YSI$date < 45,]
# adding NEP
YSI$period <- as.factor(YSI$period)
YSI$TP <- as.factor(YSI$date)
YSI <- YSI %>% group_by(TP, site) %>%
  mutate(NEP = DO.mg.L[2] - DO.mg.L[1], pH.diff = pH[2] - pH[1], pH.mean = mean(pH[1],pH[2]), temp.mean = mean(temp.C[1],temp.C[2]), DO.mean = mean(DO.mg.L[1],DO.mg.L[2]), SPC.mean = mean(SPC.uS.cm[1],SPC.uS.cm[2])) %>%
  ungroup
YSI <- YSI[YSI$period != 'dusk',]
#correcting for offset in 5th time point due to calibration problems.
YSI$NEP[YSI$date == 35] <- YSI$NEP[YSI$date == 35] - 1
YSI <- YSI %>% select(date,site,NEP:SPC.mean)

#### fluoroprobe data ####

FP <- read.csv('~/Google Drive/Recherche/LEAP Postdoc/2016/raw data/fluoroprobe.csv',header=T,stringsAsFactors=F) %>% select(-X) %>% as.data.frame
FP$date <- as.Date(FP$date, format = '%d/%m/%Y')
FP$date <- format(FP$date, '%j')
FP$date <- as.numeric(FP$date)
FP$date <- FP$date - 229 #day 1 of exp is Julian day 230
FP[FP$site == 'F1','site'] <- 'H3'
FP <- filter(FP, site %!in% c('LAKE','WELL'))
FP[FP$date == 2,'date'] <- 1
FP[FP$date == 8,'date'] <- 7
FP <- FP[order(FP$date),]
FP <- select(FP, date,site,greens:cryptos,total) %>%
  filter(date %in% Sampling.dates)

#### flow cytometry data ####

FC <- read_xlsx('~/Google Drive/Recherche/LEAP Postdoc/2016/raw data/bacterial_counts.xlsx')
FC[FC$pond == 'F1','pond'] <- 'H3'
FC <- select(FC, day, pond, mean.cells.ul) %>% 
  filter(!is.na(mean.cells.ul)) %>%
  mutate('day' = as.numeric(day)) %>% 
  rename('date' = day, 'site' = pond, 'BA' = mean.cells.ul) %>%
  filter(date %in% Sampling.dates) %>% arrange(date,site)

#### ecoplates data ####

EP <- read_xlsx('~/Google Drive/Recherche/LEAP Postdoc/2016/raw data/Ecoplates/Ecoplates_clean.xlsx', sheet = 'by_substrate')
EP2 <- read_xlsx('~/Google Drive/Recherche/LEAP Postdoc/2016/raw data/Ecoplates/Ecoplates_clean.xlsx', sheet = 'by_guild') %>% select(date:Amines_amides)
EP <- left_join(EP, EP2, by = c('date','site')) %>%
  filter(date %in% Sampling.dates) %>% arrange(date,site)
rm(EP2)

#### bind and clean ####

data <- inner_join(FP,YSI, by = c('date','site')) %>%
  inner_join(FC, by = c('date','site')) %>%
  inner_join(EP, by = c('date','site')) %>%
  left_join(treat, by = c('site')) %>% 
  select(-gly.target.ppb,-imi.target.ppb,-water) %>%
  mutate(nut = as.numeric(factor(nut, levels=c('low','high')))) %>%
  select(date, site, gly:pond.id, NEP:SPC.mean, greens:total, BA, AWCD:Amines_amides, everything())

#adding ordered factor for GAMs
data$o.nut <- as.ordered(data$nut.fac)

#rescaling pesticides gradients from 0 to 1 to compare effect with nutrient factor
data$sc.gly <- rescale(data$gly, c(0,1))
data$sc.imi <- rescale(data$imi, c(0,1))

#adding site factor
data$site.f <- as.factor(data$site)

#adding date factor
data$date.f <- as.factor(data$date)

#reordering
data <- select(data, date:pond.id,o.nut:date.f,everything())

#### some GAMMs ####

ba.model <- gam(log10(BA) ~ o.nut + s(date,k=6) + ti(date,sc.gly, k=6) + ti(date,sc.gly, by = o.nut, k=6) + ti(date,sc.imi, k=6) + ti(date,sc.imi, by = o.nut, k=6) + ti(date,sc.gly,sc.imi, k=6) + ti(date,sc.gly,sc.imi, by = o.nut, k=6) + s(date, site.f, bs='fs',k=5, m=2), data=data, method = 'REML')

ba.model <- gam(log10(BA) ~ o.nut + s(date,k=4) + ti(date,sc.gly, k=5) + ti(date,sc.gly, by = o.nut, k=5) + ti(date,sc.imi, k=4) + ti(date,sc.imi, by = o.nut, k=4) + ti(date,sc.gly,sc.imi, k=4) + s(date, site.f, bs='fs',k=4, m=2), data=data, method = 'REML')
summary(ba.model)
gam.check(ba.model)
plot(ba.model)

ba.model <- gam(log10(BA) ~ o.nut + te(date,sc.gly,sc.imi, by = o.nut, k=4) + s(date, site.f, bs='fs',k=4, m=2), data=data, method = 'REML')
summary(ba.model)
gam.check(ba.model)
plot(ba.model)
fit1 <- fitted(ba.model)

ba.model2 <- gam(log10(BA) ~ o.nut + s(date,k=4) + s(sc.gly, k = 4) + s(sc.imi, k = 4) + ti(date,sc.gly, k=4) + ti(date,sc.gly, by = o.nut, k=4) + ti(date,sc.imi, k=4) + ti(date,sc.imi, by = o.nut, k=4) + ti(date,sc.gly,sc.imi, k=4) + ti(date,sc.gly,sc.imi, by = o.nut, k=4) + s(date, site.f, bs='fs',k=3, m=2), data=data, method = 'REML')
fit2 <- fitted(ba.model2)
summary(ba.model2)
gam.check(ba.model2)

plot(fit2~fit1)
abline(a=0,b=1)

ba.model3 <- gam(log10(BA) ~ o.nut + s(date,k=4) + s(sc.gly, k = 4) + s(sc.imi, k = 4) + ti(date,sc.gly, k=4) + ti(date,sc.gly, by = o.nut, k=4) + ti(date,sc.imi, k=4) + ti(date,sc.imi, by = o.nut, k=4) + ti(date,sc.gly,sc.imi, k=4) + s(date, site.f, bs='fs',k=3, m=2), data=data, method = 'REML')
fit3 <- fitted(ba.model3)
summary(ba.model3)
gam.check(ba.model3)

plot(fit2~fit3)

fitted <- as.data.frame(predict(ba.model3, se.fit = T,exclude='s(date,site)'))
fitted$lwr <- fitted$fit - 1.96*fitted$se.fit
fitted$upr <- fitted$fit + 1.96*fitted$se.fit
ylims <- range(c(fitted$lwr,fitted$upr))

plot_smooth(ba.model3, view="date", cond=list('sc.gly'=0,'sc.imi'=0,'o.nut'='low'), rm.ranef=TRUE, rug=F,col=gly.cols[1],print.summary = F,hide.label = T,yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l',legend_plot_all = F, h0=NA,ylim = ylims)
plot_smooth(ba.model3, view="date", cond=list('sc.gly'=0.5,'sc.imi'=0,'o.nut'='low'), rm.ranef=TRUE, rug=F,col=gly.cols[4],add=T)
plot_smooth(ba.model3, view="date", cond=list('sc.gly'=1,'sc.imi'=0,'o.nut'='low'), rm.ranef=TRUE, rug=F,col=gly.cols[8],add=T)

plot_smooth(ba.model3, view="date", cond=list('sc.gly'=0,'sc.imi'=0,'o.nut'='high'), rm.ranef=TRUE, rug=F,col=gly.cols[1],print.summary = F,hide.label = T,yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l',legend_plot_all = F, h0=NA,ylim = ylims)
plot_smooth(ba.model3, view="date", cond=list('sc.gly'=0.5,'sc.imi'=0,'o.nut'='high'), rm.ranef=TRUE, rug=F,col=gly.cols[4],add=T)
plot_smooth(ba.model3, view="date", cond=list('sc.gly'=1,'sc.imi'=0,'o.nut'='high'), rm.ranef=TRUE, rug=F,col=gly.cols[8],add=T)

## contour plots

zlims <- range(fitted(ba.model3))
par(mfrow=c(2,6),mar=c(2,2,2,2),oma=c(2,2,0,0),cex=1)
fvisgam(ba.model3, view = c('sc.gly','sc.imi'), cond = list('date' = 1, 'o.nut' = 'low'), zlim=zlims, add.color.legend=F,hide.label=T,plot.type = 'contour', lwd=1.5,color = inferno(100), main = NULL,rm.ranef = T,dec=1,print.summary = F)
mtext('LN day 1',3,cex=0.7)
fvisgam(ba.model3, view = c('sc.gly','sc.imi'), cond = list('date' = 7, 'o.nut' = 'low'), zlim=zlims, add.color.legend=F,hide.label=T,plot.type = 'contour', lwd=1.5,color = inferno(100), main = NULL,rm.ranef = T,dec=1,print.summary = F)
mtext('LN day 7',3,cex=0.7)
fvisgam(ba.model3, view = c('sc.gly','sc.imi'), cond = list('date' = 15, 'o.nut' = 'low'), zlim=zlims, add.color.legend=F,hide.label=T,plot.type = 'contour', lwd=1.5,color = inferno(100), main = NULL,rm.ranef = T,dec=1,print.summary = F)
mtext('LN day 15',3,cex=0.7)
fvisgam(ba.model3, view = c('sc.gly','sc.imi'), cond = list('date' = 30, 'o.nut' = 'low'), zlim=zlims, add.color.legend=F,hide.label=T,plot.type = 'contour', lwd=1.5,color = inferno(100), main = NULL,rm.ranef = T,dec=1,print.summary = F)
mtext('LN day 30',3,cex=0.7)
fvisgam(ba.model3, view = c('sc.gly','sc.imi'), cond = list('date' = 35, 'o.nut' = 'low'), zlim=zlims, add.color.legend=F,hide.label=T,plot.type = 'contour', lwd=1.5,color = inferno(100), main = NULL,rm.ranef = T,dec=1,print.summary = F)
mtext('LN day 35',3,cex=0.7)
fvisgam(ba.model3, view = c('sc.gly','sc.imi'), cond = list('date' = 43, 'o.nut' = 'low'), zlim=zlims, add.color.legend=T,hide.label=T,plot.type = 'contour', lwd=1.5,color = inferno(100), main = NULL,rm.ranef = T,dec=1,print.summary = F)
mtext('LN day 43',3,cex=0.7)
fvisgam(ba.model3, view = c('sc.gly','sc.imi'), cond = list('date' = 1, 'o.nut' = 'high'), zlim=zlims, add.color.legend=F,hide.label=T,plot.type = 'contour', lwd=1.5,color = inferno(100), main = NULL,rm.ranef = T,dec=1,print.summary = F)
mtext('HN day 1',3,cex=0.7)
fvisgam(ba.model3, view = c('sc.gly','sc.imi'), cond = list('date' = 7, 'o.nut' = 'high'), zlim=zlims, add.color.legend=F,hide.label=T,plot.type = 'contour', lwd=1.5,color = inferno(100), main = NULL,rm.ranef = T,dec=1,print.summary = F)
mtext('HN day 7',3,cex=0.7)
fvisgam(ba.model3, view = c('sc.gly','sc.imi'), cond = list('date' = 15, 'o.nut' = 'high'), zlim=zlims, add.color.legend=F,hide.label=T,plot.type = 'contour', lwd=1.5,color = inferno(100), main = NULL,rm.ranef = T,dec=1,print.summary = F)
mtext('HN day 15',3,cex=0.7)
fvisgam(ba.model3, view = c('sc.gly','sc.imi'), cond = list('date' = 30, 'o.nut' = 'high'), zlim=zlims, add.color.legend=F,hide.label=T,plot.type = 'contour', lwd=1.5,color = inferno(100), main = NULL,rm.ranef = T,dec=1,print.summary = F)
mtext('HN day 30',3,cex=0.7)
fvisgam(ba.model3, view = c('sc.gly','sc.imi'), cond = list('date' = 35, 'o.nut' = 'high'), zlim=zlims, add.color.legend=F,hide.label=T,plot.type = 'contour', lwd=1.5,color = inferno(100), main = NULL,rm.ranef = T,dec=1,print.summary = F)
mtext('HN day 35',3,cex=0.7)
fvisgam(ba.model3, view = c('sc.gly','sc.imi'), cond = list('date' = 43, 'o.nut' = 'high'), zlim=zlims, add.color.legend=T,hide.label=T,plot.type = 'contour', lwd=1.5,color = inferno(100), main = NULL,rm.ranef = T,dec=1,print.summary = F)
mtext('HN day 43',3,cex=0.7)
mtext('glyphosate',1,outer=T,line=1)
mtext('imidacloprid',2,outer=T)

## linear model

ba.mod.lin <- lmer(log10(BA) ~ date.f + nut.fac + date.f:nut.fac + date.f:sc.gly + date.f:sc.imi + date.f:sc.gly:sc.imi + date.f:sc.gly:nut.fac +
                     date.f:sc.imi:nut.fac + date.f:sc.gly:sc.imi:nut.fac + (1|site.f),data)
ba.mod.lin <- lmer(log10(BA) ~ o.nut + date.f:o.nut + date.f:sc.gly + date.f:sc.imi + date.f:sc.gly:sc.imi + date.f:sc.gly:o.nut +
                     date.f:sc.imi:o.nut + date.f:sc.gly:sc.imi:o.nut + (1|site.f),data)
plot_model(ba.mod.lin)
performance(ba.mod.lin)
