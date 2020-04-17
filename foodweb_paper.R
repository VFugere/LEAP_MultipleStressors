## Code to analyze LEAP 2016 MS data
## Vincent Fugère 2019-2020

rm(list=ls())
library(tidyverse)
library(scales)
library(shape)
library(readxl)
library(RColorBrewer)
library(plotrix)
library(randomcoloR)
library(skimr)
library(magrittr)
library(corrgram)
library(readxl)
library(sjPlot)
library(sjlabelled)
library(sjmisc)

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

library(wesanderson)
#three.cols <- c('#058A51','#130CA6','#EB0C00')
three.cols <- wes_palette('Zissou1',5)[c(3,1,5)]
glycolfunc <- colorRampPalette(c("gray80", three.cols[1]))
gly.cols <- glycolfunc(8)
imicolfunc <- colorRampPalette(c("gray80", three.cols[2]))
imi.cols <- imicolfunc(8)
bothcolfunc <- colorRampPalette(c("gray80", three.cols[3]))
both.cols <- bothcolfunc(8)

allcols <- c(gly.cols,gly.cols,imi.cols,imi.cols,both.cols,both.cols)
pchs <- c(1,0)

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
  select(date, site, nut, gly, imi, NEP:SPC.mean, greens:total, BA, AWCD:Amines_amides, everything())

#### Figure: did treatments work? ####

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

######

## nutrients

pdf('~/Desktop/figure2.pdf',width=11.5,height = 8.5,pointsize = 12)
layout(rbind(c(1,1,2,2),c(3,3,4,5),c(6,6,7,8)))
par(cex=1,mar=c(4,4,1,1))

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

#### data exploration: correlations in the dataset ####

corrgram(data, order=F, lower.panel=panel.shade,
         upper.panel=panel.pie, text.panel=panel.txt,
         main=NULL)

#hard to see because of non-linear effects. Splitting by sampling date
pdf('~/Desktop/corrgrams.pdf',height = 12,width = 12,onefile = T,pointsize = 8)
for(d in Sampling.dates){
  data.sub <- filter(data, date == d) %>% select(-date,-site)
  data.sub <- data.sub[,apply(data.sub, 2, sd) != 0]
  corrgram(data.sub, order=F, lower.panel=panel.shade,
           upper.panel=panel.pie, text.panel=panel.txt,
           main=paste0('Day ',d))
}
dev.off()


#for talk
vars <- c('total','BA','NEP','use','Glycogen')
var.names <- c(expression(chlorophyll~italic(a)~(mu*g~L^-1)),
               expression(bacterial~abundance~(cells~mu*l^-1)),
               expression(NEP~(Delta*mg~DO~L^-1)),
               expression(C~substrates~used),
               expression(glycogen~respiration~(OD~units~above~blanks)))
var.log <- c(T,T,F,F,F)

# #ecoplates exploratory analysis
# vars <- names(data[13:52])
# var.names <- names(data[13:52])
# var.log <- rep(F,40)

lwds <- seq(from = 1, to = 2.5, length.out = 8)
ln.alpha <- seq(from = 0.6, to = 0.9, length.out = 8)

#pdf('~/Desktop/MSplots.pdf',width=8,height = 6,pointsize = 12,onefile = T)
par(mfrow=c(3,2),mar=c(2,2,1,1),oma=c(2.5,2.5,1,1),cex=1)

for(v in 1:length(vars)){

tmp <- data[,c(colnames(data)[1:5],vars[v])]

for(letter in c('C','D','E','H','J','K')){
  if(var.log[v] == T){
    plotfunctions::emptyPlot(xlim=c(1,43),ylim=range(tmp[,6]),yaxt='n',xaxt='n',ann=F, bty='l', log = 'y')
  }else{
    plotfunctions::emptyPlot(xlim=c(1,43),ylim=range(tmp[,6]),yaxt='n',xaxt='n',ann=F, bty='l')
  }
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
  tmp.sub <- filter(tmp, str_detect(site, letter))
  if(letter %in% c('C','D')){
    cols <- gly.cols
  }else if(letter %in% c('E','H')){
    cols <- imi.cols
  }else{
    cols <- both.cols
  }
  for(i in 1:8){
    pond <- unique(tmp.sub$site)[i]
    sub.sub <- filter(tmp.sub, site == pond)
    points(x=sub.sub$date,y=sub.sub[,6],type='l',lwd=lwds[i],col=alpha(cols[i],ln.alpha[i]))
  }
  abline(v=pulse.dates[1:2],lty=3)
}
mtext('day of experiment',side=1,outer=T,line=1)
mtext(var.names[v],side=2,outer=T,line=1)

}

#dev.off()

#for Jesse

ba <- select(data, date:imi, BA) %>% mutate_at(vars('site':'imi'), as.factor)
hist(ba$BA) #gamma
plot(BA~date,ba)
library(lme4)
f1 <- formula(BA ~ date+gly+imi+nut)
f1 <- update(f1, . ~ . + .:date)
f1 <- update(f1, . ~ . + (1+date|site))
# m1 <- glmer(f1, ba, family=Gamma(link='log'))
# plot(m1)
# anova(m1)
# f2 <- update(f1, log(BA) ~ .)
# m2 <- lmer(f1, ba)
# summary(m2)
# plot(m2) #bad

library(brms)



m3 <- brm(f1, ba, family=Gamma(link='log'), cores=4, iter = 10000)
summary(m3)
plot(m3)
pp_check(m3)
marginal_effects(m3)
m3 -> m.fact

ba <- ba %>% mutate_at(vars('gly','imi'), as.numeric)
m3 <- brm(f1, ba, family=Gamma(link='log'), cores=8, iter = 10000)
summary(m3)
plot(m3)
pp_check(m3)
marginal_effects(m3)

#### Ecoplates ordination (in progress) ####

#com2 <- select(data, date:imi,Polymers:Amines_amides)
com2 <- select(data, date:imi,`Pyruvic acid methyl ester`:Putrescine)
com2$fill.alpha <- 0.1
com2$fill.alpha[com2$date == 43] <- 0.8
com2$pt.cex <- 1
com2$pt.cex[com2$date > 1 & com2$date < 43] <- 0.3
com2$col <- rep(c(gly.cols,gly.cols,imi.cols,imi.cols,both.cols,both.cols),6)
com2$pch <- rep(rep(c(21,22,21,22,22,21),each=8),6)
com2$site <- as.factor(com2$site)

#spe <- com2 %>% select(Polymers:Amines_amides) %>% as.matrix
spe <- com2 %>% select(`Pyruvic acid methyl ester`:Putrescine) %>% as.matrix
spe <- spe-min(spe)+0.01
spe <- log1p(spe)
row.names(spe) <- com2$site

library(vegan)
ordi <- metaMDS(spe, distance = 'bray', k = 2, autotransform = FALSE, trymax = 1000)

pdf('~/Desktop/ordi.pdf',height = 5,width = 6,pointsize = 12)
par(mfrow=c(1,1),mar=c(4,4,1,1),oma=c(0,0,0,0))
g<-ordi$points[,1:2]
plot(g[,2] ~ g[,1], type = "n",yaxt='n',xaxt='n',ann=F,xlim=c(-0.05,0.05),ylim=c(-0.05,0.05))
title(xlab='NMDS dimension 1',cex.lab=1,line = 2.5)
axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
title(ylab='NMDS dimension 2',cex.lab=1,line = 2.5)
axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
for(i in 1:nlevels(com2$site)){
  subdat <- g[row.names(g) == levels(com2$site)[i],]
  subdat <- cbind(subdat, com2[com2$site == levels(com2$site)[i],c('pch','col')])
  lncol <- com2$col[com2$site == levels(com2$site)[i]][1]
  #arrows(x0 = subdat[1,1], x1 = subdat[5,1], y0 = subdat[1,2], y1 = subdat[5,2], lty=1, lwd=1, col = alpha(lncol,1))
  lines(x = subdat[,1], y = subdat[,2], lty=1, lwd=1, col = alpha(lncol,0.5))
  #segments(x0=subdat[1,1],x1=subdat[5,1],y0=subdat[1,2],y1=subdat[5,2], col = alpha(cols[subdat$col],0.8))
  #subdat <- subdat[c(1,5),]
  #subdat$pch[1] <- subdat$pch[1] - 15
  points(subdat[6,2]~subdat[6,1], pch = subdat$pch, col=1, bg = subdat$col, cex = 1.3)
}
stress.val <- round(ordi$stress,2)
legend('topright',bty='n',legend='Stress = 0.11')
dev.off()
