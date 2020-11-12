## Code to analyze LEAP 2016 MS data
## Vincent Fugère 2017-2020
## This is a data exploration script for the ecoplates data

rm(list=ls())
library(tidyverse)
library(readxl)
library(skimr)
library(magrittr)
library(scales)
library(shape)
library(RColorBrewer)
library(plotrix)
library(randomcoloR)
library(vegan)
library(corrgram)

devtools::source_url('https://raw.githubusercontent.com/VFugere/Rfuncs/master/utils.R')

#### basic data ####

load('~/Google Drive/Recherche/LEAP Postdoc/2016/MSdata.RData')

pchs <- c(1,0)

gly.cols <- c('gray90',brewer.pal(9, 'Reds')[2:8])
imi.cols <- c('gray90',brewer.pal(9, 'Blues')[2:8])
both.cols <- c('gray90',brewer.pal(9, 'Greens')[2:8])

glycolfunc <- colorRampPalette(gly.cols)
imicolfunc <- colorRampPalette(imi.cols)
bothcolfunc <- colorRampPalette(both.cols)

allcols <- c(gly.cols,gly.cols,imi.cols,imi.cols,both.cols,both.cols)

#### plots ####

corrgram(merged.data[,c(1,3,4,15:24,71:101)], order=TRUE, lower.panel=panel.shade,
         upper.panel=panel.pie, text.panel=panel.txt,
         main=NULL)

#for talk
vars <- c('BA','use','Glycogen')
var.names <- c(expression(bacterial~abundance~(cells~L^-1)),
               expression(C~substrates~used),
               expression(glycogen~respiration~(OD~units~above~blanks)))
var.log <- c(T,F,F)

lwds <- seq(from = 1, to = 2.5, length.out = 8)
ln.alpha <- seq(from = 0.6, to = 0.9, length.out = 8)

pdf('~/Desktop/MSplots.pdf',width=8,height = 6,pointsize = 12,onefile = T)
par(mfrow=c(3,2),mar=c(2,2,1,1),oma=c(2.5,2.5,1,1),cex=1)

for(v in 1:length(vars)){
  
  tmp <- merged.data[,c(colnames(merged.data)[1:5],vars[v])]
  
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
      sub.sub <- filter(tmp.sub, site == pond) %>% as.data.frame
      points(x=sub.sub$date,y=sub.sub[,6],type='l',lwd=lwds[i],col=alpha(cols[i],ln.alpha[i]))
    }
    abline(v=pulse.dates[1:2],lty=3)
  }
  mtext('day of experiment',side=1,outer=T,line=1)
  mtext(var.names[v],side=2,outer=T,line=1)
  
}

dev.off()

#### correlation between FC, Ecoplates, 16S read number and DNA concentration

dna.conc <- read_xlsx('/Users/vincentfugere/Google Drive/Recherche/LEAP Postdoc/2016/raw data/Sequencing/reads_DNAconc_vincent.xlsx', sheet = 'DNA_conc') %>%
  select(Sample,DNA.ng.ml,ecoplates.date) %>% rename('site' = Sample, 'date' = ecoplates.date) %>%
  filter(!is.na(date))

reads <- read_xlsx('/Users/vincentfugere/Google Drive/Recherche/LEAP Postdoc/2016/raw data/Sequencing/reads_DNAconc_vincent.xlsx', sheet = 'reads_clean') %>%
  select(site,ecoplates.date,nonchim,run) %>% rename('date' = ecoplates.date) %>%
  filter(!is.na(date))

library(mgcv)
library(itsadug)

dat <- merged.data %>% select(date,site,gly:imi,NEP:Amines_amides) %>% left_join(dna.conc) %>%
  left_join(reads)

corrgram(dat, order=TRUE, lower.panel=panel.shade,
         upper.panel=panel.pie, text.panel=panel.txt,
         main=NULL)

boxplot(nonchim~run,dat)

dat$site <- as.factor(dat$site)

# m1 <- gam(BA~s(max,k=7)+s(max,site,bs='fs',k=3),data=dat)
# plot_smooth(m1,view='max')
# summary(m1)

pdf('~/Desktop/BAmethods.pdf',width=6,height=6,pointsize = 12)
par(mfrow=c(2,2),mar=c(4,4,1,1),oma=c(0,0,0,0),cex=1)
plot(BA~max,dat,pch=16,col='gray70',xlab='max AWCD',ylab='cell count (FC)',log='y')
plot(BA~day.of.reading,dat,pch=16,col='gray70',xlab='days before 0.5 AWCD',ylab='cell count (FC)',log='y')
plot(BA~DNA.ng.ml,dat,pch=16,col='gray70',xlab='DNA concentration',ylab='cell count (FC)',log='xy')
plot(BA~nonchim,dat,pch=16,col=c('gray70','red')[dat$run],xlab='non-chimeric 16S reads',ylab='cell count (FC)',log='xy')
legend('bottomleft',legend=c('run 1','run 2'),pch=16,col=c('gray70','red'),bty='n')
dev.off()

#### for Jesse ####

ba <- select(merged.data, date:imi, BA) %>% mutate_at(vars('site':'imi'), as.factor)
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

#com2 <- select(merged.data, date:imi,Polymers:Amines_amides)
com2 <- select(merged.data, date:imi,`Pyruvic acid methyl ester`:Putrescine)
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


# OLD AND USELESS
# # ploting the data exlcuding 3rd pulse
# 
# glycol <- brewer.pal(8,'Reds')
# imicol <- brewer.pal(8,'Blues')
# mscol <- brewer.pal(8,'Purples')
# 
# MS.dat <- epdat[epdat$date < pulse.dates[3],]
# MS.dat <- droplevels(MS.dat)
# 
# MS.dat$array <- treat$array[match(MS.dat$site,treat$pond)]
# MS.dat$col <- imicol[MS.dat$imi.col.idx]
# MS.dat$col[MS.dat$array == 'C'] <- glycol[MS.dat$gly.col.idx[MS.dat$array == 'C']]
# MS.dat$col[MS.dat$array == 'D'] <- glycol[MS.dat$gly.col.idx[MS.dat$array == 'D']]
# MS.dat$col[MS.dat$array == 'J'] <- mscol[MS.dat$gly.col.idx[MS.dat$array == 'J']]
# MS.dat$col[MS.dat$array == 'K'] <- mscol[MS.dat$gly.col.idx[MS.dat$array == 'J']]
# MS.dat$lty <- 1
# MS.dat[MS.dat$site == 'C1' | MS.dat$site == 'D1' | MS.dat$site == 'E1' | MS.dat$site == 'H1' | MS.dat$site == 'J1' | MS.dat$site == 'K1', 'col'] <- 'black'
# MS.dat[MS.dat$site == 'C1' | MS.dat$site == 'D1' | MS.dat$site == 'E1' | MS.dat$site == 'H1' | MS.dat$site == 'J1' | MS.dat$site == 'K1', 'lty'] <- 2
# 
# layout(rbind(c(1,2),c(3,5),c(3,4)),heights=c(1,0.2,0.8))
# par(mar = c(4.2,4.2,1,1))
# 
# plot(max~date,data = subset(MS.dat, array == 'C' | array == 'D'), type='n',yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l',ylim=range(MS.dat$max),log='y')
# title(ylab=max~AWCD~(density), cex.lab=1)
# axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
# title(xlab="day", cex.lab=1)
# axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
# for (i in 1:16){
#   tmp <- subset(MS.dat, array == 'C' | array == 'D')
#   tmp <- subset(tmp, site == levels(tmp$site)[i])
#   points(max~date,tmp,type='l',col=alpha(tmp$col[1],0.8),lty=tmp$lty,lwd=seq(1,1.7,0.1)[tmp$gly.col.idx[1]])
#   points(max~date,tmp,type='p',pch=tmp$pch[1],col=alpha(tmp$col[1],0.8))
# }
# abline(v=pulse.dates[1:2], lty=3)
# text(x=pulse.dates[1],y=range(MS.dat$max)[1],'glyphosate:pulse 1',pos=4)
# text(x=pulse.dates[2],y=range(MS.dat$max)[1],'pulse 2',pos=4)
# 
# plot(max~date,data = subset(MS.dat, array == 'E' | array == 'H'), type='n',yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l',ylim=range(MS.dat$max),log='y')
# title(ylab=max~AWCD~(density), cex.lab=1)
# axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
# title(xlab="day", cex.lab=1)
# axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
# for (i in 17:32){
#   tmp <- subset(MS.dat, array == 'E' | array == 'H')
#   tmp <- subset(tmp, site == levels(tmp$site)[i])
#   points(max~date,tmp,type='l',col=alpha(tmp$col[1],0.8),lty=tmp$lty,lwd=seq(1,1.7,0.1)[tmp$imi.col.idx[1]])
#   points(max~date,tmp,type='p',pch=tmp$pch[1],col=alpha(tmp$col[1],0.8))
# }
# abline(v=pulse.dates[1:2], lty=3)
# text(x=pulse.dates[1],y=range(MS.dat$max)[1],'imidacloprid:pulse 1',pos=4)
# text(x=pulse.dates[2],y=range(MS.dat$max)[1],'pulse 2',pos=4)
# 
# plot(max~date,data = subset(MS.dat, array == 'J' | array == 'K'), type='n',yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l',ylim=range(MS.dat$max),log='y')
# title(ylab=max~AWCD~(density), cex.lab=1)
# axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
# title(xlab="day", cex.lab=1)
# axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
# for (i in 33:48){
#   tmp <- subset(MS.dat, array == 'J' | array == 'K')
#   tmp <- subset(tmp, site == levels(tmp$site)[i])
#   points(max~date,tmp,type='l',col=alpha(tmp$col[1],0.8),lty=tmp$lty,lwd=seq(1,1.7,0.1)[tmp$gly.col.idx[1]])
#   points(max~date,tmp,type='p',pch=tmp$pch[1],col=alpha(tmp$col[1],0.8))
# }
# abline(v=pulse.dates[1:2], lty=3)
# text(x=pulse.dates[1],y=range(MS.dat$max)[1],'both pesticides:pulse 1',pos=4)
# text(x=pulse.dates[2],y=range(MS.dat$max)[1],'pulse 2',pos=4)
# 
# #modeling variable
# 
# MS.mod <- lmer(log(max) ~ 1 + date + fTime:gly.conc + nut + fTime:imi.conc + (1|site), data = MS.dat)
# summary(MS.mod)
# 
# #model coefficients
# MS.coefs <- as.data.frame(summary(MS.mod)$coefficients[-c(1,2),1:2])
# nb.tp <- nlevels(MS.dat$fTime)
# MS.coefs$ci <- MS.coefs$'Std. Error' * 1.96
# MS.coefs$lwr <- MS.coefs$Estimate - MS.coefs$ci
# MS.coefs$upr <- MS.coefs$Estimate + MS.coefs$ci
# MS.coefs$date <- c(as.numeric(levels(MS.dat$fTime))[1],rep(as.numeric(levels(MS.dat$fTime)),2))
# MS.coefs$time.stamp <- MS.coefs$date + c(0,rep(-0.4,nb.tp),rep(0.4,nb.tp))
# coef.col <- alpha(c('darkgreen',glycol[8],imicol[8]),0.8)
# MS.coefs$col <- c(coef.col[1],rep(coef.col[2],nb.tp),rep(coef.col[3],nb.tp))
# MS.coefs <- MS.coefs[-c(2,nb.tp+2),]
# 
# plot(Estimate~date,data=MS.coefs,type='n',yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l',xlim=range(MS.dat$date),ylim=range(c(MS.coefs$lwr,MS.coefs$upr)))
# title(ylab='model coefficient', cex.lab=1)
# title(xlab="day", cex.lab=1)
# axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
# axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
# abline(h=0,lty=1,col='gray')
# arrows(x0=MS.coefs$time.stamp,y0=MS.coefs$lwr,y1=MS.coefs$upr,length=0,lwd=1.5,col=MS.coefs$col)
# points(x=MS.coefs$time.stamp,y=MS.coefs$Estimate,col=MS.coefs$col,pch=16)
# abline(v=pulse.dates[1:2], lty=3)
# legend(x=10,y=0.1,bty='n',legend=c('nutrient addition (press treatment)','glyphosate:time','imidacloprid:time'),pch=16,col=coef.col,cex=1,y.intersp=1.2)
# 
# ### nb substrates used
# 
# layout(rbind(c(1,2),c(3,5),c(3,4)),heights=c(1,0.2,0.8))
# par(mar = c(4.2,4.2,1,1))
# 
# plot(use~date,data = subset(MS.dat, array == 'C' | array == 'D'), type='n',yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l',ylim=range(MS.dat$use))
# title(ylab=substrates~used~(diversity), cex.lab=1)
# axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
# title(xlab="day", cex.lab=1)
# axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
# for (i in 1:16){
#   tmp <- subset(MS.dat, array == 'C' | array == 'D')
#   tmp <- subset(tmp, site == levels(tmp$site)[i])
#   points(use~date,tmp,type='l',col=alpha(tmp$col[1],0.8),lty=tmp$lty,lwd=seq(1,1.7,0.1)[tmp$gly.col.idx[1]])
#   points(use~date,tmp,type='p',pch=tmp$pch[1],col=alpha(tmp$col[1],0.8))
# }
# abline(v=pulse.dates[1:2], lty=3)
# text(x=pulse.dates[1],y=range(MS.dat$use)[1],'glyphosate:pulse 1',pos=4)
# text(x=pulse.dates[2],y=range(MS.dat$use)[1],'pulse 2',pos=4)
# 
# plot(use~date,data = subset(MS.dat, array == 'E' | array == 'H'), type='n',yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l',ylim=range(MS.dat$use))
# title(ylab=substrates~used~(diversity), cex.lab=1)
# axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
# title(xlab="day", cex.lab=1)
# axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
# for (i in 17:32){
#   tmp <- subset(MS.dat, array == 'E' | array == 'H')
#   tmp <- subset(tmp, site == levels(tmp$site)[i])
#   points(use~date,tmp,type='l',col=alpha(tmp$col[1],0.8),lty=tmp$lty,lwd=seq(1,1.7,0.1)[tmp$imi.col.idx[1]])
#   points(use~date,tmp,type='p',pch=tmp$pch[1],col=alpha(tmp$col[1],0.8))
# }
# abline(v=pulse.dates[1:2], lty=3)
# text(x=pulse.dates[1],y=range(MS.dat$use)[1],'imidacloprid:pulse 1',pos=4)
# text(x=pulse.dates[2],y=range(MS.dat$use)[1],'pulse 2',pos=4)
# 
# plot(use~date,data = subset(MS.dat, array == 'J' | array == 'K'), type='n',yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l',ylim=range(MS.dat$use))
# title(ylab=substrates~used~(diversity), cex.lab=1)
# axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
# title(xlab="day", cex.lab=1)
# axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
# for (i in 33:48){
#   tmp <- subset(MS.dat, array == 'J' | array == 'K')
#   tmp <- subset(tmp, site == levels(tmp$site)[i])
#   points(use~date,tmp,type='l',col=alpha(tmp$col[1],0.8),lty=tmp$lty,lwd=seq(1,1.7,0.1)[tmp$gly.col.idx[1]])
#   points(use~date,tmp,type='p',pch=tmp$pch[1],col=alpha(tmp$col[1],0.8))
# }
# abline(v=pulse.dates[1:2], lty=3)
# text(x=pulse.dates[1],y=range(MS.dat$use)[1],'both pesticides:pulse 1',pos=4)
# text(x=pulse.dates[2],y=range(MS.dat$use)[1],'pulse 2',pos=4)
# 
# #modeling variable
# 
# MS.mod <- lmer(use ~ 1 + date + fTime:gly.conc + nut + fTime:imi.conc + (1|site), data = MS.dat)
# #MS.mod <- lmer(use ~ fTime*gly.conc*nut*imi.conc + (1|site), data = MS.dat)
# #MS.mod <- lmer(use ~ date*gly.conc*nut*imi.conc + (1|site), data = MS.dat)
# summary(MS.mod)
# 
# #model coefficients
# MS.coefs <- as.data.frame(summary(MS.mod)$coefficients[-c(1,2),1:2])
# nb.tp <- nlevels(MS.dat$fTime)
# MS.coefs$ci <- MS.coefs$'Std. Error' * 1.96
# MS.coefs$lwr <- MS.coefs$Estimate - MS.coefs$ci
# MS.coefs$upr <- MS.coefs$Estimate + MS.coefs$ci
# MS.coefs$date <- c(as.numeric(levels(MS.dat$fTime))[1],rep(as.numeric(levels(MS.dat$fTime)),2))
# MS.coefs$time.stamp <- MS.coefs$date + c(0,rep(-0.4,nb.tp),rep(0.4,nb.tp))
# coef.col <- alpha(c('darkgreen',glycol[8],imicol[8]),0.8)
# MS.coefs$col <- c(coef.col[1],rep(coef.col[2],nb.tp),rep(coef.col[3],nb.tp))
# MS.coefs <- MS.coefs[-c(2,nb.tp+2),]
# 
# plot(Estimate~date,data=MS.coefs,type='n',yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l',xlim=range(MS.dat$date),ylim=range(c(MS.coefs$lwr,MS.coefs$upr)))
# title(ylab='model coefficient', cex.lab=1)
# title(xlab="day", cex.lab=1)
# axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
# axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
# abline(h=0,lty=1,col='gray')
# arrows(x0=MS.coefs$time.stamp,y0=MS.coefs$lwr,y1=MS.coefs$upr,length=0,lwd=1.5,col=MS.coefs$col)
# points(x=MS.coefs$time.stamp,y=MS.coefs$Estimate,col=MS.coefs$col,pch=16)
# abline(v=pulse.dates[1:2], lty=3)
# legend(x=6,y=-0.5,bty='n',legend=c('nutrient addition (press treatment)','glyphosate:time','imidacloprid:time'),pch=16,col=coef.col,cex=1,y.intersp=1.2)
# 
# par(mfrow=c(1,1))
# 
# ### gly-treated ponds including 3rd pulse
# 
# sites2rm <- c(18:25,27:32)
# gly.dat <- epdat[!(as.numeric(epdat$site) %in% sites2rm),]
# gly.dat <- droplevels(gly.dat)
# gly.dat$lty <- 1
# 
# cols <- brewer.pal(8,'Reds')
# cols[9] <- 'black'
# 
# gly.dat$lty[gly.dat$site == 'E1' | gly.dat$site == 'H1'] <- 2
# gly.dat$gly.col.idx[gly.dat$site == 'E1' | gly.dat$site == 'H1'] <- 9
# 
# par(mfrow = c(2,1))
# 
# plot(max~date,data = gly.dat, type='n',yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l',ylim=range(gly.dat$max),log='y')
# title(ylab=max~AWCD~(density), cex.lab=0.7)
# axis(2,cex.axis=1,lwd=0,lwd.ticks=1,cex.axis=0.7)
# title(xlab="day", cex.lab=1,cex.lab=0.7)
# axis(1,cex.axis=1,lwd=0,lwd.ticks=1,cex.axis=0.7)
# 
# for (i in 1:34){
#   tmp <- subset(gly.dat, site == levels(gly.dat$site)[i])
#   points(max~date,tmp,type='l',col=alpha(cols[tmp$gly.col.idx[1]],0.8),lwd=seq(1,1.7,0.1)[tmp$gly.col.idx[1]],lty=tmp$lty)
#   points(max~date,tmp,type='p',pch=tmp$pch[1],col=alpha(cols[tmp$gly.col.idx[1]],0.8))
# }
# 
# abline(v=pulse.dates, lty=3)
# text(x=pulse.dates[1],y=range(gly.dat$max)[2]-0.3,'pulse 1',pos=4,cex=0.7)
# text(x=pulse.dates[2],y=range(gly.dat$max)[1]+0.05,'pulse 2',pos=4,cex=0.7)
# text(x=pulse.dates[3],y=range(gly.dat$max)[1]+0.05,'pulse 3',pos=4,cex=0.7)
# 
# plot(use~date,data = gly.dat, type='n',yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l',ylim=range(gly.dat$use))
# title(ylab=substrates~used~(diversity), cex.lab=0.7)
# axis(2,cex.axis=1,lwd=0,lwd.ticks=1,cex.axis=0.7)
# title(xlab="day", cex.lab=1,cex.lab=0.7)
# axis(1,cex.axis=1,lwd=0,lwd.ticks=1,cex.axis=0.7)
# 
# for (i in 1:34){
#   tmp <- subset(gly.dat, site == levels(gly.dat$site)[i])
#   points(use~date,tmp,type='l',col=alpha(cols[tmp$gly.col.idx[1]],0.8),lwd=seq(1,1.7,0.1)[tmp$gly.col.idx[1]],lty=tmp$lty)
#   points(use~date,tmp,type='p',pch=tmp$pch[1],col=alpha(cols[tmp$gly.col.idx[1]],0.8))
# }
# 
# abline(v=pulse.dates, lty=3)
# text(x=pulse.dates[1],y=7,'pulse 1',pos=4,cex=0.7)
# text(x=pulse.dates[2],y=range(gly.dat$use)[1]+0.5,'pulse 2',pos=4,cex=0.7)
# text(x=pulse.dates[3],y=range(gly.dat$use)[1]+0.5,'pulse 3',pos=4,cex=0.7)
# 
# ### ER plots
# 
# sites2rm <- c(17:32)
# ER.dat <- epdat[!(as.numeric(epdat$site) %in% sites2rm),]
# ER.dat <- droplevels(ER.dat)
# p3dat <- ER.dat[ER.dat$date == 49,]
# p3dat$bm.pre.stress <- ER.dat$max[ER.dat$date == 43]
# p3dat$div.pre.stress <- ER.dat$use[ER.dat$date == 43]
# p3dat$gly.conc[p3dat$gly.conc == 0] <- p3dat$gly.conc[p3dat$gly.conc == 0]+3
# p3dat$bstart <- p3dat$bm.pre.stress
# p3dat$dstart <- p3dat$div.pre.stress
# 
# par(mfrow=c(3,2))
# par(mar = c(4,4.2,1,1))
# 
# plot(max~bstart,data = p3dat,type='n',yaxt='n',xaxt='n',ann=F,bty='l',log='xy')
# title(ylab=density~after~phase~2, cex.lab=0.7)
# axis(2,cex.axis=1,lwd=0,lwd.ticks=1, cex.axis=0.7)
# title(xlab=density~before~phase~2,cex.lab=0.7)
# axis(1,cex.axis=1,lwd=0,lwd.ticks=1,cex.axis=0.7)
# points(max~bstart,data = p3dat,pch=p3dat$pch,col=alpha('goldenrod4',0.5),cex=1.5)
# 
# plot(use~bstart,data = p3dat,type='n',yaxt='n',xaxt='n',ann=F,bty='l',log='x')
# title(ylab=diversity~after~phase~2, cex.lab=0.7)
# axis(2,cex.axis=1,lwd=0,lwd.ticks=1, cex.axis=0.7)
# title(xlab=density~before~phase~2,cex.lab=0.7)
# axis(1,cex.axis=1,lwd=0,lwd.ticks=1,cex.axis=0.7)
# points(use~bstart,data = p3dat,pch=p3dat$pch,col=alpha('goldenrod4',0.5),cex=1.5)
# 
# plot(max~gly.conc,data = p3dat,type='n',yaxt='n',xaxt='n',ann=F,bty='l',xlim=c(3,10),log='y')
# title(ylab=density~after~phase~2, cex.lab=0.7)
# axis(2,cex.axis=1,lwd=0,lwd.ticks=1, cex.axis=0.7)
# title(xlab=phase~1~glyphosate~treatment~(mu*g~L^-1~pulse^-1),cex.lab=0.7)
# axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=3:10,labels=c('0',parse(text= paste("italic(e^", 4:10, ")", sep=""))),cex.axis=0.7)
# axis.break(axis=1,breakpos=3.5,bgcol="white",breakcol="black",style="slash",brw=0.02)
# points(max~gly.conc,data = p3dat,pch=p3dat$pch,col=alpha('goldenrod4',0.5),cex=1.5)
# 
# plot(use~gly.conc,data = p3dat,type='n',yaxt='n',xaxt='n',ann=F,bty='l',xlim=c(3,10))
# title(ylab=diversity~after~phase~2, cex.lab=0.7)
# axis(2,cex.axis=1,lwd=0,lwd.ticks=1, cex.axis=0.7)
# title(xlab=phase~1~glyphosate~treatment~(mu*g~L^-1~pulse^-1),cex.lab=0.7)
# axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=3:10,labels=c('0',parse(text= paste("italic(e^", 4:10, ")", sep=""))),cex.axis=0.7)
# axis.break(axis=1,breakpos=3.5,bgcol="white",breakcol="black",style="slash",brw=0.02)
# points(use~gly.conc,data = p3dat,pch=p3dat$pch,col=alpha('goldenrod4',0.5),cex=1.5)
# 
# plot(max~dstart,data = p3dat,type='n',yaxt='n',xaxt='n',ann=F,bty='l',log='y')
# title(ylab=density~after~phase~2, cex.lab=0.7)
# axis(2,cex.axis=1,lwd=0,lwd.ticks=1, cex.axis=0.7)
# title(xlab=diversity~before~phase~2,cex.lab=0.7)
# axis(1,cex.axis=1,lwd=0,lwd.ticks=1,cex.axis=0.7)
# points(max~dstart,data = p3dat,pch=p3dat$pch,col=alpha('goldenrod4',0.5),cex=1.5)
# 
# plot(use~dstart,data = p3dat,type='n',yaxt='n',xaxt='n',ann=F,bty='l')
# title(ylab=diversity~after~phase~2, cex.lab=0.7)
# axis(2,cex.axis=1,lwd=0,lwd.ticks=1, cex.axis=0.7)
# title(xlab=diversity~before~phase~2,cex.lab=0.7)
# axis(1,cex.axis=1,lwd=0,lwd.ticks=1,cex.axis=0.7)
# points(use~dstart,data = p3dat,pch=p3dat$pch,col=alpha('goldenrod4',0.5),cex=1.5)
# 
# 
