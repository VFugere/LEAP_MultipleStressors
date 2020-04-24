## Code to analyze LEAP 2016 MS data
## Vincent Fugère 2019-2020
## This is a raw/rough data exploration script:
## Older code for UQAM LEAP talk in April 2019
## Plus some new raw code for data exploration
## Clean analyses are in script foodweb_paper.R

rm(list=ls())

library(tidyverse)
library(skimr)
library(magrittr)

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
library(mgcViz)

devtools::source_url('https://raw.githubusercontent.com/VFugere/Rfuncs/master/vif.R')
devtools::source_url('https://raw.githubusercontent.com/VFugere/Rfuncs/master/utils.R')

#### load data ####

load('~/Google Drive/Recherche/LEAP Postdoc/2016/MSdata.RData')

#### colour palette & plotting parameters ####

# mauve <- brewer.pal(9,'BuPu')[8]
# jaune <- brewer.pal(11,'RdYlBu')[6]
# vert <- brewer.pal(11,'RdYlGn')[11]
# rouge <- brewer.pal(11,'RdYlBu')[1]
# four.cols <- c(jaune,rouge,mauve,vert)

#three.cols <- c('#058A51','#130CA6','#EB0C00')
#library(wesanderson)
#three.cols <- wes_palette('Zissou1',5)[c(3,1,5)]
#three.cols <- brewer.pal(3,'Dark2')[c(2,3,1)]
#four.cols <- c('gray80',three.cols)

# glycolfunc <- colorRampPalette(c(four.cols[1], four.cols[2]))
# gly.cols <- glycolfunc(8)
# imicolfunc <- colorRampPalette(c(four.cols[1], four.cols[3]))
# imi.cols <- imicolfunc(8)
# bothcolfunc <- colorRampPalette(c(four.cols[1], four.cols[4]))
# both.cols <- bothcolfunc(8)

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

#### data exploration: correlations in the dataset ####

corrgram(merged.data, order=F, lower.panel=panel.shade,
         upper.panel=panel.pie, text.panel=panel.txt,
         main=NULL)

#hard to see because of non-linear effects. Splitting by sampling date
pdf('~/Desktop/corrgrams.pdf',height = 12,width = 12,onefile = T,pointsize = 8)
for(d in Sampling.dates){
  data.sub <- filter(merged.data, date == d) %>% select(-(date:o.nut), -site.f, -date.f)
  data.sub <- data.sub[,apply(data.sub, 2, sd) != 0]
  corrgram(data.sub, order=F, lower.panel=panel.shade,
           upper.panel=panel.pie, text.panel=panel.txt,
           main=paste0('Day ',d))
}
dev.off()

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

##### Time series plot ####

# vars <- c('total','BA','NEP','use','Glycogen')
# var.names <- c(expression(chlorophyll~italic(a)~(mu*g~L^-1)),
#                expression(bacterial~abundance~(cells~mu*L^-1)),
#                expression(net~change~'in'~DO~(Delta*mg~DO~L^-1)),
#                expression(C~substrates~used),
#                expression(glycogen~respiration~(OD~units~above~blanks)))
# var.log <- c(T,T,F,F,F)

vars <- c('BA','greens','NEP','use')
var.names <- c(expression(bacterial~abundance~(cells~mu*L^-1)),
               expression(chlorophyll~italic(a)~(mu*g~L^-1)),
               expression(daytime~Delta*DO~(mg~L^-1)),
               expression(C~substrates~used))
var.log <- c(T,T,F,F)

#pdf('~/Desktop/MSplots.pdf',width=9,height = 4.5,pointsize = 12,onefile = T)
pdf('~/Desktop/MSplots.pdf',width=14,height = 10,pointsize = 12,onefile = T)
par(mfrow=c(4,3),mar=c(2,2,1,1),oma=c(2.5,2.5,1,1),cex=1)
#par(mfrow=c(2,3),mar=c(2,2,1,1),oma=c(2.5,2.5,1,1),cex=1)
pchs.2 <- c(16,15)
#pchs.2 <- c(16,15)

for(v in 1:length(vars)){
  
  tmp <- merged.data[,c(colnames(merged.data)[1:13],vars[v])]
  #tmp$date.idx <- tmp$date
  tmp$date.idx <- tmp$date - 0.6
  tmp$date.idx[tmp$nut.num == 2] <- tmp$date[tmp$nut.num == 2] + 0.6
  for(letter in c('C|D','E|H','J|K')){
  #for(letter in c('C','E','K','D','H','J')){
    if(var.log[v] == T){
      plotfunctions::emptyPlot(xlim=c(1,43),ylim=range(tmp[,14]),yaxt='n',xaxt='n',ann=F, bty='l', log = 'y')
    }else{
      plotfunctions::emptyPlot(xlim=c(1,43),ylim=range(tmp[,14]),yaxt='n',xaxt='n',ann=F, bty='l')
    }
    axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
    axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
    tmp.sub <- filter(tmp, str_detect(site, letter))
    for(i in 1:16){
    #for(i in 1:8){
      pond <- unique(tmp.sub$site)[i]
      sub.sub <- filter(tmp.sub, site == pond)
      points(x=sub.sub$date.idx,y=sub.sub[,14],type='o',lwd=1.2, pch=pchs.2[sub.sub$nut.num], col=alpha(allcols[sub.sub$pond.id],1))
    }
    abline(v=pulse.dates[1:2],lty=3)
    if(letter == 'C|D'){mtext(var.names[v],side=2,line=3)}
  }
  #mtext(var.names[v],side=2,outer=T,line=1)
  #mtext('day of experiment',side=1,outer=T,line=1)
}

mtext('day of experiment',side=1,outer=T,line=1)

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
